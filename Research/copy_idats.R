#!/usr/bin/env R
## Script name: pullRedcap_manual.R
## Purpose: source global functions for copying idat using REDCap and save CSV
## Author: Jonathan Serrano
## Date Created: March 17, 2022
## Copyright (c) NYULH Jonathan Serrano, 2026

# =============================================================================
# REDCap IDAT Puller
#
# Pulls methylation sample records from REDCap for a set of RD-numbers, writes a
# minfi-compatible sample sheet, copies the matching .idat files from the
# research and clinical share drives, and records each sample's array type.
#
# You can provide input in two ways:
#   1. Command-line arguments: either a vector of RD-numbers (RD-12-345 ...) or a
#      single .csv/.xlsx path whose first column lists RD-numbers.
#   2. Edit DEFAULT_INPUT variable below, if not running from command line
#
# Outputs (written into WORK_DIR):
#   samplesheet_og.csv            - the generated sample sheet
#   array_types_sample_sheet.csv  - detected array type per sample with no PHI
#   duplicated_samples.csv        - samples flagged DUPLICATE (if any)
#   samples_missing_sentrix.csv   - records with no Sentrix ID (if any)
#   missing_idats_log.csv / *.txt - copy/lookup failures (if any)
# =============================================================================


# -----------------------------------------------------------------------------
# Configuration / global variables
# -----------------------------------------------------------------------------

# REDCap API access
REDCAP_API_URL   <- "https://redcap.nyumc.org/apps/redcap/api/"
REDCAP_API_TOKEN <- "XXXXXXXXXXXXXXXXXXXX"
COPY_IDATS <- TRUE
SAVE_ARRAY_CSV <- TRUE

# Directory where the sample sheet and copied idats are written
WORK_DIR <- '/Volumes/CBioinformatics/Methylation/Clinical_Runs/pull_redcap_idats'
IDAT_OUTPUT <- file.path(WORK_DIR, "idats")
CURR_DATE <- format(Sys.Date(), "%B_%d_%Y")

# Manual input, used only when no command-line arguments are supplied.
# May be a vector of RD-numbers or a path to a .csv/.xlsx file.
DEFAULT_INPUT <- file.path(WORK_DIR, "IDH_Val_Cases.csv")
# c( "RD-20-123", "RD-21-123", "RD-21-456")

# .idat search locations
RESEARCH_IDAT_DIR <- "/Volumes/snudem01labspace/idats"
CLINICAL_IDAT_DIR <- "/Volumes/molecular/MOLECULAR/iScan"

# Network mount that must be available before running
CLINICAL_SMB <- "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular"

# REDCap fields required to build the sample sheet
REDCAP_FIELDS <- c(
    "record_id", "b_number", "primary_tech", "second_tech", "run_number",
    "barcode_and_row_column", "accession_number", "tm_number", "arrived"
)

# Output file names
OUTPUT_CSV  <- "samplesheet_og.csv"
ARRAY_TYPES_FILE  <- "array_types_sample_sheet.csv"

# Package management
CRAN_PACKAGES <- c(
    "data.table", "openxlsx", "jsonlite", "RCurl", "readxl",
    "stringr", "dplyr", "crayon", "fs", "cli", "httr"
)


# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------

#' Test whether a package is not installed
not_installed <- function(pkg) {!pkg %in% rownames(installed.packages())}

#' Test whether a value is usable (non-null, non-empty, non-NA, non-blank)
is_valid <- function(x) {!is.null(x) && length(x) > 0 && !any(is.na(x)) && all(nzchar(x))}

#' Print a data frame as a single multi-line message
msg_df <- function(dat) message(paste0(capture.output(as.data.frame(dat)), collapse = "\n"))

save_csv <- function(df, file_name) {
    csv_out_path <- file.path(WORK_DIR, file_name)
    message("Saving CSV file:\n", csv_out_path)
    msg_df(df)
    utils::write.csv(df, file = csv_out_path, quote = FALSE, row.names = FALSE)
}


#' Install (if needed) and load the packages required by this script
load_packages <- function() {
    repos <- getOption("repos")
    repos["CRAN"] <- "http://cran.us.r-project.org"
    options(repos = repos)
    
    for (pkg in CRAN_PACKAGES) {
        if (not_installed(pkg)) {
            install.packages(pkg, dependencies = TRUE, ask = FALSE)
        }
    }
    
    # minfi (Bioconductor) is required to read the array type from idats
    if (not_installed("minfi")) {
        if (not_installed("BiocManager")) install.packages("BiocManager")
        BiocManager::install("minfi", update = FALSE, ask = FALSE)
    }
    
    suppressPackageStartupMessages({
        lapply(CRAN_PACKAGES, library, character.only = TRUE, logical.return = TRUE)
    })
    suppressPackageStartupMessages("minfi")
}


#' Verify the required network mount is accessible
check_mounts <- function() {
    if (!dir.exists(CLINICAL_IDAT_DIR)) {
        message("\nPATH does not exist, ensure path is mounted:\n")
        message(crayon::white$bgRed$bold(CLINICAL_IDAT_DIR))
        message("\nYou must mount the network Z-drive path:\n")
        message(crayon::white$bgRed$bold(CLINICAL_SMB), "\n")
        stop("Required network mount is not accessible.")
    }
    message("\n", crayon::bgGreen("Z-drive path is accessible"), "\n")
}


#' Switch into a directory, stopping if it does not exist
set_work_dir <- function(path) {
    if (!dir.exists(path)) stop("Location not found: ", path)
    setwd(path)
}


write_log <- function(info, log_file) {
    message(crayon::bgBlue("~~~Message logged~~~"), "\n", info)
    message(crayon::bgGreen("To file:"), " ", log_file)
    LOG_OUT <-  file.path(WORK_DIR, log_file)
    utils::write.table(
        info, file = LOG_OUT, append = TRUE, quote = FALSE,
        sep = ",", row.names = FALSE, col.names = FALSE
    )
}


# -----------------------------------------------------------------------------
# Input resolution
# -----------------------------------------------------------------------------

#' Read RD-numbers from the first column of a .csv, .tsv, or .xlsx file
#' @return Character vector of values (NAs removed).
parse_input_file <- function(raw_input) {
    message("Input file: ", raw_input)
    file_type <- tolower(tools::file_ext(raw_input))
    message("file_type is: ", file_type)
    values <- NULL
    if (file_type == "xlsx") {
        message("Reading sheet 1 with readxl::read_excel...")
        values <- suppressMessages(
            readxl::read_excel(raw_input, col_names = FALSE, sheet = 1)[[1]]
        )
    }
    if (file_type == "csv") {
        message("Reading with read.delim...")
        values <- read.delim(raw_input, header = FALSE, sep = ",", colClasses = "character",
                             row.names = NULL)[[1]]
    }
    if (file_type == "tsv") {
        message("Reading with read.delim...")
        values <- read.delim(raw_input, header = FALSE, colClasses = "character",
                             row.names = NULL)[[1]]
    }
    stopifnot(!is.null(values))
    values <- as.character(values)
    values <- values[!is.na(values)]
    return(values)
}


#' Get RD-numbers from raw input (a vector, or a single file path)
#'
#' A multi-element input is treated as an explicit list of RD-numbers.
#' A single element is parsed as a file if it exists
#' Anything else is treated as a single RD-number
#'
#' @param raw_input Character vector of RD-numbers, or a single file path.
#' @return Character vector of RD-numbers
resolve_rd_numbers <- function(raw_input) {
    rd_numbers <- NULL
    if (length(raw_input) == 1) {
        if (file.exists(raw_input)) {
            rd_numbers <- parse_input_file(raw_input)
        } else {
            rd_numbers <- raw_input
        }
    }
    
    rd_numbers <- rd_numbers[grepl("^RD-", rd_numbers)]
    rd_numbers <- stringr::str_trim(trimws(rd_numbers))
    
    if (length(rd_numbers) == 0) {
        stop(
            "Your RD-numbers input is not valid!\n",
            "Check that RD-numbers (e.g. RD-26-123) are in the first column of your ",
            "input sheet or passed as arguments:\n",
            paste(raw_input, collapse = ", ")
        )
    }
    message("Input RD-number(s):")
    msg_df(data.frame(rd_numbers))
    return(rd_numbers)
}


# -----------------------------------------------------------------------------
# REDCap
# -----------------------------------------------------------------------------

#' Query REDCap for the given RD-numbers
#'
#' @param rd_numbers Character vector of RD-numbers (record_id values).
#' @param token REDCap API token.
#' @param fields Character vector of fields to export
#' @return A data frame of the exported records
search_redcap <- function(rd_numbers, token, fields = REDCAP_FIELDS) {
    if (!is_valid(token)) stop("You must provide a REDCap API token!")
    
    result <- jsonlite::fromJSON(httr::content(
        httr::POST(
            REDCAP_API_URL,
            body = c(
                list(
                    token = REDCAP_API_TOKEN, content = "record", action = "export",
                    format = "json", type = "flat", rawOrLabel = "raw",
                    exportDataAccessGroups = "false", returnFormat = "json"
                ),
                stats::setNames(as.list(rd_numbers), sprintf("records[%d]", seq_along(rd_numbers) - 1L)),
                stats::setNames(as.list(REDCAP_FIELDS), sprintf("fields[%d]", seq_along(REDCAP_FIELDS) - 1L))
            ),
            encode = "form"
        ),
        as = "text", encoding = "UTF-8"
    ))
    
    missing <- !rd_numbers %in% result$record_id
    
    if (any(missing)) {
        message("Some RD-numbers were not found in REDCap!")
        message(paste0(capture.output(rd_numbers[missing]), collapse = "\n"))
        MISSING_CSV <- paste(CURR_DATE, "input_not_in_redcap_log.csv", sep = "_")
        missing_df <- data.frame(Not_Found = rd_numbers[missing])
        write_log(missing_df, MISSING_CSV)
        message("Check which cases were not found in REDCap in:\n", "redcap_not_found_log.csv")
    }
    
    result_df <- as.data.frame(result)
    return(result_df)
}


# -----------------------------------------------------------------------------
# Sample sheet
# -----------------------------------------------------------------------------

#' Build and write a minfi sample sheet from REDCap records
#'
#' Rows flagged DUPLICATE are split out into duplicated_samples.csv and excluded
#' from the written sample sheet.
#'
#' @param df Data frame of REDCap records
#' @param sentrix_id Two-column data frame: split barcode (ID) and position
write_samplesheet <- function(df, sentrix_id) {
    message(crayon::bgCyan("~~~Writing samplesheet from REDCap records to:"),
            "\n", OUTPUT_CSV)
    
    basenames <- file.path(IDAT_OUTPUT, df$barcode_and_row_column)
    df <- df[!is.na(df[, "barcode_and_row_column"]), , drop = FALSE]
    
    samplesheet <- data.frame(
        Sample_Name      = df[, "record_id"],
        DNA_Number       = df[, "b_number"],
        Sentrix_ID       = sentrix_id[, 1],
        Sentrix_Position = sentrix_id[, 2],
        SentrixID_Pos    = df[, "barcode_and_row_column"],
        Basename         = basenames,
        RunID            = df$run_number,
        MP_num           = df$tm_number,
        Date             = df$arrived
    )
    
    samplesheet <- samplesheet[!is.na(samplesheet$SentrixID_Pos), , drop = FALSE]
    
    is_duplicate <- stringr::str_detect(samplesheet$SentrixID_Pos, "DUPLICATE")
    is_duplicate[is.na(is_duplicate)] <- FALSE
    
    if (any(is_duplicate)) {
        message("Dropping duplicated samples!!")
        duplicated_csv <- samplesheet[is_duplicate, ]
        save_csv(duplicated_csv, "duplicated_samples.csv")
    }
    
    samplesheet <- samplesheet[!is_duplicate, ]
    save_csv(samplesheet, OUTPUT_CSV)
}


# -----------------------------------------------------------------------------
# IDAT copying
# -----------------------------------------------------------------------------

#' Copy idat files into the current working directory
#' @param files Character vector of full paths to .idat files.
copy_idat_files <- function(idat_paths) {
    
    report_copy_status <- function(paths) {
        copied <- basename(paths)
        copied <- copied[copied != ""]
        success <- file.exists(copied)
        message(".idat idat_paths that failed to copy:")
        if (all(success)) cat("none", "\n") else print(copied[!success])
        invisible(all(success))
    }
    
    readable <- fs::file_access(idat_paths, mode = "read")
    
    if (any(!readable)) {
        info <- paste("Cannot read idat file:", idat_paths[!readable], collapse = "\n")
        write_log(info, "read_error_idat.csv")
        idat_paths <- idat_paths[readable]
    }
    
    if (length(idat_paths) == 0) {
        report_copy_status(idat_paths)
        return(invisible(NULL))
    }
    
    cli::cli_progress_bar("Copying idat_paths", total = length(idat_paths))
    for (fi in idat_paths) {
        tryCatch(
            fs::file_copy(fi, file.path(IDAT_OUTPUT, basename(fi)), overwrite = TRUE),
            error = function(e) {
                info <- paste("Failed to copy:", fi)
                cli::cli_alert_danger(info)
                write_log(info, "missing_idat_files.csv")
            }
        )
        cli::cli_progress_update()
    }
    cli::cli_progress_done()
    
    report_copy_status(idat_paths)
}


require_mount <- function(idat_dir) stopifnot(dir.exists(idat_dir))

idat_bases_from_files <- function(idat_paths) {
    if (length(idat_paths) == 0) return(character(0))
    
    parts <- stringr::str_split_fixed(basename(idat_paths), "_", 3)
    return(unique(paste0(parts[, 1], "_", parts[, 2])))
}

idats_complete <- function(idat_paths, bases_needed) {
    expected <- length(unique(bases_needed)) * 2
    actual <- length(unique(basename(idat_paths)))
    return(expected == actual)
}

log_missing_idats <- function(idat_paths, bases_needed) {
    if (idats_complete(idat_paths, bases_needed)) return(invisible(NULL))
    
    message(crayon::bgRed("Still missing idat files not in External folder:"))
    
    bases_found <- idat_bases_from_files(idat_paths)
    missing_samples <- bases_needed[!(bases_needed %in% bases_found)]
    
    message("The following samples are missing:")
    msg_df(missing_samples)
    
    save_csv(data.frame(Missing_Samples = missing_samples), "missing_idats_log.csv")
    message(
        crayon::bgRed("Check the log file to see which idats were not found:"),
        " missing_idats_log.csv"
    )
    
    return(invisible(NULL))
}

find_external_idats <- function(missing_bases) {
    external_idat_dir <- file.path(RESEARCH_IDAT_DIR, "External")
    
    message(crayon::bgRed("The following idats are missing:"))
    msg_df(missing_bases)
    message(crayon::bgGreen("Searching the External folder for more idats..."))
    
    red_green_files <- paste0(
        rep(missing_bases, each = 2), c("_Grn.idat", "_Red.idat")
    )
    direct_idats <- file.path(external_idat_dir, red_green_files)
    
    if (all(file.exists(direct_idats))) return(direct_idats)
    
    other_idats <- dir(
        external_idat_dir, pattern = ".idat",
        full.names = TRUE, recursive = TRUE
    )
    found <- stringr::str_detect(
        other_idats, pattern = paste(missing_bases, collapse = "|")
    )
    
    if (!any(found)) {
        message(crayon::bgRed("Still missing idat files not in External folder:"))
        msg_df(missing_bases)
        return(NULL)
    }
    
    message(
        crayon::bgGreen("Found extra idats in External folder:"),
        " ", external_idat_dir
    )
    
    idats_to_add <- other_idats[found]
    msg_df(idats_to_add)
    
    return(idats_to_add)
}

add_external_idats <- function(idat_paths, ssheet, bases_needed) {
    external_idat_dir <- file.path(RESEARCH_IDAT_DIR, "External")
    
    if (idats_complete(idat_paths, bases_needed)) {
        message("All idats detected in folders!")
        return(idat_paths)
    }
    
    message(
        crayon::bgRed("Still missing some idats! Checking External Folder:"),
        " ", external_idat_dir
    )
    
    bases_found <- idat_bases_from_files(idat_paths)
    still_missing <- !(bases_needed %in% bases_found)
    
    if (!any(still_missing)) return(idat_paths)
    
    message("Missing idats:")
    msg_df(ssheet[still_missing, , drop = FALSE])
    
    idats_to_add <- find_external_idats(bases_needed[still_missing])
    
    if (length(idat_paths) > 0) {
        idat_paths <- unique(c(idat_paths, setdiff(idats_to_add, idat_paths)))
        log_missing_idats(idat_paths, bases_needed)
    } else {
        idat_paths <- idats_to_add
    }
    
    return(idat_paths)
}


#' Locate idat files for a sample sheet and copy any that are missing locally
#'
#' Searches the research and clinical idat drives for each sample's Grn/Red
#' files, falls back to the research "External" folder for anything still
#' missing, logs unresolved samples, and copies whatever was found into the run
#' directory.
#'
#' @param samplesheet_file Path to the sample sheet CSV.
#' @return Invisibly the vector of resolved idat paths.
resolve_and_copy_idats <- function(samplesheet_file = OUTPUT_CSV) {
    require_mount(RESEARCH_IDAT_DIR)
    require_mount(CLINICAL_IDAT_DIR)
    
    if (!file.exists(samplesheet_file)) {
        message("Cannot find your sheet named:", samplesheet_file)
        stopifnot(file.exists(samplesheet_file))
    }
    
    ssheet <- utils::read.csv(samplesheet_file, strip.white = TRUE)
    barcode <- as.vector(ssheet$Sentrix_ID)
    sentrix_pos <- ssheet$SentrixID_Pos
    bases_needed <- as.vector(ssheet$SentrixID_Pos)
    
    all_fi <- character(0)
    
    for (idat_dir in c(RESEARCH_IDAT_DIR, CLINICAL_IDAT_DIR)) {
        dir_names <- file.path(idat_dir, barcode)
        green_files <- file.path(dir_names, paste0(sentrix_pos, "_Grn.idat"))
        red_files <- file.path(dir_names, paste0(sentrix_pos, "_Red.idat"))
        all_fi <- c(all_fi, green_files, red_files)
    }
    
    all_fi <- all_fi[fs::file_exists(all_fi)]
    
    if (length(all_fi) == 0) {
        all_fi <- add_external_idats(all_fi, ssheet, bases_needed)
    }
    
    if (length(all_fi) == 0) {
        warning(crayon::bgRed("No .idat files found!"))
        message(
            "Check worksheet for typos and if the barcode folder exists in the search path(s):"
        )
        message(RESEARCH_IDAT_DIR, "\nor\n", CLINICAL_IDAT_DIR)
        stop(crayon::bgRed(paste(
            "No .idat files found for these sample(s)!",
            "The case(s) may have not been run yet."
        )))
    }
    
    message("Files found: ")
    msg_df(all_fi)
    
    all_fi <- add_external_idats(all_fi, ssheet, bases_needed)
    
    message("Checking if idats exist in run directory...")
    
    current_idats <- basename(
        dir(IDAT_OUTPUT, pattern = "\\.idat$", recursive = FALSE)
    )
    idats_found <- basename(all_fi) %in% current_idats
    
    if (all(idats_found)) {
        message(".idat files already copied to run directory")
    } else {
        copy_idat_files(all_fi[!idats_found])
    }
    
    return(invisible(all_fi))
}


# -----------------------------------------------------------------------------
# Array type detection
# -----------------------------------------------------------------------------

#' Detect and record the array type for each sample
#'
#' Reads each sample's idats with minfi and writes the array annotation to
#' array_types_sample_sheet.csv. Rows marked "NO IDAT FILE" are skipped.
#'
#' @param targets Data frame read from the sample sheet.
save_array_types <- function(targets) {
    targets$ArrayType <- ""
    targets <- targets[!grepl("NO IDAT FILE", targets$SentrixID_Pos), , drop = FALSE]
    
    for (idx in seq_len(nrow(targets))) {
        current_sam <- targets[idx, ]
        rg_set <- minfi::read.metharray.exp(targets = current_sam,
                                            force = TRUE, verbose = TRUE)
        targets$ArrayType[idx] <- rg_set@annotation[["array"]]
    }
    
    targets_df <- targets[, c("Sample_Name", "SentrixID_Pos", "ArrayType")]
    save_csv(targets_df, ARRAY_TYPES_FILE)
}


# -----------------------------------------------------------------------------
# Orchestration
# -----------------------------------------------------------------------------

#' Pull records, write the sample sheet, copy idats, and record array types
#' @param rd_numbers Character vector of RD-numbers.
#' @param token REDCap API token.
pull_redcap_idats <- function(rd_numbers, token) {
    
    request <- list(token = token, content = "version")
    
    is_token_valid <- tryCatch(
        httr::POST(REDCAP_API_URL, body = request, encode = "form"),
        error = function(e) NULL
    )
    
    if (is.null(is_token_valid) || httr::http_error(is_token_valid)) {
        stop("Your REDCap API Token is invalid: ", token)
    }
    
    if (!dir.exists(IDAT_OUTPUT)) {dir.create(IDAT_OUTPUT, recursive = TRUE)}
    stopifnot(length(rd_numbers) > 0)
    
    records_found <- search_redcap(rd_numbers, token)
    
    # Records without a Sentrix ID cannot be processed; log and drop them
    missing_sentrix <- is.na(records_found$barcode_and_row_column) | records_found$barcode_and_row_column == ""
    if (any(missing_sentrix)) {
        message("Some samples have no SentrixID and will be dropped!")
        dropped <- records_found[missing_sentrix, 1]
        save_csv(dropped, "samples_missing_sentrix.csv")
        records_found <- records_found[!missing_sentrix, , drop = FALSE]
    }
    
    sentrix_id <- as.data.frame(
        stringr::str_split_fixed(records_found[, "barcode_and_row_column"], "_", 2)
    )
    
    if (nrow(sentrix_id) == 0) {
        message("Input cases have not been run or do not have Sentrix ID in REDCap:")
        message(paste(capture.output(records_found), collapse = "\n"))
        stopifnot(nrow(sentrix_id) > 0)
    }
    
    write_samplesheet(df = records_found, sentrix_id = sentrix_id)
    
    if (COPY_IDATS == TRUE) {
        Sys.sleep(5)
        resolve_and_copy_idats(samplesheet_file = OUTPUT_CSV)
    }
    
    if (SAVE_ARRAY_CSV == TRUE) {
        targets <- utils::read.csv(OUTPUT_CSV, strip.white = TRUE, row.names = NULL)
        save_array_types(targets)
    }
}


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

main <- function(){
    load_packages()
    check_mounts()
    set_work_dir(WORK_DIR)
    
    # Command-line arguments take priority over DEFAULT_INPUT
    cli_args  <- commandArgs(trailingOnly = TRUE)
    raw_input <- if (length(cli_args) > 0) cli_args else DEFAULT_INPUT
    rd_numbers <- resolve_rd_numbers(raw_input)
    pull_redcap_idats(rd_numbers = rd_numbers, token = REDCAP_API_TOKEN)
}

main()
