#!/usr/bin/env R
## Script name: copy_idats.R
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
#   2. The DEFAULT_INPUT variable below, used when no arguments are passed.
#
# Outputs (written into OUTPUT_DIR):
#   samplesheet_og.csv            - the generated sample sheet
#   array_types_sample_sheet.csv  - detected array type per sample
#   duplicated_samples.csv        - samples flagged DUPLICATE (if any)
#   samples_missing_sentrix.csv   - records with no Sentrix ID (if any)
#   missing_idats_log.csv / *.txt - copy/lookup failures (if any)
# =============================================================================


# -----------------------------------------------------------------------------
# Configuration / global variables
# -----------------------------------------------------------------------------

# REDCap API access
REDCAP_API_TOKEN <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXX"
REDCAP_API_URL   <- "https://redcap.nyumc.org/apps/redcap/api/"


# Default Directory where the sample sheet and copied idats are written
OUTPUT_DIR <- "/Volumes/CBioinformatics/Methylation/Clinical_Runs/pull_redcap_idats"

# Manual input, used only when no command-line arguments are supplied.
# May be a vector of RD-numbers or a path to a .csv/.xlsx file.
DEFAULT_INPUT <- "/path/to/your/sample_list.xlsx"

# .idat search locations
RESEARCH_IDAT_DIR <- "/Volumes/snudem01labspace/idats"
CLINICAL_IDAT_DIR <- "/Volumes/molecular/MOLECULAR/iScan"

# Network mount that must be available before running
REQUIRED_MOUNT     <- "/Volumes/molecular/MOLECULAR LAB ONLY"
REQUIRED_MOUNT_SMB <- "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular"

# REDCap fields required to build the sample sheet
REDCAP_FIELDS <- c(
    "record_id", "b_number", "primary_tech", "second_tech", "run_number",
    "barcode_and_row_column", "accession_number", "tm_number", "arrived"
)

# Output file names
SAMPLESHEET_FILE  <- "samplesheet_og.csv"
ARRAY_TYPES_FILE  <- "array_types_sample_sheet.csv"

# Package management
CRAN_PACKAGES <- c(
    "data.table", "openxlsx", "jsonlite", "RCurl", "readxl",
    "stringr", "dplyr", "crayon", "fs", "cli"
)
MIN_REDCAPAPI_VERSION <- "2.7.4"

# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------

#' Test whether a package is not installed
#'
#' @param pkg Character scalar package name.
#' @return TRUE if the package is not installed, FALSE otherwise.
not_installed <- function(pkg) {!pkg %in% rownames(installed.packages())}


#' Test whether a value is usable (non-null, non-empty, non-NA, non-blank)
#'
#' @param x Value to validate.
#' @return TRUE if x is a usable, non-empty value; FALSE otherwise.
is_valid <- function(x) {!is.null(x) && length(x) > 0 && !any(is.na(x)) && all(nzchar(x))}


#' Print a data frame as a single multi-line message
#'
#' @param dat Object coercible to a data frame.
#' @return Invisibly NULL (called for the side effect of messaging).
msg_table <- function(dat) message(paste0(capture.output(as.data.frame(dat)), collapse = "\n"))


# -----------------------------------------------------------------------------
# Environment setup
# -----------------------------------------------------------------------------

#' Install (if needed) and load the packages required by this script
#'
#' Installs missing CRAN packages, ensures a recent enough redcapAPI, installs
#' minfi from Bioconductor if absent, and attaches the packages used directly.
#'
#' @return Invisibly NULL.
load_packages <- function() {
    repos <- getOption("repos")
    repos["CRAN"] <- "http://cran.us.r-project.org"
    options(repos = repos)
    
    for (pkg in CRAN_PACKAGES) {
        if (not_installed(pkg)) {
            install.packages(pkg, dependencies = TRUE)
        }
    }
    
    # redcapAPI must be recent enough to provide exportRecordsTyped()
    if (not_installed("redcapAPI") ||
        utils::packageVersion("redcapAPI") < MIN_REDCAPAPI_VERSION) {
        install.packages("redcapAPI", dependencies = TRUE)
    }
    
    # minfi (Bioconductor) is required to read the array type from idats
    if (not_installed("minfi")) {
        if (not_installed("BiocManager")) install.packages("BiocManager")
        BiocManager::install("minfi", update = FALSE, ask = FALSE)
    }
    
    suppressPackageStartupMessages({
        library("redcapAPI")
        library("dplyr")
    })
}


#' Verify the required network mount is accessible
#'
#' Stops execution with a clear message if the Z-drive mount is missing.
#'
#' @return Invisibly NULL.
check_mounts <- function() {
    if (!dir.exists(REQUIRED_MOUNT)) {
        cat("\nPATH does not exist, ensure path is mounted:\n")
        cat(crayon::white$bgRed$bold(REQUIRED_MOUNT))
        cat("\nYou must mount the network Z-drive path:\n")
        cat(crayon::white$bgRed$bold(REQUIRED_MOUNT_SMB), "\n")
        stop("Required network mount is not accessible.")
    }
    message("\n", crayon::bgGreen("Z-drive path is accessible"), "\n")
}


#' Switch into a directory, stopping if it does not exist
#'
#' @param path Directory to switch into.
#' @return Invisibly the path.
set_working_dir <- function(path) {
    if (!dir.exists(path)) stop("Location not found: ", path)
    setwd(path)
}


# -----------------------------------------------------------------------------
# Input resolution
# -----------------------------------------------------------------------------

#' Read RD-numbers from the first column of a .csv or .xlsx file
#'
#' @param path Path to a .csv or .xlsx file.
#' @return Character vector of values (NAs removed).
read_rd_numbers <- function(path) {
    message("input file: ", path)
    if (endsWith(path, ".xlsx")) {
        message("FileType is .xlsx, reading with readxl::read_excel...")
        values <- suppressMessages(
            readxl::read_excel(path, col_names = FALSE, sheet = 1)[[1]]
        )
    } else {
        message("FileType is .csv, reading with read.delim...")
        values <- read.delim(
            path, header = FALSE, sep = ",",
            colClasses = "character", row.names = NULL
        )[[1]]
    }
    values <- as.character(values)
    values[!is.na(values)]
}


#' Resolve RD-numbers from raw input (a vector, or a single file path)
#'
#' A multi-element input is treated as an explicit list of RD-numbers. A single
#' element ending in .csv/.xlsx is read from disk; anything else is treated as a
#' single RD-number.
#'
#' @param input Character vector of RD-numbers, or a single file path.
#' @return Character vector of RD-numbers.
resolve_rd_numbers <- function(input) {
    rd_numbers <- NULL
    if (length(input) > 1) {
        rd_numbers <- as.character(input)
    }
    if (endsWith(input, ".csv") || endsWith(input, ".xlsx")) {
        rd_numbers <- read_rd_numbers(input)
    }
    
    if (is.null(rd_numbers)) {
        rd_numbers <- as.character(input)
    }
    
    rd_numbers <- rd_numbers[grepl("^RD-", rd_numbers)]
    rd_numbers <- trimws(rd_numbers)
    rd_numbers <- stringr::str_trim(rd_numbers)
    
    if (length(rd_numbers) == 0) {
        stop(
            "Your RD-numbers input is not valid!\n",
            "Check that RD-numbers (e.g. RD-26-123) are in the first column of your ",
            "input sheet or passed as arguments:\n",
            paste(input, collapse = ", ")
        )
    }
    message("Input RD-number(s):\n", paste0(capture.output(rd_numbers), collapse = "\n"))
    return(rd_numbers)
}


# -----------------------------------------------------------------------------
# REDCap
# -----------------------------------------------------------------------------

#' Query REDCap for the given RD-numbers
#'
#' @param rd_numbers Character vector of RD-numbers (record_id values).
#' @param token REDCap API token.
#' @param fields Character vector of fields to export.
#' @return A data frame of the exported records.
search_redcap <- function(rd_numbers, token, fields = REDCAP_FIELDS) {
    if (!is_valid(token)) {
        stop("You must provide a REDCap API token!")
    }
    
    rcon <- redcapAPI::redcapConnection(REDCAP_API_URL, token)
    result <- redcapAPI::exportRecordsTyped(
        rcon,
        records            = rd_numbers,
        fields             = fields,
        dag                = FALSE,
        factors            = FALSE,
        form_complete_auto = FALSE,
        format             = "csv"
    )
    
    missing <- !rd_numbers %in% result$record_id
    if (any(missing)) {
        message("Some RD-numbers were not found in REDCap!")
        message(paste0(capture.output(rd_numbers[missing]), collapse = "\n"))
    }
    
    as.data.frame(result)
}


# -----------------------------------------------------------------------------
# Sample sheet
# -----------------------------------------------------------------------------

#' Build and write a minfi sample sheet from REDCap records
#'
#' Rows flagged DUPLICATE are split out into duplicated_samples.csv and excluded
#' from the written sample sheet.
#'
#' @param df Data frame of REDCap records.
#' @param sentrix_id Two-column data frame: split barcode (ID) and position.
#' @param output_file Path for the generated sample sheet.
#' @return Invisibly NULL.
write_samplesheet <- function(df, sentrix_id, output_file = SAMPLESHEET_FILE) {
    message(crayon::bgCyan("~~~Writing sample sheet from REDCap records to:"), "\n", output_file)
    
    basenames <- file.path(getwd(), df$barcode_and_row_column)
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
        duplicated_csv <- samplesheet[is_duplicate, , drop = FALSE]
        message(paste0(capture.output(duplicated_csv), collapse = "\n"))
        utils::write.csv(
            duplicated_csv, file = "duplicated_samples.csv",
            quote = FALSE, row.names = FALSE
        )
    }
    
    samplesheet <- samplesheet[!is_duplicate, , drop = FALSE]
    utils::write.csv(samplesheet, file = output_file, quote = FALSE, row.names = FALSE)
}


# -----------------------------------------------------------------------------
# IDAT copying
# -----------------------------------------------------------------------------

#' Copy idat files into the current working directory
#'
#' Unreadable files are logged to read_error_idat.txt and skipped; copy failures
#' are logged to missing_idat_files.txt.
#'
#' @param files Character vector of full paths to .idat files.
#' @return Invisibly NULL.
copy_idat_files <- function(files) {
    
    write_log <- function(info, log_file) {
        message(crayon::bgBlue("~~~Message logged~~~"), "\n", info)
        message(crayon::bgGreen("To file:"), " ", log_file)
        utils::write.table(
            info, file = log_file, append = TRUE, quote = FALSE,
            sep = "\t", row.names = FALSE, col.names = FALSE
        )
    }
    
    report_copy_status <- function(paths) {
        copied <- basename(paths)
        copied <- copied[copied != ""]
        success <- file.exists(copied)
        message(".idat files that failed to copy:")
        if (all(success)) cat("none", "\n") else print(copied[!success])
        invisible(all(success))
    }
    
    readable <- fs::file_access(files, mode = "read")
    if (any(!readable)) {
        info <- paste("Cannot read idat file:", files[!readable], collapse = "\n")
        write_log(info, "read_error_idat.txt")
        files <- files[readable]
    }
    
    if (length(files) == 0) {
        report_copy_status(files)
        return(invisible(NULL))
    }
    
    cli::cli_progress_bar("Copying files", total = length(files))
    for (fi in files) {
        tryCatch(
            fs::file_copy(fi, file.path(getwd(), basename(fi)), overwrite = TRUE),
            error = function(e) {
                info <- paste("Failed to copy:", fi)
                cli::cli_alert_danger(info)
                write_log(info, "missing_idat_files.txt")
            }
        )
        cli::cli_progress_update()
    }
    cli::cli_progress_done()
    
    report_copy_status(files)
}


#' Locate idat files for a sample sheet and copy any that are missing locally
#'
#' Searches the research and clinical idat drives for each sample's Grn/Red
#' files, falls back to the research "External" folder for anything still
#' missing, logs unresolved samples, and copies whatever was found into the run
#' directory.
#'
#' @param samplesheet_file Path to the sample sheet CSV.
#' @param run_dir Directory checked for already-present idats (defaults to wd).
#' @return Invisibly the vector of resolved idat paths.
resolve_and_copy_idats <- function(samplesheet_file = SAMPLESHEET_FILE,
                                   run_dir = NULL) {
    external_idat_dir <- file.path(RESEARCH_IDAT_DIR, "External")
    
    require_mount <- function(idat_dir) {
        if (!dir.exists(idat_dir)) {
            message(
                crayon::bgRed("Share drive not found, ensure path is accessible:"),
                "\n", idat_dir
            )
            stopifnot(dir.exists(idat_dir))
        }
    }
    
    # Derive the unique "barcode_position" base names from full idat paths
    idat_bases_from_files <- function(files) {
        if (length(files) == 0) return(character(0))
        parts <- stringr::str_split_fixed(basename(files), "_", 3)
        unique(paste0(parts[, 1], "_", parts[, 2]))
    }
    
    log_missing_idats <- function(files, bases_needed) {
        expected <- length(unique(bases_needed)) * 2
        actual   <- length(unique(basename(files)))
        if (expected == actual) return(invisible(NULL))
        
        message(crayon::bgRed("Still missing idat files not in External folder:"))
        bases_found <- idat_bases_from_files(files)
        missing_samples <- bases_needed[!(bases_needed %in% bases_found)]
        
        message("The following samples are missing:")
        msg_table(missing_samples)
        
        utils::write.csv(
            data.frame(Missing_Samples = missing_samples),
            "missing_idats_log.csv", row.names = FALSE, quote = FALSE
        )
        message(
            crayon::bgRed("Check the log file to see which idats were not found:"),
            " ", "missing_idats_log.csv"
        )
    }
    
    find_external_idats <- function(missing_bases) {
        message(crayon::bgRed("The following idats are missing:"))
        msg_table(missing_bases)
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
        search_pattern <- paste(missing_bases, collapse = "|")
        found <- stringr::str_detect(other_idats, pattern = search_pattern)
        
        if (!any(found)) {
            message(crayon::bgRed("Still missing idat files not in External folder:"))
            msg_table(missing_bases)
            return(NULL)
        }
        
        message(
            crayon::bgGreen("Found extra idats in External folder:"),
            " ", external_idat_dir
        )
        idats_to_add <- other_idats[found]
        msg_table(idats_to_add)
        idats_to_add
    }
    
    add_external_idats <- function(files, ssheet, bases_needed) {
        expected <- length(unique(bases_needed)) * 2
        actual   <- length(unique(basename(files)))
        
        if (expected == actual) {
            message("All idats detected in folders!")
            return(files)
        }
        
        message(
            crayon::bgRed("Still missing some idats! Checking External Folder:"),
            " ", external_idat_dir
        )
        
        if (length(files) > 0) {
            bases_found  <- idat_bases_from_files(files)
            still_missing <- !(bases_needed %in% bases_found)
        } else {
            still_missing <- rep(TRUE, length(bases_needed))
            files <- NULL
        }
        
        if (any(still_missing)) {
            message("Missing idats:")
            msg_table(ssheet[still_missing, , drop = FALSE])
            
            idats_to_add <- find_external_idats(bases_needed[still_missing])
            
            if (!is.null(files)) {
                files <- unique(c(files, setdiff(idats_to_add, files)))
                log_missing_idats(files, bases_needed)
            } else {
                files <- idats_to_add
            }
        }
        
        files
    }
    
    require_mount(RESEARCH_IDAT_DIR)
    require_mount(CLINICAL_IDAT_DIR)
    
    if (is.null(run_dir)) run_dir <- getwd()
    
    if (!file.exists(samplesheet_file)) {
        message("Cannot find your sheet named:", samplesheet_file)
        stopifnot(file.exists(samplesheet_file))
    }
    
    ssheet      <- utils::read.csv(samplesheet_file, strip.white = TRUE)
    barcode     <- as.vector(ssheet$Sentrix_ID)
    sentrix_pos <- ssheet$SentrixID_Pos
    
    # Build candidate Grn/Red paths across both idat drives
    all_fi <- character(0)
    for (idat_dir in c(RESEARCH_IDAT_DIR, CLINICAL_IDAT_DIR)) {
        dir_names   <- file.path(idat_dir, barcode)
        green_files <- file.path(dir_names, paste0(sentrix_pos, "_Grn.idat"))
        red_files   <- file.path(dir_names, paste0(sentrix_pos, "_Red.idat"))
        all_fi <- c(all_fi, green_files, red_files)
    }
    
    all_fi       <- all_fi[fs::file_exists(all_fi)]
    bases_needed <- as.vector(ssheet$SentrixID_Pos)
    
    if (length(all_fi) == 0) {
        all_fi <- add_external_idats(all_fi, ssheet, bases_needed)
    }
    
    if (length(all_fi) == 0) {
        warning(crayon::bgRed("No .idat files found!"))
        message("Check worksheet for typos and if the barcode folder exists in the search path(s):")
        message(RESEARCH_IDAT_DIR, "\nor\n", CLINICAL_IDAT_DIR)
        stop(crayon::bgRed(paste0(
            "No .idat files found for these sample(s)!  ",
            "The case(s) may have not been run yet."
        )))
    }
    
    message("Files found: ")
    msg_table(all_fi)
    
    all_fi <- add_external_idats(all_fi, ssheet, bases_needed)
    
    message("Checking if idats exist in run directory...")
    current_idats <- basename(dir(path = run_dir, pattern = "\\.idat$", recursive = FALSE))
    idats_found   <- basename(all_fi) %in% current_idats
    
    if (all(idats_found)) {
        message(".idat files already copied to run directory")
    } else {
        copy_idat_files(all_fi[!idats_found])
    }
    
    invisible(all_fi)
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
#' @return Invisibly NULL.
save_array_types <- function(targets) {
    targets$ArrayType <- ""
    targets <- targets[!grepl("NO IDAT FILE", targets$SentrixID_Pos), , drop = FALSE]
    
    for (idx in seq_len(nrow(targets))) {
        rg_set <- minfi::read.metharray.exp(
            targets = targets[idx, , drop = FALSE],
            force = TRUE, verbose = TRUE
        )
        targets$ArrayType[idx] <- rg_set@annotation[["array"]]
    }
    
    out_file <- file.path(getwd(), ARRAY_TYPES_FILE)
    message("Saving file:\n", out_file)
    utils::write.csv(
        targets[, c("Sample_Name", "SentrixID_Pos", "ArrayType")],
        file = out_file, quote = FALSE, row.names = FALSE
    )
}


# -----------------------------------------------------------------------------
# Orchestration
# -----------------------------------------------------------------------------

#' Pull records, write the sample sheet, copy idats, and record array types
#'
#' @param rd_numbers Character vector of RD-numbers.
#' @param token REDCap API token.
#' @param copy_idats Whether to locate and copy idat files.
#' @param output_file Path for the generated sample sheet.
#' @param run_dir Directory checked for already-present idats (defaults to wd).
#' @return Invisibly NULL.
pull_redcap_idats <- function(rd_numbers,
                              token,
                              copy_idats  = TRUE,
                              output_file = SAMPLESHEET_FILE,
                              run_dir     = NULL) {
    if (is.null(run_dir)) run_dir <- file.path(getwd(), "idats")
    stopifnot(length(rd_numbers) > 0)
    
    records <- search_redcap(rd_numbers, token)
    
    # Records without a Sentrix ID cannot be processed; log and drop them
    missing_sentrix <- is.na(records$barcode_and_row_column)
    if (any(missing_sentrix)) {
        message("Some samples have no SentrixID and will be dropped:")
        dropped <- records[missing_sentrix, 1]
        message(paste0(capture.output(dropped), collapse = "\n"))
        message("Saving list to file: \"samples_missing_sentrix.csv\"")
        utils::write.csv(
            dropped, "samples_missing_sentrix.csv",
            quote = FALSE, row.names = FALSE
        )
    }
    records <- records[!missing_sentrix, , drop = FALSE]
    
    sentrix_id <- as.data.frame(
        stringr::str_split_fixed(records[, "barcode_and_row_column"], "_", 2)
    )
    if (nrow(sentrix_id) == 0) {
        message("Input cases have not been run or do not have Sentrix ID in REDCap:")
        message(paste(capture.output(records), collapse = "\n"))
        stopifnot(nrow(sentrix_id) > 0)
    }
    
    write_samplesheet(df = records, sentrix_id = sentrix_id, output_file = output_file)
    
    if (copy_idats) {
        Sys.sleep(5)
        resolve_and_copy_idats(samplesheet_file = output_file, run_dir = run_dir)
    }
    
    targets <- utils::read.csv(output_file, strip.white = TRUE, row.names = NULL)
    save_array_types(targets)
}


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

load_packages()
check_mounts()
set_working_dir(OUTPUT_DIR)

# Command-line arguments take priority over DEFAULT_INPUT
cli_args  <- commandArgs(trailingOnly = TRUE)
raw_input <- if (length(cli_args) > 0) cli_args else DEFAULT_INPUT
rd_numbers <- resolve_rd_numbers(raw_input)
pull_redcap_idats(rd_numbers = rd_numbers, token = REDCAP_API_TOKEN)

