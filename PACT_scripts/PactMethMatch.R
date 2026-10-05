#!/usr/bin/env Rscript
## Script name: PactMethMatch.R
## Purpose: search REDCap for PACT samples with methylation & generate cnv PNG
## Date Created: September 2, 2021
## Version: 1.1.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2026

options(stringsAsFactors = FALSE, repos = c(CRAN = "https://cran.r-project.org"))

# Input ------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || anyNA(args[1]) || any(!nzchar(args[1]))) {
    stop("Usage: PactMethMatch.R <PACT ID or input file>")
}

PACT_INPUT <- args[[1]]

TOKEN_PATH <- "/Volumes/CBioinformatics/scripts/METH_DB_API.txt"
stopifnot(file.exists(TOKEN_PATH))
REDCAP_API_TOKEN <- trimws(readLines(TOKEN_PATH, n = 1, warn = FALSE))

read_flag <- grepl("\\.csv$", PACT_INPUT, ignore.case = TRUE)
is_sophia <- grepl("^[0-9]{2}", PACT_INPUT)
is_file_path <- grepl("/", PACT_INPUT, fixed = TRUE)

# Configuration ----------------------------------------------------------------
REDCAP_URL <- "https://redcap.nyumc.org/apps/redcap/api/"

classifier_install <- "/Volumes/CBioinformatics/Methylation/Rscripts/install_epic_v2_classifier.R"

# Returns volume path, or macOS duplicate mount such as "/Volumes/molecular-1"
fix_volume <- function(vol) {
    if (dir.exists(vol)) return(vol)
    alts <- Sys.glob(paste0(vol, "-[0-9]*"))
    if (length(alts) > 0) alts[1] else vol
}

mol_drive <- fix_volume("/Volumes/molecular")
research_vol <- fix_volume("/Volumes/snudem01labspace")

research_idat_dir <- file.path(research_vol,"idats")
clinical_idat_dir <- "/Volumes/molecular/MOLECULAR/iScan"

lab_drive <- file.path(mol_drive,"MOLECULAR LAB ONLY")

cnv_out_dir <- file.path(mol_drive, "Molecular/MethylationClassifier/CNV_PNG")
pact_data_dir <- file.path(lab_drive, "NYU PACT Patient Data")


smb_share <- "smb://shares-cifs.nyumc.org/apps/acc_pathology"
report_share <- file.path(smb_share, "molecular/Molecular/MethylationClassifier")
desktop <- path.expand("~/Desktop")
match_tsv <- file.path(desktop, paste0(basename(PACT_INPUT), "_match_log.tsv"))

main_pkgs <- c(
    "data.table", "openxlsx", "jsonlite", "readxl", "stringr",
    "tidyverse", "crayon", "tinytex", "systemfonts", "remotes",
    "dplyr", "fs", "httr", "cli"
)
brew_pkgs <- c("gcc", "llvm", "lld", "open-mpi", "pkgconf", "gdal", "proj", "apache-arrow")
redcap_fields <- c(
    "record_id", "b_number", "tm_number", "accession_number", "block",
    "diagnosis", "organ", "tissue_comments", "run_number", "nyu_mrn",
    "qc_passed", "arrived"
)
cnv_fields <- c(
    "record_id", "b_number", "primary_tech", "second_tech", "run_number",
    "barcode_and_row_column", "accession_number", "tm_number", "arrived"
)
pact_columns <- c(
    "Tumor Specimen ID", "Normal Specimen ID", "Tumor DNA/RNA Number", "MRN", "Test Number"
)

# Homebrew and compiler setup --------------------------------------------------
install_brew <- function() {
    message("Installing Homebrew...")
    system('/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"')
}

# Adds brew directory to session PATH and to ~/.Renviron
fix_brew_path <- function() {
    brews <- c("/opt/homebrew/bin/brew", "/usr/local/bin/brew")
    if (!any(file.exists(brews))) install_brew()
    brew_dir <- dirname(if (file.exists(brews[1])) brews[1] else brews[2])
    path <- unique(c(brew_dir, strsplit(Sys.getenv("PATH"), ":", fixed = TRUE)[[1]]))
    Sys.setenv(PATH = paste(path, collapse = ":"))

    renviron <- file.path(Sys.getenv("HOME"), ".Renviron")
    entry <- paste0('PATH="', paste(path, collapse = ":"), '"')
    lines <- if (file.exists(renviron)) readLines(renviron, warn = FALSE) else character()
    if (!any(grepl(entry, lines, fixed = TRUE))) {
        path_lines <- grep("^PATH=", lines)
        if (length(path_lines) > 0) lines[path_lines] <- entry else lines <- c(lines, entry)
        writeLines(lines, renviron)
    }
}

# Installs Homebrew and missing brew formulae
setup_homebrew <- function() {
    fix_brew_path()
    if (!nzchar(Sys.which("brew"))) {
        install_brew()
        fix_brew_path()
    }
    installed <- system2("brew", c("list", "--formula"), stdout = TRUE, stderr = FALSE)
    for (pkg in setdiff(brew_pkgs, installed)) system2("brew", c("install", pkg))
}

brew_prefix <- function(pkg = "") {
    system2("brew", c("--prefix", pkg), stdout = TRUE, stderr = FALSE)
}

# Points compilers and linker flags to brew llvm and apache-arrow
set_env_vars <- function() {
    Sys.unsetenv(c(
        "CC", "CXX", "OBJC", "LDFLAGS", "CPPFLAGS", "PKG_CFLAGS",
        "PKG_LIBS", "LD_LIBRARY_PATH", "R_LD_LIBRARY_PATH"
    ))
    brew <- brew_prefix()
    llvm <- brew_prefix("llvm")
    arrow <- brew_prefix("apache-arrow")
    llvm_libs <- file.path(llvm, c("lib", "lib/c++", "lib/unwind"))
    Sys.setenv(
        CC = file.path(llvm, "bin/clang"),
        CXX = file.path(llvm, "bin/clang++"),
        OBJC = file.path(llvm, "bin/clang"),
        LDFLAGS = paste(
            c(paste0("-L", llvm_libs), paste0("-Wl,-rpath,", llvm_libs[-1]), "-lunwind"),
            collapse = " "
        ),
        CPPFLAGS = paste0("-I", llvm, "/include"),
        PKG_CFLAGS = paste0("-I", c(brew, arrow), "/include", collapse = " "),
        PKG_LIBS = paste(paste0("-L", c(brew, llvm, arrow), "/lib", collapse = " "), "-larrow"),
        LD_LIBRARY_PATH = file.path(brew, "lib"),
        R_LD_LIBRARY_PATH = paste(file.path(brew, "lib"), llvm_libs[2], sep = ":"),
        DYLD_LIBRARY_PATH = file.path(arrow, "lib")
    )
    if (nzchar(Sys.which("gfortran"))) Sys.setenv(FC = Sys.which("gfortran"))
}

# R package setup --------------------------------------------------------------
install_pak <- function() {
    tryCatch(
        install.packages("pak", repos = sprintf(
            "https://r-lib.github.io/p/pak/stable/%s/%s/%s",
            .Platform$pkgType, R.Version()$os, R.Version()$arch
        )),
        error = function(e) {
            install.packages(
                "pak", ask = FALSE, dependencies = TRUE,
                repos = "https://packagemanager.rstudio.com/all/latest"
            )
        }
    )
}

# Installs missing packages with pak, then attaches all
ensure_packages <- function(pkgs) {
    installed <- rownames(installed.packages())
    missing_pkgs <- setdiff(pkgs, installed)
    if (length(missing_pkgs) > 0) {
        message("Installing missing packages: ", toString(missing_pkgs))
        if (!"pak" %in% installed) install_pak()
        for (pkg in missing_pkgs) {
            tryCatch(
                pak::pkg_install(pkg, ask = FALSE),
                error = function(e) install.packages(pkg, ask = FALSE, dependencies = TRUE)
            )
        }
    }
    attached <- vapply(pkgs, function(pkg) {
        suppressWarnings(suppressPackageStartupMessages(library(
            pkg, mask.ok = TRUE, character.only = TRUE, logical.return = TRUE
        )))
    }, logical(1))
    if (!all(attached)) stop("Failed to load packages: ", toString(pkgs[!attached]))
}

has_pkg <- function(pkg, version = NULL) {
    requireNamespace(pkg, quietly = TRUE) &&
        (is.null(version) || utils::packageVersion(pkg) == version)
}

# Installs pinned minfi, EPICv2 manifest, classifier, and conumee when missing
ensure_cnv_packages <- function() {
    manifest <- "IlluminaHumanMethylationEPICv2manifest"
    if (!has_pkg("minfi", "1.43.1") || !has_pkg(manifest, "0.1.0")) {
        Sys.setenv(R_COMPILE_AND_INSTALL_PACKAGES = "always")
        for (repo in file.path("mwsill", c("minfi", manifest))) {
            devtools::install_github(
                repo, upgrade = "always", force = TRUE, dependencies = TRUE,
                type = "source", auth_token = NULL
            )
        }
    }
    if (!has_pkg("mnp.v12epicv2") || !has_pkg("conumee2.0")) source(classifier_install)

    needed <- c("conumee2.0", "minfi", manifest, "mnp.v12epicv2")
    available <- vapply(needed, has_pkg, logical(1))
    if (!all(available)) stop("Packages not installed: ", toString(needed[!available]))
}

# REDCap API -------------------------------------------------------------------
# Posts form to REDCap and returns response text
redcap_post <- function(body) {
    response <- httr::POST(REDCAP_URL, body = c(list(token = REDCAP_API_TOKEN), body), encode = "form")
    httr::stop_for_status(response)
    httr::content(response, as = "text", encoding = "UTF-8")
}

# Builds indexed API parameters such as fields[0], fields[1]
redcap_array <- function(name, values) {
    stats::setNames(as.list(values), sprintf("%s[%d]", name, seq_along(values) - 1L))
}

# Exports records as data frame of character columns
redcap_export <- function(fields, records = NULL) {
    text <- redcap_post(c(
        list(
            content = "record", action = "export", format = "csv", type = "flat",
            csvDelimiter = "", rawOrLabel = "raw", rawOrLabelHeaders = "raw",
            exportCheckboxLabel = "false", exportSurveyFields = "false",
            exportDataAccessGroups = "false", returnFormat = "json"
        ),
        redcap_array("records", records),
        redcap_array("fields", fields)
    ))
    if (!nzchar(trimws(text))) return(data.frame())
    if (startsWith(trimws(text), "{")) stop("REDCap API returned an error: ", text)
    utils::read.csv(
        text = text, check.names = FALSE, colClasses = "character", na.strings = c("", "NA")
    )
}

# Imports one record given as named list of field values
redcap_import <- function(record) {
    invisible(redcap_post(list(
        content = "record", format = "json", type = "flat",
        data = jsonlite::toJSON(list(record), auto_unbox = TRUE, na = "null"),
        returnContent = "nothing", returnFormat = "json"
    )))
}


# PACT run inputs --------------------------------------------------------------
# Returns path to run workbook that holds Beaker or Philips export tab
get_excel_path <- function() {
    if (is_file_path) return(PACT_INPUT)
    run_year <- stringr::str_split_fixed(PACT_INPUT, "-", 3)[, if (is_sophia) 1 else 2]
    workbook <- file.path(
        pact_data_dir, "Workbook", paste0("20", run_year), PACT_INPUT, paste0(PACT_INPUT, ".xlsm")
    )
    if (is_sophia) message("Run type is Sophia, looking for workbook in:")
    message(workbook)
    workbook
}

# Picks first other workbook in run folder when expected .xlsm is absent
find_alt_workbook <- function(pact_sheet) {
    message(crayon::bgRed("PACT run worksheet not found:"), "\n", pact_sheet)
    message("Checking other files in PACT folder: ", basename(dirname(pact_sheet)))
    files <- list.files(dirname(pact_sheet), full.names = TRUE)
    workbooks <- files[grepl("\\.xlsm$|book", basename(files))]
    if (!any(grepl("\\.xlsm$", workbooks))) {
        message("\nNo .xlsm worksheet found. Checking .xlsx files and others...\n")
    }
    workbooks <- workbooks[!grepl("$", basename(workbooks), fixed = TRUE)]
    if (length(workbooks) == 0) stop("No alternative file found.")
    message(crayon::bgGreen("Using this workbook instead:"), basename(workbooks[1]))
    workbooks[1]
}

# Returns run sheet path: .xlsm workbook, or demux samplesheet for Results folder runs
get_pact_sheet <- function() {
    id_parts <- stringr::str_split_fixed(PACT_INPUT, "-", 3)
    run_dir <- file.path(pact_data_dir, "Workbook", paste0("20", id_parts[2]), PACT_INPUT)
    if (!dir.exists(run_dir)) {
        run_dir <- file.path(
            pact_data_dir, "Results", "Bioinformatics", paste0("20", id_parts[1]), PACT_INPUT
        )
        if (!dir.exists(run_dir)) stop("PACT run folder not found: ", run_dir)
    }
    if (grepl("Results", run_dir)) {
        pact_sheet <- file.path(run_dir, "demux-samplesheet.csv")
    } else {
        pact_sheet <- file.path(run_dir, paste0(PACT_INPUT, ".xlsm"))
        if (!file.exists(pact_sheet)) pact_sheet <- find_alt_workbook(pact_sheet)
    }
    message("Using the following PACT sheet file:\n", pact_sheet)
    pact_sheet
}

# Reads sample identifiers from PhilipsExport tab
parse_worksheet <- function(pact_sheet) {
    message("Reading the file:\n", pact_sheet)
    sheets <- readxl::excel_sheets(pact_sheet)
    message("Excel sheet names:\n", paste(sheets, collapse = "\n"))
    stopifnot(length(sheets) > 2)
    philips_tab <- grep("PhilipsExport", sheets, ignore.case = TRUE, value = TRUE)[1]
    vals <- suppressMessages(as.data.frame(
        readxl::read_excel(pact_sheet, sheet = philips_tab, skip = 3, col_types = "text")
    ))[, pact_columns]
    vals[!is.na(vals[, 1]), ]
}

# Reads sample identifiers from demux samplesheet
parse_demux_csv <- function(pact_sheet) {
    demux <- read.csv(pact_sheet, skip = 19)
    demux <- demux[demux$Tumor_Content != 0, ]
    vals <- data.frame(
        demux$TUMOR_CASE_ID_BLOCK,
        sub("-[^-]+$", "", demux$TUMOR_CASE_ID_BLOCK),
        demux$Tumor_DNA,
        stringr::str_split_fixed(demux$Sample_ID, "_", 3)[, 1],
        demux$TM_Number
    )
    stats::setNames(vals, pact_columns)
}

# Returns sample identifiers for csv file, workbook path, or PACT run ID input
get_case_values <- function() {
    if (read_flag && grepl("-SampleSheet", PACT_INPUT)) {
        message("Parsing Data from Demux SampleSheet.csv file...")
        vals <- utils::read.csv(PACT_INPUT, skip = 19)[, c(6, 7, 9)]
        return(as.data.frame(vals[!grepl("H20|SERACARE|HAPMAP", vals[, 2]), ]))
    }
    if (read_flag) {
        message("Parsing Data from .csv file that is not a Demux SampleSheet...")
        vals <- unlist(read.csv(PACT_INPUT))
        return(as.data.frame(unique(vals[vals != ""])))
    }
    if (is_file_path) {
        message("Parsing Data from .xlsm/.xlsx file path...")
        return(parse_worksheet(PACT_INPUT))
    }
    message("Parsing Data from PACT RUN ID: ", PACT_INPUT, " finding run worksheet...")
    pact_sheet <- get_pact_sheet()
    if (grepl("\\.csv$", pact_sheet)) parse_demux_csv(pact_sheet) else parse_worksheet(pact_sheet)
}

# Returns run ID used to name output worksheet and REDCap record
get_pact_id <- function() {
    if (is_sophia) return(read.csv(get_pact_sheet(), skip = 19)$Sample_Project[1])
    if (grepl(".xls", PACT_INPUT)) {
        return(substr(basename(PACT_INPUT), 1, nchar(basename(PACT_INPUT)) - 5))
    }
    if (read_flag) substr(PACT_INPUT, 1, nchar(PACT_INPUT) - 4) else PACT_INPUT
}

# Matching PACT samples to REDCap ----------------------------------------------
# Pairs every sample identifier with Test Number of its row
build_query_table <- function(vals) {
    stopifnot("Test Number" %in% names(vals))
    queries <- data.frame(
        Test_Number = rep(as.character(vals[["Test Number"]]), times = ncol(vals)),
        query_value = trimws(unlist(lapply(vals, as.character), use.names = FALSE))
    )
    queries <- queries[!is.na(queries$query_value) & !queries$query_value %in% c("", "0"), ]

    # Adds shortened TS, TB, and TC identifiers
    short <- queries[grepl("^(TS|TB|TC)-[^-]+-", queries$query_value), ]
    short$query_value <- sub("^((TS|TB|TC)-[^-]+)-.*$", "\\1", short$query_value)
    unique(rbind(queries, short))
}

# Returns REDCap rows where any field contains sample identifier
query_cases <- function(vals, db) {
    queries <- build_query_table(vals)
    db_text <- do.call(cbind, lapply(db, as.character))
    db_text[is.na(db_text)] <- ""

    matches <- dplyr::bind_rows(lapply(seq_len(nrow(queries)), function(i) {
        value <- queries$query_value[i]
        hits <- which(matrix(grepl(value, db_text, fixed = TRUE), nrow(db_text)), arr.ind = TRUE)
        if (nrow(hits) == 0) return(NULL)
        data.frame(
            db_row = hits[, "row"],
            Test_Number = queries$Test_Number[i],
            query_value = value,
            matched_column = names(db)[hits[, "col"]],
            matched_value = db_text[hits]
        )
    }))
    if (nrow(matches) == 0) {
        return(cbind(db[0, , drop = FALSE], Test_Number = character(0)))
    }

    match_lines <- sprintf(
        "Match found for '%s' (%s) for %s in: \"%s\" column",
        matches$query_value, matches$matched_value, matches$Test_Number, matches$matched_column
    )
    message(paste(match_lines, collapse = "\n"))
    cat(match_lines, file = match_tsv, append = TRUE, sep = "\n")

    output <- db[matches$db_row, , drop = FALSE]
    output$Test_Number <- matches$Test_Number
    output <- unique(output)
    rownames(output) <- NULL
    output
}

# Fills missing Test Numbers from PACT sheet by accession number
fill_test_numbers <- function(output, vals) {
    if (all(output$Test_Number %in% vals$`Test Number`)) {
        message("All NGS Test Numbers Found in Methylation Database")
        return(output)
    }
    message("Not all NGS do not have methylation")
    unfilled <- is.na(output$Test_Number)
    accession_row <- match(
        output$accession_number[unfilled], vals$`Tumor Specimen ID`, incomparables = NA
    )
    output$Test_Number[unfilled] <- vals$`Test Number`[accession_row]
    if (anyNA(output$Test_Number)) {
        warning(
            "Some samples still missing NGS Numbers:\n",
            paste(output$record_id[is.na(output$Test_Number)], collapse = "\n")
        )
    }
    return(output)
}

# Adds report status and smb report links for MGDM runs
add_report_links <- function(output) {
    run <- output$run_number
    has_report <- grepl("MGDM", run)
    year <- gsub("MC", "", sub("-.*$", "", run))
    year <- paste0("20", ifelse(nchar(year) > 2, substring(year, 3), year))
    link <- file.path(report_share, year, run, paste0(output$record_id, ".html"))
    output$report_complete <- ifelse(has_report, "YES", "NOT_YET_RUN")
    output$`Report Link` <- ifelse(has_report, link, "")
    output$`Report Path` <- output$`Report Link`
    output
}

# Report paths -----------------------------------------------------------------
# Returns mounted volume paths of MGDM reports that do not exist
missing_reports <- function(report_paths) {
    volume_paths <- sub(smb_share, "/Volumes", report_paths, fixed = TRUE)
    volume_paths <- volume_paths[grepl("MGDM", volume_paths)]
    volume_paths[!file.exists(volume_paths)]
}

# Points paths that end in missing html name to matching file found in run_dir
swap_report_file <- function(report_paths, html, run_dir, new_template = FALSE) {
    found <- dir(run_dir, pattern = sub("\\.html$", "", html), full.names = TRUE)
    if (length(found) > 1) {
        found <- found[paste0(sub("_.*", "", basename(found)), ".html") == html]
        to_swap <- which(sub("_.*", "", basename(report_paths)) == html)
    } else {
        to_swap <- grep(html, report_paths)
    }
    if (length(found) == 0) return(report_paths)

    new_paths <- stringr::str_replace(report_paths[to_swap], html, basename(found[1]))
    if (new_template) {
        new_paths <- file.path(paste0(dirname(new_paths), "-new-template"), basename(found[1]))
    }
    message("Updating file path:\n", new_paths)
    report_paths[to_swap] <- new_paths
    report_paths
}

# Corrects report path year, then repairs paths to reports that do not exist
check_report_paths <- function(report_paths) {
    # Sets year folder from run number prefix
    parts <- stringr::str_split_fixed(report_paths, "/", 11)
    mgdm <- grepl("MGDM", parts[, 10])
    parts[mgdm, 9] <- paste0("20", sub("-.*$", "", parts[mgdm, 10]))
    report_paths[mgdm] <- apply(parts[mgdm, , drop = FALSE], 1, paste, collapse = "/")

    missing <- missing_reports(report_paths)
    if (length(missing) == 0) return(report_paths)
    message("Fixing broken file paths...")

    # Substitutes similarly named run folders for run folders that do not exist
    run_dirs <- unique(dirname(missing))
    for (i in which(!dir.exists(run_dirs))) {
        similar <- dir(dirname(run_dirs[i]), pattern = basename(run_dirs[i]), full.names = TRUE)
        if (length(similar) > 0) run_dirs[i] <- similar[1]
    }
    for (html in basename(missing)) {
        message("Fixing path for missing report: ", html)
        for (run_dir in run_dirs) report_paths <- swap_report_file(report_paths, html, run_dir)
    }

    # Checks "-new-template" run folders for reports still missing
    for (report in missing_reports(report_paths)) {
        message("Fixing path for missing report: ", basename(report))
        report_paths <- swap_report_file(
            report_paths, basename(report), paste0(dirname(report), "-new-template"),
            new_template = TRUE
        )
    }

    missing <- missing_reports(report_paths)
    if (length(missing) > 0) {
        message(
            crayon::bgRed("Some paths to html reports need editing in MethylMatch.xlsx sheet!"), "\n",
            crayon::bgRed("Fix the following paths in worksheet 'Report Path' column that do not exist:"), "\n"
        )
        message(paste(missing, collapse = "\n"), "\n")
    }
    report_paths
}

# Writes match workbook with report hyperlinks to Desktop
write_match_xlsx <- function(output, pact_id) {
    output$`Report Path` <- check_report_paths(output$`Report Path`)
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, pact_id)
    openxlsx::writeData(wb, sheet = pact_id, x = output)
    link_col <- which(names(output) == "Report Link")
    for (i in which(output$`Report Link` != "")) {
        link <- structure(
            output$`Report Link`[i],
            names = paste0(output$record_id[i], ".html"),
            class = "hyperlink"
        )
        openxlsx::writeData(wb, sheet = pact_id, x = link, startCol = link_col, startRow = i + 1)
    }
    xlsx_path <- file.path(desktop, paste0(pact_id, "_MethylMatch.xlsx"))
    openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
    return(xlsx_path)
}

# IDAT files -------------------------------------------------------------------
idat_names <- function(sentrix_ids) {
    c(paste0(sentrix_ids, "_Grn.idat"), paste0(sentrix_ids, "_Red.idat"))
}

stop_if_unmounted <- function(path) {
    if (!dir.exists(path)) {
        stop(crayon::bgRed("Share drive not found, ensure path is accessible:"), "\n", path)
    }
}

# Writes minfi samplesheet for REDCap records that have Sentrix barcode
write_samplesheet <- function(rds) {
    records <- redcap_export(cnv_fields, records = rds)
    has_barcode <- !is.na(records$barcode_and_row_column) & nzchar(records$barcode_and_row_column)
    records <- records[has_barcode, , drop = FALSE]
    if (nrow(records) == 0) stop("No REDCap records with barcode_and_row_column were returned")

    sentrix <- stringr::str_split_fixed(records$barcode_and_row_column, "_", 2)
    samplesheet <- data.frame(
        Sample_Name = records$record_id,
        DNA_Number = records$b_number,
        Sentrix_ID = sentrix[, 1],
        Sentrix_Position = sentrix[, 2],
        SentrixID_Pos = records$barcode_and_row_column,
        Basename = file.path(getwd(), records$barcode_and_row_column),
        RunID = records$run_number,
        MP_num = records$tm_number,
        tech = records$primary_tech,
        tech2 = records$second_tech,
        Date = records$arrived
    )
    message("Writing REDCap data to samplesheet.csv")
    message(paste(capture.output(samplesheet), collapse = "\n"))
    utils::write.csv(samplesheet, "samplesheet.csv", quote = FALSE, row.names = FALSE)
}

# Copies readable idat files that are not yet in working directory
copy_idats <- function(idat_files) {
    readable <- fs::file_access(idat_files, mode = "read")
    if (any(!readable)) {
        write(
            paste("Cannot read idat file:", idat_files[!readable]),
            "read_error_idat.txt", append = TRUE
        )
    }
    already_copied <- basename(idat_files) %in% list.files(pattern = "\\.idat$")
    to_copy <- idat_files[readable & !already_copied]
    if (length(to_copy) == 0) return(message(".idat files already copied to run directory"))

    cli::cli_progress_bar("Copying files", total = length(to_copy))
    for (idat in to_copy) {
        tryCatch(
            fs::file_copy(idat, file.path(getwd(), basename(idat)), overwrite = TRUE),
            error = function(e) {
                cli::cli_alert_danger(paste("Failed to copy:", idat))
                write(paste("Failed to copy:", idat), "missing_idat_files.txt", append = TRUE)
            }
        )
        cli::cli_progress_update()
    }
    cli::cli_progress_done()
}

# Finds idat pairs for samplesheet.csv on idat drives and copies them to working directory
get_idats <- function() {
    stop_if_unmounted(research_idat_dir)
    stop_if_unmounted(clinical_idat_dir)
    samplesheet <- utils::read.csv("samplesheet.csv", strip.white = TRUE)
    expected <- idat_names(unique(samplesheet$SentrixID_Pos))

    idats <- unlist(lapply(
        c(research_idat_dir, clinical_idat_dir), file.path,
        samplesheet$Sentrix_ID, idat_names(samplesheet$SentrixID_Pos)
    ))
    idats <- idats[file.exists(idats)]

    absent <- setdiff(expected, basename(idats))
    if (length(absent) > 0) {
        external_dir <- file.path(research_idat_dir, "External")
        message(
            crayon::bgRed("Still missing some idats! Checking External Folder:"), "\n", external_dir
        )
        wanted <- idat_names(unique(sub("_(Grn|Red)\\.idat$", "", absent)))
        external <- file.path(external_dir, wanted)
        external <- external[file.exists(external)]
        if (!all(wanted %in% basename(external)) && dir.exists(external_dir)) {
            nested <- list.files(external_dir, pattern = "\\.idat$", full.names = TRUE, recursive = TRUE)
            nested <- nested[basename(nested) %in% setdiff(wanted, basename(external))]
            external <- unique(c(external, nested))
        }
        idats <- c(idats, external)
    }

    idats <- idats[!duplicated(basename(idats))]
    if (length(idats) == 0) {
        stop(
            "No .idat files found for these samples. Checked:\n",
            research_idat_dir, "\n", clinical_idat_dir
        )
    }

    absent <- setdiff(expected, basename(idats))
    if (length(absent) > 0) {
        utils::write.csv(
            data.frame(Missing_Samples = unique(sub("_(Grn|Red)\\.idat$", "", absent))),
            "missing_idats_log.csv", row.names = FALSE, quote = FALSE
        )
        warning(
            "Some IDAT files were not found. See missing_idats_log.csv:\n",
            paste(absent, collapse = "\n")
        )
    }
    message("Files found:\n", paste(idats, collapse = "\n"))
    copy_idats(idats)
}

# CNV PNG creation -------------------------------------------------------------
# Drops RD-numbers that already have CNV PNG in output folder
skip_existing_pngs <- function(rds) {
    out_pngs <- file.path(cnv_out_dir, paste0(rds, "_cnv.png"))
    finished <- file.exists(out_pngs)
    if (any(finished)) {
        message(crayon::bgGreen("Cases below already exist and will be skipped:"))
        message(paste(out_pngs[finished], collapse = "\n"))
        rds <- rds[!finished]
        if (length(rds) > 0) {
            message(crayon::bgGreen("The following CNV PNG will be generated:"))
            message(paste(rds, collapse = "\n"))
        }
    }
    rds
}

# Plots CNV PNG for one sample to Desktop
make_cnv_png <- function(sample_name, sentrix_id) {
    png_path <- file.path(desktop, paste0(sample_name, "_cnv.png"))
    rgset <- minfi::read.metharray(file.path(getwd(), sentrix_id), verbose = TRUE, force = TRUE)
    mset <- minfi::preprocessIllumina(rgset, bg.correct = TRUE, normalize = "controls")
    cnv <- mnp.v12epicv2::MNPcnv(mset, sex = NULL, main = sample_name)
    message("Saving file to:\n", png_path)
    png(filename = png_path, width = 1820, height = 1040, res = 150)
    conumee2.0::CNV.genomeplot(
        cnv,
        chr = paste0("chr", 1:22),
        main = sample_name,
        bins_cex = "sample_level",
        cols = c("darkred", "salmon", "lightgrey", "lightgreen", "darkgreen")
    )
    invisible(dev.off())
}

# Gets idats for RD-numbers and plots CNV PNG for each sample to Desktop
make_cnv_pngs <- function(rds) {
    message("\nRD-numbers with idats:\n", paste(rds, collapse = "\n"))
    message("--------", crayon::bgMagenta("Starting CNV PNG Creation"), "--------")
    ensure_cnv_packages()
    write_samplesheet(rds)
    get_idats()
    targets <- read.csv("samplesheet.csv")
    targets <- targets[targets$Sample_Name %in% rds & grepl("_R0", targets$SentrixID_Pos), ]

    tryCatch(
        {
            if (nrow(targets) > 0) {
                mapply(
                    make_cnv_png,
                    as.character(targets$Sample_Name), as.character(targets$SentrixID_Pos)
                )
            } else {
                message("The RD-number(s) do not have idat files in REDCap:\n")
                print(targets)
            }
            pngs <- file.path(desktop, paste0(targets$Sample_Name, "_cnv.png"))
            created <- file.exists(pngs)
            if (any(!created)) {
                message("The following failed to be created:")
                print(pngs[!created])
                message("Try running again or check GitHub troubleshooting")
            }
            if (any(created)) {
                message("The following were created successfully:")
                print(pngs[created])
            }
            while (!is.null(dev.list())) dev.off()
        },
        error = function(e) {
            message("The following error occured:\n", conditionMessage(e))
            message("\nTry checking the troubleshooting section on GitHub")
            stop(crayon::bgRed("CNV PNG generation failed"))
        }
    )
}

# Copies CNV PNGs created today from Desktop to output folder
copy_cnv_pngs <- function() {
    pngs <- dir(desktop, "_cnv.png", full.names = TRUE)
    pngs <- pngs[as.Date(file.info(pngs)$ctime) == Sys.Date()]
    if (length(pngs) == 0) return(message("No CNV files found on Desktop to copy"))

    copies <- file.path(cnv_out_dir, basename(pngs))
    message("\nCopying png files to Molecular folder:\n", cnv_out_dir, "\n")
    message(paste(pngs, collapse = "\n"))
    fs::file_copy(path = pngs, new_path = copies)
    if (any(!file.exists(copies))) {
        message("The following failed to copy from the desktop:\n")
        print(basename(copies[!file.exists(copies)]))
        message(crayon::bgRed("Manually copy any methylation CNV PNGs that failed to copy"))
    }
}

# Main -------------------------------------------------------------------------
message("\n================ Parameters input ================\n")
message("token: [provided]\nPACT_INPUT: ", PACT_INPUT, "\n")

# Compilers and packages
if (!"devtools" %in% rownames(installed.packages())) {
    install.packages("devtools", ask = FALSE, dependencies = TRUE)
}
suppressPackageStartupMessages(library("devtools"))
setup_homebrew()
set_env_vars()
ensure_packages(main_pkgs)
if (!"mnp.v12epicv2" %in% rownames(installed.packages())) {
    stop("You need to install classifier package and pre-reqs first use all_installer.R")
}
suppressWarnings(suppressPackageStartupMessages(library("mnp.v12epicv2")))

if (!dir.exists(lab_drive)) {
    message("Network share is not mounted:\n", crayon::bgRed(file.path(smb_share, "molecular")))
    stop("Molecular shared drive is not mounted")
}

# PACT sample identifiers, with MRN from Beaker or Philips export tab
excel_path <- get_excel_path()
sheets <- readxl::excel_sheets(excel_path)
export_tab <- sheets[grepl(if (is_sophia) "Beaker" else "Philips", sheets)]
tab_mrns <- readxl::read_xlsx(excel_path, sheet = export_tab)[, c("Specimen ID", "MRN")]

vals <- get_case_values()
if (nrow(vals) == 0) stop("No samples found in PACT input: ", PACT_INPUT)
mrn_row <- match(vals$`Test Number`, tab_mrns$`Specimen ID`, incomparables = NA)
if (anyNA(mrn_row)) {
    stop(
        "Test Number not found in '", export_tab, "' tab: ",
        toString(vals$`Test Number`[is.na(mrn_row)])
    )
}
vals$beaker_mrn <- as.character(tab_mrns$MRN[mrn_row])
message("Values from PACT demux csv used for matching to METH REDCap DB:")
message(paste(capture.output(vals), collapse = "\n"))

# REDCap matches
message("Pulling REDCap data...")
db <- redcap_export(redcap_fields)
if (nrow(db) == 0) stop("REDCap API returned no records")
record_rows <- grepl(pattern = "RD-", db$record_id)
db <- db[record_rows,]
rownames(db) <- NULL

output <- query_NULLoutput <- query_cases(vals, db)
if (nrow(output) == 0) {
    warning("No Methylation Cases on this PACT run, generating blank file")
    output[1, ] <- "NONE"
}
pact_id <- get_pact_id()
tm_numbers_df <- fill_test_numbers(output, vals)
output <- add_report_links(tm_numbers_df)

if (read_flag) {
    write.csv(output, file = "meth_sample_data.csv",
              quote = FALSE, row.names = FALSE)
}

if (anyNA(output$Test_Number)) {
    message("Potential NGS Number in 'tm_number' field of REDCap!")
    has_ngs <- grepl("NGS", output$tm_number)
    output$Test_Number[has_ngs] <- output$tm_number[has_ngs]
}

meth_xlsx <- write_match_xlsx(output, pact_id)

# CNV PNGs for cases with completed methylation report
rds <- grep("^RD-", output$record_id[output$report_complete == "YES"], value = TRUE)
if (length(rds) == 0) {
    message(crayon::bgGreen("The PACT run has no cases with methylation."))
    quit(status = 0)
}

rds <- skip_existing_pngs(rds)

if (length(rds) == 0) {
    message(crayon::bgGreen("No CNV png images to generate. Check the output directory."))
    quit(status = 0)
}

make_cnv_pngs(rds)
copy_cnv_pngs()
