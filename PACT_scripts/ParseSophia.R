#!/usr/bin/env Rscript
## Script name: pactParse.R
## Purpose: source of global scripts and generate PACT -SampleSheet.csv file
## Author: Jonathan Serrano
## Date Created: February 10, 2026
## Version: 1.1.0
## Copyright (c) NYULH Jonathan Serrano, 2026

args <- commandArgs(TRUE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

key_run_ids <- "NB551709|NB501073|VH01471"
DEBUG <- FALSE

# Function to assign parameter values based on args ---------------------------
assign_param <- function(args, idx) {
    return(if (!is.na(args[idx]) && nzchar(args[idx])) args[idx] else NULL)
}

# Assign values to variables & print ------------------------------------------
INPUT <- assign_param(args, 1)
RUNID <- assign_param(args, 2)

API_URL <- "https://redcap.nyumc.org/apps/redcap/api/"

TOKEN_PATH <- "/Volumes/CBioinformatics/scripts/PACT_DB_API.txt"
stopifnot(file.exists(TOKEN_PATH))
API_TOKEN <- trimws(readLines(TOKEN_PATH, n = 1, warn = FALSE))

# Functions to print parameters  ) --------------------------------------------
brick_message <- function(phrase) {
    dsh <- paste(rep("-", nchar(phrase)), collapse = "")
    message(paste(dsh, phrase, dsh, sep = "\n"))
}

print_parameters <- function(INPUT, RUNID, API_TOKEN) {
    brick_message("Parameters INPUT")
    message(sprintf(
        "\nINPUT:  %s\nRUNID:  %s\nAPI_TOKEN:  %s\n",
        c(INPUT, "NULL")[1], c(RUNID, "NULL")[1], c(API_TOKEN, "NULL")[1]))
}

print_parameters(INPUT, RUNID, API_TOKEN)

# Helper functions to message as data frame and red background ----------------
msg_df <- function(...) {
    message(paste0(capture.output(as.data.frame(...)), collapse = "\n"), "\n")
}

mkRed <- function(...) {return(crayon::white$bgRed$bold(...))}
mkBlu <- function(...) {return(crayon::white$bgBlue$bold(...))}

# FUN: Checks if z-drive is accessible to the RScript -------------------------
checkMounts <- function() {
    if (!dir.exists("/Volumes/molecular/MOLECULAR LAB ONLY")) {
        zDrive = "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular"
        stop(paste("Path is not accessible:",
                   mkRed("/Volumes/molecular/MOLECULAR LAB ONLY"),
                   "Ensure the network Z-drive is mounted:",
                   mkRed(zDrive), sep = "\n"), call. = FALSE)
    }
}

options(repos = c(getOption("repos"), CRAN = "http://cran.us.r-project.org"))

# Function to silently load library without conflict warnings -----------------
libLoad <- function(pkg) {
    if (!pkg %in% loadedNamespaces())
    suppressWarnings(suppressPackageStartupMessages(
        library(pkg, character.only = TRUE, warn.conflicts = FALSE, quietly = TRUE)
    ))
}

# Installs package(s) if not already installed then loads the library ---------
SilentInstall <- function(pkgs) {
    for (pkg in pkgs) {
        if (!pkg %in% rownames(installed.packages()))
            install.packages(pkg, dependencies = TRUE, ask = FALSE)
        libLoad(pkg)
    }
}


# Function to load or install missing required packages -----------------------
loadPacks <- function() {
    SilentInstall("devtools")
    pkgs <- c("RCurl", "jsonlite", "digest", "readxl", "stringr",
              "httr", "fs", "crayon", "dplyr", "tibble", "data.table")
    SilentInstall(pkgs)
}


# Retuns message of any columns and characters removed ------------------------
MsgChangesMade <- function(mainSheet, badChars = c("[\r\n]", ",", " ")) {
    for (patt in badChars) {
        for (i in seq_along(mainSheet)) {
            col <- ifelse(is.na(mainSheet[[i]]), "", mainSheet[[i]])
            idxMatch <- which(stringr::str_detect(col, pattern = patt))
            if (length(idxMatch) > 0) {
                replaced_df <- as.data.frame(mainSheet[idxMatch, i])
                replace_col <- names(mainSheet)[i]
                colnames(replaced_df) <- replace_col
                if (patt == "[\r\n]") {
                    patt <- "[\\r\\n]"
                }
                message(mkBlu(sprintf(
                    'The following "%s" will be removed from the column: %s',
                                patt, replace_col)))
                msg_df(replaced_df)
            }
        }
    }
}


CountControls <- function(mainSheet, controlNames = "NTC_H20|SC_SERACARE|NC_HAPMAP") {
    controlSamples <- grepl(pattern = controlNames, mainSheet$Sample_ID)
    if (table(controlSamples)[['TRUE']] != 3) {
        message(mkRed("There are not 3 Control samples"))
        message("Either NTC_H20, SC_SERACARE, or NC_HAPMAP is missing")
        message("OR\nThere are extra controls added in this run:")
        msg_df(mainSheet$Sample_Name[controlSamples])
    } else{
        controlIndexes <- which(controlSamples == TRUE)
        mainSheet[controlIndexes, 'Paired_Normal'] <- ""
    }
    return(mainSheet)
}


# Removes and fixes newlines, commas, and blanks from samplesheet -------------
sanitizeSheet <- function(mainSheet) {
    MsgChangesMade(data.frame(Tumor_Type = mainSheet$Tumor_Type), " ")
    mainSheet$Tumor_Type <- gsub(" ", "-", mainSheet$Tumor_Type)
    MsgChangesMade(mainSheet)
    mainSheet <- data.frame(lapply(mainSheet, function(col) gsub("[\r\n, ]", "", col)))

    pn_blank <- mainSheet$Paired_Normal == 0 | is.na(mainSheet$Paired_Normal)
    mainSheet$Paired_Normal[pn_blank] <- ""
    mainSheet$Tumor_Type[mainSheet$Tumor_Type == 0] <- "NA"
    mainSheet[, 1:ncol(mainSheet)] <- sapply(mainSheet[, 1:ncol(mainSheet)],
                                             function(x) {gsub("\\\\", "-", x)})
    mainSheet <- CountControls(mainSheet)

    mainSheet$Date_of_Procedure <- gsub("/", "-", mainSheet$Date_of_Procedure)
    mainSheet$Date_of_Procedure <- gsub(" ", "", mainSheet$Date_of_Procedure)


    mainSheet$Tumor_Specimen_Site <- gsub(" ", "-", mainSheet$Tumor_Specimen_Site)
    mainSheet$Tumor_Specimen_Site <- gsub(",", "-", mainSheet$Tumor_Specimen_Site)

    if (any(duplicated(mainSheet$Sample_ID))) {
        message(mkRed("There are duplicated Sample_ID in the SampleSheet"))
        msg_df(mainSheet$Sample_ID[duplicated(mainSheet$Sample_ID)])
    }

    mainSheet[] <- lapply(mainSheet, function(col) {
        col <- gsub(" ", "-", as.character(col))
        return(col)
    })

    if (any( stringr::str_detect(pattern = "...1", colnames(mainSheet)))) {
        toDrop <- min(which(stringr::str_detect(pattern = "...1", colnames(mainSheet)))) - 1
        mainSheet <- mainSheet[, c(1:toDrop)]
    }

    raw_dates <- mainSheet$Date_of_Procedure
    raw_dates <- head(raw_dates, -3)

    are_missing <- which(nchar(raw_dates) != 15)

    if (length(are_missing) > 0) {
        message(crayon::bgRed("The following cases have wrong date format:"), "\n")
        msg_df(mainSheet[are_missing,])
    }

    mainSheet[] <- lapply(mainSheet, function(x) {
        if (is.character(x)) {
            x <- gsub("--", "-", x, fixed = TRUE)
        }
        return(x)
    })

    return(mainSheet)
}


GetExcelData <- function(worksheetPath, sheetNum, shRange, toSkip = 0, include_head = FALSE) {
    sheetData <- suppressMessages(as.data.frame(readxl::read_excel(
        worksheetPath, sheet = sheetNum, na = "", range = shRange,
        col_types = "text", col_names = include_head, skip = toSkip)))
    sheetData[is.na(sheetData)] <- ""
    return(sheetData)
}


MakeFileHeader <- function(PACT_ID) {
    td <- as.character(paste0(format(Sys.time(), "%m/%d/%y")))
    header1 <- c(
        "[Header]",
        "IEMFileVersion",
        "Investigator",
        "Project Name",
        "Experiment Name",
        "Date",
        "Workflow",
        "Application",
        "Assay",
        "Description",
        "Chemistry",
        "[Reads]",
        "151",
        "151",
        "[Settings]",
        "AdapterSequenceRead1",
        "AdapterSequenceRead2"
    )

    header2 <- c(
        "",
        "4",
        "Name",
        PACT_ID,
        PACT_ID,
        td,
        "GenerateFASTQ",
        "FASTQ Only",
        "KAPA-IDT",
        "Description",
        "Capture Enrichment",
        "",
        "",
        "",
        "",
        "AGATCGGAAGAGCACACGTC",
        "AGATCGGAAGAGCGTCGTGT"
    )
    theHeader <- cbind(header1, header2)
    sheetHead <- rbind(theHeader, c("", ""), c("[Data]", ""))
    return(sheetHead)

}


check_extra_cases <- function(sheet1_data, mainSheet) {
    beaker_bNumT <- mainSheet$Tumor_DNA[mainSheet$Tumor_DNA != "0" & !is.na(mainSheet$Tumor_DNA)]
    beaker_bNumN <- mainSheet$Normal_DNA[mainSheet$Normal_DNA != "0" & !is.na(mainSheet$Normal_DNA)]
    beaker_b_numbs <- c(beaker_bNumT, beaker_bNumN)
    first_page_b <- sheet1_data$`DNA #`[
        !grepl(sheet1_data$`DNA #`, pattern = "H20|SERA|HAP") &
            sheet1_data$`DNA #` != "0" &
            sheet1_data$`DNA #` != "" & !is.na(sheet1_data$`DNA #`)
    ]
    missed <- first_page_b %in% beaker_b_numbs
    sheet_1_cols <- c(4, 5, 6, 14, 15)
    sheet1_names <- c("TM", "Surgical", "DNA_Number", "Type", "Tissue")
    if (any(missed == FALSE)) {
        message(mkRed("!!!Potential Fillers:"))
        missing_b <- first_page_b[missed == FALSE]
        missed_idx <- which(sheet1_data$`DNA #` %in% missing_b)
        missed_rows_b <- sheet1_data[missed_idx, sheet_1_cols]
        colnames(missed_rows_b) <- sheet1_names
        missed_rows_b[missed_rows_b == ""] <- '"BLANK"'
        msg_df(missed_rows_b)
    }
    norm_rows <- which(mainSheet$Paired_Normal == "")
    first_page_TM <- sheet1_data$`TM #`[sheet1_data$`TM #` != "0" &
                                            sheet1_data$`TM #` != "" &
                                            !is.na(sheet1_data$`TM #`)
    ]
    beaker_tm_Tumor <- unique(mainSheet$TM_Number[mainSheet$TM_Number != "" & mainSheet$TM_Number != "0" &
                                                      !is.na(mainSheet$TM_Number)])
    beaker_tm_Norm <- mainSheet$Case_ID[norm_rows]
    control_tm <- grepl(beaker_tm_Norm, pattern = "HAPMAP|NTC|SERA")
    beaker_tm_Norm <- beaker_tm_Norm[!control_tm]
    beaker_tm_Norm <- unique(stringr::str_split_fixed(beaker_tm_Norm, "_", 2)[,1])
    beaker_tm <- c(beaker_tm_Tumor, beaker_tm_Norm)
    missed_tm <- first_page_TM %in% beaker_tm
    if (any(missed_tm == FALSE)) {
        message(mkRed("!!!Potential Fillers:"))
        missing_tm <- first_page_TM[missed_tm == FALSE]
        missed_idx <- which(sheet1_data$`TM #` %in% missing_tm)
        missed_rows_tm <- sheet1_data[missed_idx, sheet_1_cols]
        colnames(missed_rows_tm) <- sheet1_names
        missed_rows_tm[missed_rows_tm == ""] <- '"BLANK"'
        msg_df(missed_rows_tm)
    }
}



fix_filler_zeros <- function(sheet1_data, mainSheet) {
    is_filler <- stringr::str_detect(sheet1_data$`Type & Tissue`,
                                     pattern = "(?i)fill|val|qc")

    check_extra_cases(sheet1_data, mainSheet)

    if (any(is_filler)) {
        message("Run has filler cases to be output with 0_ below:")
        msg_df(sheet1_data[is_filler, 1:5])
        tm_vals <- sheet1_data$`TM #`[is_filler]
        dna_vals <- sheet1_data$`DNA #`[is_filler]

        matches <- logical(length(mainSheet$Sample_ID))

        for (i in seq_along(mainSheet$Sample_ID)) {
            sample_id <- mainSheet$Sample_ID[i]

            tm_match <- stringr::str_detect(sample_id, stringr::fixed(tm_vals))
            dna_match <- stringr::str_detect(sample_id, stringr::fixed(dna_vals))

            matches[i] <- any(tm_match & dna_match, na.rm = TRUE)
        }

        replace0 <- which(matches)

        for (idx in replace0) {
            curr_vals <- mainSheet[idx, "Sample_ID"]
            new_end <- stringr::str_split_fixed(curr_vals, "_", 2)[, 2]
            new_zeros <- paste0("0_", new_end)
            mainSheet[idx, "Sample_ID"] <- new_zeros
            normals2fix <- grepl(curr_vals[[1]], mainSheet[, "Paired_Normal"] )
            if (any(normals2fix)) {
                curr_norm <- mainSheet[normals2fix, "Paired_Normal"]
                new_end <- stringr::str_split_fixed(curr_norm, "_", 2)[, 2]
                new_norm <- paste0("0_", new_end)
                mainSheet[normals2fix, "Paired_Normal"] <- new_norm
            }
        }
    }
    return(mainSheet)
}


make_zero <- function(mainSheet, col_name) {
    to_zero <- which(mainSheet[, col_name] == "")
    mainSheet[to_zero, col_name] <- 0
    return(mainSheet)
}


build_sam_rows <- function(zero_sam_rows, rawData) {
    empty_matrix <- matrix(nrow = nrow(zero_sam_rows), ncol = ncol(rawData))
    empty_df <- as.data.frame(empty_matrix)
    colnames(empty_df) <- colnames(rawData)
    con_sams <- paste(0, RUNID, zero_sam_rows$`Surgical #`, zero_sam_rows$`DNA #`, sep = "_")

    empty_df$Sample_ID <- make.unique(con_sams, sep = "-")
    empty_df$I7_Index_ID <- zero_sam_rows$I7_Index_ID
    empty_df$index <- zero_sam_rows$index
    empty_df$Run_Number <- rawData$Run_Number[1:nrow(zero_sam_rows)]
    empty_df$Case_ID <- paste(zero_sam_rows$`Surgical #`, zero_sam_rows$`DNA #`, sep = "_")
    empty_df$TUMOR_CASE_ID_BLOCK <- 0
    empty_df$Tumor_Content <- 0
    to_fix <- is.na(empty_df)
    empty_df[to_fix] <- ""
    return(empty_df)
}


strip_newlines <- function(x) {
    if (is.data.frame(x)) {
        x[] <- lapply(x, function(col) {
            if (is.factor(col)) {
                col <- as.character(col)
            }
            if (is.character(col)) return(gsub("[\r\n]+", "", col, perl = TRUE))
            return(col)
        })
        return(x)
    }
    x[] <- gsub("[\r\n]+", "", x, perl = TRUE)
    return(x)
}


check_data_format <- function(mainSheet) {
    control <- mainSheet$Case_ID %in% c("SC_SERACARE", "NTC_H20", "NC_HAPMAP")
    toCheck <- mainSheet[!control,]
    total_hyphens <- stringr::str_count(toCheck$Date_of_Procedure, pattern = "-") == 3
    all_dates <- all(nchar(toCheck$Date_of_Procedure) == 15) & all(total_hyphens)
    if (!all_dates == TRUE) {
        message(crayon::bgRed("DATES FORMAT SEEMS INCORRECT:"))
        msg_df(toCheck$Date_of_Procedure)
    }

    not_zero <- toCheck$TUMOR_CASE_ID_BLOCK != "0"
    tumors <- toCheck[not_zero,]
    total_dash <- all(stringr::str_count(tumors$TUMOR_CASE_ID_BLOCK, pattern = "-") == 2)

    if (!total_dash == TRUE) {
        message(crayon::bgRed("TUMOR_CASE_ID_BLOCK FORMAT SEEMS INCORRECT:"))
        msg_df(tumors$TUMOR_CASE_ID_BLOCK)
    }

    total_dash <- all(stringr::str_count(tumors$Case_ID, pattern = "-") == 3)
    if (!total_dash == TRUE) {
        message(crayon::bgRed("Case_ID FORMAT SEEMS INCORRECT:"))
        msg_df(tumors$Case_ID)
    }

    total_underscores <- all(stringr::str_count(tumors$Sample_ID, pattern = "_") == 3)
    if (!total_underscores == TRUE) {
        message(crayon::bgRed("Sample_ID tumors UNDERSCORES FORMAT SEEMS INCORRECT:"))
        msg_df(tumors$Sample_ID)
    }

    normals <- toCheck[!not_zero,]

    total_underscores <- all(stringr::str_count(normals$Sample_ID, pattern = "_") == 3)
    if (!total_underscores == TRUE) {
        message(crayon::bgRed("Sample_ID normals UNDERSCORES FORMAT SEEMS INCORRECT:"))
        msg_df(normals$Sample_ID)
    }
}


WriteMainSheet <- function(mainSheet, sheetHead, out_end = "SampleSheet.csv") {
    csv_file <- paste(mainSheet[1, "Run_Number"], out_end, sep = "-")
    out_path <- file.path(fs::path_home(), "Desktop", csv_file)

    message(crayon::bgGreen("Writing output, check samplesheet here:"), "\n", out_path)

    if (any(mainSheet$Tumor_Type == "NA")) {
        mainSheet[mainSheet$Tumor_Type == "NA", "Tumor_Type"] <- ""
    }

    sheetHead <- strip_newlines(sheetHead)
    mainSheet <- strip_newlines(mainSheet)
    check_data_format(mainSheet)

    utils::write.table(sheetHead, sep = ",", file = out_path, row.names = FALSE,
    col.names = FALSE, quote = FALSE, na = "", eol = "\n")

    suppressWarnings(utils::write.table(
        mainSheet, sep = ",", file = out_path, row.names = FALSE,
        col.names = TRUE, append = TRUE, quote = FALSE, na = "", eol = "\n"))
    return(out_path)
}


BuildMainSheet <- function(sheet1_df, rawData, RUNID, pact_run) {
    toKeep <- grepl("cont|Cont", sheet1_df$`Type & Tissue`)
    stopifnot(any(toKeep))
    control_rows <- sheet1_df[toKeep, ]
    type_idx <- which(colnames(sheet1_df) == "Type & Tissue")
    tissue_idx <-  type_idx + 1

        fillers <- grepl("fill|val|qc", sheet1_df$`Type & Tissue`,
                     ignore.case = TRUE) |
        grepl("fill|val|qc", sheet1_df[, tissue_idx], ignore.case = TRUE)

        controls_df <- build_sam_rows(control_rows, rawData)
    if (any(fillers)) {
        message("FILLERS Detected for on following rows:")
        filler_rows <- sheet1_df[fillers, c(1:5, type_idx, tissue_idx)]

        msg_df(filler_rows)
        row.names(filler_rows) <- NULL

        fill_sams_df <- build_sam_rows(filler_rows, rawData)
        tumor_fillers <- grepl("Tumor|tumor", filler_rows$`Type & Tissue`)
        tumor_fill_dna <- filler_rows$`DNA #`[tumor_fillers]
        normal_fill_dna <- filler_rows$`DNA #`[!tumor_fillers]
        fill_sams_df$Tumor_DNA[tumor_fillers] <- tumor_fill_dna
        fill_sams_df$Normal_DNA[!tumor_fillers] <- normal_fill_dna
        fill_sams_df$Date_of_Procedure <- 0
        fill_sams_df$Tumor_Specimen_Site <- filler_rows[,ncol(filler_rows)]
        fill_sams_df$TM_Number <- filler_rows$`TM #`
        fill_sams_df$TUMOR_CASE_ID_BLOCK <- filler_rows$`Surgical #`
        fill_sams_df$Paired_Normal[!tumor_fillers] <- 0
        fill_sams_df$Paired_Normal[tumor_fillers] <- fill_sams_df$Sample_ID[!tumor_fillers]
        fillers_main <- rbind(fill_sams_df, controls_df)
        fillers_main[] <- lapply(fillers_main, function(x) gsub(",", "", x))
        fillers_main <- make_zero(fillers_main, "TUMOR_CASE_ID_BLOCK")
        fillers_main <- make_zero(fillers_main, "Tumor_DNA")
        fillers_main <- make_zero(fillers_main, "Normal_DNA")
        fillers_main$Sample_Project <- paste0(rawData$Sample_Project[1])
        csv_file <- paste("SampleSheet", "FILLERS.csv", sep = "-")
        out_path <- file.path(fs::path_home(), "Desktop", csv_file)
        sheetHead <- MakeFileHeader(pact_run)
        out_path <- WriteMainSheet(fillers_main, sheetHead, out_end = csv_file)
    }

    if (DEBUG == TRUE) {
        message("\n", crayon::bgBlue("Appending control rows:"))
        msg_df(controls_df)
    }

    mainSheet <- rbind(rawData, controls_df)
    mainSheet[] <- lapply(mainSheet, function(x) gsub(",", "", x))

    mainSheet <- make_zero(mainSheet, "TUMOR_CASE_ID_BLOCK")
    mainSheet <- make_zero(mainSheet, "Tumor_DNA")
    mainSheet <- make_zero(mainSheet, "Normal_DNA")
    mainSheet$Sample_Project <- paste0(rawData$Sample_Project[1])
    return(mainSheet)
}


GrabFirstSheet <- function(worksheetPath, section = "bottom") {
    shNames <- readxl::excel_sheets(worksheetPath)
    sh <- which(grepl("PACT-", shNames, ignore.case = TRUE))[1]
    sh <- ifelse(is.na(sh), 1, sh)

    rawData <- GetExcelData(worksheetPath, sh, shRange = "A6:O100", include_head = TRUE)

    lastRow <- which(rawData$`DNA #` == "Lot #")[1]
    if (is.na(lastRow)) {
        stop(mkRed("Cannot find ending of samplesheet!"))
    }
    if (section == "bottom") {
        sheetBottom <- rawData[c(lastRow:nrow(rawData)), ]
        return(sheetBottom)
    } else{
        sheetTop <- rawData[c(1:lastRow), ]
        toKeep <- sheetTop$`ul H2O to 35uL` != ""
        sheetTop <- sheetTop[toKeep,]
        rownames(sheetTop) <- NULL
        return(sheetTop)
    }
}


Err_runID <- function(RUNID) {
    err_msg <- 'Keyword "Run ID:" not found in SampleSheet "DNA #" column'
    new_id <- paste('\nDefaulting to script INPUT RUNID:', RUNID)
    message(mkRed(err_msg), new_id)
}


GrabRunNumber <- function(worksheetPath, RUNID) {
    sheetBottom <- GrabFirstSheet(worksheetPath, section = "bottom")
    runID_row <- which(stringr::str_detect(sheetBottom$`DNA #`, "Run ID:"))

    # If 'Run ID:' keyword not found, use provided RUNID as default
    if (length(runID_row) != 1) {
        message(mkRed("Run ID: not found in column 'DNA #, trying column 8..."))
        runID_row <- which(stringr::str_detect(sheetBottom[, 8], "Run ID:"))
        if (length(runID_row) != 1) {
            Err_runID(RUNID)
            return(RUNID)
        } else{
            message("runID_row found in column 8! row: ", runID_row)
        }
    }

    runID_find <- paste(sheetBottom[runID_row, ])
    runID_col <- which(stringr::str_detect(runID_find, key_run_ids))

    if (length(runID_col) == 0) {
        Err_runID(RUNID)
        return(RUNID)
    }

    run_number <- paste(sheetBottom[runID_row, runID_col])
    run_number <-
        stringr::str_replace_all(run_number, c("Run ID:" = "", " " = ""))

    return(run_number)
}


GetRawSamplesheet <- function(worksheetPath, sheet_kwd = "Bioinformatics") {
    shNames <- readxl::excel_sheets(worksheetPath)
    if (DEBUG == TRUE) {
        msg_df(data.frame(`Sheet names in Workbook` = shNames))
    }

    sh <- which(grepl(sheet_kwd, shNames, ignore.case = TRUE))[1]
    if (is.na(sh)) {
        sh = length(shNames)
        message(mkRed(paste0('"', sheet_kwd, '"'),
            'not in Excel sheet names, defaulting to last sheet:', sh))
    }
    msgRd <- paste0('Reading Excel Sheet named \"', shNames[sh], '\" from file:')
    message(msgRd,'\n', worksheetPath, '\n')
    rawData <- GetExcelData(worksheetPath, sh, shRange = "A1:O200", include_head = TRUE)
    toDrop <- which(rawData[, 1] == "")[1] - 1
    rawData <- rawData[1:toDrop, ]

    return(rawData)
}


# Checks if the user input RunID matches worksheet RunID ----------------------
CheckRunIDMatch <- function(sheetRunID, RUNID) {
    if (sheetRunID != RUNID) {
        message(mkRed("SampleSheet and Input RUNID do not match!"))
        message(sprintf("Samplesheet RunID: %s\nUser input RunID:  %s",
                        crayon::bgBlue(sheetRunID), crayon::bgGreen(RUNID)))
        message("Using RUNID from inside worksheet: ", sheetRunID)
        return(sheetRunID)
    }
    return(RUNID)
}


# Parses File when no samplesheet.csv tab is detected -------------------------
parse_workbook <- function(worksheetPath, RUNID) {
    rawData <- GetRawSamplesheet(worksheetPath)
    rawData$Sample_Project <- INPUT
    sheetRunID <- GrabRunNumber(worksheetPath, RUNID)
    beakerData <- GetExcelData(worksheetPath, "BeakerExport", shRange = NULL, include_head = TRUE)
    RUNID <- CheckRunIDMatch(sheetRunID, RUNID)
    if (!"Sample_Project" %in% colnames(rawData)) {
        stop("Input worksheet is outdated or missing Sample_Project header:\n", worksheetPath)
    }

    sheet_pact_run <- paste0(beakerData$`PACT Batch Number`[1])

    if (sheet_pact_run != INPUT) {
        warning("Sheet PACT_ID is: ", sheet_pact_run, "\nINPUT PACT_ID is: ", INPUT)
    }

    pact_run <- INPUT

    rawData$Run_Number <- RUNID
    message(crayon::bgBlue(paste("PACT ID is:", pact_run)), "\n")
    sheet1_df <- GetExcelData(worksheetPath, 1, shRange = "B6:O100", include_head = TRUE)
    mainSheet <- BuildMainSheet(sheet1_df, rawData, sheetRunID, pact_run)
    sheet1_data <- GrabFirstSheet(worksheetPath, section = "Top")
    mainSheet <- fix_filler_zeros(sheet1_data, mainSheet)
    sheet_pact_id <- mainSheet$Sample_Project[1]
    if (sheet_pact_id != INPUT) {
        warning("Sheet PACT_ID is: ", sheet_pact_id, "\nINPUT PACT_ID is: ", INPUT)
    }
    mainSheet$Sample_Project <- INPUT
    mainSheet$Sample_Project <- mainSheet$Sample_Project[1]
    mainSheet$Date_of_Procedure <- stringr::str_replace_all(
        mainSheet$Date_of_Procedure, pattern = " ", replacement = "-")
    mainSheet$Tumor_Specimen_Site <- stringr::str_replace_all(
        mainSheet$Tumor_Specimen_Site, pattern = " ", replacement = "-")
    mainSheet$Tumor_Type <- stringr::str_replace_all(
        mainSheet$Tumor_Type, pattern = " ", replacement = "-")

    return(mainSheet)
}

# Uploads the output CSV samplesheet to REDCap and sends email notification with the attachment ---
pushToRedcap <- function(outVals, API_TOKEN = NULL) {
    stopifnot(!is.null(API_TOKEN) & !is.null(outVals))
    runID <- outVals[["RUNID"]]
    outFile <- outVals[["outFile"]]
    pactID <- stringr::str_to_upper(outVals[["PACT_ID"]])

    stopifnot(is.character(runID), is.character(outFile), length(outFile) == 1L, file.exists(outFile))

    record <- c(record_id = pactID, pact_run_number = runID)
    data <- jsonlite::toJSON(list(as.list(record)), auto_unbox = TRUE)

    # Create the record entry first
   result <- RCurl::postForm(
       API_URL, token = API_TOKEN, content = "record", format = "json",
       type = "flat", data = data)

    csv_import_status <- "400"

    # Import the CSV file after creating the record_id
    if (file.exists(outFile)) {
        csv_import_res <- httr::POST(
            url = API_URL, body = list(
                token = API_TOKEN, content = "file", action = "import", record = pactID,
                field = "pact_csv_sheet", returnFormat = "json",
                file = httr::upload_file(path = outFile, type = "text/csv")), encode = "multipart")
        httr::warn_for_status(csv_import_res)
        csv_import_status <- httr::status_code(csv_import_res)
    }

    # Import the record data to trigger email with uploaded csv attached if successful
    if (csv_import_status != "200") return(message("httr POST of csv file failed"))
    record <- c(record_id = pactID, pact_run_number = runID, pact_csv_email = "pact_csv_email")
    data <- jsonlite::toJSON(list(as.list(record)), auto_unbox = TRUE)
    notify_response <- httr::POST(url = API_URL, body = list(
        token = API_TOKEN, content = "record", format = "json", type = "flat", data = data
        ), encode = "form"
    )

    httr::warn_for_status(notify_response)

}


# Parses xlsx file and writes as csv file -------------------------------------
parseExcelFile <- function(worksheetPath, RUNID = NULL) {
    mainSheet <- parse_workbook(worksheetPath, RUNID)
    mainSheet <- sanitizeSheet(mainSheet)
    PACT_ID <- INPUT #mainSheet[1, "Sample_Project"]
    sheetHead <- MakeFileHeader(PACT_ID)
    outFile <- WriteMainSheet(mainSheet, sheetHead)
    outVals <- c(RUNID = mainSheet[1, "Run_Number"],
                 outFile = outFile, PACT_ID = mainSheet[1, "Sample_Project"])
    return(outVals)
}


# Filters list of possible files in the directory for worksheet ---------------
filterFiles <- function(potentialFi) {
    wbFiles <- stringr::str_which(basename(potentialFi), pattern = "xlsm")
    if (length(wbFiles) == 0) {
        message("No .xlsm worksheet found. Checking if .xlsx file exists:\n")
        wbFiles <- stringr::str_which(basename(potentialFi), pattern = "book")
    }
    potentialFi <- potentialFi[wbFiles]
    potentialFi <- potentialFi[!stringr::str_detect(potentialFi, "\\$")]
    return(potentialFi)
}


CheckOtherFiles <- function(worksheetPath, RUNID) {
    notFoundMsg <- crayon::bgRed("The PACT run worksheet was not found:")
    parentFolder <- dirname(worksheetPath)
    potentialFi <- list.files(path = parentFolder, full.names = TRUE)
    if (length(potentialFi) > 1) {
        message(mkRed("Checking other existing files:"), "\n")
        message(paste(potentialFi, sep = "\n", collapse = "\n"))
        potentialFi <- filterFiles(potentialFi)
    } else{
        stop(paste("No PACT worksheet was found in the directory:", parentFolder))
    }
    if (file.exists(potentialFi[1])) {
        message(crayon::bgGreen("Now trying to read:"), "\n", potentialFi[1], "\n")
        pfile <- potentialFi[1]
    } else{
        nextMsg <-
            paste(notFoundMsg, potentialFi[1], "Trying:", potentialFi[2], sep = "\n")
        message(nextMsg)
        pfile <- potentialFi[2]
    }
    stopifnot(file.exists(pfile))
    outVals <- parseExcelFile(pfile, RUNID)
    return(outVals)
}


# FUN: Returns run type based on the input pact_id or worksheet name ----------
determineRunType <- function(INPUT) {

    runType <- dplyr::case_when(
        stringr::str_detect(INPUT, "2000|NextSeq550") ~ "NextSeq2000",
        stringr::str_detect(INPUT, "TEST")            ~ "test",
        stringr::str_detect(INPUT, "ILC")             ~ "Illumina",
        stringr::str_detect(INPUT, "TMB")             ~ "TMB",
        stringr::str_detect(INPUT, "^\\d{2}")         ~ "Sophia",
        TRUE                                          ~ "regular"
    )

    return(runType)
}


# FUN: Returns Path to xlsx file ----------------------------------------------
getExcelPath <- function(INPUT, runType) {
    if (stringr::str_detect(INPUT, .Platform$file.sep)) return(INPUT)

    lab_only = file.path("", "Volumes", "molecular", "MOLECULAR LAB ONLY")
    folder = file.path("NYU PACT Patient Data", "Workbook")

    run_year <- stringr::str_split_fixed(INPUT, "-", 3)[, 2]
    if (runType == "Sophia") {
        run_year <- stringr::str_split_fixed(INPUT, "-", 3)[, 1]
    }
    yearDir <- paste0("20", run_year)
    xlFi <- paste0(INPUT, ".xlsm")
    worksheetPath <- file.path(lab_only, folder, yearDir, INPUT, xlFi)

    if (runType == "NextSeq2000") {
        folder <- "Validations/PACT new i7-NextSeq2000/Wet Lab/Workbook"
        worksheetPath <- file.path(lab_only, folder, yearDir, xlFi)
    }
    if (runType == "test") {
        zdrive = "/Volumes/molecular/Molecular/Validation/PACT/Test_Sheets"
        worksheetPath <- file.path(zdrive, xlFi)
    }
    if (runType == "Illumina") {
        zdrive = "/Volumes/molecular/Molecular/Validation/PACT/Test_Sheets"
        worksheetPath <- file.path(zdrive, xlFi)
    }
    if (runType == "TMB") {
        zdrive = file.path(lab_only, "Validations/TMB")
        worksheetPath <- file.path(zdrive, xlFi)
    }
    if (runType == "Sophia") {
        worksheetPath <- file.path(lab_only, folder, yearDir, INPUT, xlFi)
        message("Run type is Sophia, looking for workbook in:")
    }
    message(worksheetPath, "\n")
    return(worksheetPath)
}


# Gets dataframe and saves as CSV file ----------------------------------------
writeSampleSheet <- function(INPUT, token = NULL, RUNID = NULL) {
    outVals <- NULL
    runType <- determineRunType(INPUT)
    worksheetPath <- getExcelPath(INPUT, runType)
    if (file.exists(worksheetPath)) {
        outVals <- parseExcelFile(worksheetPath, RUNID)
    } else {
        outVals <- CheckOtherFiles(worksheetPath, RUNID)
    }
    if (!is.null(outVals) && !is.null(token)) {
        pushToRedcap(outVals, token)
    }
}

# MAIN Execution start --------------------------------------------------------
loadPacks()
checkMounts()

writeSampleSheet(INPUT, API_TOKEN, RUNID)
