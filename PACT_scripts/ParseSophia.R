#!/usr/bin/env Rscript
## Script name: pactParse.R
## Purpose: source of global scripts and generate PACT -SampleSheet.csv file
## Author: Jonathan Serrano
## Date Created: February 10, 2026
## Version: 1.0.1
## Copyright (c) NYULH Jonathan Serrano, 2026

library("base"); args <- commandArgs(TRUE); gb <- globalenv(); assign("gb", gb)

key_run_ids <- "NB551709|NB501073|VH01471"

# Function to assign parameter values based on args ---------------------------
assign_param <- function(args, idx) {
    if (!is.na(args[[idx]]) && !identical(args[[idx]], "")) return(args[[idx]])
    return(NULL)
}

# Assign values to variables & print ------------------------------------------
INPUT <- assign_param(args, 1)
RUNID <- assign_param(args, 2)
token <- assign_param(args, 3)

# Functions to print parameters  ) --------------------------------------------
brick_message <- function(phrase) {
    dsh <- paste(rep("-", nchar(phrase)), collapse = "")
    message(paste(dsh, phrase, dsh, sep = "\n"))
}

print_parameters <- function(INPUT, RUNID, token) {
    brick_message("Parameters INPUT")
    message("\nINPUT:  ", c(INPUT, "NULL")[1], "\nRUNID:  ",
            c(RUNID, "NULL")[1], "\n", "token:  ", c(token, "NULL")[1], "\n")
}

print_parameters(INPUT, RUNID, token)

# Helper functions to message as data frame and red background ----------------
MsgDF <- function(...) {
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
    message("Z-drive path is accessible")
}

options(repos = c(getOption("repos"), CRAN = "http://cran.us.r-project.org"))

# Function to silently load library without conflict warnings -----------------
libLoad <- function(pkg) {
    if (pkg %in% loadedNamespaces()) return(message(pkg, "...loaded successful"))
    suppressWarnings(suppressPackageStartupMessages(
        library(pkg, character.only = TRUE, warn.conflicts = FALSE, quietly = TRUE)
    ))
    if (pkg %in% loadedNamespaces()) message(pkg, "...loaded successful")
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
    pkgs <- c("foreach", "RCurl", "jsonlite", "digest", "readxl", "stringr", "httr", "fs",
              "crayon", "dplyr", "tibble", "data.table")
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
                message(mkBlu(sprintf(
                    'The following "%s" will be removed from the column: %s',
                                patt, replace_col)))
                MsgDF(replaced_df)
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
        MsgDF(mainSheet$Sample_Name[controlSamples])
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
        MsgDF(mainSheet$Sample_ID[duplicated(mainSheet$Sample_ID)])
    }

    mainSheet[] <- lapply(mainSheet, function(col) {
        col <- gsub(" ", "-", as.character(col))
        return(col)
    })

    if (any( stringr::str_detect(pattern = "...1", colnames(mainSheet)))) {
        toDrop <- min(which(stringr::str_detect(pattern = "...1", colnames(mainSheet)))) - 1
        mainSheet <- mainSheet[, c(1:toDrop)]
    }

    return(mainSheet)
}


GetExcelData <- function(worksheetPath, sheetNum, shRange, toSkip = 0, include_head = FALSE) {
    sheetData <- suppressMessages(as.data.frame(
        readxl::read_excel(
            worksheetPath,
            sheet = sheetNum,
            na = "",
            range = shRange,
            col_types = "text",
            col_names = include_head,
            skip = toSkip
        )))
    sheetData[is.na(sheetData)] <- ""
    return(sheetData)
}


GetPactName <- function(worksheetPath) {
    sheet_nam <- readxl::excel_sheets(worksheetPath)
    beaker_sh <- which(grepl("Beaker", sheet_nam))
    beaker_df <- suppressWarnings(as.data.frame(
        readxl::read_excel(worksheetPath, sheet = beaker_sh, col_names = TRUE)))
    pact_col <- which(grepl("batch", colnames(beaker_df), ignore.case = TRUE))
    pact_name <- beaker_df[1, pact_col]
    message("pact_name: ", pact_name)
    return(pact_name)
}


WriteFileHeader <- function(worksheetPath) {
    pact_name <- GetPactName(worksheetPath)
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
        pact_name,
        pact_name,
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
    return(theHeader)

}


fix_filler_zeros <- function(sheet1_data, mainSheet) {
    is_filler <- stringr::str_detect(sheet1_data$`Type & Tissue`,
                                     pattern = "(?i)filler|val")
    if (any(is_filler)) {
        message("Run has filler cases to be output with 0_ below:")
        MsgDF(sheet1_data[is_filler, 1:5])
        tm_vals <- sheet1_data$`TM #`[is_filler]
        dna_vals <- sheet1_data$`DNA #`[is_filler]

        replace0 <- which(vapply(
            mainSheet$Sample_ID,
            function(x) {
                any(
                    stringr::str_detect(x, stringr::fixed(tm_vals)) &
                        stringr::str_detect(x, stringr::fixed(dna_vals))
                )
            },
            logical(1)
        ))
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


BuildMainSheet <- function(sheet1_df, rawData, RUNID, pact_run) {
    toKeep <- grepl("cont|Cont", sheet1_df$`Type & Tissue`)
    control_rows <- sheet1_df[toKeep,]

    empty_matrix <- matrix(nrow = nrow(control_rows), ncol = ncol(rawData))
    empty_df <- as.data.frame(empty_matrix)
    colnames(empty_df) <- colnames(rawData)
    con_sams <- paste(0, RUNID, control_rows$`Surgical #`, control_rows$`DNA #`, sep = "_")

    empty_df$Sample_ID <- make.unique(con_sams, sep = "-")
    empty_df$I7_Index_ID <- control_rows$I7_Index_ID
    empty_df$index <- control_rows$index
    empty_df$Run_Number <- rawData$Run_Number[1:nrow(control_rows)]
    empty_df$Case_ID <- paste(control_rows$`Surgical #`, control_rows$`DNA #`, sep = "_")
    empty_df$TUMOR_CASE_ID_BLOCK <- 0
    empty_df$Tumor_Content <- 0
    to_fix <- is.na(empty_df)
    empty_df[to_fix] <- ""
    message(crayon::bgBlue("Appending control rows:"))
    message(paste0(capture.output(empty_df), collapse = "\n"))
    mainSheet <- rbind(rawData, empty_df)
    mainSheet[] <- lapply(mainSheet, function(x) gsub(",", "", x))

    mainSheet <- make_zero(mainSheet, "TUMOR_CASE_ID_BLOCK")
    mainSheet <- make_zero(mainSheet, "Tumor_DNA")
    mainSheet <- make_zero(mainSheet, "Normal_DNA")
    mainSheet$Sample_Project <- paste0(rawData$Sample_Project[1])
    return(mainSheet)
}


GetSheetHeading <- function(worksheetPath) {
    sheetHead <- WriteFileHeader(worksheetPath)
    sheetHead <- rbind(sheetHead, c("", ""), c("[Data]", ""))
    return(sheetHead)
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
    sheetBottom <- GrabFirstSheet(worksheetPath)
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


GetRawSamplesheet <- function(worksheetPath) {
    shNames <- readxl::excel_sheets(worksheetPath)
    sh <- which(grepl("Bioinformatics", shNames, ignore.case = TRUE))[1]
    if (is.na(sh)) {
        message(mkRed(
            '"Bioinformatics" not in Excel sheet names, defaulting to last sheet'))
        sh = length(shNames)
    }
    msgRd <- paste0('Reading Excel Sheet named \"', shNames[sh], '\" from file:')
    message(crayon::bgGreen(msgRd),'\n', worksheetPath)
    rawData <- GetExcelData(worksheetPath, sh, shRange = "A1:O200", include_head = TRUE)
    toDrop <- which(rawData[, 1] == "")[1] - 1
    rawData <- rawData[1:toDrop, ]

    return(rawData)
}


WriteMainSheet <- function(mainSheet, sheetHead) {
    csv_file <- paste(mainSheet[1, "Run_Number"], "SampleSheet.csv", sep = "-")
    out_path <- file.path(fs::path_home(), "Desktop", csv_file)

    message(crayon::bgGreen("Writing output, check your samplesheet here:"),
            "\n", out_path)

    if (any(mainSheet$Tumor_Type == "NA")) {
        mainSheet[mainSheet$Tumor_Type == "NA", "Tumor_Type"] <- ""
    }

    write.table(sheetHead, sep = ",", file = out_path, row.names = FALSE,
        col.names = FALSE, quote = FALSE)

    suppressWarnings(
        write.table(
            mainSheet,
            sep = ",",
            file = out_path,
            row.names = FALSE,
            col.names = TRUE,
            append = TRUE,
            quote = FALSE
        )
    )
    return(out_path)
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
    sheetRunID <- GrabRunNumber(worksheetPath, RUNID)
    beakerData <- GetExcelData(worksheetPath, "BeakerExport", shRange = NULL, include_head = TRUE)
    RUNID <- CheckRunIDMatch(sheetRunID, RUNID)
    if (!"Sample_Project" %in% colnames(rawData)) {
        stop("Input worksheet is outdated or missing Sample_Project header:\n", worksheetPath)
    }
    pact_run <- paste0(beakerData$`PACT Batch Number`[1])
    rawData$Run_Number <- RUNID
    message(crayon::bgBlue(paste("PACT ID is:", pact_run)))
    sheet1_df <- GetExcelData(worksheetPath, 1, shRange = "B6:O100", include_head = TRUE)
    mainSheet <- BuildMainSheet(sheet1_df, rawData, sheetRunID, pact_run)
    sheet1_data <- GrabFirstSheet(worksheetPath, section = "Top")
    mainSheet <- fix_filler_zeros(sheet1_data, mainSheet)
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
pushToRedcap <- function(outVals, token = NULL) {
    stopifnot(!is.null(token) & !is.null(outVals))

    runID <- outVals[["RUNID"]]
    outFile <- outVals[["outFile"]]
    pactID <- outVals[["PACT_ID"]]

    stopifnot(is.character(runID), is.character(outFile), length(outFile) == 1L, file.exists(outFile))
    api_url <- "https://redcap.nyumc.org/apps/redcap/api/"

    record <- c(record_id = pactID, pact_run_number = runID, pact_csv_email = "pact_csv_email")
    data <- jsonlite::toJSON(list(as.list(record)), auto_unbox = TRUE)

    # Import the CSV file to REDCap first
    csv_import_res <- httr::POST(
        url = api_url,
        body = list(
            token = token,
            content = "file",
            action = "import",
            record = pactID,
            field = "pact_csv_sheet",
            returnFormat = "json",
            file = httr::upload_file(path = outFile, type = "text/csv")
        ), encode = "multipart"
    )

    csv_import_status <- httr::status_code(csv_import_res)
    print(csv_import_status)

    if (httr::http_error(csv_import_res)) {
        stop("REDCap file import failed with HTTP status ", csv_import_status)
    }

    # Import the record data to trigger email with uploaded csv attached
    notify_response <- httr::POST(
        url = api_url,
        body = list(
            token = token,
            content = "record",
            format = "json",
            type = "flat",
            data = data
        ),
        encode = "form"
    )

    notify_status <- httr::status_code(notify_response)
    print(notify_status)
}


# Parses xlsx file and writes as csv file -------------------------------------
parseExcelFile <- function(worksheetPath, RUNID = NULL) {
    shNames <- readxl::excel_sheets(worksheetPath)
    MsgDF(data.frame(`Sheet names in Workbook` = shNames))
    sheetHead <- GetSheetHeading(worksheetPath)
    mainSheet <- parse_workbook(worksheetPath, RUNID)
    mainSheet <- sanitizeSheet(mainSheet)
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
    message(worksheetPath)
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

writeSampleSheet(INPUT, token, RUNID)
