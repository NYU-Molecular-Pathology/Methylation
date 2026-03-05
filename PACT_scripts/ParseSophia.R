#!/usr/bin/env Rscript
## Script name: pactParse.R
## Purpose: source of global scripts and generate PACT -SampleSheet.csv file
## Author: Jonathan Serrano
## Date Created: February 10, 2026
## Version: 1.0.1
## Copyright (c) NYULH Jonathan Serrano, 2026

library("base"); args <- commandArgs(TRUE); gb <- globalenv(); assign("gb", gb)

# Function to assign parameter values based on args ---------------------------
assign_param <- function(args, idx) {
    if (idx == 3) {
        if (length(args) >= idx) {
            return(args[idx])
        } else {
            return("NONE")
        }
    }
    if (length(args) >= idx) args[idx] else stop("arg ", idx, " is NULL")
}

# Functions to print parameters  ) --------------------------------------------
brick_message <- function(phrase) {
    dsh <- paste(rep("-", nchar(phrase)), collapse = "")
    message(paste(dsh, phrase, dsh, sep = "\n"))
}

print_parameters <- function(token, input, runID) {
    brick_message("Parameters input")
    message("token:   ", token)
    message("PACT ID: ", input)
    message("Run ID:  ", runID, "\n")
}

# Assign values to variables & print ------------------------------------------
token <- assign_param(args, 1)
input <- assign_param(args, 2)
runID <- assign_param(args, 3)
print_parameters(token, input, runID)

# Helper functions to message as data frame and red background ----------------
MsgDF <- function(...) {
    return(message(paste0(
        capture.output(as.data.frame(...)), collapse = "\n"
    )))
}


mkRed <- function(...) {return(crayon::white$bgRed$bold(...))}

# FUN: Checks if z-drive is accessible to the RScript -------------------------
checkMounts <- function() {
    molecDrive = "/Volumes/molecular/MOLECULAR LAB ONLY"
    zDrive = "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular"
    failMount <- ifelse(dir.exists(molecDrive), T, F)
    if (failMount != TRUE) {
        mountMsg <- paste(
            "PATH does not exist, ensure path is mounted:", mkRed(molecDrive),
            "You must mount the network Z-drive path:", mkRed(zDrive),
            sep = "\n"
        )
        message("\n", mountMsg, "\n")
        stopifnot(!any(failMount == TRUE))
    } else {
        message("Z-drive path is accessible")
    }
}

options(repos = c(getOption("repos"), CRAN = "http://cran.us.r-project.org"))

# Function to silently load library without conflict warnings -----------------
libLoad <- function(pkg) {
    suppressWarnings(suppressPackageStartupMessages(
        library(pkg, character.only = TRUE, warn.conflicts = FALSE, quietly = TRUE)
    ))
    message(pkg, "...loaded successful")
}

# Installs package(s) if not already installed then loads the library ---------
SilentInstall <- function(pkgs) {
    for (pkg in pkgs) {
        if (!pkg %in% rownames(installed.packages()))
            install.packages(pkg, dependencies = TRUE, ask = FALSE)
        libLoad(pkg)
    }
}


# Functions to load or install missing required packages ----------------------
loadPacks <- function() {
    SilentInstall("devtools")
    pkgs <- c("foreach", "jsonlite", "RCurl", "readxl", "stringr", "httr", "fs",
              "crayon", "dplyr", "tibble", "data.table")
    SilentInstall(pkgs)
}


MsgChangesMade <- function(mainSheet, badChars = c("[\r\n]", ",", " ")) {
    for (ptrn in badChars) {
        for (i in seq_along(mainSheet)) {
            col <- ifelse(is.na(mainSheet[[i]]), "", mainSheet[[i]])
            idxMatch <- which(stringr::str_detect(col, pattern = ptrn))
            if (length(idxMatch) > 0) {
                message(sprintf('The following "%s" will be removed from the column: %s',
                                ptrn, names(mainSheet)[i]))
                MsgDF(mainSheet[idxMatch, i])
            }
        }
    }
}


CountControls <- function(mainSheet,
                          controlNames = "NTC_H20|SC_SERACARE|NC_HAPMAP") {
    controlSamples <- grepl(pattern = controlNames, mainSheet$Sample_ID)
    if (table(controlSamples)[['TRUE']] != 3) {
        message(crayon::bgRed("There are not 3 Control samples"))
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

    if (any(duplicated(mainSheet$Sample_ID))) {
        message(crayon::bgRed("There are duplicated Sample_ID in the SampleSheet"))
        MsgDF(mainSheet$Sample_ID[duplicated(mainSheet$Sample_ID)])
    }

    mainSheet[] <- lapply(mainSheet, function(col) {
        col <- gsub(" ", "-", as.character(col))
        return(col)
    })

    return(mainSheet)
}


GetExcelData <- function(worksheetPath, sheetNum, shRange, toSkip = 0, cm = FALSE) {
    sheetData <- suppressMessages(as.data.frame(
        readxl::read_excel(
            worksheetPath,
            sheet = sheetNum,
            na = "",
            range = shRange,
            col_types = "text",
            col_names = cm,
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


fix_filler_zeros <- function(rawData, mainSheet) {
    is_filler <- stringr::str_detect(rawData$`Type & Tissue`,
                                     pattern = "(?i)filler")
    is_val <- is_filler <- stringr::str_detect(rawData$`Type & Tissue`,
                                               pattern = "(?i)val")
    is_filler <- is_filler | is_val
    if (any(is_filler)) {
        message("Run has filler cases to be output with 0_ below:")
        MsgDF(rawData[is_filler, 1:5])
        to_be_zero <- rawData$Sample_ID[is_filler]
        replace0 <- which(mainSheet$Sample_ID %in% to_be_zero)
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


BuildMainSheet <- function(sheet1_df, rawData, runID, pact_run) {
    toKeep <- grepl("cont|Cont", sheet1_df$`Type & Tissue`)
    control_rows <- sheet1_df[toKeep,]

    empty_matrix <- matrix(nrow = nrow(control_rows), ncol = ncol(rawData))
    empty_df <- as.data.frame(empty_matrix)
    colnames(empty_df) <- colnames(rawData)
    con_sams <- paste(0, runID, control_rows$`Surgical #`, control_rows$`DNA #`, sep = "_")

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
    return(mainSheet)
}


GetSheetHeading <- function(worksheetPath) {
    sheetHead <- WriteFileHeader(worksheetPath)
    sheetHead <- rbind(sheetHead, c("", ""), c("[Data]", ""))
    return(sheetHead)
}


GrabSheetBottom <- function(worksheetPath) {
    shNames <- readxl::excel_sheets(worksheetPath)
    sh <- which(grepl("PACT-", shNames, ignore.case = TRUE))[1]
    sh <- ifelse(is.na(sh), 1, sh)

    rawData <- GetExcelData(worksheetPath, sh, shRange = "A6:O100", cm = TRUE)

    lastRow <- which(rawData$`DNA #` == "Lot #")[1]
    if (is.na(lastRow)) {
        stop(mkRed("Cannot find ending of samplesheet!"))
    }
    sheetBottom <- rawData[c(lastRow:nrow(rawData)), ]
    return(sheetBottom)
}


Err_runID <- function(runID) {
    err_msg <- crayon::bgRed(
        'Keyword "Run ID:" not found in SampleSheet "DNA #" column'
    )
    new_id <- paste('\nDefaulting to script input RUNID:', runID)
    message(err_msg, new_id)
}


GrabRunNumber <- function(worksheetPath, runID) {
    key_run_ids <- "NB551709|NB501073|VH01471"
    sheetBottom <- GrabSheetBottom(worksheetPath)
    runID_row <- which(stringr::str_detect(sheetBottom$`DNA #`, "Run ID:"))

    # If 'Run ID:' keyword not found, use provided runID as default
    if (length(runID_row) != 1) {
        message(mkRed("Run ID: not found in column 'DNA #, trying column 8..."))
        runID_row <- which(stringr::str_detect(sheetBottom[, 8], "Run ID:"))
        if (length(runID_row) != 1) {
            Err_runID(runID)
            return(runID)
        } else{
            message("runID_row found in column 8! row: ", runID_row)
        }
    }

    runID_find <- paste(sheetBottom[runID_row, ])
    runID_col <- which(stringr::str_detect(runID_find, key_run_ids))

    if (length(runID_col) == 0) {
        Err_runID(runID)
        return(runID)
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
        miss_msg <-
            '"Bioinformatics" not detected in Excel sheetnames, defaulting to last sheet'
        message(crayon::bgRed(miss_msg))
        sh = length(shNames)
    }
    msgRd <-
        paste0('Reading Excel Sheet named \"', shNames[sh], '\" from file:')
    message(crayon::bgGreen(msgRd),'\n',worksheetPath)
    rawData <- GetExcelData(worksheetPath, sh, shRange = "A1:L100", cm = TRUE)
    toDrop <- which(rawData[, 1] == "")[1] - 1
    rawData <- rawData[1:toDrop, ]


    return(rawData)
}


WriteMainSheet <- function(mainSheet, sheetHead) {
    csv_file <- paste(mainSheet[1, "Run_Number"], "SampleSheet.csv", sep = "-")
    out_path <- file.path(fs::path_home(), "Desktop", csv_file)

    message("Writing output, check your samplesheet here:\n", out_path)

    if (any(mainSheet$Tumor_Type == "NA")) {
        mainSheet[mainSheet$Tumor_Type == "NA", "Tumor_Type"] <- ""
    }

    write.table(
        sheetHead,
        sep = ",",
        file = out_path,
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )

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
CheckRunIDMatch <- function(sheetRunID, runID) {
    if (sheetRunID != runID) {
        message(mkRed("SampleSheet and Input RUNID do not match!"))
        message(sprintf("Samplesheet RunID: %s\nUser input RunID:  %s",
                        crayon::bgBlue(sheetRunID), crayon::bgGreen(runID)))
        message("Using runID from inside worksheet: ", sheetRunID)
    }
}


# Parses File when no samplesheet.csv tab is detected -------------------------
parse_workbook <- function(worksheetPath, runID) {
    rawData <- GetRawSamplesheet(worksheetPath)
    sheetRunID <- GrabRunNumber(worksheetPath, runID)
    CheckRunIDMatch(sheetRunID, runID)
    pact_run <-  paste0(GetExcelData(worksheetPath, 1, shRange = "G3:G3", cm = FALSE)[1,1])
    message(crayon::bgBlue(paste("PACT ID is:", pact_run)))
    sheet1_df <- GetExcelData(worksheetPath, 1, shRange = "B6:O100", cm = TRUE)
    mainSheet <-
            BuildMainSheet(sheet1_df, rawData, sheetRunID, pact_run)
    mainSheet <- fix_filler_zeros(rawData, mainSheet)

    return(mainSheet)
}


pushToRedcap <- function(outVals, token = NULL) {
    stopifnot(!is.null(token) & !is.null(outVals))

    runID <- outVals[[1]]
    outFile <- outVals[[2]]
    stopifnot(is.character(runID), is.character(outFile), length(outFile) == 1L, file.exists(outFile))

    api_url <- "https://redcap.nyumc.org/apps/redcap/api/"
    to_json_data <- function(df) {
        jsonlite::toJSON(unname(split(df, seq_len(nrow(df)))), auto_unbox = TRUE)
    }

    # 1) create/update record + set pact_run_number
    rec_json <- to_json_data(
        data.frame(record_id = runID, pact_run_number = runID, stringsAsFactors = FALSE))

    httr::POST(
        api_url,
        body = list(
            token = token, content = "record", action = "import",
            format = "json", type = "flat",
            overwriteBehavior = "overwrite",
            returnContent = "nothing", returnFormat = "json",
            data = rec_json
        ),
        encode = "form", httr::timeout(120)
    )

    # 2) upload file to pact_csv_sheet
    httr::POST(
        api_url,
        body = list(
            token = token, content = "file", action = "import",
            record = runID, field = "pact_csv_sheet",
            file = httr::upload_file(outFile)
        ),
        encode = "multipart", httr::timeout(120)
    )

    # 3) set email flag pact_csv_email = "pact_csv_email"
    email_json <- to_json_data(data.frame(
        record_id = runID, pact_csv_email = "pact_csv_email", stringsAsFactors = FALSE))

    httr::POST(
        api_url,
        body = list(
            token = token, content = "record", action = "import",
            format = "json", type = "flat", overwriteBehavior = "overwrite",
            returnContent = "nothing", returnFormat = "json", data = email_json
        ),
        encode = "form", httr::timeout(120)
    )
}


# Parses xlsx file and writes as csv file -------------------------------------
parseExcelFile <- function(worksheetPath, runID = NULL) {
    shNames <- readxl::excel_sheets(worksheetPath)
    MsgDF(data.frame(`Sheet names in Workbook` = shNames))
    sheetHead <- GetSheetHeading(worksheetPath)
    mainSheet <- parse_workbook(worksheetPath, runID)
    mainSheet <- sanitizeSheet(mainSheet)
    outFile <- WriteMainSheet(mainSheet, sheetHead)
    outVals <- c(runID = mainSheet[1, "Run_Number"], outFile = outFile)
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


CheckOtherFiles <- function(worksheetPath, runID) {
    notFoundMsg <- crayon::bgRed("The PACT run worksheet was not found:")
    message(notFoundMsg,"\n", worksheetPath, "\n")
    parentFolder <- dirname(worksheetPath)
    potentialFi <- list.files(path = parentFolder, full.names = TRUE)
    if (length(potentialFi) > 1) {
        message(crayon::bgRed("Checking other existing files:"), "\n")
        message(paste(potentialFi, sep = "\n"))
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
    outVals <- parseExcelFile(pfile, runID)
    return(outVals)
}


# FUN: Returns run type based on the input pact_id or worksheet name ----------
determineRunType <- function(input) {

    runType <- dplyr::case_when(
        stringr::str_detect(input, "2000|NextSeq550") ~ "NextSeq2000",
        stringr::str_detect(input, "TEST")            ~ "test",
        stringr::str_detect(input, "ILC")             ~ "Illumina",
        stringr::str_detect(input, "TMB")             ~ "TMB",
        TRUE                                          ~ "regular"
    )

    return(runType)
}


# FUN: Returns Path to xlsx file ----------------------------------------------
getExcelPath <- function(input, runType) {
    if (stringr::str_detect(input, .Platform$file.sep)) {
        return(input)
    }
    drive = file.path("", "Volumes", "molecular", "MOLECULAR LAB ONLY")
    folder = file.path("NYU PACT Patient Data", "Workbook")
    runyr <- stringr::str_split_fixed(input, "-", 3)[, 1]
    yearDir <- paste0("20", runyr)
    xlFi <- paste0(input, ".xlsm")
    worksheetPath <- file.path(drive, folder, yearDir, input, xlFi)

    if (runType == "NextSeq2000") {
        folder <- "Validations/PACT new i7-NextSeq2000/Wet Lab/Workbook"
        worksheetPath <- file.path(drive, folder, yearDir, xlFi)
    }
    if (runType == "test") {
        drive = "/Volumes/molecular/Molecular/Validation/PACT/Test_Sheets"
        worksheetPath <- file.path(drive, xlFi)
    }
    if (runType == "Illumina") {
        drive = "/Volumes/molecular/Molecular/Validation/PACT/Test_Sheets"
        worksheetPath <- file.path(drive, xlFi)
    }
    if (runType == "TMB") {
        drive = "/Volumes/molecular/MOLECULAR LAB ONLY/Validations/TMB"
        worksheetPath <- file.path(drive, xlFi)
    }
    return(worksheetPath)
}


# Gets dataframe and saves as CSV file ----------------------------------------
writeSampleSheet <- function(input, token = NULL, runID = NULL) {
    outVals <- NULL
    runType <- determineRunType(input)
    worksheetPath <- getExcelPath(input, runType)
    if (file.exists(worksheetPath)) {
        outVals <- parseExcelFile(worksheetPath, runID)
    } else {
        outVals <- CheckOtherFiles(worksheetPath, runID)
    }
    if (!is.null(outVals) & !is.null(token)) {
        pushToRedcap(outVals, token)
    }
}

# MAIN Execution start ---------------------------------------- ---------------

loadPacks()
checkMounts()

writeSampleSheet(input, token, runID)
