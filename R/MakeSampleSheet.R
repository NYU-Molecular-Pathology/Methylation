#!/usr/bin/env Rscript
## Script name: MakeSampleSheet.R
## Purpose: Functions to copy Methylation xlsm worksheet and write to csv
## Date Created: August 13, 2022
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025

gb <- globalenv(); assign("gb", gb)
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
cpInLnk2 = "https://github.com/NYU-Molecular-Pathology/Methylation/main/R/MakeSampleSheet.R"

# Global Paths ----------------------------------------------------------------
rschSheets = "/Volumes/snudem01labspace/Methylation_Worksheets"
worksheetDirPath = "/Volumes/molecular/MOLECULAR LAB ONLY/NYU-METHYLATION/WORKSHEETS/"
researchWsPath = "/Volumes/snudem01labspace/Methylation_Worksheets/"
validation_dir = "/Volumes/molecular/MOLECULAR LAB ONLY/Validations/MethylationEPIC V2.0 Kit/WORK SHEET"
epicv1_val_dir = "/Volumes/CBioinformatics/Validations/Methylation/EPIC_V1_RUNS"
clinDrv = "/Volumes/molecular/MOLECULAR LAB ONLY/NYU-METHYLATION"
methDir = "/Volumes/CBioinformatics/Methylation/Clinical_Runs"


# Helper Functions ------------------------------------------------------------

# Prints messages in blue background
bkBlue <- function(...) {
    return(crayon::bgBlue$white$bold(paste0(...)))
}

# Prints messages in red background
bkRed <- function(...) {
    return(crayon::bgRed$white$bold(paste(..., collapse = " ")))
}

# Prints message of the current function executed in yellow background
msgFunName <- function(pthLnk, funNam) {
    message("Executing function: ", crayon::black$bgYellow(funNam),
            " from RScript in:\n", pthLnk, "\n")
}

# Messages what parameter(s) or variables were passed to a function
msgParams <- function(...) {
    vars <- data.frame(...)
    message("Param passed:\n",
            paste0(crayon::bgGreen(names(vars)), "=", vars, " "),
            "\n")
}

# Returns the captured ouput of arguments evaluated as a character string on newlines
msg_df <- function(obj_dat) {
    if (class(obj_dat) %in% c("character", "data.frame")) {
        message(paste0(capture.output(obj_dat), collapse = "\n"))
    } else {
        message(paste0(capture.output(as.data.frame(obj_dat)), collapse = "\n"))
    }
}


# FUNC: Returns run year from the methylation runID for the workbook year path
grabYear <- function(runID) {
    msgFunName(cpInLnk2, "grabYear")
    msgParams(runID)

    runYear <- NULL
    runID <- gsub(" ", "", runID)
    yr <- stringr::str_split_fixed(runID, "-", 2)[, 1]
    runYear <- yr
    if (nchar(yr) > 2) {
        runYear <- substring(yr, 3)
    }
    if (nchar(yr) > 0) {
        runYear <- paste0("20", runYear)
    }
    return(runYear)
}


# FUNC: Returns the newest xlsm file in the lab worksheet directory
GetNewestFile <- function(prevMC) {
    newestFile <- which.max(file.info(prevMC)$mtime)
    prevMC <- sub(".xlsm", "", basename(prevMC))
    newestRun <- paste0(prevMC[newestFile])
    message(crayon::bgCyan("List of Runs Found:\n"))
    msg_df(prevMC)
    return(prevMC)
}


# Returns a text string of the latest modified run name
# If non-clinical run, then a research directory is returned
listMolecularSheets <- function(runID = gb$runID, getAll = F) {
    msgFunName(cpInLnk2, "listMolecularSheets")
    msgParams(getAll, gb$runID)

    isMC <- sjmisc::str_contains(runID, "MGDM") | sjmisc::str_contains(runID, "MC")
    is_validation <- sjmisc::str_contains(runID, "VAL")
    runYear <- grabYear(runID)

    if (isMC) {
        wsPath <- file.path(gb$clinDrv, "WORKSHEETS", runYear)
        prevMC <- dir(path = wsPath, pattern = "MGDM", full.names = T)
    } else {
        wsPath <- file.path(gb$rschSheets, runYear)
        prevMC <- dir(path = wsPath, pattern = "MR", full.names = T)
    }
    if (is_validation) {
        wsPath <- file.path(gb$validation_dir, runYear)
        prevMC <- dir(path = wsPath, pattern = "MGDM", full.names = T)
    }
    if (sjmisc::str_contains(runID, "EPICV1")) {
        message(gb$epicv1_val_dir)
        wsPath <- gb$epicv1_val_dir
        prevMC <- dir(path = wsPath, pattern = "MGDM", full.names = T)
    }
    if (getAll) {
        return(prevMC)
    } else{
        return(GetNewestFile(prevMC))
    }
}


# Verifies if input runID is valid by comparing to available worksheet names
checkValidRun <- function(runID) {
    msgFunName(cpInLnk2, "checkValidRun")

    is_validation <- sjmisc::str_contains(runID, "VAL")
    ws.list <- listMolecularSheets(runID, getAll = T)
    found <- paste0(runID, ".xlsm") %in% basename(ws.list)

    message("Is the runID valid? ", found)
    if (!found) {
        message(bkRed("runID", gb$runID, "is not valid"))
        message(bkBlue(gb$runID, ".xlsm", " not found in worksheets folders:"))
        message(worksheetDirPath, "\n", researchWsPath)
        message(bkBlue("Ensure file exists in the worksheet directory!"))
        stop(
            paste("Excel Worksheet could not be found in either directory:\n"),
            worksheetDirPath,
            "\n",
            researchWsPath
        )
    }
    return(found)
}


# Copies the excel worksheet to the current directory using fs copy
FsCopyFile <- function(fileLoc) {
    if (!file.exists(fileLoc)) {
        return(message("File does not exist:\n", fileLoc))
    }

    message("\n\nCopying worksheet from Worksheets Folder:")
    if (!file.exists(file.path(getwd(), basename(fileLoc)))) {
        fs::file_copy(fileLoc, getwd(), overwrite = F)
    }
}


# FUNC: Copies lab xlsm worksheet from network share to the working directory
copyWorksheetFile <- function(runID = NULL, runYear = NULL) {
    msgFunName(cpInLnk2, "copyWorksheetFile")

    if (is.null(runID)) {
        runID = paste0(basename(getwd()))
    }

    runYear <- grabYear(runID)
    msgParams(runID, runYear)
    isMC = sjmisc::str_contains(runID, "MGDM") | sjmisc::str_contains(runID, "MC")
    assign("isMC", isMC)
    message("\nIs methylation run Clinical? ", isMC)
    is_validation <- sjmisc::str_contains(runID, "VAL")
    mountLoc <- ifelse(isMC, file.path(gb$clinDrv, "WORKSHEETS"), gb$rschSheets)
    if (is_validation) {
        mountLoc <- gb$validation_dir
    }
    if (dir.exists(mountLoc)) {
        fiPath <- file.path(mountLoc, runYear, paste0(runID, ".xlsm"))
        if (sjmisc::str_contains(runID, "EPICV1")) {
            fiPath <- file.path(gb$epicv1_val_dir, paste0(runID, ".xlsm"))
        }
        message("\nCopying file from:\n", fiPath)
        FsCopyFile(fiPath)
    } else {
        message("\n", bkRed("Folder not found:"), mountLoc)
        message("\nEnsure the path to this network Volume is correct\n")
    }
}


#FUNC: Returns the full path to the xlsm workbook in the working directory
GrabSampleSheet <- function() {
    msgFunName(cpInLnk2, "GrabSampleSheet")
    samSh <- dir(path = getwd(), full.names = T, ".xlsm")
    if (!(length(samSh) > 0)) {
        return(NULL)
    }
    if (length(samSh) > 1) {
        message("Multiple samplesheets found:\n")
        msg_df(samSh)
        samSh <- samSh[stringr::str_detect(samSh, pattern = "\\$", negate = T)]
        samSh <- samSh[!stringi::stri_detect_fixed(samSh, "~$")]
    }
    message("Using following samplesheet:\n", samSh[1])
    stopifnot(length(samSh) > 0)
    return(samSh[1])
}


# FUNC: Returns message of the total number of samples in workbook or if error
MsgSampleCount <- function(worksheet, thisSh) {
    msgFunName(cpInLnk2, "MsgSampleCount")
    templateDir = "Clinical_Methylation/methylation_run_TEMPLATE.xlsm"
    if (length(thisSh) == 0) {
        warning("No .xlsm sheet found")
    }
    if (length(worksheet) == 0) {
        warning("Samplesheet ", thisSh[1],
                " has an invalid format or is missing an integer in Cell B4")
        message("You may have to manually edit samplesheet")
        message("Try copying data into a template file:\n", templateDir)
        stopifnot(length(worksheet) > 0)
    } else {
        message("Total sample count found is: ", worksheet[1])
    }
}


# FUNC: Returns the total number of samples in the run by reading
# Cell B4 in Sheet1 of the XLSM lab worksheet
getTotalSamples <- function(thisSh = NULL) {
    msgFunName(cpInLnk2, "getTotalSamples")
    thisSh <-  ifelse(is.null(thisSh), GrabSampleSheet(), thisSh)
    thisSh <- thisSh[!stringi::stri_detect_fixed(thisSh, "~$")]
    worksheet <- suppressMessages(
        readxl::read_excel(thisSh[1], col_names = "Total",  range = "B4:B4")
    )
    MsgSampleCount(worksheet, thisSh)
    return(as.integer(paste0(worksheet[1])))
}


# FUNC: Reads the workbook cell F4 and returns the run date
ReadSheetDate <- function(sampleSheet) {
    msgFunName(cpInLnk2, "ReadSheetDate")
    msgParams(sampleSheet)
    wsDate <- suppressMessages(as.data.frame(
        readxl::read_excel(
            sampleSheet,
            sheet = 1,
            col_names = F,
            range = "F4:F4",
            trim_ws = T
        ))[1])
    names(wsDate) = "Date"
    message("DATE: ", wsDate$Date)
    return(wsDate)
}

#' Reads the raw_input tab of an xlsm workbook as a samplesheet for minfi
#' Reads the first worksheet containing "raw" in its name, processes it as a samplesheet, 
#' and returns the run ID, total sample count, or the worksheet.
#'
#' @param run_ID Logical. If TRUE, return the run ID.
#' @param totalSam Logical. If TRUE, return the total number of samples.
#' @param wks Logical. If TRUE, return the raw_input sheet as data frame.
#'
#' @return Character (run ID), integer (total samples), or data.frame (worksheet), depending on input
readSampleSheet <- function(run_ID = FALSE, totalSam = FALSE, wks = FALSE) {

    msgFunName(cpInLnk2, "readSampleSheet")
    msgParams(run_ID, totalSam, wks)

    sam_col <- "Sample_Name"
    xlsm_wb <- GrabSampleSheet()
    message(paste("Reading worksheet named:", xlsm_wb))

    # Returns the index(es) of sheet names containing "raw"
    raw_index <- grep("raw", readxl::excel_sheets(xlsm_wb), ignore.case = TRUE)[1]

    worksheet <- suppressMessages(readxl::read_excel(xlsm_wb,
        sheet = raw_index, col_names = TRUE, col_types = "text", trim_ws = TRUE)
    )

    to_keep <- worksheet[, sam_col] != 0 & !is.na( worksheet[, sam_col])
    worksheet <- worksheet[to_keep, ]

    wsDate <- ReadSheetDate(xlsm_wb)
    worksheet$Date <- paste0(wsDate$Date[1])

    if (is.null(gb$runID)) {
        gb$runID <- paste0(stringr::str_split_fixed(basename(xlsm_wb), pattern = ".xlsm", 2)[1,1])
    }

    if (gb$runID != worksheet$Project[1]) {
        message("The Batch ID in the samplesheet: ", worksheet$Project[1])
        message("Does not Match the input Run ID: ", gb$runID)
        stopifnot(gb$runID == worksheet$Project[1])
    }

    if (run_ID == TRUE) return(worksheet$Project[1])
    if (totalSam == TRUE) return(getTotalSamples(sampleSheet))
    if (wks == TRUE) return(worksheet)

    return(worksheet)
}


GenerateSheetWarning <- function(columnName, df) {
    defaultVal <- if (columnName == "Tech") "NA" else if (columnName == "MP_number") "none" else "blank"
    warnMsg <- sprintf("No %s column found in samplesheet: Assigning as %s", columnName, defaultVal)
    if (is.null(df[[columnName]])) {
        message(crayon::bgRed(warnMsg))
        df[[columnName]] <- defaultVal
    }
    return(df)
}


checkForIssues <- function(condition, warningMessage, colsToPrint) {
    if (any(condition)) {
        warning(warningMessage)
        message(bkRed(warningMessage))
        msg_df(colsToPrint)
        stop(warningMessage)
    }
}


# FUNC: Validates the raw_input tab of the xlsm worksheet and header names
checkSampleSheet <- function(df) {
    msgFunName(cpInLnk2, "checkSampleSheet")
    message("Params: df = \n", paste0(capture.output(df), collapse = "\n"))

    df <- GenerateSheetWarning("Tech", df)
    df <- GenerateSheetWarning("MP_number", df)
    df <- GenerateSheetWarning("b_number", df)

    missingName <- df$Sample_Name == "0" | is.na(df$Sample_Name)
    if (any(missingName)) {
        df <- df[!missingName,]
    }
    missMsg <-
        "Some samples are missing RD-Numbers or are 0! Check samplesheet.csv"

    checkForIssues(missingName, missMsg, df[, c(1, 3, 8:11)])

    dupes <- duplicated(df$Sample_Name)
    dupeMsg <-
        "Duplicated sample name found: check df$Sample_Name in samplesheet.csv"
    checkForIssues(dupes, dupeMsg, df[, c(1, 3, 8:11)])
    if (df$Sample_Name[1] == "Control") {
        stop("Error: The first item should be 'control', not 'Control'.")
    }

    missingControls <- stringr::str_count(df$Sample_Name, "control|RD-") == 0
    rd_msg <-
        "Some samples are missing RD-numbers or a control! Check Sample_Name in samplesheet.csv"
    checkForIssues(missingControls, rd_msg, df[, c(1, 3, 8:11)])

    if (nrow(df) < 8) {
        message(bkRed("Warning: samplesheet.csv contains less than 8 samples!"))
        message("Check .xslm sheet tab 'raw_labels' for any errors")
        msg_df(df[, 1:11])
    }

    return(df)
}


# FUNC: Writes a .csv samplesheet for redcap & minfi with defined input headers
writeSampleSheet <- function(df, bn = NULL, sampleName, dnaNumber, Sentrix) {
    msgFunName(cpInLnk2, "writeSampleSheet")
    if (is.null(bn)) {bn = file.path(gb$methDir, df$Batch, df$Sentrix)}
    toKeep <- df[, sampleName] != 0 & !is.na(df[, sampleName])
    df <- df[toKeep,]
    splitSentrix = as.data.frame(stringr::str_split_fixed(df[, "Sentrix_ID"], "_", 2))
    samplesheet_csv = data.frame(
        Sample_Name =  df[, sampleName],
        DNA_Number = df[, dnaNumber],
        Sentrix_ID = splitSentrix[, 1],
        Sentrix_Position = splitSentrix[, 2],
        SentrixID_Pos = df[, Sentrix],
        Basename = paste0(bn),
        RunID = df$Batch,
        MP_num = df$MP_number,
        tech = df$Tech,
        tech2 = df$Tech2,
        Date = df$Date
    )
    runID <- paste0(df$Batch)[1]
    is_validation <- sjmisc::str_contains(runID, "VAL")

    if (is_validation) {
        samplesheet_csv$Sample_Name <-
            paste(samplesheet_csv$Sample_Name, "VAL", sep = "_")
    }
    write.csv(
        samplesheet_csv,
        file = "samplesheet.csv",
        quote = F,
        row.names = F
    )
}


# FUNC: Checks header names of .xlsm sheet for mismatch
checkHeaders <- function(worksheet) {
    msgFunName(cpInLnk2, "checkHeaders")
    hdrs <- dimnames(worksheet)[[2]]
    msg_df(data.frame(HEADERS = hdrs))
    Var.names = c(
        sampleName = paste(hdrs[9]), # "Sample_Name"
        Sentrix = paste(hdrs[1]),    # "Sentrix_ID"
        dnaNumber = paste(hdrs[10]), # "b_number"
        mpnum = paste(hdrs[8])       # "MP_number"
    )
    msg_df(data.frame(VarNames = Var.names))
    var.default = c("Sample_Name", "Sentrix_ID", "b_number", "MP_number")
    if (!(all(var.default %in% hdrs))) {
        missing_head <- paste(var.default[!var.default %in% hdrs], collapse = "\n")
        err_msg <-
            "\nThe following header(s) in the 'raw_input' tab are missing or in the wrong column:\n"
        stop(err_msg, missing_head)
    }
    return(Var.names)
}

# FUNC: returns the latest RD in the input list of numbers
get_newest_rd <- function(rd_num){
    parsed <- strcapture(
        "^RD-(\\d+)-(\\d+)$",
        rd_num,
        proto = list(year = integer(), num = integer())
    )
    rd_max <- rd_num[which.max(parsed$year * 1e6 + parsed$num)]
    return(rd_max)
}


get_next_rd <- function(rd_nums) {
    suffix_num <- suppressWarnings(as.integer(sub(".*-", "", rd_nums)))
    
    if (all(is.na(suffix_num))) {
        return(NA_character_)
    } else {
        last_rd <- rd_nums[which.max(suffix_num)]
        last_num <- as.integer(sub(".*-", "", last_rd))
        prefix <- sub("-[^-]+$", "", last_rd)
        
        return(paste0(prefix, "-", last_num + 1))
    }
}


is_present <- function(x) {
    !is.null(x) && !is.na(x) && as.character(x) != "" && as.character(x) != "0"
}

# FUNC: Checks RC DB for any record_id matching FILLER annotated rows in samplesheet using b & mp numbers
check_rd_fillers <- function(df) {
    has_filler <- stringr::str_detect(df$Sample_Name, stringr::regex("fill", ignore_case = TRUE))
    if (any(has_filler)) {
        apiUrl = "https://redcap.nyumc.org/apps/redcap/api/"
        rcon <- redcapAPI::redcapConnection(apiUrl, gb$token)
        flds <- c("record_id", "b_number", "tm_number", "accession_number", "block",
                  "barcode_and_row_column", "tissue_comments", "nyu_mrn")
        params = list(rcon, fields = flds, survey = F, dag = F, factors = F, form_complete_auto = F)
        message("Checking REDCap Database for fillers...")
        dbCols <- do.call(redcapAPI::exportRecordsTyped, c(params))
        db <- as.data.frame(dbCols)
        db <- db[grepl(pattern = "RD-", db$record_id), ]
        db <- db[!grepl(pattern = "_R0", db$barcode_and_row_column), ]
        db <- db[!grepl(pattern = "NO IDAT FILE", db$barcode_and_row_column), ]
        db <- db[!grepl(pattern = "DUPLICATE", db$barcode_and_row_column), ]
        rownames(db) <- NULL
        
        for (filler in which(has_filler)) {
            curr_row <- df[filler,]
            curr_bnum <- curr_row$b_number
            rd_num <- ""
            if (is_present(curr_bnum)) {
                match_b <- grepl(curr_bnum, db$b_number) | grepl(curr_bnum, db$block)
                match_b[is.na(match_b)] <- FALSE
                if (any(match_b)) {
                    rd_num <- db$record_id[match_b]
                    if (length(rd_num) > 1) {
                        rd_num <- get_newest_rd(rd_num)
                    }
                    message("Filler B-number found: ", curr_bnum, " ", rd_num)
                    df[filler, "Sample_Name"] <- rd_num
                    next
                }
                curr_mp <- curr_row$MP_number
                if (is_present(curr_mp)) {
                    match_mp <- apply(db[c("accession_number", "tm_number", "b_number")], 1,
                                      function(x) {any(grepl(curr_mp, x, fixed = TRUE), na.rm = TRUE)}
                    )
                    if (any(match_mp)) {
                        matched_rows <- db[match_mp, , drop = FALSE]
                        if (nrow(matched_rows) > 1) {
                            rd_num <- get_newest_rd(matched_rows$record_id)
                        } else {
                            rd_num <- matched_rows$record_id
                        }
                        message("Filler MP match found: ", curr_mp, " ", rd_num)
                        df[filler, "Sample_Name"] <- rd_num
                        next
                    }
                }
                if (rd_num == "") {
                    message("No matching RD number found for filler:\n",
                            paste(curr_row[1, c("b_number", "MP_number")], collapse = " "))
                    message("Creating new RD-number...")
                    mgdms <- grepl(pattern = paste0("RD-", gb$runID, "-"), df$Sample_Name)
                    all_rd <- df$Sample_Name[mgdms]
                    rd_num <- get_next_rd(all_rd)
                    df[filler, "Sample_Name"] <- rd_num
                }
            }
        }
    }
    return(df)
}


# FUNC: Validates RD-numbers and workbook formatting
FormatSampleData <- function(worksheet, runID, sampleNumb) {
    msgFunName(cpInLnk2, "FormatSampleData")
    msgParams(runID, sampleNumb)

    hdrs <- checkHeaders(worksheet)
    df = as.data.frame(worksheet)[1:sampleNumb,]
    rd_missing <- is.na(df[, 1]) | df[, 1] == 0

    if (any(rd_missing)) {
        df <- df[!rd_missing, ]
    }

    df$Notes <- paste(df$Notes[1])
    df <- gb$checkSampleSheet(df)
    df <- check_rd_fillers(df)
    return(df)
}


# FUNC: Reads .xlsm lab worksheet and outputs a CSV file formatted for minfi
readSheetWrite <- function(sampleNumb = NULL, runID = NULL) {
    msgFunName(cpInLnk2, "readSheetWrite")

    if (file.exists("samplesheet.csv")) {
        message(crayon::white$bgGreen("samplesheet.csv already exists!"))
        return(message("Delete existing file to create a new samplesheet.csv"))
    }

    if (is.null(sampleNumb)) {
        sampleNumb <- getTotalSamples()
    }
    if (is.null(runID)) {
        runID <- paste0(basename(getwd()))
    }
    if (!is.integer(sampleNumb) | sampleNumb < 8) {
        stop("Check samplesheet .xlsm cell B4 for valid integer total arrays")
    }

    msgParams(sampleNumb, runID)

    raw_input_df <- readSampleSheet(wks = TRUE)
    run_csv_df <- FormatSampleData(raw_input_df, runID, sampleNumb)

    writeSampleSheet(
        run_csv_df,
        bn = file.path(gb$methDir, runID, run_csv_df[, "Sentrix_ID"]),
        sampleName = "Sample_Name",
        dnaNumber = "b_number",
        Sentrix = "Sentrix_ID"
    )
}
