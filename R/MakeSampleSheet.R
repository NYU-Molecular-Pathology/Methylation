#!/usr/bin/env Rscript
## Script name: MakeSampleSheet.R
## Purpose: source of global scripts to copy Methylation worksheet and write to csv samplesheet
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023

gb <- globalenv(); assign("gb", gb)
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
cpInLnk2 = "https://github.com/NYU-Molecular-Pathology/Methylation/main/R/MakeSampleSheet.R"
rschSheets = "/Volumes/snudem01labspace/Methylation_Worksheets"
worksheetDirPath = "/Volumes/molecular/MOLECULAR LAB ONLY/NYU-METHYLATION/WORKSHEETS/"
researchWsPath = "/Volumes/snudem01labspace/Methylation_Worksheets/"
validation_dir = "/Volumes/molecular/MOLECULAR LAB ONLY/Validations/MethylationEPIC V2.0 Kit/WORK SHEET"
epicv1_val_dir = "/Volumes/CBioinformatics/Validations/Methylation/EPIC_V1_RUNS"

bkBlue <- function(...) {
    return(crayon::bgBlue$white$bold(paste0(...)))
}

bkRed <- function(...) {
    return(crayon::bgRed$white$bold(paste(..., collapse = " ")))
}

msgFunName <- function(pthLnk, funNam) {
    message("Executing function: ", crayon::black$bgYellow(funNam),
            " from RScript in:\n", pthLnk, "\n")
}

msgParams <- function(...) {
    vars <- data.frame(...)
    message("Param passed:\n", paste0(crayon::bgGreen(names(vars)), "=", vars, " "), "\n")
}


grabYear <- function(runID) {
    msgFunName(cpInLnk2, "grabYear")
    msgParams(runID)
    rnum <- NULL
    yr <- stringr::str_split_fixed(runID, "-", 2)[, 1]

    if (nchar(yr) > 2) {
        rnum <- substring(yr, 3)
    } else{
        rnum <- yr
    }
    if (nchar(yr) > 0) {
        rnum <- paste0("20", rnum)
    }
    return(rnum)
}


GetNewestFile <- function(prevMC) {
    newestFile <- which.max(file.info(prevMC)$mtime)
    prevMC <- sub(".xlsm", "", basename(prevMC))
    newestRun <- paste0(prevMC[newestFile])
    message(crayon::bgCyan("List of Runs Found:\n"))
    message(paste0(capture.output(as.data.frame(prevMC)), collapse = "\n"))
    return(prevMC)
}


# Returns a text string of the latest modified Run name, if isMC=False then research directory is returned
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

# Verifies if the runID inputted is valid
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
        message(bkBlue("Check worksheet directory if the file is exists as .xlsm and not .xlsx"))
        stop("Excel Worksheet could not be found in the input directory")
    }
    return(found)
}


FsCopyFile <- function(fileLoc) {
    if (file.exists(fileLoc)) {
        message("\n\nCopying worksheet from Worksheets Folder:")
        if (!file.exists(file.path(getwd(), basename(fileLoc)))) {
            fs::file_copy(fileLoc, getwd(), overwrite = F)
        }
    } else {
        message("\nFile not found:\n", fileLoc)
        message("\nMake sure your path is correct, try print(gb$copyWorksheetFile)\n")
    }
}


# FUN: copies the molecular or research lab Worksheet xlsm to cwd
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
        if (sjmisc::str_contains(runID, "EPICV1")){
            fiPath <- file.path(gb$epicv1_val_dir, paste0(runID, ".xlsm"))
        }
        message("\nCopying file from:\n", fiPath)
        FsCopyFile(fiPath)
    } else {
        message("\n", bkRed("Folder not found:"), mountLoc)
        message("\nEnsure the path to this network Volume is correct\n")
    }
}

GrabSampleSheet <- function() {
    msgFunName(cpInLnk2, "GrabSampleSheet")
    samSh <- dir(path=getwd(), full.names=T, ".xlsm")
    if (!(length(samSh) > 0)) {
        return(NULL)
    }
    if (length(samSh) > 1) {
        message("Multiple samplesheets found:\n")
        message(paste0(capture.output(samSh), collapse="\n"))
        samSh <- samSh[stringr::str_detect(samSh,pattern = "\\$",negate = T)]
        samSh <- samSh[!stringi::stri_detect_fixed(samSh, "~$")]
    }
    message("Using following samplesheet:\n", samSh[1])
    stopifnot(length(samSh) > 0)
    return(samSh[1])
}


MsgSampleCount <- function(worksheet, thisSh) {
    msgFunName(cpInLnk2, "MsgSampleCount")
    templateDir = "Clinical_Methylation/methylation_run_TEMPLATE.xlsm"
    if (length(thisSh)==0) {warning("No .xlsm sheet found")}
    if (length(worksheet) == 0) {
        warning("Samplesheet ", thisSh[1], " has an invalid format or is missing an integer in Cell B4")
        message("You may have to manually edit samplesheet")
        message("Try copying data into a template file:\n", templateDir)
        stopifnot(length(worksheet) > 0)
    } else {
        message("Total sample count found is: ", worksheet[1])
    }
}

# Returns Total Sample Count in the run
getTotalSamples <- function(thisSh=NULL) {
    msgFunName(cpInLnk2, "getTotalSamples")
    thisSh <-  ifelse(is.null(thisSh), GrabSampleSheet(), thisSh)
    thisSh <- thisSh[!stringi::stri_detect_fixed(thisSh, "~$")]
    worksheet <- suppressMessages(
        readxl::read_excel(thisSh[1], col_names="Total", range="B4:B4")
    )
    MsgSampleCount(worksheet, thisSh)
    return(as.integer(paste0(worksheet[1])))
}

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


# FUN: Reads the csv samplesheet for minfi input
readSampleSheet <- function(run_ID=F, totalSam=F, wks=F) {
    msgFunName(cpInLnk2, "readSampleSheet")
    msgParams(run_ID, totalSam, wks)

    sampleSheet <- GrabSampleSheet()
    message(paste0("Reading worksheet named: ", sampleSheet))

    worksheet <- suppressMessages(readxl::read_excel(
        sampleSheet,
        sheet = 2,
        col_names = T,
        col_types = "text",
        trim_ws = T
    ))
    to_keep <- worksheet$Sample_Name != 0
    worksheet <- worksheet[to_keep, ]

    wsDate <- ReadSheetDate(sampleSheet)
    worksheet$Date <- paste0(wsDate$Date[1])

    if (is.null(gb$runID)){
        gb$runID <- paste0(stringr::str_split_fixed(basename(sampleSheet), pattern = ".xlsm", 2)[1,1])
    }

    if (gb$runID != worksheet$Project[1]){
        message("The Batch ID in the samplesheet: ", worksheet$Project[1])
        message("Does not Match the input Run ID: ", gb$runID)
        stopifnot(gb$runID == worksheet$Project[1])
    }

    if (run_ID == T) {
        return(worksheet$Project[1])
    }
    if (totalSam == T) {
        return(getTotalSamples())
    }
    if (wks == T) {
        return(worksheet)
    }
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
        message(paste0(capture.output(colsToPrint), collapse = "\n"))
        stop(warningMessage)
    }
}


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
    missMsg <- "Some samples are missing RD-Numbers or are 0! Check samplesheet.csv"
    checkForIssues(missingName, missMsg, df[, c(1, 3, 8:11)])
    dupes <- duplicated(df$Sample_Name)
    dupeMsg <- "Duplicated sample name found: check df$Sample_Name in samplesheet.csv"
    checkForIssues(dupes, dupeMsg, df[, c(1, 3, 8:11)])
    if (df$Sample_Name[1] == "Control") {
        stop("Error: The first item should be 'control', not 'Control'.")
    }
    missingControls <- stringr::str_count(df$Sample_Name, "control|RD-") == 0
    rd_msg <- "Some samples are missing RD-numbers or control is missing! Check df$Sample_Name in samplesheet.csv"
    checkForIssues(missingControls, rd_msg, df[, c(1, 3, 8:11)])

    if (nrow(df) < 8) {
        message(bkRed("Warning: samplesheet.csv contains less than 8 samples!"))
        message("Check .xslm sheet tab 'raw_labels' for any errors")
        message(paste0(capture.output(df[, 1:11]), collapse = "\n"))
    }

    return(df)
}

# FUN: translates the xlsm excel file to the .csv samplesheet for redcap and minfi
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
        samplesheet_csv$Sample_Name <- paste(samplesheet_csv$Sample_Name, "VAL", sep = "_")
    }
    write.csv(samplesheet_csv, file = "samplesheet.csv", quote = F, row.names = F)
}

# FUN: Checks header names of .xlsm sheet for mismatch
checkHeaders <- function(worksheet) {
    msgFunName(cpInLnk2, "checkHeaders")
    hdrs <- dimnames(worksheet)[[2]]
    message(paste0(capture.output(data.frame(HEADERS = hdrs)), collapse = "\n"))
    Var.names = c(
        sampleName = paste(hdrs[9]), # "Sample_Name"
        Sentrix = paste(hdrs[1]),    # "Sentrix_ID"
        dnaNumber = paste(hdrs[10]), # "b_number"
        mpnum = paste(hdrs[8])       # "MP_number"
    )
    message(paste0(capture.output(data.frame(VarNames = Var.names)), collapse="\n"))

    var.default = c("Sample_Name", "Sentrix_ID", "b_number", "MP_number")
    if (!(all(var.default %in% hdrs))) {
        missing_head <- paste(var.default[!var.default %in% hdrs], collapse = "\n")
        err_msg <-
            "\nThe following header(s) in the 'raw_input' tab are missing or in the wrong column:\n"
        stop(err_msg, missing_head)
    }
    return(Var.names)
}

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
    df <- checkSampleSheet(df)
    return(df)
}

# FUN: reads the .xlsm worksheet and outputs the .csv methyl experiment for MINFI
readSheetWrite <- function(sampleNumb = NULL, runID = NULL) {
    msgFunName(cpInLnk2, "readSheetWrite")
    if (!file.exists("samplesheet.csv")) {
        if (is.null(sampleNumb)) {sampleNumb <- getTotalSamples()}
        if (is.null(runID)) {runID <- paste0(basename(getwd()))}
        if (!is.integer(sampleNumb) | sampleNumb < 8) {
            stop("Check samplesheet .xlsm cell B4 for valid integer total arrays")
        }
        msgParams(sampleNumb, runID)

        worksheet = readSampleSheet(wks = T)
        df <- FormatSampleData(worksheet, runID, sampleNumb)

        writeSampleSheet(
            df, bn = file.path(gb$methDir, runID, df[, "Sentrix_ID"]),
            sampleName = "Sample_Name", dnaNumber = "b_number", Sentrix = "Sentrix_ID"
        )

    } else {
        message("\n", crayon::white$bgGreen("samplesheet.csv already exists!"), "\n")
        message("Create a new samplesheet.csv file by deleting the existing file")
    }

}
