#!/usr/bin/env Rscript
## Script name: CopyOutput.R
## Purpose: Functions to copy methylation output files to zdrive and REDCap
## Date Created: June 13, 2022
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025

options(stringsAsFactors = FALSE)
gb <- globalenv(); assign("gb", gb)
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
cpOutLnk = "https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/R/CopyOutput.R"
rschDrv = "/Volumes/snudem01labspace/FINAL_PDF_Reports_Brain"
metrics_dir <-
    "/Volumes/molecular/Molecular/MethylationClassifier/Methylation_QC_metrics"
UPLOAD_LOG_TSV = "redcap_failed_file_uploads.tsv"

# Helper function to print out a function's name and path on Github
msgFunName <- function(pthLnk, funNam) {
    message("Executing function: ", funNam, " from RScript in:\n", pthLnk, "\n")
}

# Helper functions to color text background in messages
mkBlue <- function(strVar) {return(crayon::white$bgBlue(strVar))}
mkGrn <- function(strVar) {return(crayon::white$bgGreen$bold(strVar))}
mkRed <- function(strVar) {return(crayon::white$bgRed$bold(strVar))}

CheckDirMake <- function(newFolder) {
    if (!dir.exists(newFolder)) {
        fs::dir_create(newFolder)
    }
}

MsgDF <- function(datObj) {
    if (class(datObj) != "data.frame") {
        datObj <- as.data.frame(datObj)
    }
    message(paste0(capture.output(datObj), collapse = "\n"))
}


# Helper archive function: creates a "previous" folder when reports already exists
SavePrevDir <- function(newFolder) {
    msgFunName(cpOutLnk, "SavePrevDir")
    message("\n", mkRed("Output folder already exists!"), "\n")
    message(mkRed('Copying existing reports to a new folder named "Previous"'))
    oldFi = dir(path = newFolder, full.names = T)
    prevs = file.path(newFolder, "previous")
    dir.create(prevs)
    message(prevs)
    fs::file_copy(path = oldFi, new_path = prevs)
}


Copy2TempDir <- function(fi2copy, runID) {
    msgFunName(cpOutLnk, "Copy2TempDir")
    tempDir <- file.path(fs::path_home(), runID)
    if (!dir.exists(tempDir)) {dir.create(tempDir)}
    fs::file_copy(fi2copy, tempDir, overwrite = T)
    fi2copy <- dir(tempDir, pattern = ".html", full.names = T)
    return(fi2copy)
}


CopyHtmlFiles <- function(newFolder, runID) {
    msgFunName(cpOutLnk, "CopyHtmlFiles")
    message(mkBlue("Copying Reports to output folder:"), "\n", newFolder)
    fi2copy <- dir(getwd(), pattern = ".html", full.names = T)
    if (length(fi2copy) >= 1) {
        message(mkGrn("Found the following html reports..."), "\n")
        MsgDF(fi2copy)
    }else{
        return(message(mkRed("No html reports to copy found!"), "\n"))
    }
    fi2copy <- Copy2TempDir(fi2copy, runID)
    #fs::file_copy(fi2copy, newFolder, overwrite = T)
}


# FUN: Copies Reports to Z drive
copy.to.clinical <- function(clinOut, runID, runYear) {
    msgFunName(cpOutLnk, "copy.to.clinical")
    newFolder <- file.path(clinOut, runYear, runID)
    CheckDirMake(newFolder)
    if (!dir.exists(newFolder)) {
        message("Output Folder could not be created, htmls not copied to:\n", newFolder)
    }
    if (length(dir(path = newFolder, pattern = ".html", full.names = T)) > 0) {
        SavePrevDir(newFolder)
    }else{
        CopyHtmlFiles(newFolder, runID)
    }
}


# Checks REDCap record to see last edit in log
redcap_last_edit <- function(rcon, record_id) {
    redcap_logs <- redcapAPI::exportLogging(rcon, logtype = "record", record = record_id)
    last_edit <- redcap_logs[1, "timestamp"]
    last_action <- redcap_logs[1, "action"]
    last_edit_date <- base::as.Date(last_edit, format = "%Y-%m-%d")
    is_today <- last_edit_date == Sys.Date()
    return(data.frame(is_today = is_today, last_action = redcap_logs[1, "action"]))
}


# Checks if field is already filled in REDCap returns boolean
checkRedcapRecord <- function(recordName, fieldName = 'classifier_pdf') {
    msgFunName(cpOutLnk, "checkRedcapRecord")
    is_val_control <- sjmisc::str_contains(recordName, "_control_")
    if (is_val_control) return(FALSE)

    formData <- list(
        "token" = gb$ApiToken,
        content = 'record',
        action = 'export',
        format = 'json',
        type = 'flat',
        csvDelimiter = '',
        'records[0]' = recordName,
        rawOrLabel = 'raw',
        rawOrLabelHeaders = 'raw',
        exportCheckboxLabel = 'false',
        exportSurveyFields = 'false',
        exportDataAccessGroups = 'false',
        returnFormat = 'json'
    )

    MsgDF(as.data.frame(formData))
    response <- httr::POST(url = gb$apiLink, body = formData, encode = "form")
    result <- httr::content(response)[[1]]
    return(result[[fieldName]])
}


MakeLogFile <- function(infoData, logFile) {
    message("Writing log:\n", infoData, "To file: ", logFile)
    write.table(infoData, file = logFile, append = T, quote = F,
                sep = '\t', row.names = F, col.names = F)
}


writeLogFi <- function(recordName, isHtml = T, logFile = "upload_log.tsv") {
    message("Check ", logFile, "\n", mkRed(recordName[1]), " already has data in REDCap")
    if (isHtml == T) {
        i = paste(recordName[1], "already has an html file in REDCap\n")
        MakeLogFile(i, logFile)
    } else{
        i = paste(recordName[1], "already has a classifier_value filled in REDCap\n")
        MakeLogFile(i, logFile)
    }
}


RedcapRcurl <- function(datarecord = NULL) {
    msgFunName(cpOutLnk, "RedcapRcurl")
    if (!is.null(datarecord)) {
        rcon <-
            redcapAPI::redcapConnection(apiLink, token = gb$ApiToken)
        redcapAPI::importRecords(rcon, data = datarecord,
                                 overwriteBehavior = "normal")
        message(crayon::bgBlue("Record Data Uploaded:"), "\n",
                paste(datarecord, collapse = "\n"))
    }
}


WarnSentrix <- function(record, isEmpty) {
    msgFunName(cpOutLnk, "WarnSentrix")

    if (record$barcode_and_row_column == isEmpty) {
        txtMsg <- " already has the same Sentrix ID on REDCap: "
        message(mkRed(record$record_id), txtMsg, isEmpty)
    } else{
        warning(mkRed("Different Sentrix ID is already on REDCap!"))
        message("Samplesheet Sentrix ID is: ", record$barcode_and_row_column)
        message("The Sentrix ID on REDCap is: ", isEmpty)
        logFi_out <- paste(gb$runID, "import_log.tsv", sep = "_")
        writeLogFi(record$record_id, isHtml = F, logFile = logFi_out)
    }
}


ValidateRedImport <- function(record) {
    is_validation <- sjmisc::str_contains(gb$runID, "VAL")
    is_val <- sjmisc::str_contains(record$record_id, "VAL")
    if (is_validation == T & is_val == F) {
        record$record_id <- paste0(record$record_id, "_VAL")
    }
    isEmpty <- checkRedcapRecord(record$record_id, "barcode_and_row_column")
    if (isEmpty == "") {
        #json_data <- jsonlite::toJSON(list(as.list(record)), auto_unbox = T)
        datarecord <- as.data.frame(record)
        RedcapRcurl(datarecord)
    } else{
        WarnSentrix(record, isEmpty)
    }
}


loopRedcapImport <- function(data) {
    msgFunName(cpOutLnk, "loopRedcapImport")
    if (!is.null(data)) {
        for (n in 1:nrow(data)) {
            record = c(data[n,])
            is_validation <- sjmisc::str_contains(gb$runID, "VAL")
            is_val <- sjmisc::str_contains(record$record_id, "VAL")
            if (is_validation == T & is_val == F) {
                record$record_id <- paste0(record$record_id, "_VAL")
            }
            ValidateRedImport(record)
        }
    } else{
        message("Nothing imported to REDCap. Dataframe is null!")
    }
}


GetRedcapCsv <- function(samsheet) {
    msgFunName(cpOutLnk, "GetRedcapCsv")
    if (is.null(samsheet)) {
        samsheet <- dir(path = getwd(), full.names = T, pattern = "_Redcap.csv", recursive = F)
    }
    if (length(samsheet) == 1) {
        data <- read.csv(samsheet, stringsAsFactors = F)
        if (any(duplicated(data$record_id))) {
            message(mkRed("Remove duplicate rows in the REDCap csv dataframe:"))
            MsgDF(data$record_id[duplicated(data$record_id)])
            data = data[!duplicated(data$record_id),]
        }
        return(data)
    } else {
        message("no _Redcap.csv file found or multiple files exist:/n")
        MsgDF(samsheet)
        stopifnot(length(samsheet) == 1)
    }
}


CheckImportData <- function(rawCsv) {
    msgFunName(cpOutLnk, "CheckImportData")
    message("Checking REDCap for existing data:")
    toImport <- unlist(lapply(rawCsv$record_id, FUN = function(rd) {
        return(checkRedcapRecord(rd, "classifier_value") == "")}))
    if (any(!toImport)) {
        message("Records have existing classifier_value and will not be over-written in REDCap:")
        toSkip <- rawCsv[!toImport, ]
        MsgDF(toSkip)
        log_fi_out <- paste(gb$runID, "import_log.tsv", sep = "_")
        invisible(lapply(1:nrow(toSkip), function(x) {
            writeLogFi(as.data.frame(toSkip[x, ]), isHtml = F, logFile = log_fi_out)
        }))
        rawCsv <- rawCsv[toImport, ]
    }
    return(rawCsv)
}


# REDCap: API call & Upload ---------------------------------------------------
# uploads the redcap classifier values must convert to JSON first
importDesktopCsv <- function(rcon, samsheet = NULL) {
    msgFunName(cpOutLnk, "importDesktopCsv")
    rawCsv <- GetRedcapCsv(samsheet)
    data <- rawCsv
    #data <- CheckImportData(rawCsv)
    if (nrow(data > 0)) {
        res <- redcapAPI::importRecords(rcon, data, "normal", "ids", logfile = "REDCapImportLog.txt")
        message("REDCap Response:\n", res)
    } else{
        message("No new data to import to REDCap")
    }
}


AddPngFilePath <- function(sh_Dat) {
    msgFunName(cpOutLnk, "AddPngFilePath")
    nfldr = file.path(stringr::str_split_fixed(gb$clinDrv, " ", 2)[1],
                      "MethylationClassifier")
    pathNam = file.path(nfldr, paste0(gb$runID, "_CNVs"), paste0(sh_Dat$record_id, "_cnv.png"))
    rms = paste(c("control_", "low_"), collapse = '|')
    pathNam <- stringr::str_replace(pathNam, rms, "")
    pathNam <- stringr::str_replace(pathNam, "//", "/")
    sh_Dat$cnv_file_path <- pathNam
    return(sh_Dat)
}


# Uploads any created cnv png files to redcap database
uploadCnPng <- function() {
    msgFunName(cpOutLnk, "uploadCnPng")
    rcon <- redcapAPI::redcapConnection(apiLink, gb$ApiToken)
    samSh <- gb$GrabSampleSheet()
    sampleNumb <- gb$getTotalSamples()
    sh_Dat <- suppressMessages(as.data.frame(readxl::read_excel(
        samSh, sheet = 3, range = "A1:N97", col_types = c("text")))[1:sampleNumb, 1:13])
    missing_rd <- sh_Dat$record_id == 0 | is.na(sh_Dat$record_id)
    if (any(missing_rd)) {
        sh_Dat <- sh_Dat[!missing_rd,]
    }
    sh_Dat <- AddPngFilePath(sh_Dat = sh_Dat)
    records <- sh_Dat$record_id
    for (idx in 1:length(records)) {
        pth = sh_Dat$cnv_file_path[idx]
        recordName = paste0(records[idx])
        message(mkBlue("Importing CNV Record:"), "\n", recordName, " ", pth)
        redcapAPI::importFiles(rcon, pth, recordName, field = "methyl_cn", overwrite = F, repeat_instance = 1)
    }
}


# Copy Output cnv Files if generated
copy.cnv.files <- function(newFolder, runID, runYear = NULL) {
    msgFunName(cpOutLnk, "copy.cnv.files")
    if (is.null(runYear)) {
        runYear = paste0(format(Sys.Date(), "%Y"))
    }
    cnv_folder <- file.path(newFolder, paste0(runID, "_CNVs/"))
    cnvNames <- dir(path = getwd(), full.names = T, "*_cnv.png")
    if (length(cnvNames) > 2) {
        message(paste0("Copying PNG to: ", cnv_folder))
        MsgDF(cnvNames)
        CheckDirMake(cnv_folder)
        if (dir.exists(cnv_folder)) {
            fs::file_copy(cnvNames, cnv_folder)
        }
    }
    uploadCnPng()
}

# Returns Total Sample Count in the run
getTotalSamples <- function(thisSh = NULL) {
    msgFunName(cpOutLnk, "getTotalSamples")

    templateDir = "Clinical_Methylation/methylation_run_TEMPLATE.xlsm"
    thisSh <-  ifelse(is.null(thisSh), gb$GrabSampleSheet(), thisSh)
    thisSh <- thisSh[!stringi::stri_detect_fixed(thisSh, "~$")]
    if (length(thisSh) == 0) {
        print("No .xlsm sheet, defaulting to 16 total samples")
        return(16)
    }
    worksheet <- suppressMessages(readxl::read_excel(thisSh[1], col_names = "Total", range = "B4:B4"))
    if (length(worksheet) == 0) {
        warning("Samplesheet ", thisSh[1], " is invalid format or no integer in Cell B4 found.")
        message("Manually edit samplesheet to fix-Try copying data into the template file:\n", templateDir)
        stopifnot(length(worksheet) > 0)
    } else {
        message("Total sample count found is: ", worksheet[1])
    }
    totNumb <- paste0(worksheet[1])
    return(as.integer(totNumb))
}


# Imports the xlsm sheet 3 data
importRedcapStart <- function(nfldr) {
    msgFunName(cpOutLnk, "importRedcapStart")

    samSh <- gb$GrabSampleSheet()
    sampleNumb <- getTotalSamples(samSh)
    sh_Dat <- suppressMessages(as.data.frame(readxl::read_excel(
        samSh, sheet = 3, range = "A1:N97", col_types = c("text")))[1:sampleNumb, 1:13])
    missing_rd <- sh_Dat$record_id == 0 | is.na(sh_Dat$record_id)
    if (any(missing_rd)) {
        sh_Dat <- sh_Dat[!missing_rd,]
    }
    sh_Dat <- AddPngFilePath(sh_Dat)
    runID <- paste0(sh_Dat$run_number[1])
    sh_Dat <- gb$NameControl(sh_Dat, runId = runID)
    sh_Dat <- CheckImportData(sh_Dat)
    if (nrow(sh_Dat) > 0) {
        loopRedcapImport(sh_Dat)
    } else{
        message(mkGrn("No new data to import from SampleSheet Data:"), "\n", samSh)
    }
}


DoRedcapApi <- function(rcon, recordName, runID) {
    msgFunName(cpOutLnk, "DoRedcapApi")
    message(mkBlue("Importing Record:"))
    data = data.frame(record_id = recordName, run_number = runID)
    logfi = paste0(recordName, "_redcapLog.txt")
    tryCatch(
        expr = {
            redcapAPI::importRecords(
                rcon, data, overwriteBehavior = "normal",
                returnContent = "ids", logfile = logfi)
        },
        error = function(e) {
            rdMsg <- paste(data$record_id, "failed import data to REDCap:")
            message(mkRed("DATA:"),"\n")
            MsgDF(data)
            message(mkRed(rdMsg), "\n", e$message)
        }
    )
}


CheckSarcRDnumber <- function(record) {
    msgFunName(cpOutLnk, "CheckSarcRDnumber")
    isSarc <- ifelse(stringr::str_detect(record, pattern = "sarc"),
                     yes = T,
                     no = F)
    if (isSarc == T) {
        record <- stringr::str_split_fixed(record, pattern = "_", 2)[1,1]
    }
    return(record)
}

callApiImport <- function(rcon, recordName, runID) {
    msgFunName(cpOutLnk, "callApiImport")
    is_validation <- sjmisc::str_contains(runID, "VAL")
    has_val <- sjmisc::str_contains(recordName, "VAL")
    if (is_validation == T & has_val == F) {
        recordName <- paste0(recordName, "_VAL")
    }
    isEmpty <- checkRedcapRecord(recordName, "subgroup")
    if (isEmpty == '') {
        DoRedcapApi(rcon, recordName, runID)
    } else{
        message(mkRed(recordName), " already has an assigned subgroup in REDCap: ", isEmpty)
        log_fi_out <- paste(runID, "import_log.tsv", sep = "_")
        writeLogFi(recordName, isHtml = F, logFile = log_fi_out)
    }
}


callApiFile <- function(rcon, recordName, ovwr = F, fiPath = NULL, fld = NULL) {
    msgFunName(cpOutLnk, "callApiFile")
    if (is.null(fiPath)) {
        recordFi <- paste0(recordName, ".html")
        message("\n", gb$mkBlue("Importing Record File:"), recordFi)
        fiPath <- file.path(getwd(), recordFi)
    }
    if (ovwr == F) {
        log_fi_out <- paste(gb$runID, "import_log.tsv", sep = "_")
        writeLogFi(recordName, logFile = log_fi_out)
    } else{
        if (is.null(fld)) {
            fld <- "classifier_pdf"
        }
        message("Uploading file:\n", fiPath)
        message("To REDCap Record: ", recordName)
        if (file.exists(fiPath)) {
            redcapAPI::importFiles(rcon,
                                   file = fiPath,
                                   record = recordName,
                                   field = fld)
        } else{
            message("REDCap file upload failed: ", fiPath)
            log_fi_out <- paste(gb$runID, "import_log.tsv", sep = "_")
            writeLogFi(recordName, logFile = log_fi_out)
        }
    }
}


# Creates QC record and uploads reports to redcap
uploadToRedcap <- function(file.list, deskCSV = T, runNumb = NULL) {
    msgFunName(cpOutLnk, "uploadToRedcap")
    rcon <- redcapAPI::redcapConnection(apiLink, gb$ApiToken)
    runID <- ifelse(is.null(runNumb), gb$runID, runNumb)
    message(paste(file.list))
    htmlLi <- stringr::str_replace_all(basename(file.list), "_QC_FAILED", "")
    htmlLi <- stringr::str_replace_all(basename(htmlLi), ".html", "")
    names(file.list) <- htmlLi
    for (n_file in 1:length(file.list)) {
        recordName <- names(file.list)[n_file]
        fiPath <- file.list[[n_file]]
        is_validation <- sjmisc::str_contains(runNumb, "VAL")
        has_val <- sjmisc::str_contains(recordName, "VAL")
        if (is_validation == T & has_val == F) {
            recordName <- paste0(recordName, "_VAL")
        }
        recordName2 <- CheckSarcRDnumber(recordName)
        has_val <- sjmisc::str_contains(recordName2, "VAL")
        if (is_validation == T & has_val == F) {
            recordName2 <- paste0(recordName2, "_VAL")
        }
        callApiImport(rcon, recordName2, runID)
        isEmpty <- checkRedcapRecord(recordName2) == ''

        if (isEmpty == F) {
            message(paste(recordName, "already has a file in REDCap"))
            message("Overwrite uploading will be set to FALSE")
        }

        callApiFile(rcon, recordName, isEmpty, fiPath)
    }
    if (deskCSV == T) {
        try(importDesktopCsv(rcon), outFile = "importDesktopRedcapLog.txt")
    }
}

# Imports the xlsm sheet 3 data
importSingle <- function(sh_Dat) {
    msgFunName(cpOutLnk, "importSingle")
    sh_Dat <- AddPngFilePath(sh_Dat)
    is_validation <- sjmisc::str_contains(sh_Dat$run_number[1], "VAL")
    is_val <- sjmisc::str_contains(sh_Dat$record_id, "VAL")
    if (is_validation == T & is_val == F) {
        sh_Dat$record_id <- paste0(sh_Dat$record_id, "_VAL")
    }
    recordEmpty <- checkRedcapRecord(sh_Dat$record_id, fieldName = "well_number")
    record = sh_Dat$record_id
    if (recordEmpty == '') {
        loopRedcapImport(sh_Dat)
    } else{
        message(crayon::white$bgBlue("Record Data not Uploaded:"), "\n", record[1])
    }

    htmlEmpty <- checkRedcapRecord(paste0(record[1]), fieldName = "classifier_pdf")
    if (htmlEmpty == '') {
        uploadToRedcap(file.list = paste0(record[1], ".html"), deskCSV = F)
    } else{
        message(crayon::white$bgRed("Record already has an HTML in REDCap:"), "\n", record[1])
    }
}


MakeOutputDir <- function(runYear, clinDrv, runID, isMC) {
    msgFunName(cpOutLnk, "MakeOutputDir")
    researchOutDir = file.path(rschDrv, runID)
    clinicalOutDir = file.path(clinDrv, "Results", runYear, runID)
    runYear = ifelse(isMC, paste0("20", stringr::str_split_fixed(runID, "-", 2)[1]), runYear)
    newFolder <- ifelse(isMC == T, clinicalOutDir, researchOutDir)
    cat(crayon::white$bgCyan("Output Folder is:\n", newFolder))
    CheckDirMake(newFolder)
    return(newFolder)
}


CallApiFileForce <- function(rcon, recordName) {
    msgFunName(cpOutLnk, "CallApiFileForce")
    is_validation <- sjmisc::str_contains(gb$runID, "VAL")
    has_val <- sjmisc::str_contains(recordName, "VAL")
    if (is_validation == T & has_val == F) {
        recordName <- paste0(recordName, "_VAL")
    }

    recordFi <- paste0(recordName, ".html")
    message("\n", gb$mkBlue("Importing Record File:"), paste0(" ", recordFi))
    fiPath <- file.path(getwd(), recordFi)
    fld <- "classifier_pdf"

    body <- list(
        token = rcon$token,
        content = 'file',
        action = 'import',
        record = recordName,
        field = fld,
        file = httr::upload_file(fiPath),
        returnFormat = 'csv'
    )
    res <-
        tryCatch(
            httr::POST(url = rcon$url, body = body, config = rcon$config),
            error = function(cond) {
                list(status_code = "400")
            }
        )
    if (res$status_code == "200") {
        message("REDCap file upload successful: ", fiPath)
    }else{
        message("REDCap file upload failed: ", fiPath)
    }

}


ForceCallApiFile <- function(rcon, recordName, ovwr = T) {
    msgFunName(cpOutLnk, "ForceCallApiFile")
    uploadField <- "classifier_pdf" # Default Value
    if (stringr::str_detect(recordName, "_sarc")) {
        uploadField <- "classifier_pdf_other"
    }
    if (stringr::str_detect(recordName, "_sarcV13")) {
        uploadField <- "sarcoma_v13_report"
    }

    recordName <- CheckSarcRDnumber(recordName)
    recordFi <- dir(getwd(), ignore.case = T, full.names = T,
                    pattern = paste0(recordName, ".*\\.html$"))[1]

    message("\n", mkBlue("Uploading File:"), " ", recordFi)
    message(mkBlue("To Record:"), " ", recordName)
    message(mkBlue("In Field:"), " ", uploadField)

    if (ovwr == F) {
        log_fi_out <- paste(gb$runID, "import_log.tsv", sep = "_")
        writeLogFi(recordName, logFile = log_fi_out)
    } else{
        tryCatch(
            expr = {
                suppressWarnings(
                    redcapAPI::importFiles(
                        rcon = rcon,
                        file = recordFi,
                        record = recordName,
                        field = uploadField,
                        overwrite = ovwr,
                        repeat_instance = 1
                    )
                )
            },
            error = function(e) {
                message(recordFi, " was not imported to REDCap")
                message(mkRed(e$message))
            }
        )
    }
}


UploadToRedcapOnly <- function(file.list, token=NULL) {
    msgFunName(cpOutLnk, "UploadToRedcapOnly")
    stopifnot(!is.null(token))
    MsgDF(file.list)
    rcon <- redcapAPI::redcapConnection(gb$apiLink, token)
    htmlLi <- stringr::str_replace_all(basename(file.list), ".html", "")
    MsgDF(htmlLi)
    for (recordName in htmlLi) {
        CallApiFileForce(rcon, recordName)
    }
}


ForceUploadToRedcap <- function(file.list, token = NULL, deskCSV = T) {
    msgFunName(cpOutLnk, "ForceUploadToRedcap")
    stopifnot(!is.null(token))
    print(as.data.frame(file.list))

    rcon <- redcapAPI::redcapConnection(apiLink, token)
    htmlLi <- stringr::str_replace_all(basename(file.list), ".html", "")

    MsgDF(htmlLi)
    for (recordName in htmlLi) {
        ForceCallApiFile(rcon, recordName, T)
    }
    if (deskCSV == T) {
        try(importDesktopCsv(rcon), outFile = "importDesktopRedcapLog.txt")
    }
}


grabAllRecords <- function(flds, rcon) {
    msgFunName(cpOutLnk, "grabAllRecords")
    message("Pulling REDCap data...")
    library("dplyr")
    params = list(rcon, fields = flds, survey = F, dag = F,
                  factors = F, form_complete_auto = F)
    dbCols <- do.call(redcapAPI::exportRecordsTyped, c(params))
    rd_df <- as.data.frame(dbCols)
    rd_df <- rd_df %>% dplyr::mutate_all(~stringr::str_replace_all(., ",", ""))
    return(rd_df)
}


GrabSpecificRecords <- function(flds, rd_num, rcon) {
    msgFunName(cpOutLnk, "GrabSpecificRecords")
    message("Pulling REDCap data...")
    library("dplyr")
    params = list(rcon, records = rd_num, fields = flds, survey = F, dag = F, factors = F, form_complete_auto = F)
    dbCols <- do.call(redcapAPI::exportRecordsTyped, c(params))
    rd_df <- as.data.frame(dbCols)
    redcap_db <- rd_df %>% dplyr::mutate_all(~stringr::str_replace_all(., ",", ""))
    return(redcap_db)
}


rename_reorder_output <- function(output) {
    library("rlang")
    library("dplyr")
    rename_map <- c(
        "RD-number" = "RD.number",
        "B-number" = "B.number",
        "TM-number" = "TM.number",
        "Log2sqrt(M*U)" = "Log2sqrt.M.U.",
        "Log2(M/U)" = "Log2.M.U.",
        "log2sqrt(R*G)" = "log2sqrt.R.G.",
        "log2(R/G)" = "log2.R.G.",
        "BS_log2sqrt(R*G)" = "BS_log2sqrt.R.G.",
        "BS_log2(R/G)" = "BS_log2.R.G.",
        "log2sqrt(H*L)" = "log2sqrt.H.L.",
        "log2(H/L)" = "log2.H.L."
    )
    MsgDF(output)

    output <- output %>% dplyr::rename(!!!rename_map)

    ordered_cols <- c(
        "RunID",
        "RD-number",
        "B-number",
        "TM-number",
        "Log2sqrt(M*U)",
        "Log2(M/U)",
        "log2sqrt(R*G)",
        "log2(R/G)",
        "BS_log2sqrt(R*G)",
        "BS_log2(R/G)",
        "log2sqrt(H*L)",
        "log2(H/L)",
        "Pvalue",
        "record_id",
        "subgroup_score",
        "classifier_score",
        "classifier_value",
        "subgroup",
        "b_number",
        "primary_tech",
        "run_number",
        "accession_number",
        "block",
        "tm_number",
        "diagnosis"
    )

    if (any(is.na(output))) {
        output[is.na(output)] <- ""
    }
    message("\n>>>>Ordering columns:\n",
            paste(ordered_cols, collapse = "\n"))
    message("\n>>>>Current columns:\n",
            paste(colnames(output), collapse = "\n"))
    if (!"record_id" %in% colnames(output)) {
        output$record_id <- output$`RD-number`
    }
    missing_cols <- !(ordered_cols %in% colnames(output))

    if (any(missing_cols)) {
        missed_cols <- ordered_cols[missing_cols]
        message("The following column(s) are missing:\n",
                paste(missed_cols, collapse = "\n"))
    }

    final_output <- output[, ordered_cols]
    return(final_output)
}


CheckOutputScoresQC <- function(output, runID, redcap_db, fieldsToPull, rcon) {
    msgFunName(cpOutLnk, "CheckOutputScoresQC")
    newCols <- colnames(redcap_db)
    library("dplyr")
    for (col_id in newCols) {
        output[col_id] <- ""
    }
    rownames(redcap_db) <- redcap_db$record_id

    for (xrow in 1:nrow(output)) {
        currRow <- output$`RD.number`[xrow]
        output[xrow, newCols] <- redcap_db[currRow, newCols]
    }

    totalNA <- which(is.na(output$classifier_value))
    if (length(totalNA) > 3) {
        find_redCsv <- dir(path = getwd(), pattern = "_Redcap.csv")[1]
        if (is.na(find_redCsv)) {
            desk_path <- file.path(fs::path_home(), "Desktop", runID)
            find_redCsv <- dir(path = desk_path, pattern = "_Redcap.csv")[1]
        }
        redcap_dat <- as.data.frame(read.csv(find_redCsv))
        redcap_dat <- unique(redcap_dat)
        redcapAPI::importRecords(rcon, redcap_dat, "normal", "ids", logfile = "REDCapImportLog.txt")
        if (nrow(redcap_dat) != nrow(redcap_db)) {
            if (nrow(redcap_dat) > nrow(redcap_db)) {
                missing_sam <- which(!redcap_dat$record_id %in% redcap_db$record_id)
                message("Sample is missing from REDCap Pull:",
                        redcap_dat$record_id[missing_sam])
                redcap_dat <- redcap_dat[-missing_sam,]

            }else{
                missing_sam <- which(!redcap_db$record_id %in% redcap_dat$record_id)
                message("Sample is missing from REDCap Pull:",
                        redcap_db$record_id[missing_sam])
                redcap_db <- redcap_db[-missing_sam,]
            }
        }
        redcap_dat$block <- redcap_db$block
        redcap_dat$diagnosis <- redcap_db[redcap_dat$record_id, "diagnosis"]
        redcap_dat$tm_number <- redcap_db$tm_number
        missing_cols <- setdiff(fieldsToPull, colnames(redcap_dat))

        if (length(missing_cols) > 0) {
            for (n_col in missing_cols) {
                redcap_dat[, n_col] <- redcap_db[, n_col]
            }
        }

        data_subset <- redcap_dat[ , fieldsToPull]

        newCols <- c(
            "record_id",
            "subgroup_score",
            "classifier_score",
            "classifier_value",
            "subgroup",
            "b_number",
            "primary_tech",
            "run_number",
            "accession_number",
            "block",
            "tm_number",
            "diagnosis"
        )
        for (new_col in newCols) {
            output[, new_col] <- redcap_dat[, new_col]
        }

        # if (is.null(redcap_dat$tm_number)) {
        #     redcap_dat$tm_number <- output$tm_number
        # }

        if (all(is.na(output$subgroup_score))) {
            output$subgroup_score <- redcap_dat$subgroup_score
        }

        if (all(is.na(output$classifier_score))) {
            output$classifier_score <- redcap_dat$classifier_score
        }

        if (all(is.na(output$classifier_value))) {
            output$classifier_value <- redcap_dat$classifier_value
        }

        if (all(is.na(output$subgroup))) {
            output$subgroup <- redcap_dat$subgroup
        }

        filtered_output <- output %>% dplyr::select(-one_of(colnames(data_subset)))
        merged_data <- merge(filtered_output, data_subset, by.x = "RD.number",
                             by.y = "record_id")
        MsgDF(merged_data)
        return(merged_data)
    } else {
        return(output)
    }
}


get_QC_metric_data <- function(output_fi, runDir, runID) {
    msgFunName(cpOutLnk, "get_QC_metric_data")
    runYear <- paste0(20, stringr::str_split_fixed(runID, "-", 2)[1, 1])
    if (is.null(output_fi)) {
        xlsx_file <- paste("Meth_QC_metrics", runYear,"runs.xlsx", sep = "_")
        file_to_read <- file.path(metrics_dir, runYear, xlsx_file)
        output <- as.data.frame(readxl::read_excel(file_to_read))
    } else{
        currDir <- file.path(runDir, output_fi)
        message("Reading file: ", currDir)
        output <- as.data.frame(read.csv(currDir))
    }
    return(output)
}


CombineClassAndQC <- function(output_fi = NULL, token, runDir = NULL, runID = NULL) {
    msgFunName(cpOutLnk, "CombineClassAndQC")

    apiUrl = "https://redcap.nyumc.org/apps/redcap/api/"
    rcon <- redcapAPI::redcapConnection(apiUrl, token)

    fieldsToPull <- c(
        "record_id", "run_number", "b_number", "tm_number", "block",
        "accession_number", "subgroup_score", "classifier_value", "subgroup",
        "classifier_score", "primary_tech", "diagnosis"
    )

    if (is.null(runDir)) {runDir <- getwd()}
    if (is.null(runID)) {runID <- basename(getwd())}

    output <- get_QC_metric_data(output_fi, runDir, runID)

    controlRows <- which(stringr::str_detect(output$`RD.number`, "control"))
    output[controlRows, "RD.number"] <-
        paste(output[controlRows, "RunID"], output[controlRows, "RD.number"], sep = "_")

    redcap_db <- GrabSpecificRecords(flds = fieldsToPull, rd_num = output$RD.number, rcon)

    output <- CheckOutputScoresQC(output, runID, redcap_db, fieldsToPull, rcon)
    final_output <- rename_reorder_output(output)
    MsgDF(final_output)
    outFile <- file.path(runDir, paste(runID, "QC_and_Classifier_Scores.csv", sep = "_"))
    write.csv(final_output, file = outFile, row.names = F, quote = F)

}
