#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
gb <- globalenv(); assign("gb", gb)
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
cpOutLnk = "https://github.com/NYU-Molecular-Pathology/Methylation/edit/main/CopyOutput.R"
rschDrv = "/Volumes/snudem01labspace/FINAL_PDF_Reports_Brain"

msgFunName <- function(pthLnk, funNam) {
    message("\nExecuting function: ", funNam, " from RScript in:\n", pthLnk, "\n")
}

mkBlue <- function(strVar) {return(crayon::white$bgBlue(strVar))}
mkGrn <- function(strVar) {return(crayon::white$bgGreen$bold(strVar))}
mkRed <- function(strVar) {return(crayon::white$bgRed$bold(strVar))}

CheckDirMake <- function(newFolder) {if (!dir.exists(newFolder)) {fs::dir_create(newFolder)}}

MsgDF <- function(datObj) {
    if(class(datObj) != "data.frame"){datObj <- as.data.frame(datObj)}
    message(paste0(capture.output(datObj), collapse="\n"))
}


# FUN: Generate CNV image
saveCNVhtml <- function(data) {
    msgFunName(cpOutLnk, "saveCNVhtml")
    # get RGSet for sample
    RGsetEpic <- minfi::read.metharray(file.path(getwd(), 'idats', unlist(data[4])), verbose = F, force = T)
    # Save the CNV as a png file
    MsetEpic <- minfi::preprocessRaw(RGsetEpic)
    sample_id <- unlist(data[1])
    xx <- mnp.v11b6::MNPcnv(MsetEpic, sex = unlist(data[6]), main = sample_id)
    hg <- mnp.v11b6::MNPcnvggplotly(xx, getTables = F)
    htmlwidgets::saveWidget(hg, paste(sample_id, "_cnv.html", sep = ""))
}


# Helper archive function: creates a "previous" folder when reports already exists
SavePrevDir <- function(newFolder) {
    msgFunName(cpOutLnk, "SavePrevDir")
    message("\n", mkRed( 'Output folder already exists, moving existing reports to new folder named "Previous"'), "\n")
    oldFi = dir(path = newFolder, full.names = T)
    prevs = file.path(newFolder, "previous")
    dir.create(prevs)
    fs::file_copy(path = oldFi, new_path = prevs)
}


Copy2TempDir <- function(fi2copy, runID){
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
    if(length(fi2copy)>=1){
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
    if (length(dir(path = newFolder, pattern=".html", full.names = T)) > 0) {
        SavePrevDir(newFolder)
    }else{
        CopyHtmlFiles(newFolder, runID)
    }
}


# Checks if field is already filled in REDCap returns boolean
checkRedcapRecord <- function(recordName, fieldName = 'classifier_pdf') {
      msgFunName(cpOutLnk, "checkRedcapRecord")    
    url = gb$apiLink
        formData <- list(
            "token" = gb$ApiToken,
            content = 'record',
            action = 'export',
            format = 'json',
            type = 'flat',
            csvDelimiter = '',
            'records[0]' = recordName,
            #'fields[0]'='classifier_pdf',
            rawOrLabel = 'raw',
            rawOrLabelHeaders = 'raw',
            exportCheckboxLabel = 'false',
            exportSurveyFields = 'false',
            exportDataAccessGroups = 'false',
            returnFormat = 'json'
        )
        response <- httr::POST(url, body = formData, encode = "form")
        result <- httr::content(response)[[1]]
        return(result[[fieldName]])
    }


MakeLogFile <- function(i, logFile){
    write.table(i, file = logFile, append = TRUE, quote = F, sep = '\t', row.names = F, col.names = F)
}


writeLogFi <- function(recordName, isHtml = T, logFile = "upload_log.tsv") {
    message("Check ", logFile, "\n", mkRed(recordName[1]), " already has data in REDCap")
    if (isHtml == T) {
        i = paste(recordName[1], "already has an html file in REDCap\n")
        MakeLogFile(i, logFile)
    } else{
        i = "\n\nThe following data fields already have classifier_value filled in REDCap:\n"
        MakeLogFile(i, logFile)
        i = paste(names(recordName), recordName, sep = ":", collapse = " ")
        MakeLogFile(i, logFile)
    }
}

RedcapRcurl <- function(datarecord) {
    message(crayon::bgBlue("Record Data Uploaded:"), datarecord)
    RCurl::postForm(
        apiLink,
        token = gb$ApiToken,
        content = 'record',
        format = 'json',
        type = 'flat',
        data = datarecord,
        overwriteBehavior = 'normal'
    )
}


WarnSentrix <- function(record, isEmpty) {
    if (record$barcode_and_row_column == isEmpty) {
        message(mkRed(record$record_id), " already has the same Sentrix ID on REDCap: ", isEmpty)
    } else{
        warning(mkRed("Sentrix ID is already on REDCap!"))
        message("Samplesheet Sentrix ID is: ", record$barcode_and_row_column)
        message("The Sentrix ID on REDCap is: ", isEmpty)
        writeLogFi(record$record_id, isHtml = F, logFile = "import_log.tsv")
    }
}


ValidateRedImport <- function(record) {
    isEmpty <- checkRedcapRecord(record$record_id, "barcode_and_row_column")
    if (isEmpty == "") {
        RedcapRcurl(datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox = T))
    } else{
        WarnSentrix(record, isEmpty)
    }
}


loopRedcapImport <- function(data) {
    if (!is.null(data)) {
        for (n in 1:nrow(data)) {
            record = c(data[n,])
            ValidateRedImport(record)
        }
    } else{
        message("Nothing imported to REDCap. Dataframe is null!")
    }
}


GetRedcapCsv <- function(samsheet) {
    msgFunName(cpOutLnk, "GetRedcapCsv")
    if (is.null(samsheet)) {
        samsheet <- dir(path = getwd(), full.names = T, "_Redcap.csv", recursive = F)
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
        invisible(lapply(1:nrow(toSkip), function(x) {
            writeLogFi(as.data.frame(toSkip[x, ]), isHtml = F)
        }))
        rawCsv <- rawCsv[toImport, ]
    }
    return(rawCsv)
}


# REDCap: API call & Upload --------------------------------------------------------------------------------
# uploads the redcap classifier values must convert to JSON first
importDesktopCsv <- function(rcon, samsheet = NULL) {
    msgFunName(cpOutLnk, "importDesktopCsv")
    rawCsv <- GetRedcapCsv(samsheet)
    data <- rawCsv
    #data <- CheckImportData(rawCsv)
    if (nrow(data > 0)) {
        cat(redcapAPI::importRecords(rcon, data, "normal", "ids", logfile = "REDCapImportLog.txt"))
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
            message(paste0(capture.output(data), collapse="\n"))
            message(mkRed(rdMsg), "\n", e$message)
            }
        )
}


CheckSarcRDnumber <- function(record){
    isSarc <- ifelse(stringr::str_detect(record, pattern="_sarc"), T, F)
    if(isSarc==T){
        record <- stringr::str_split_fixed(record, pattern = "_", 2)[1,1]
    }
    return(record)
}

callApiImport <- function(rcon, recordName, runID) {
    msgFunName(cpOutLnk, "callApiImport")
    isEmpty <- checkRedcapRecord(recordName, "subgroup")
    if (isEmpty == '') {
        DoRedcapApi(rcon, recordName, runID)
    } else{
        message(mkRed(recordName), " already has an assigned subgroup in REDCap: ", isEmpty)
        writeLogFi(recordName, isHtml = F, logFile = "import_log.tsv")
    }
}


callApiFile <- function(rcon, recordName, ovwr = T) {
    recordFi <- paste0(recordName, ".html")
    message("\n", gb$mkBlue("Importing Record File:"), paste0(" ", recordFi))
    fiPath <- file.path(getwd(), recordFi)
    if (ovwr == F) {
        writeLogFi(recordName)
    } else{
        
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
                error = function(cond){
                    list(status_code = "200")
                }
            )
        if(res$status_code=="200"){
            message("REDCap file upload successful: ", fiPath)
        }else{
            message("REDCap file upload failed: ", fiPath)
        }
        
    }
}


# Creates QC record and uploads reports to redcap
uploadToRedcap <- function(file.list, deskCSV = T, runNumb = NULL) {
    msgFunName(cpOutLnk, "uploadToRedcap")
    rcon <- redcapAPI::redcapConnection(apiLink, gb$ApiToken)
    runID <- ifelse(is.null(runNumb), gb$runID, runNumb)
    message(paste(file.list))
    htmlLi <- stringr::str_replace_all(basename(file.list), ".html", "")
    message(paste(htmlLi))
    for (recordName in htmlLi) {
        recordName2 <- CheckSarcRDnumber(recordName)
        callApiImport(rcon, recordName2, runID)
        isEmpty <- checkRedcapRecord(recordName2) == ''
        callApiFile(rcon, recordName, isEmpty)
    }
    if (deskCSV == T) {
        try(importDesktopCsv(rcon), outFile = "importDesktopRedcapLog.txt")
    }
}

# Imports the xlsm sheet 3 data
importSingle <- function(sh_Dat) {
    msgFunName(cpOutLnk, "importSingle")
    sh_Dat <- AddPngFilePath(sh_Dat)
    recordEmpty <- checkRedcapRecord(sh_Dat$record_id, fieldName = "well_number")
    record = sh_Dat$record_id
    if (recordEmpty == '') {
        loopRedcapImport(sh_Dat)
    } else{
        message(crayon::white$bgBlue("Record Data not Uploaded:"), "\n", record[1])
    }
    uploadToRedcap(file.list = paste0(record[1], ".html"), deskCSV = F)
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


RsyncCopyFiles <- function(file.list, newFolder) {
    msgFunName(cpOutLnk, "CopyFilesOut")
    if (any(!file.exists(file.list))) {
        missed <- paste(file.list[!file.exists(file.list)], "\n")
        message("File(s) do not exist:\n", missed)
    } else{
        missed <- !file.exists(file.path(newFolder, basename(file.list)))
        if (any(missed)) {
            file.list <- file.list[missed]
            cmnd = paste("rsync -a", file.list, file.path(newFolder, basename(file.list)))
            message(cmnd)
            for (foo in cmnd) {
                system(foo)
            }
        }
    }
}

CopyFilesOut <- function(file.list, newFolder, runID) {
    msgFunName(cpOutLnk, "CopyFilesOut")
    message("\nCopying Existing Reports to Folder...\n",
            newFolder, "\n", mkBlue("Files to copy:"), "\n")
    MsgDF(file.list)
    file.list <- Copy2TempDir(file.list, runID)
    tryCatch(
        expr = {
            fs::file_copy(file.list, newFolder, overwrite = F)
            },
        error = function(e) {
            message(e, "\n", mkRed("Trying other file copy method:"),"\n")
            RsyncCopyFiles(file.list, newFolder)
            }
        )
}


# FUN: Copies Reports to Z drive
copy2outFolder <- function(clinDrv = NULL, runID, runYear = NULL) {
    msgFunName(cpOutLnk, "copy2outFolder")

    runYear <- ifelse(is.null(runYear), paste0(format(Sys.Date(), "%Y")), runYear)
    clinDrv <- ifelse(is.null(clinDrv), gb$clinDrv, clinDrv)
    isMC = sjmisc::str_contains(runID, "MGDM") | sjmisc::str_contains(runID, "MC")
    newFolder <- MakeOutputDir(runYear, clinDrv, runID, isMC)
    oldFi = dir(path = newFolder, full.names = T)
    message(mkGrn("Clinical Drive output folder:"), "\n", newFolder)
    if (length(oldFi) > 0) {
        SavePrevDir(newFolder) # saves any old files
    }
    file.list <- dir(path = getwd(), ".html", full.names = T)
    CopyFilesOut(file.list, newFolder, runID)
    if (isMC==T) {
        clinOut = file.path(stringr::str_split_fixed(clinDrv, " ", 2)[1], "MethylationClassifier")
        message("Clinical Drive output folder:\n", clinOut)
        importRedcapStart(clinOut)
        copy.to.clinical(clinOut, runID, runYear)
    }
    if (length(dir(getwd(), "_cnv.png", recursive = F)) > 2) {
        copy.cnv.files(newFolder, runID, runYear)
    }
    return(file.list)
}


CallApiFileForce <- function(rcon, recordName) {
    msgFunName(cpOutLnk, "CallApiFileForce")
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
            error = function(cond){
                list(status_code = "200")
            }
        )
    if(res$status_code=="200"){
        message("REDCap file upload successful: ", fiPath)
    }else{
        message("REDCap file upload failed: ", fiPath)
    }
    
}


ForceCallApiFile <- function(rcon, recordName, ovwr = T) {
  msgFunName(cpOutLnk, "ForceCallApiFile")
  uploadField = "classifier_pdf"
  if(stringr::str_detect(recordName, "_sarc")){
    uploadField = "classifier_pdf_other"
  }
  recordFi <- paste0(recordName, ".html")
  message("\n", mkBlue("Importing Record File:"), paste0(" ", recordFi))
  recordName <- CheckSarcRDnumber(recordName)
  if (ovwr == F) {
    writeLogFi(recordName)
  } else{
    tryCatch(
      expr = {
        suppressWarnings(
          redcapAPI::importFiles(
            rcon = rcon,
            file = file.path(getwd(), recordFi),
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
    message(paste0(capture.output(file.list), collapse="\n"))
    rcon <- redcapAPI::redcapConnection(gb$apiLink, token)
    htmlLi <- stringr::str_replace_all(basename(file.list), ".html", "")
    message(paste0(capture.output(htmlLi), collapse="\n"))
    for (recordName in htmlLi) {CallApiFileForce(rcon, recordName)}
}


ForceUploadToRedcap <- function(file.list, token=NULL, deskCSV = T) {
    msgFunName(cpOutLnk, "ForceUploadToRedcap")
    stopifnot(!is.null(token))
    print(as.data.frame(file.list))
    msgFunName(cpOutLnk, "uploadToRedcap"); message(paste0(capture.output(file.list), collapse="\n"))
    rcon <- redcapAPI::redcapConnection(apiLink, token)
    htmlLi <- stringr::str_replace_all(basename(file.list), ".html", "")
    message(paste0(capture.output(htmlLi), collapse="\n"))
    for (recordName in htmlLi) {ForceCallApiFile(rcon, recordName, T)}
    if (deskCSV == T) {try(importDesktopCsv(rcon), outFile = "importDesktopRedcapLog.txt")}
}

CombineClassAndQC <- function(output = NULL, token, runDir = NULL, runID = NULL) {

    if(is.null(output)){
        currDir <- "/Volumes/molecular/Molecular/MethylationClassifier/Methylation_QC_metrics/2023"
        output <- as.data.frame(readxl::read_excel(file.path(currDir, "Meth_QC_metrics_2023_runs.xlsx")))
    }else{
        output <- as.data.frame(read.csv(file.path(runDir, output)))
    }

    fieldsToPull <- c("record_id", "run_number", "b_number", "tm_number", "block", "accession_number",
                      "subgroup_score", "classifier_value", "subgroup", "classifier_score", "primary_tech", "diagnosis")

    apiUrl = "https://redcap.nyumc.org/apps/redcap/api/"
    rcon <- redcapAPI::redcapConnection(apiUrl, token)
    db <- grabAllRecords(fieldsToPull, rcon)
    newCols <- colnames(db)

    for (col_id in newCols) {
        output[col_id] <- ""
    }

    controlRows <- which(output$`RD-number` == "control")
    output[controlRows, "RD-number"] <- paste(output[controlRows, "RunID"], output[controlRows, "RD-number"], sep = "_")

    for (xrow in 1:nrow(output)) {
        currRow <- which(db$record_id == output$`RD-number`[xrow])
        stopifnot(length(currRow) == 1)
        redcap_row <- db[currRow, newCols]
        output[xrow, newCols] <- redcap_row
    }
    outFile <- file.path(runDir, paste(runID,"QC_and_Classifier_Scores.csv", sep = "_"))
    write.csv(output, file = outFile, row.names = F, quote = F)

}

