#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE);gb <- globalenv(); assign("gb", gb)
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
cpOutLnk = "https://github.com/NYU-Molecular-Pathology/Methylation/edit/main/CopyOutput.R"
rschDrv = "/Volumes/snudem01labspace/FINAL_PDF_Reports_Brain"
msgFunName <- function(pthLnk, funNam){
    message("\nExecuting function: ", funNam, " from RScript in:\n", pthLnk,"\n")
}
mkBlue <- function(strVar){return(crayon::white$bgBlue(strVar))}
mkGrn <- function(strVar){return(crayon::white$bgGreen$bold(strVar))}

# FUN: Generate CNV image
saveCNVhtml <- function(data) {
    msgFunName(cpOutLnk, "saveCNVhtml")
    # get RGSet for sample
    RGsetEpic <- minfi::read.metharray(file.path(getwd(),'idats',unlist(data[4])), verbose = F, force = T)
    # Save the CNV as a png file
    MsetEpic <- minfi::preprocessRaw(RGsetEpic)
    sample_id <- unlist(data[1])
    xx <- mnp.v11b6::MNPcnv(MsetEpic,sex = unlist(data[6]),main = sample_id)
    hg <- mnp.v11b6::MNPcnvggplotly(xx,getTables = F)
    htmlwidgets::saveWidget(hg, paste(sample_id, "_cnv.html", sep = ""))
}

# Helper archive function: creates a "previous" folder when reports already exists
save.prev.folder <- function(prevs, oldFi){
    msgFunName(cpOutLnk, "save.prev.folder")
    message("\n",crayon::white$bgRed('Output folder already exists, moving existing reports to new folder named "Previous"'),"\n")
    dir.create(prevs)
    fs::file_copy(path=oldFi, new_path=prevs)
}

CopyHtmlFiles <- function(newFolder) {
    msgFunName(cpOutLnk, "CopyHtmlFiles")
    message(mkBlue("Copying Reports to output folder:"), "\n", newFolder)
    message(mkGrn("Now copying html reports..."),"\n")
    fi2copy <- dir(getwd(), pattern = "*.html", full.names = T)
    print(fi2copy)
    file.copy(fi2copy, newFolder, overwrite=F, copy.mode = F, copy.date = T)
    #fs::file_copy(fi2copy, newFolder, overwrite=F)
}

CheckDirMake <- function(newFolder){
    if (!dir.exists(newFolder)) {dir.create(newFolder)}
}

# FUN: Copies Reports to Z drive
copy.to.clinical <- function(clinOut, runID, runYear) {
    msgFunName(cpOutLnk, "copy.to.clinical")

    newFolder <- file.path(clinOut, runYear, runID)

    CheckDirMake(newFolder)
    if (dir.exists(newFolder)) {
        oldFi=dir(path=newFolder, full.names=T)
        prevs=file.path(newFolder,"previous")
        if (length(oldFi) > 0) {save.prev.folder(prevs,oldFi)}
        CopyHtmlFiles(newFolder)
    }
}

# Checks if field is already filled in REDCap returns boolean
checkRedcapRecord <- function(recordName, fieldName='classifier_pdf'){
    url = gb$apiLink
    formData <- list(
        "token"=gb$ApiToken,
        content='record',
        action='export',
        format='json',
        type='flat',
        csvDelimiter='',
        'records[0]'=recordName,
        #'fields[0]'='classifier_pdf',
        rawOrLabel='raw',
        rawOrLabelHeaders='raw',
        exportCheckboxLabel='false',
        exportSurveyFields='false',
        exportDataAccessGroups='false',
        returnFormat='json'
    )
    response <- httr::POST(url, body = formData, encode = "form")
    result <- httr::content(response)[[1]]
    return(result[[fieldName]]=="")
}

writeLogFi <- function(recordName, isHtml=T){
    logFile = "upload_log.tsv"
    message("Check ",logFile,"\n", crayon::white$bgRed(recordName[1]), " already has an data in REDCap:")

    if(isHtml==T){
    i = paste(recordName[1], "already has an html file in REDCap\n")
    write.table(i, file = logFile, append = TRUE, quote = F, sep = '\t', row.names = F, col.names = F)
    }else{
        i = "\n\nThe following data fields already have classifier_value filled in REDCap:\n"
        write.table(i, file = logFile, append = TRUE, quote = F, sep = '\t', row.names = F, col.names = F)
        i = paste(names(recordName), recordName, sep=":", collapse=" ")
        write.table(i, file = logFile, append = TRUE, quote = F, sep = ' ', row.names = F, col.names = F)
    }


}

loopRedcapImport <-function(data){
    for (n in 1:nrow(data)) {
        datarecord = jsonlite::toJSON((as.list(data[n,])), auto_unbox=T)
        message("~~",crayon::bgBlue("Record Uploaded:"))
        print(datarecord)
        RCurl::postForm(
            apiLink,
            token = gb$ApiToken,
            content = 'record',
            format = 'csv',
            type = 'flat',
            data = datarecord,
            returnFormat = 'csv',
            overwriteBehavior = 'normal'
        )
    }
}

# REDCap: API call & Upload
# uploads the redcap classifier values must convert to JSON first
importDesktopCsv <- function(rcon,samsheet=NULL) {
    msgFunName(cpOutLnk, "importDesktopCsv")
    rcon <- redcapAPI::redcapConnection(apiLink, gb$ApiToken)
    ur=paste0(rcon$url);tk=rcon$token
    if(is.null(samsheet)){
        samsheet=dir(path=getwd(),full.names=T,"_Redcap.csv",recursive=F)} else {samsheet=samsheet}
    if (length(samsheet) < 1) {message("Redcap headers csv file not found in:\n", paste(samsheet))}
    if (length(samsheet) == 1) {
        data<-read.csv(samsheet, stringsAsFactors=F)
        if(any(duplicated(data$record_id))){
            warning("Remove duplicate rows in the REDCap csv dataframe:")
            print(data$record_id[duplicated(data$record_id)])
            data = data[!duplicated(data$record_id),]
        }
        message("Checking REDCap for existing data:")
        toImport <- unlist(lapply(data$record_id, FUN=function(rd){checkRedcapRecord(rd, "classifier_value")}))
        if(any(!toImport)){
        message("The following records already have an existing classifier_value and will not be over-written in REDCap:")
            toSkip <- data[!toImport,]
            print(toSkip)
            invisible(lapply(1:length(toSkip), function(x){writeLogFi(as.data.frame(toSkip[x,]), isHtml=F)}))
            data <- data[toImport,]
        }
        if (nrow(data > 0)) {
            cat(redcapAPI::importRecords(rcon, data, "normal", "ids", logfile = "REDCapImportLog.txt"))
        }else{message("No new data to import to REDCap")}
    } else {message("no _Redcap.csv file found")}
}

# Copy Output cnv Files if generated
copy.cnv.files <- function(newFolder, runID, runYear=NULL) {
    msgFunName(cpOutLnk, "copy.cnv.files")
    if (is.null(runYear)){runYear=paste0(format(Sys.Date(), "%Y"))}
    cnv_folder <- file.path(newFolder,paste0(runID,"_CNVs/"))
    cnvNames <- dir(path=getwd(), full.names=T, "*_cnv.png")
    if(length(cnvNames)>2){
        message(paste0("Copying PNG to: ", cnv_folder));print(as.data.frame(cnvNames))
        CheckDirMake(cnv_folder)
        if (dir.exists(cnv_folder)){fs::file_copy(cnvNames,cnv_folder)}
    }
}

GrabSampleSheet <- function(){
    samSh <- dir(path=getwd(), full.names=T, ".xlsm")
    if(length(samSh)>1){
        message("Multiple samplesheets found:\n")
        print(samSh)
        samSh <- samSh[stringr::str_detect(samSh,pattern = "\\$",negate = T)]
    }
    message("Using following samplesheet:\n", samSh[1])
    stopifnot(length(samSh)>0)
    return(samSh[1])
}

# Returns Total Sample Count in the run
getTotalSamples <- function(thisSh=NULL){
    msgFunName(cpOutLnk, "getTotalSamples")

    templateDir = "Clinical_Methylation/methylation_run_TEMPLATE.xlsm"
    thisSh <-  ifelse(is.null(thisSh), GrabSampleSheet(), thisSh)
    thisSh <- thisSh[!stringi::stri_detect_fixed(thisSh, "~$")]
    if(length(thisSh)==0){print("No .xlsm sheet, defaulting to 16 total samples");return(16)}
    worksheet <- suppressMessages(readxl::read_excel(thisSh[1], col_names="Total", range="B4:B4"))
    if (length(worksheet) == 0) {
        warning("Samplesheet ", thisSh[1]," is invalid format or no integer in Cell B4 found.")
        message("Manually edit samplesheet to fix-Try copying data into the template file:\n", templateDir)
        stopifnot(length(worksheet)>0)
    } else {message("Total sample count found is: ", worksheet[1])}
    totNumb <- paste0(worksheet[1])
    return(as.integer(totNumb))
}


AddPngFilePath <- function(nfldr,sh_Dat){
    pathNam = file.path(nfldr,paste0(gb$runID,"_CNVs"), paste0(sh_Dat$record_id, "_cnv.png"))
    rms=paste(c("control_","low_"),collapse ='|')
    pathNam <- stringr::str_replace(pathNam, rms, "")
    pathNam <- stringr::str_replace(pathNam, "//", "/")
    sh_Dat$cnv_file_path <- pathNam
    return(sh_Dat)
}

# Uploads any created cnv png files to redcap database
uploadCnPng <- function() {
    msgFunName(cpOutLnk, "uploadCnPng")
    rcon <- redcapAPI::redcapConnection(apiLink, gb$ApiToken)
    samSh <- GrabSampleSheet()
    sampleNumb <- getTotalSamples()
    nfldr = file.path(stringr::str_split_fixed(gb$clinDrv, " ", 2)[1],"MethylationClassifier")
    sh_Dat <-suppressMessages(
        as.data.frame(readxl::read_excel(samSh,sheet=3,range="A1:N97", col_types=c("text")))[1:sampleNumb,1:13])
    sh_Dat <- AddPngFilePath(nfldr,sh_Dat)
    records <- sh_Dat$record_id
    for (idx in 1:length(records)) {
        pth = sh_Dat$cnv_file_path[idx]
        recordName = paste0(records[idx])
        message(mkBlue("Importing CNV Record:"),"\n", recordName, " ", pth)
        redcapAPI::importFiles(
            rcon,
            pth,
            recordName,
            field = "methyl_cn",
            overwrite = T,
            repeat_instance = 1
        )
    }
}


# Imports the xlsm sheet 3 data
importRedcapStart <- function(nfldr){
    msgFunName(cpOutLnk, "importRedcapStart")
    samSh <- GrabSampleSheet()
    sampleNumb <- getTotalSamples(samSh)
    sh_Dat <-
        suppressMessages(as.data.frame(readxl::read_excel(samSh,sheet=3,range="A1:N97", col_types=c("text")))[1:sampleNumb,1:13])

    sh_Dat <- AddPngFilePath(nfldr,sh_Dat)

    if (!is.null(sh_Dat)){
        for(n in 1:nrow(sh_Dat)){
            record=c(sh_Dat[n,])
            datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
            res <-
                RCurl::postForm(
                    uri= apiLink,
                    token = gb$ApiToken,
                    content = 'record',
                    format = 'json',
                    type = 'flat',
                    data = datarecord,
                    overwriteBehavior = 'normal'
                )
            message("Record Uploaded")
            print(res)
        }
    }
}

callApiImport <- function(rcon,
                          recordName,
                          runID){
    message(mkBlue("Importing Record Data:"))
    data = data.frame(record_id = recordName, run_number = runID)
    logfi = paste0(recordName,"_redcapLog.txt")
    tryCatch(
        expr = {
            cat(redcapAPI::importRecords(
                rcon,
                data,
                overwriteBehavior = "normal",
                returnContent = "ids",
                logfile=logfi
            ),sep = "\n\n")
        },
        error = function(e) {message(crayon::white$bgRed(paste(data$record_id, "failed import data to REDCap:")), "\n", e$message)}
    )
}

callApiFile <- function(rcon,
                        recordName,
                        ovwr=T)
{
    recordFi <- paste0(recordName, ".html")
    message("\n", mkBlue("Importing Record File:"),paste0(" ", recordFi))

    if(ovwr==F){writeLogFi(recordName)}
    tryCatch(
        expr = {
            suppressWarnings(redcapAPI::importFiles(
                rcon = rcon,
                file = file.path(getwd(), recordFi),
                record = recordName,
                field = "classifier_pdf",
                overwrite = ovwr,
                repeat_instance=1
            ))
        },
        error = function(e) {
            message(recordFi, " was not imported to REDCap")
            message(crayon::white$bgRed(e$message))
        }
    )
}


# Creates QC record and uploads reports to redcap
uploadToRedcap <- function(file.list,
                           deskCSV = T,
                           runNumb = NULL)
{
    msgFunName(cpOutLnk, "uploadToRedcap")
    rcon <- redcapAPI::redcapConnection(apiLink, gb$ApiToken)
    runID <- ifelse(is.null(runNumb), gb$runID, runNumb)
    if (deskCSV == F) {
        htmlLi <- stringr::str_replace_all(basename(file.list), ".html", "")
        for (recordName in htmlLi) {
            callApiImport(rcon, recordName, runID)
            isEmpty <- checkRedcapRecord(recordName)
            callApiFile(rcon, recordName, isEmpty)
        }
    }
    if (deskCSV == T) {importDesktopCsv(rcon)}
}

# Imports the xlsm sheet 3 data
importSingle <- function(sh_Dat) {
    msgFunName(cpOutLnk, "importSingle")
    nfldr = file.path(stringr::str_split_fixed(gb$clinDrv, " ", 2)[1],"MethylationClassifier")
    sh_Dat <- AddPngFilePath(nfldr,sh_Dat)
    record = c(sh_Dat[1,])
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox = T)
    print(sh_Dat)
    res <- RCurl::postForm(
        uri= apiLink,
        token = gb$ApiToken,
        content = 'record',
        format = 'json',
        type = 'flat',
        data = datarecord,
        overwriteBehavior='normal'
        )
    message(crayon::white$bgBlue("Record Uploaded"))
    invisible(res)
    uploadToRedcap(file.list = paste0(record[1], ".html"), deskCSV = F)
}

MakeOutputDir <- function(runYear, clinDrv, runID, isMC){
    researchOutDir = file.path(rschDrv,runID)
    clinicalOutDir = file.path(clinDrv, "Results", runYear, runID)
    runYear = ifelse(isMC, paste0("20", stringr::str_split_fixed(runID, "-", 2)[1]), runYear)
    newFolder <- ifelse(isMC==T, clinicalOutDir, researchOutDir)
    cat(crayon::white$bgCyan("Output Folder is:\n", newFolder))
    CheckDirMake(newFolder)
    return(newFolder)
}

CopyFilesOut <- function(file.list, newFolder){
    message("\nCopying Existing Reports to Folder...\n", newFolder,"\n",mkBlue("Files to copy:"), "\n")
    print(file.list)
    #lapply(file.list, function(foo) {
        #destDir = file.path(newFolder, basename(foo))
        tryCatch(
            #expr = {fs::file_copy(foo, destDir, overwrite = F)},
            expr = {file.copy(file.list, newFolder, overwrite = F, copy.mode = F)},
            error = function(e) {message(e,"\nTrying other file copy method:\n")
                if(!file.exists(file.list)) {message("File ", foo, " does not exist")} else{
                    cmnd = paste("cp -p", file.list, file.path(newFolder, basename(file.list)))
                    print(cmnd)
                    for (foo in cmnd) {
                        system(foo)
                    }
                    }}
            )
   # })
}


# FUN: Copies Reports to Z drive
copy2outFolder <-function(clinDrv = NULL, runID, runYear = NULL) {
    msgFunName(cpOutLnk, "copy2outFolder")

    runYear <- ifelse(is.null(runYear), paste0(format(Sys.Date(), "%Y")), runYear)
    clinDrv <- ifelse(is.null(clinDrv), gb$clinDrv, clinDrv)
    isMC = sjmisc::str_contains(runID, "MGDM") | sjmisc::str_contains(runID, "MC")

    newFolder <- MakeOutputDir(runYear, clinDrv, runID, isMC)
    oldFi = dir(path = newFolder, full.names = T)

    if (length(oldFi) > 0) {
        save.prev.folder(file.path(newFolder, "previous"), oldFi) # saves any old files
    }

    file.list = dir(path = getwd(), "*.html", full.names = T)
    CopyFilesOut(file.list, newFolder)

    if (isMC) {
        cnList <- dir(getwd(), "_cnv.png", recursive = F)
        if (length(cnList) > 2) {
            copy.cnv.files(newFolder, runID, runYear)
            uploadCnPng()
        }
        clinOut = file.path(stringr::str_split_fixed(clinDrv, " ", 2)[1],"MethylationClassifier")
        importRedcapStart(clinOut)
        copy.to.clinical(clinOut, runID, runYear)
    }
    return(file.list)
}
