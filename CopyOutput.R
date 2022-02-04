#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE);gb <- globalenv(); assign("gb", gb)
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
cpOutLnk = "https://github.com/NYU-Molecular-Pathology/Methylation/edit/main/CopyOutput.R"

msgFunName <- function(pthLnk, funNam){
message("\nExecuting function: ", funNam, " from RScript in:\n", pthLnk,"\n")
}

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
save.prev.folder <- function(prevs,oldFi){
        msgFunName(cpOutLnk, "save.prev.folder")
    cat("\n")
    cat(crayon::white$bgRed('Output folder already exists, moving existing reports to new folder named "Previous"'))
    cat("\n")
    dir.create(prevs) 
    fs::file_copy(path=oldFi, new_path=prevs)
}

# FUN: Copies Reports to Z drive
copy.to.clinical <- function(clinOut, runID, runYear) {
        msgFunName(cpOutLnk, "copy.to.clinical")
    newFolder <- file.path(clinOut, runYear,runID); message(newFolder)
    if (!dir.exists(newFolder)) {dir.create(newFolder)}
    if (dir.exists(newFolder)) {
        oldFi=dir(path=newFolder, full.names=T)
        prevs=file.path(newFolder,"previous")
        if (length(oldFi) > 0) {save.prev.folder(prevs,oldFi)}
        mmm1 <- "Copying Reports to output folder:"
        message(crayon::white$bgBlue$bold(mmm1),"\n",newFolder)
        fi2copy <- dir(getwd(), pattern="*.html", full.names = T)
        message(crayon::white$bgGreen$bold("Now copying html reports..."),"\n")
        fs::file_copy(path=fi2copy,newFolder)
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
    if (length(samsheet) < 1) {message("Redcap headers csv file not found")}
    if (length(samsheet) == 1) {
        data<-read.csv(samsheet, stringsAsFactors=F)
        if(any(duplicated(data$record_id))){
            message("Removing duplicates in redcap dataframe")
            data = data[!duplicated(data$record_id),]
        }
        cat(redcapAPI::importRecords(rcon,data,"normal","ids",returnData = F))
        for (n in 1:nrow(data)) {
            datarecord = jsonlite::toJSON((as.list(data[n,])), auto_unbox=T)
            print(datarecord)
            message("~~",crayon::bgBlue("Record Uploaded:"))
            RCurl::postForm(ur, token = tk, content = 'record', format = 'csv', type = 'flat', data = datarecord, returnFormat='csv')
        }
    } else {message("no redcap file found")}
}

# Copy Output cnv Files if generated
copy.cnv.files <- function(newFolder, runID, runYear=NULL) {
        msgFunName(cpOutLnk, "copy.cnv.files")
    if (is.null(runYear)){runYear=paste0(format(Sys.Date(), "%Y"))}
    cnv_folder <- file.path(newFolder,paste0(runID,"_CNVs/"))
    cnvNames <- dir(path=getwd(), full.names=T, "*_cnv.png")
    if(length(cnvNames)>2){
        message(paste0("Copying PNG to: ", cnv_folder));print(as.data.frame(cnvNames))
        if (!dir.exists(cnv_folder)){dir.create(cnv_folder)}
        if (dir.exists(cnv_folder)){fs::file_copy(cnvNames,cnv_folder)}
    }
}

# Uploads any created cnv png files to redcap database
uploadCnPng <- function() {
        msgFunName(cpOutLnk, "uploadCnPng")
    rcon <- redcapAPI::redcapConnection(apiLink, gb$ApiToken)
    samSh <- dir(path=getwd(), full.names=T, ".xlsm")
    sampleNumb <- gb$getTotalSamples()
    sh_Dat <-as.data.frame(readxl::read_excel(samSh,sheet=3,range="A1:N97", col_types=c("text")))
    sampleNumb=as.integer(sampleNumb)
    sh_Dat = sh_Dat[1:sampleNumb,]
    filnm =paste0(sh_Dat$tm_number,"_", sh_Dat$b_number, "_cnv.png")
    pathNam=file.path(getwd(),filnm)
    rms=paste(c("control_","low_"),collapse ='|')
    pathNam <- stringr::str_replace(pathNam, rms, "")
    pathNam <- stringr::str_replace(pathNam, "//", "/")
    records <- sh_Dat$record_id
    for (idx in 1:length(records)) {
        pth = pathNam[idx]
        recordName = paste0(records[idx])
        message("Importing Record:");cat(recordName, sep = "\n");cat(pth, sep = "\n")
        redcapAPI::importFiles(rcon, pth, recordName, field="methyl_cn", repeat_instance=1)
    }
}


# Imports the xlsm sheet 3 data
importRedcapStart <- function(nfldr){
        msgFunName(cpOutLnk, "importRedcapStart")
    rcon <- redcapAPI::redcapConnection(apiLink, gb$ApiToken)
    uri=paste0(rcon$url); tk=rcon$token
    samSh <- dir(path=getwd(), full.names=T, ".xlsm")
    sampleNumb <- getTotalSamples()
    sh_Dat <-as.data.frame(readxl::read_excel(samSh,sheet=3,range="A1:N97", col_types=c("text")))[,1:13]
    sampleNumb=as.integer(sampleNumb);sh_Dat = sh_Dat[1:sampleNumb,]
    runName <- gb$runID
    pathNam = file.path(nfldr,paste0(runName,"_CNVs"), paste0(sh_Dat$record_id, "_cnv.png"))
    rms=paste(c("control_","low_"),collapse ='|')
    pathNam <- stringr::str_replace(pathNam, rms, "")
    pathNam <- stringr::str_replace(pathNam, "//", "/")
    sh_Dat$cnv_file_path <- pathNam
    if (!is.null(sh_Dat)){
        for(n in 1:nrow(sh_Dat)){
            record=c(sh_Dat[n,])
            datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
            res<-RCurl::postForm(uri,token=tk,content='record',format='json',type='flat',data=datarecord)
            message("Record Uploaded")
            print(res)
        }
    }
}

# Creates QC record and uploads reports to redcap
uploadToRedcap <- function(file.list, deskCSV=T) {
        msgFunName(cpOutLnk, "uploadToRedcap")
    rcon <- redcapAPI::redcapConnection(apiLink, gb$ApiToken)
    message("\nFiles to Import:\n"); print(file.list)
    if(deskCSV==T){importDesktopCsv(rcon)}else{
    recordName <-stringr::str_replace_all(string = paste0(file.list), ".html", "")
    runIDs <- rep(gb$runID, length(recordName))
    records <- basename(recordName)
    for (idx in 1:length(records)) {
        pth = file.list[idx]
        recordName = paste0(records[idx])
        runID = paste0(runIDs[idx])
        message(crayon::white$bgBlue("Importing Record Report:"))
        data = data.frame(record_id = recordName, run_number = runID)
        redcapAPI::importRecords(rcon, data, overwriteBehavior = "normal",returnContent = "ids", returnData = F)
        redcapAPI::importFiles(rcon = rcon,file = pth,record = recordName,field = "classifier_pdf",repeat_instance = 1)
    }
    }
}

# Imports the xlsm sheet 3 data
importSingle <- function(sh_Dat) {
        msgFunName(cpOutLnk, "importSingle")
    nfldr = file.path(stringr::str_split_fixed(gb$clinDrv, " ", 2)[1],"MethylationClassifier")
    rcon <- redcapAPI::redcapConnection(apiLink, gb$ApiToken)
    uri = paste0(rcon$url)
    tk = rcon$token
    filnm = paste0(sh_Dat$record_id, "_cnv.png")
    pathNam = file.path(nfldr, paste0(gb$runID, "_CNVs"), filnm)
    rms = paste(c("control_", "low_"), collapse = '|')
    pathNam <- stringr::str_replace(pathNam, rms, "")
    pathNam <- stringr::str_replace(pathNam, "//", "/")
    sh_Dat$cnv_file_path <- pathNam
    record = c(sh_Dat[1,])
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox = T)
    print(sh_Dat)
    res <-
        RCurl::postForm(
            uri,
            token = tk,
            content = 'record',
            format = 'json',
            type = 'flat',
            data = datarecord
        )
    message(crayon::white$bgBlue("Record Uploaded"))
    invisible(res)
    uploadToRedcap(file.list = paste0(record[1], ".html"), deskCSV = F)
}

# FUN: Copies Reports to Z drive
copy2outFolder <-function(clinDrv = NULL, runID, runYear = NULL) {
        msgFunName(cpOutLnk, "copy2outFolder")
    if (is.null(runYear)) {runYear = paste0(format(Sys.Date(), "%Y"))}
    if (is.null(clinDrv)) {clinDrv <- gb$clinDrv}
    isMC = sjmisc::str_contains(runID, "MGDM") | sjmisc::str_contains(runID, "MC")
    
    researchOutDir = file.path("/Volumes/snudem01labspace/FINAL_PDF_Reports_Brain",runID)
    clinicalOutDir = file.path(clinDrv, "Results", runYear, runID)
    runYear = ifelse(isMC, paste0("20", stringr::str_split_fixed(runID, "-", 2)[1]), runYear)
    newFolder <- ifelse(isMC==T, clinicalOutDir, researchOutDir)
    
    cat(crayon::white$bgCyan("Output Folder is:\n", newFolder))
    if (!dir.exists(newFolder)) {dir.create(newFolder)}
    oldFi = dir(path = newFolder, full.names = T)
    prevs = file.path(newFolder, "previous")
    
    if (length(oldFi) > 0) {save.prev.folder(prevs, oldFi)} # saves any old files
    file.list = dir(path = getwd(), "*.html", full.names = T)
    #newFolder <- base::rep(newFolder,length(file.list))
    
    message("\nCopying Existing Reports to Folder...\n")
    message(crayon::bgBlue(newFolder),"\n")
    message("Files to copy:\n")
    print(file.list)
    #fs::file_copy(file.list, newFolder, overwrite = T)
    #lapply(file.list, function(foo){fs::file_copy(foo, newFolder, overwrite=T)})
    lapply(file.list,
           function(foo) {
               destDir = file.path(newFolder, basename(foo))
               tryCatch(
                   expr = {
                       if(file.exists(destDir)){unlink(destDir)}
                       fs::file_copy(foo, destDir, overwrite = T)
                   },
                   error = function(cond) {
                       message(cond,"\n")
                       message("Trying other file copy method:")
                       if(!file.exists(foo)){print(paste(foo, "does not exist"))}else{
                           cmnd = paste("cp", foo, newFolder)
                           system(cmnd)}
                   })  
           })
    
    if (isMC) {
        cnList <- dir(getwd(), "_cnv.png",recursive = F)
        hasCn <- length(cnList) > 2
        if (hasCn) {
            copy.cnv.files(newFolder, runID, runYear)
            uploadCnPng()
        }
        clinOut = file.path(stringr::str_split_fixed(clinDrv, " ", 2)[1],"MethylationClassifier")
        importRedcapStart(clinOut)
        copy.to.clinical(clinOut, runID, runYear)
    }
    return(file.list)
}
