#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library("base")
#token=NULL;inputSheet=NULL
dsh="\n================"

# Main arguments input in comandline
token<-args[1]; inputSheet<-args[2]

# Displays the Input args -----
message(dsh,"\nParameters input",dsh)
message("token: ",token)
message("Worksheet: ", inputSheet,".xlsx\n")

# Check Input Params -----
stopifnot(exists("token")|!is.null(token))
stopifnot(exists("inputSheet")|!is.null(inputSheet))

# FUN: Checks if z-drive is accessible to the Rscript
checkMounts <- function(){
    molecDrive = "/Volumes/molecular/MOLECULAR LAB ONLY"
    failMount <- ifelse(dir.exists(molecDrive),T,F)
    if(failMount!=T){
        cat("PATH does not exist, ensure path is mounted:")
        cat(crayon::white$bgRed$bold(molecDrive))
        cat("You must mount the network Z-drive path:")
        cat("smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular\n")
        stopifnot(!any(failMount==T))
    } else {message("Z-drive path is accessible")}
}

# Function to silently load library without conflict warnings
libLoad <- function(libName) {
    lib.opts <- list(package = libName, character.only = T, verbose = F, warn.conflicts = F)
    suppressWarnings(suppressPackageStartupMessages(do.call(library, c(lib.opts))))
    message(libName, " ...loaded successful")
}

# Functions to load or install missing required packages -----
loadPacks <- function(){
    if(suppressWarnings(!require("redcapAPI"))){
        params=list('nutterb/redcapAPI', dependencies=T, upgrade="always", type="source")
        do.call(devtools::install_github,c(params))
    } else{libLoad("redcapAPI")}
    pkgs = c("foreach","jsonlite","RCurl","readxl","stringr")
    rlis = getOption("repos");rlis["CRAN"] = "http://cran.us.r-project.org"
    options(repos = rlis)
    invisible(lapply(pkgs, function(pk){
        if(suppressWarnings(!require(pk, character.only=T))){
            install.packages(pk,dependencies=T, verbose=T, repos="http://cran.us.r-project.org", type="both")
            libLoad(pk)
        }else{libLoad(pk)}}))
}

# Returns Path to xlsx file -----
getExcelPath <- function(inputSheet){
    drive = file.path("", "Volumes", "molecular", "MOLECULAR LAB ONLY")
    folder = file.path("NYU PACT Patient Data", "Workbook")
    runyr <- stringr::str_split_fixed(inputSheet, "-", 3)[, 2]
    return(file.path(drive, folder, paste0("20", runyr),inputSheet,paste0(inputSheet, ".xlsx")))
}

# Parses xlsx file and writes as csv file -----
parseExcelFile <- function(inputFi){
    shNames <- readxl::excel_sheets(inputFi)
    sh <- which(grepl("SampleSheet", shNames, ignore.case = T))[1]
    runId <- paste0(head(readxl::read_excel(inputFi, sheet = shNames[sh]))[3, 2])
    sheetHead <- as.data.frame(readxl::read_excel(inputFi,sheet = shNames[sh], na="", range="A1:B17", col_types = "text", col_names=F))
    sheetHead[is.na(sheetHead)] <- ""
    sheetHead <- rbind(sheetHead,c("",""),c("[Data]",""))
    sheetVals <- as.data.frame(readxl::read_excel(inputFi,sheet = shNames[sh],skip = 19,col_types = "text"))
    mainSheet <- sheetVals[!is.na(sheetVals[,1]),]
    for(i in 1:ncol(mainSheet)){
    mainSheet[,i] <- sapply(mainSheet[,i], function(x) { gsub("[\r\n]", "", x) })
    }
    batchID <- mainSheet[1,"Run_Number"]
    outFile <- file.path("~","Desktop",paste(batchID,"SampleSheet.csv",sep="-"))
    write.table(sheetHead,sep=",", file=outFile, row.names=F, col.names=F)
    suppressWarnings(write.table(mainSheet,sep=",", file=outFile, row.names=F, col.names=T, append=T))
    return(c(runId=runId, outFile=outFile))
}

# Generate Email notification and attach csv file
emailNotify <- function(record,rcon){
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
    res<-RCurl::postForm(rcon$url, token=rcon$token, content='record',format='json',type='flat',
                         data = datarecord, returnContent = 'ids', returnFormat = 'csv')
    cat(res);message(dsh,"\nEmail Notification Created",dsh)
}

# Connect to REDCap and send email attachments of csv file ----
pushToRedcap <- function(runId,outFile,token){
    rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", token)
    record = data.frame(record_id = runId, pact_run_number = runId)
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
    res <- RCurl::postForm(rcon$url,token=rcon$token,content='record',format='json',
                           type='flat',data=datarecord, returnContent = 'count', returnFormat = 'csv')
    cat(res,sep="\n")
    redcapAPI::importFiles(rcon=rcon, file=outFile, record=runId, field="pact_csv_sheet", repeat_instance=1)
    record$pact_csv_email<-"pact_csv_email"
    emailNotify(record, rcon)
    #unlink(outFile)
}

# Gets dataframe and saves as CSV file -----
writeSampleSheet <- function(inputSheet, token){
    inputFi <- getExcelPath(inputSheet)
    if (file.exists(inputFi)) {
        outputVals <- suppressMessages(parseExcelFile(inputFi))
        pushToRedcap(runId=outputVals[[1]], outFile=outputVals[[2]], token)
    } else {
        message(crayon::bgRed("The PACT run worksheet was not found:"),"\n", inputFi, dsh)
        stopifnot(file.exists(inputFi))
        }
}

loadPacks()
checkMounts()
writeSampleSheet(inputSheet, token)
