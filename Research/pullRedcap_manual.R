#!/usr/bin/env R
library("base"); gb <- globalenv(); assign("gb", gb)
args <- commandArgs(TRUE)
if(!require("devtools")){install.packages("devtools", quiet=T)}

# Input Arguments ------------------------------------------
args[1] -> token
args[2] -> inputSheet

dsh="\n================"
dsh2="================\n"

# Default Paths
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
moVol = "/Volumes/molecular"
rsVol = "/Volumes/snudem01labspace"
rsch.idat = paste0(file.path(rsVol,"idats"))
clin.idat = paste0(file.path(moVol, "MOLECULAR/iScan"))

# REDcap Heading Fields to pull for SampleSheet -----
flds = c("record_id","b_number","tm_number","accession_number","block","diagnosis",
         "organ","tissue_comments","run_number", "nyu_mrn")

# Load redcapAPI Package -----
if(suppressWarnings(!require("redcapAPI"))){
    params=list('nutterb/redcapAPI', dependencies=T, upgrade="always", type="source")
    do.call(devtools::install_github,c(params))
}

supM <- function(sobj){return(suppressMessages(suppressWarnings(sobj)))}

# FUN: Checks if z-drive is accessible to the Rscript
checkMounts <- function(){
    molecDrive = "/Volumes/molecular/MOLECULAR LAB ONLY"
    zDrive = "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular"
    failMount <- ifelse(dir.exists(molecDrive),T,F)
    if(failMount!=T){
        cat("\nPATH does not exist, ensure path is mounted:\n")
        cat(crayon::white$bgRed$bold(molecDrive))
        cat("\nYou must mount the network Z-drive path:\n")
        cat(crayon::white$bgRed$bold(zDrive),"\n")
        stopifnot(!any(failMount==T))
    } else {message("\n",crayon::bgGreen("Z-drive path is accessible"),"\n")}
}

# Functions to load packages and get redcap info -----
loadPacks <- function(){
    pkgs = c("data.table", "foreach", "openxlsx","jsonlite","RCurl","readxl","stringr","tidyverse","crayon")
    rlis = getOption("repos")
    rlis["CRAN"] = "http://cran.us.r-project.org"
    options(repos = rlis)
    invisible(lapply(pkgs, function(pk){
        if(suppressWarnings(!require(pk, character.only=T))){
            install.packages(pk,dependencies=T, verbose=T, repos="http://cran.us.r-project.org", type="both")
        }}))
    if(!require("redcapAPI")){install.packages("redcapAPI", dependencies = T, type="both",ask=F)}
    if(!require("remotes")){install.packages("remotes", dependencies=T)}
    library("redcapAPI")
    library(dplyr)
    require('foreach')
}

# API Call functions -----
grabAllRecords <- function(flds, rcon){
    params = list(rcon, fields=flds, labels=F, dates = F, survey = F, dag = F, factors=F, form_complete_auto=F)
    dbCols <- do.call(redcapAPI::exportRecords, c(params))
    return(as.data.frame(dbCols))
}

search.redcap <- function(rd_numbers, token=NULL, flds=NULL) {
    if(is.null(token)){message("You must provide an ApiToken!")};stopifnot(!is.null(token))
    rcon <- redcapAPI::redcapConnection(gb$apiLink, token)
    if (is.null(flds)){flds = c("record_id","b_number","primary_tech","second_tech","run_number",
                                "barcode_and_row_column","accession_number","arrived")}
    result <- redcapAPI::exportRecords(rcon,records = rd_numbers,fields = flds, dag = F,factors = F,
                                       labels = F,dates = F, form_complete_auto = F,format = 'csv')
    return(as.data.frame(result))
}

# FUN: Sets your directory and sources the helper functions
sourceFuns <- function(workingPath = NULL) {
    mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/"
    script.list <- c("R/SetRunParams.R","R/CopyInputs.R","PACT_scripts/generateCNV.R")
    if (is.null(workingPath)) {workingPath = getwd()}
    scripts <- paste0(mainHub, script.list)
    invisible(lapply(scripts, function(i){devtools::source_url(i)}))
    gb$setDirectory(workingPath)
    return(gb$defineParams(loadClassifier=F))
}

readInfo <- function(inputSheet) {
    # Detect if file is xlsx or csv
    readFlag <- endsWith(inputSheet,".csv")==T
    
    if (readFlag == T) {
        message("FileType is .csv, executing read.delim...")
        rds <- read.delim(inputSheet, sep=",", colClasses=character(), row.names=NULL)[,1]
    } else{
        message("FileType is .xlsx, executing readxl::read_excel...")
        rds <- readxl::read_excel(inputSheet, sheet = 1)[, 1]
    }
    if(typeof(rds)!="character"){
        warning('Converting RD-numbers to "Character" [1], readxl or read.delim should not output typeof=="list" update version of tibble')
        rds <- as.data.frame(rds)[,1]
    }
    return(rds)
}

makeSampleSheet <- function(df, samplesheet_ID, bn = NULL, outputFi="samplesheet_og.csv") {
    if (is.null(bn)) {bn = file.path(getwd(), df$barcode_and_row_column)}
    message(crayon::bgCyan("~~~Writing from redcap samplesheet.csv using dataframe:"))
    # Drop Null/Missing sentrix IDs
    df<- df[!is.na(df[, "barcode_and_row_column"]),]
    df<- df[!is.null(df[, "barcode_and_row_column"]),]
    # Make Data Frame
    samplesheet_csv = data.frame(
        Sample_Name = df[, "record_id"],
        DNA_Number = df[,"b_number"],
        Sentrix_ID = samplesheet_ID[, 1],
        Sentrix_Position = samplesheet_ID[, 2],
        SentrixID_Pos = df[, "barcode_and_row_column"],
        Basename = paste0(bn),
        RunID = df$run_number,
        MP_num = df$accession_number,
        Date = df$arrived
    )
    samplesheet_csv <- samplesheet_csv[!is.na(samplesheet_csv$SentrixID_Pos),]
    write.csv(samplesheet_csv, file = outputFi, quote = F,row.names = F)
}

grabRDCopyIdat <-
function(rd_numbers, token, copyIdats=T, outputFi="samplesheet_og.csv"){
    ApiToken <- token
    result_raw <- gb$search.redcap(rd_numbers, token)

    result <- result_raw[!is.na(result_raw$barcode_and_row_column),]
    samplesheet_ID = as.data.frame(stringr::str_split_fixed(result[,"barcode_and_row_column"],"_",2))
    if(nrow(samplesheet_ID)>0) {
        message(
            "The RD-numbers you entered have not been run yet or do not have idat files in REDCap:\n\n",
            paste(capture.output(result_raw), collapse = "\n")
        )
        stopifnot(nrow(samplesheet_ID)>0)
    }
    gb$makeSampleSheet(result, samplesheet_ID, bn = NULL, outputFi=outputFi) # writes API export as minfi dataframe sheet
    # copies idat files from return to current directory
    if(copyIdats==T){
        supM(gb$get.idats(csvNam = outputFi))
    }
}

fillMissingDat <- function(targets, col_samNames="Sample_Name", originalFi="samplesheet_og.csv"){
  newTarg <- read.csv(originalFi, strip.white=T, row.names=NULL)
  targets <- merge(newTarg,targets, by=col_samNames, all=F, suffixes = c("",".xyzq"))
  dupeDrop <- grepl(".xyzq", colnames(targets))==F
  targets <- targets[,dupeDrop]
  write.csv(targets, file="samplesheet.csv", quote=F, row.names=F)
  targets <- read.csv("samplesheet.csv", strip.white=T, row.names=NULL)
  if(class(targets)!="data.frame"){targets <- as.data.frame(targets)}
  return(targets)
}

fixBaseName <- function(targets, runDir, col_sentrix) {
  if(class(targets)!="data.frame"){targets <- as.data.frame(targets)}
  senCol <- min(which(grepl(col_sentrix, colnames(targets)) == T))
  targets$Basename <- file.path(runDir, targets[, senCol]) # writes path to idat files
  return(targets)
}

# Search REDCap Worksheets for MRN Match for output -------------------------------------
loadPacks()
checkMounts()
sourceFuns()

# Example Use
rds <- readInfo(inputSheet)
grabRDCopyIdat(rd_numbers=rds, token=token)
