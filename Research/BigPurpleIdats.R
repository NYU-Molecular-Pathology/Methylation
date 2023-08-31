#!/usr/bin/env R
library("base"); gb <- globalenv(); assign("gb", gb)
args <- commandArgs(TRUE)
if(!require("devtools")){install.packages("devtools", quiet=T)}

# Input Arguments ------------------------------------------
args[1] -> token
args[2] -> inputSheet

# Default Paths --------------------------------------------
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
userNam <- paste0(Sys.info()[["user"]])

moVol = file.path("/mnt",userNam,"molecular")
rsVol = file.path("/mnt",userNam,"snudem01labspace")

clin.idat = file.path(moVol, "MOLECULAR","iScan")
rsch.idat = file.path(rsVol,"idats")

# REDcap Heading Fields for API pull to create SampleSheet -----
flds = c("record_id","b_number","tm_number","accession_number","block","diagnosis",
         "organ","tissue_comments","run_number", "nyu_mrn")

# Load redcapAPI Package -----
if(suppressWarnings(!require("redcapAPI"))){
    params=list('nutterb/redcapAPI', dependencies=T, upgrade="always", type="source")
    do.call(devtools::install_github,c(params))
}

# Variables for messages
dsh="\n================"
dsh2="================\n"

supM <- function(sobj){return(suppressMessages(suppressWarnings(sobj)))}

#  Copy idats and Worksheets creation
writeFromRedcap <- function(df, samplesheet_ID, bn = NULL) {
    if (is.null(bn)) {bn = file.path(getwd(), df$barcode_and_row_column)}
    message(crayon::bgCyan("~~~Writing from redcap samplesheet.csv using dataframe:"))

    df<- df[!is.na(df[, "barcode_and_row_column"]),]
    df<- df[!is.null(df[, "barcode_and_row_column"]),]
	    print(df)
    samplesheet_csv = data.frame(
        Sample_Name = df[, "record_id"],
        DNA_Number = df[,"b_number"],
        Sentrix_ID = samplesheet_ID[, 1],
        Sentrix_Position = samplesheet_ID[, 2],
        SentrixID_Pos = df[, "barcode_and_row_column"],
        Basename = paste0(bn),
        RunID = df$run_number,
        MP_num = df$accession_number,
        tech = df$primary_tech,
        tech2 = df$second_tech,
        Date = df$arrived
    )
    samplesheet_csv <- samplesheet_csv[!is.na(samplesheet_csv$SentrixID_Pos),]
    print(samplesheet_csv)
    write.csv(samplesheet_csv, file = "samplesheet.csv", quote = F,row.names = F)
}

search.redcap <- function(rd_numbers, ApiToken=NULL) {
    if(is.null(ApiToken)){message("You must provide an ApiToken!")}
    rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", ApiToken)
    flds = c("record_id","b_number","primary_tech","second_tech","run_number","barcode_and_row_column","accession_number","arrived")
    result <- redcapAPI::exportRecords(
	    rcon, records = rd_numbers, fields = flds, dag = F, factors = F, labels = F, dates = F, form_complete_auto = F, format = 'csv')
    result <- as.data.frame(result)
    return(result)
}

get.idats2<-function(csvNam = "samplesheet.csv"){
    rsch.idat <- gb$rsch.idat;clin.idat <- gb$clin.idat
    if(!dir.exists(rsch.idat)){warnMount(rsch.idat)}; if(!dir.exists(clin.idat)){warnMount(clin.idat)}
    stopifnot(dir.exists(rsch.idat)|dir.exists(clin.idat))
    if (file.exists(csvNam)) {
        allFi <- gb$getAllFiles(idatDir = c(rsch.idat, clin.idat), csvNam = csvNam)
        allFi = allFi[file.exists(allFi)]
        if (length(allFi) > 0) {
            message("Files found: "); print(allFi)
            cur.idat <- dir(pattern = "*.idat$")
            bcds <- paste0(basename(allFi))
            if (all(bcds %in% cur.idat)) {message(".idat files already copied")}
            if (!all(bcds %in% cur.idat)) {gb$copyBaseIdats(allFi[!(bcds %in% cur.idat)])}
        } else {message("No .idat files found! Check worksheet and input folder path")}
    } else {message(paste("Cannot find your sheet named:", csvNam))}
}

# FUN: Copies .idat files to your directory and saves samplesheet.csv
get.rd.info <- function(rd_numbers=NULL, token=NULL, sh_name=NULL){
    if (is.null(rd_numbers)){
	    message("No RD-numbers found, Input RD-numbers using get.rd.info(rd_numbers)")
	    return(NULL)
    }
	print(rd_numbers)
    if (is.null(sh_name)) {sh_name = "samplesheet.csv"}
    result <- gb$search.redcap(rd_numbers, token)
    samplesheet_ID = as.data.frame(stringr::str_split_fixed(result[, "barcode_and_row_column"], "_", 2))
    writeFromRedcap(result, samplesheet_ID) # writes API export as minfi dataframe sheet
    gb$get.idats2(csvNam = sh_name)  # copies idat files from return to current directory
    return(result)
}

# FUN: Sets your directory and sources the helper functions
sourceFuns <- function(workingPath = NULL) {
    mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/"
    script.list <- c("R/SetRunParams.R","R/CopyInputs.R","PACT_scripts/generateCNV.R")
    if (is.null(workingPath)) {workingPath = getwd()}
    scripts <- paste0(mainHub, script.list)
    invisible(lapply(scripts, function(i){devtools::source_url(i)}))
    gb$setDirectory(workingPath)
    return(gb$defineParams(loadClassifier=F, clin.idat=clin.idat, rsch.idat=rsch.idat))
}

# FUN: Checks if z-drive is accessible to the Rscript
checkMounts <- function(){
    mntAccess <- dir.exists(clin.idat) && dir.exists(rsch.idat)
    failMount <- ifelse(mntAccess, T, F)
    if(failMount!=T){
        cat("\nidat drive paths are not accessible, ensure both paths are mounted:\n")
        cat("Drive mounted:", dir.exists(clin.idat))
        cat(crayon::white$bgRed$bold(moVol))
        cat("Drive mounted:", dir.exists(rsch.idat))
        cat(crayon::white$bgRed$bold(rsVol))
        stopifnot(!any(failMount==T))
    } else {message("\n",crayon::bgGreen("idat drive paths are accessible"),"\n")}
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
    message(paste0(capture.output(samplesheet_csv), collapse="\n"))
}

grabRDCopyIdat <- function(rd_numbers, token, copyIdats=T, outputFi="samplesheet_og.csv"){
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
            supM(get.idats(csvNam = outputFi))
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

# Search REDCap Worksheets for MRN Match for output -------------------------------------
loadPacks()
checkMounts()
sourceFuns()

# Example Use
if(!is.na(token) & !is.na(inputSheet)){
    rds <- readInfo(inputSheet)
    grabRDCopyIdat(rd_numbers=rds, token=token, copyIdats=T)
}

