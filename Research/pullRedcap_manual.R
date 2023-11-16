#!/usr/bin/env R
## ---------------------------
## Script name: pullRedcap_manual.R
## Purpose: source of global scripts imported for copying idat files to cwd and create minfi samplesheet.csv
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------

library("base"); gb <- globalenv(); assign("gb", gb)
formals(library)$quietly <- T
formals(library)$warn.conflicts <- F
formals(require)$warn.conflicts <- F
formals(install.packages)$dependencies <- T
formals(install.packages)$verbose <- T
formals(install.packages)$ask <- F

args <- commandArgs(T)
if(!requireNamespace("devtools")){install.packages("devtools", quiet=T)}

# Input Arguments ------------------------------------------
try(args[1] -> token)
try(args[2] -> inputSheet)
try(args[3] -> copyToFolder)

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
if(suppressWarnings(!requireNamespace("redcapAPI"))){
    params = list('nutterb/redcapAPI', dependencies = T, upgrade = "always", type = "source")
    do.call(devtools::install_github,c(params))
}

if(!(utils::packageVersion("redcapAPI") >= "2.7.4")){
    install.packages("redcapAPI", ask=F, update=T, dependencies=T)
}

are_valid <- function(...) {
    args <- list(...)
    return(all(sapply(args, function(x) !is.null(x) && !is.na(x) && length(x) > 0 && nzchar(x))))
}

if(!are_valid(copyToFolder)){
    copyToFolder <- getwd()
}

setwd(copyToFolder)

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
    if (Sys.info()[['sysname']]=="Darwin") {
        bothType <- "both"
    }else{
        bothType <- "source"
    }
    invisible(lapply(pkgs, function(pk){
        if(suppressWarnings(!requireNamespace(pk))){
            install.packages(pk,dependencies=T, verbose=T, repos="http://cran.us.r-project.org", type=bothType)
        }}))
    if(!requireNamespace("redcapAPI")){install.packages("redcapAPI", dependencies = T, type=bothType,ask=F)}
    if(!requireNamespace("remotes")){install.packages("remotes", dependencies=T)}
    library("redcapAPI")
    library(dplyr)
    requireNamespace('foreach')
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
    result <- redcapAPI::exportRecordsTyped(rcon,records = rd_numbers,fields = flds, dag = F,factors = F,
                                       labels = F,dates = F, form_complete_auto = F,format = 'csv')
    allFound <- rd_numbers %in% result$record_id
    missing <- which(allFound == FALSE)

    if(length(missing) > 0){
        message("Some RD-numbers were not found in REDCap!")
        message(paste0(capture.output(rd_numbers[missing]), collapse="\n"))
    }

    return(as.data.frame(result))
}


# Sets default variable paths/names
getDefaults <- function() {
    cbVol <- switch (Sys.info()[['sysname']],
                     "Darwin" = "/Volumes/CBioinformatics/Methylation",
                     "Linux" = "~/molecpathlab/production/Methylation"
    )
    moVol = "/Volumes/molecular"
    rsVol = "/Volumes/snudem01labspace"
    defaultParams <- data.frame(
        mnp.pk.loc = paste0(file.path(cbVol, "classifiers/mnp.v11b6")),
        ApiToken = "",
        methDir = paste0(file.path(cbVol, "Clinical_Runs")),
        clinDrv = paste0(file.path(moVol, "MOLECULAR LAB ONLY/NYU-METHYLATION")),
        rschOut = paste0(file.path(rsVol, "FINAL_PDF_Reports_Brain")),
        clinOut = paste0(file.path(moVol, "MOLECULAR/MethylationClassifier")),
        rsch.idat = paste0(file.path(rsVol,"idats")),
        clin.idat = paste0(file.path(moVol, "MOLECULAR/iScan")),
        QC_file = paste0("~/Methyl_QC.Rmd"),
        baseDir = paste0(file.path(cbVol, "Clinical_Runs")),
        stringsAsFactors=F
    )
    return(defaultParams)
}

setVar <- function(valueName,val){
    return(assign(valueName, val, envir=.GlobalEnv))
}

assignVar <- function(varStr, assignedVal){
    return(
        tryCatch(
            expr = {if(!is.null(base::get(varStr))){cat("")}},
            error = {setVar(varStr, assignedVal)},
            finally = message(varStr," = ",assignedVal)
        )
    )
}

# returns assigned global variable values
getSetvars <- function() {
    assignedVars <- data.frame(
        mnp.pk.loc = gb$mnp.pk.loc,
        ApiToken = gb$ApiToken,
        methDir = gb$methDir,
        clinDrv = gb$clinDrv,
        rschOut = gb$rschOut,
        clinOut = gb$clinOut,
        rsch.idat = gb$rsch.idat,
        clin.idat = gb$clin.idat,
        QC_file = gb$QC_file,
        baseDir = gb$baseDir,
        stringsAsFactors = F
    )
    return(assignedVars)
}

# Sets the default parameters for a methylation Run
defineParams <- function(mnp.pk.loc = NULL,
                         ApiToken = NULL,
                         methDir = NULL,
                         clinDrv = NULL,
                         rschOut = NULL,
                         clinOut = NULL,
                         rsch.idat = NULL,
                         clin.idat = NULL,
                         QC_file = NULL,
                         isMC = T,
                         baseDir = NULL,
                         runID = NULL,
                         loadClassifier = T) {
    defVars <- getDefaults()
    inVars <- list(
        mnp.pk.loc,
        ApiToken,
        methDir,
        clinDrv,
        rschOut,
        clinOut,
        rsch.idat,
        clin.idat,
        QC_file,
        baseDir
    )
    i=1:length(inVars)
    for (x in 1:length(defVars)) {
        if(!is.null(inVars[[x]])){
            defVars[1,x] <- inVars[[x]]
        }
    }
    invisible(lapply(i,function(x){if(!is.null(inVars[[x]])){setVar(names(defVars[x]), inVars[[x]])}}))
    i=1:length(defVars)
    invisible(sapply(i,FUN=function(i){assignVar(names((defVars[i])), paste0(defVars[,i]))}))
    cbVol <- switch (Sys.info()[['sysname']],
                     "Darwin" = "/Volumes/CBioinformatics/",
                     "Linux" = "~/molecpathlab/production/")
    if (!isMC) {
        methDir = rschOut
        assign("workDir", cbVol)
    }
}



setDirectory <- function(foldr) {
    bsDir = paste("cd", foldr)
    if (dir.exists(foldr)) {
        system(bsDir)
        setwd(foldr)
        assign("workDir", foldr)
    } else{
        warning(paste("Location Not Found:", foldr))
    }
}


# FUN: Sets your directory and sources the helper functions
sourceFuns <- function(workingPath = NULL) {
    mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/"
    script.list <- c("R/CopyInputs.R")
    if (is.null(workingPath)) {workingPath = getwd()}
    scripts <- paste0(mainHub, script.list)
    invisible(lapply(scripts, function(i){devtools::source_url(i)}))
    setDirectory(workingPath)
    return(defineParams(loadClassifier=F))
}

readInfo <- function(inputSheet) {
    # Detect if file is xlsx or csv
    readFlag <- endsWith(inputSheet,".csv")==T
    if (readFlag == T) {
        message("FileType is .csv, executing read.delim...")
        rds <- read.delim(inputSheet, header = F, sep=",", colClasses=character(), row.names=NULL)[,1]
    } else{
        message("FileType is .xlsx, executing readxl::read_excel...")
        rds <- readxl::read_excel(inputSheet, col_names = F, sheet = 1)[, 1]
    }
    if ( typeof(rds) != "character" ) {
        warning('Converting RD-numbers to class "character", your readxl should not output typeof == "list". Update your version of the `tibble` package in the future!')
        rds <- as.data.frame(rds)[,1]
    }
    rds <- rds[!is.na(rds)]
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
    toDrop <- stringr::str_detect(samplesheet_csv$SentrixID_Pos, "DUPLICATE")
    if(any(toDrop)==T){
        message("Dropping duplicated samples!!")
        otherCsv <- samplesheet_csv[toDrop,]
        message(paste0(capture.output(otherCsv), collapse="\n"))
        write.csv(otherCsv, file = "duplicated_samples.csv", quote = F,row.names = F)
    }
    toKeep <- stringr::str_detect(samplesheet_csv$SentrixID_Pos, "DUPLICATE", T)
    samplesheet_csv <- samplesheet_csv[toKeep,]
    write.csv(samplesheet_csv, file = outputFi, quote = F,row.names = F)
}

grabRDCopyIdat <- function(rd_numbers, token, copyIdats=T, outputFi="samplesheet_og.csv", idatPath=NULL){
    if(is.null(idatPath)){idatPath<- file.path(getwd(),"idats")}
    ApiToken <- token
    result_raw <- gb$search.redcap(rd_numbers, token)
    result <- result_raw[!is.na(result_raw$barcode_and_row_column),]
    samplesheet_ID = as.data.frame(stringr::str_split_fixed(result[,"barcode_and_row_column"],"_",2))
    if(nrow(samplesheet_ID)==0) {
        message(
            "The RD-numbers you entered have not been run yet or do not have idat files in REDCap:\n\n",
            paste(capture.output(result_raw), collapse = "\n")
        )
        stopifnot(nrow(samplesheet_ID)>0)
    }
    gb$makeSampleSheet(result, samplesheet_ID, bn = NULL, outputFi=outputFi) # writes API export as minfi dataframe sheet
    # copies idat files from return to current directory
    if(copyIdats==T){
        supM(gb$get.idats(csvNam = outputFi, runDir=idatPath))
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
    stopifnot(col_sentrix %in% colnames(targets))
    senCol <- min(which(grepl(col_sentrix, colnames(targets)) == T))
    targets$Basename <- file.path(runDir, targets[, senCol]) # writes path to idat files
    return(targets)
}

FixBaseName <- function(targets, runDir=NULL, col_sentrix) {
    if(is.null(runDir)){runDir <- file.path(getwd(), "idats")}
    if(class(targets)!="data.frame"){targets <- as.data.frame(targets)}
    stopifnot(col_sentrix %in% colnames(targets))
    senCol <- min(which(grepl(col_sentrix, colnames(targets)) == T))
    targets$Basename <- file.path(runDir, targets[, senCol]) # writes path to idat files
    return(targets)
}

# Search REDCap Worksheets for MRN Match for output -------------------------------------
loadPacks()
if (Sys.info()[['sysname']]=="Darwin") {checkMounts()}
sourceFuns()

if(are_valid(inputSheet, token)) {
    rd_numbers <- readInfo(inputSheet)
    grabRDCopyIdat(rd_numbers = rd_numbers, token = token)
}

# Example Use ---------------------------------------------------------------------------
#rds <- readInfo(inputSheet = "~/Desktop/MySampleSheet.xlsx")
# OR
#rds <- c("RD-22-123", "RD-21-345", "RD-20-678")
# THEN run
#grabRDCopyIdat(rd_numbers=rds, token=token)
