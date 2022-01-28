#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library("base"); gb <- globalenv(); assign("gb", gb)
dsh="\n================"
dsh2="================\n"

# Main arguments input in comandline
token <- args[1]
inputSheet <- args[2]

# Default Paths
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
moVol = "/Volumes/molecular"
rsVol = "/Volumes/snudem01labspace"
rsch.idat = paste0(file.path(rsVol,"idats"))
clin.idat = paste0(file.path(moVol, "MOLECULAR/iScan")
                   
# Displays the Input args -----
message(dsh,"Parameters input",dsh2)
message("token: ",token)
message("inputSheet: ", inputSheet,"\n")

stopifnot(!is.na(token))
stopifnot(!is.na(inputSheet))

readFlag <- endsWith(inputSheet,".csv")==T
stopifnot(rlang::is_bool(readFlag))

# REDcap Heading Fields -----
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
}

# API Call functions -----
grabAllRecords <- function(flds, rcon){
    params = list(rcon, fields=flds, labels=F, dates = F, survey = F, dag = F, factors=F, 
                  form_complete_auto=F)
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

grabRDs <- function(rd_numbers, token){
    result <- gb$search.redcap(rd_numbers, token)
    result <- result[!is.na(result$barcode_and_row_column),]
    samplesheet_ID = as.data.frame(stringr::str_split_fixed(result[,"barcode_and_row_column"],"_",2))
    gb$writeFromRedcap(result, samplesheet_ID) # writes API export as minfi dataframe sheet
    gb$get.idats()  # copies idat files from return to current directory
}

readInfo <- function(inputSheet, readFlag) {
    if (readFlag == T) {
        rds <- as.data.frame(read.csv(inputSheet))[, 1]
        
    } else{
        rds <- readxl::read_excel(inputSheet, sheet = 1)[, 1]
    }
    return(rds)
}

# Search REDCap Worksheets for MRN Match for output -------------------------------------
loadPacks()
checkMounts()

rds <- readInfo(inputSheet,readFlag)
message("Your RDs:")
print(rds)
grabRDs(rds, token)
