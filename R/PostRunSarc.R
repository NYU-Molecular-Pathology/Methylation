#!/usr/bin/env Rscript
## ---------------------------
## Script name: PostRunSarc.R
## Purpose: Source global scripts for executing Clinical sarcoma pipeline
## Date Last Modified: January 12, 2024
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------
gb <- globalenv(); assign("gb", gb)
if (!requireNamespace("devtools")) {install.packages("devtools", quiet = T)}

args <- commandArgs(TRUE)
library("base"); gb <- globalenv(); assign("gb", gb)

# Main arguments input in comandline
token <- args[1]
runID <- args[2]

PullNeedsSarcoma <- function(rd_numbers, token){
    message("Pulling REDCap data...")
    apiUrl = "https://redcap.nyumc.org/apps/redcap/api/"
    pull_flds = c("record_id", "barcode_and_row_column", "needs_sarcoma", "classifier_pdf_other")
    rcon <- redcapAPI::redcapConnection(apiUrl, token)
    params = list(rcon, records = rd_numbers, fields = pull_flds, survey = F, dag = F, factors=F, form_complete_auto=F)
    dbCols <- do.call(redcapAPI::exportRecordsTyped, c(params))
    return(as.data.frame(dbCols))
}


StartSarcWorkflow <- function(rd_numbers, token, runID){
    runLocal <- T
    selectRDs <- baseFolder <- NULL

    ghLink <- "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R"
    devtools::source_url(file.path(ghLink, "Report-Scripts/SourceLoadGitHub.R"))

    targets = as.data.frame(read.csv("samplesheet.csv"))
    rownames(targets) <- targets[,1]
    
    gb$MakeSarcomaReport(worksheet = NULL, targets = targets[rd_numbers,])

    file.list <- dir(getwd(), pattern = "_sarc.html", full.names = T)
    if(length(file.list) > 0){
        gb$ForceUploadToRedcap(file.list, token, F)
    }
}


CheckNeedsSarcoma <- function(rd_numbers, token, runID){
    pulled_rds <- PullNeedsSarcoma(rd_numbers, token)
    if(nrow(pulled_rds) > 0){
        valid_sams <-
            pulled_rds$needs_sarcoma == "yes" &
            !is.na(pulled_rds$barcode_and_row_column) &
            is.na(pulled_rds$classifier_pdf_other)
        valid_sams[is.na(valid_sams)] <- FALSE
    }else{
        valid_sams <- FALSE
    }
    if(any(valid_sams)){
        rd_to_run <- pulled_rds$record_id[valid_sams]
        message(crayon::bgGreen("Generating Sarcoma Reports for:"), "\n", paste(rd_to_run, sep = "\n"))
        StartSarcWorkflow(rd_to_run, token, runID)
    }
}

runPath <- file.path("/Volumes/CBioinformatics/Methylation/Clinical_Runs", runID)
setwd(runPath)

data <- as.data.frame(read.csv("samplesheet.csv"))
CheckNeedsSarcoma(rd_numbers = data[-1,1], token, runID)
