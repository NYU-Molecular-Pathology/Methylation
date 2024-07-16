#!/usr/bin/env Rscript
## Script name: methylExpress_sarcoma.R
## Purpose: source global scripts to generate sarcoma reports for EPICv1 samples
## Date Created: December 19, 2023
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024

gb <- globalenv(); assign("gb", gb)

# Default Variables --------------------------------------------------
gb$token <- gb$ApiToken <- token <- "XXXXXXXXXXXX_API_TOKEN_XXXXXXXXXX"
selectRDs <- baseFolder <- NULL

if (!requireNamespace("devtools", quietly = T)) {
    install.packages("devtools", dependencies = T, ask = F, type = "binary")
}

args <- base::commandArgs(trailingOnly = T)

# Check if at least one argument is provided
get_rd_args <- function(args) {
    rd_numbers <- NULL
    if (length(args) == 0) {
        stop("No RD-numbers provided!")
    } else {
        rd_numbers <- as.vector(args)
        message("RD-number(s) input:")
        message(paste0(capture.output(rd_numbers), collapse = "\n"))
    }
    return(rd_numbers)
}

readFlag <- endsWith(args[1], ".csv")

if (readFlag == FALSE) {
    rd_numbers <- get_rd_args(args)
    rd_csv <- NULL
} else {
    rd_csv <- args[1]
    rd_numbers <- NULL
}

# rd_numbers = c("RD-23-1382")  # NULL # change to any individual numbers
# rd_csv = NULL     # a list of RD-numbers to run in .csv column 1

if (token == "XXXXXXXXXXXX_API_TOKEN_XXXXXXXXXX") {
    stop("A valid REDCap API token for 'All Samples Database' is required")
}

runID <- "MR23-rerun_sarc"
getSarcHtml <- forcedUpload <- makeNewSheet <- runLocal <- T
sarc_dir <- "/Volumes/CBioinformatics/Methylation/Clinical_Runs/Sarcoma_runs"
system(paste("cd", sarc_dir))
setwd(sarc_dir)
outputFi = "samplesheet_sarc.csv"
ghLink <-
    "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R"

# Source All Rscripts  -------------------------------------------------------
devtools::source_url(file.path(ghLink, "Report-Scripts/SourceLoadGitHub.R"))

if (makeNewSheet == T) {
    if (file.exists(outputFi)) file.remove(outputFi)
    file.list <- dir(getwd(), pattern = ".html|.idat", full.names = T)
    if (length(file.list) > 0) file.remove(file.list)
    gb$MakeLocalSampleSheet(runID, token, samSheetIn = rd_csv, rd_numbers, outputFi)
}

if (getSarcHtml == T) {
    gb$MakeSarcomaReport(targets = as.data.frame(read.csv(outputFi)))
} else {
    gb$StartCustomRun(redcapUp = F)
}

if (forcedUpload == T) {
    file.list <- dir(getwd(), pattern = ".html", full.names = T)
    gb$ForceUploadToRedcap(file.list, token, F)
}
