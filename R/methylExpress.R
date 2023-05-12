#!/usr/bin/env Rscript
## ---------------------------
## Script name: methylExpress.R
## Purpose: Source global scripts for executing Clinical methylation pipeline
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------

library("base"); gb <- globalenv(); assign("gb", gb); args <- commandArgs(TRUE)
if(!require("devtools")){install.packages("devtools", quiet=T)}

# Parameters Input trailing commandline -------------------------------------------------------
args[1] -> token      # <- NULL
args[2] -> runID      # <- "21-MGDM_TEST"
args[3] -> selectRDs  # <- NULL
args[4] -> baseFolder # <- NULL
args[5] -> redcapUp   # <- TRUE
args[6] -> runLocal   # <- FALSE
#args[7] -> forcedUpload
forcedUpload <- F
# Source and Load Functions and Packages --------------------------------------------------------------------------
LoadGitHubScripts <- function(ghRepo, scriptList){scripts <- file.path(ghRepo, scriptList)
    return(invisible(lapply(scripts, function(i){message("Sourcing: ", i); invisible(devtools::source_url(i))})))
}

# Source GitHub Scripts --------------------------------------------------------------------------
mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R"
scriptList <- c("LoadInstallPackages.R", "SetRunParams.R", "MakeSampleSheet.R", "CopyInputs.R", "CopyOutput.R", "pipelineHelper.R", "CustomRuns.R")
rmdScripts <- c("ClassTables.R", "MLH1_Functions.R", "PipeLineU.R", "RedcapOutput.R", "TsneFunctions.R", "cnvggplotly.R")
LoadGitHubScripts(mainHub, scriptList)
invisible(suppressWarnings(LoadGitHubScripts(file.path(mainHub,"Report-Scripts"), rmdScripts)))

# Message and Check Input Args ------------------------------------------------------------------
gb$CheckInputArg(token, gb); gb$CheckInputArg(runID, gb); gb$CheckInputArg(selectRDs, gb); gb$CheckInputArg(baseFolder, gb)
gb$CheckInputArg(redcapUp, gb, T); gb$CheckInputArg(runLocal, gb, F)
#gb$CheckInputArg(forcedUpload, gb, F)

# Assign Parameters if Defined -------------------------------------------------------------------
gb$ApiToken <- gb$token <- token
baseFolder <- gb$CheckBaseFolderInput(baseFolder)
selectRDs <- gb$AssignArgs(runID, baseFolder, token, selectRDs, redcapUp, gb)
gb$reportMd <- reportMd <- file.path(fs::path_home(), "report.Rmd")

# Execute Pipeline Functions ----------------------------------------------------------------------
gb$PrepareRun(token, baseFolder, runID, runLocal=runLocal) # If running local set runLocal = TRUE
gb$CheckIdatsCopied()
selectRDs <- gb$GetPriorityCases(selectRDs) # Prioritizes select RD-numbers and BN cases
gb$StartRun(selectRDs, emailNotify=T, redcapUp=redcapUp) # Can be changed to default false
# gb$MakeSarcomaReport()

if(forcedUpload == T){
    file.list <- dir(getwd(), pattern = ".html", full.names = T)
    print(as.data.frame(file.list))
    gb$ForceUploadToRedcap(file.list, token, F)
}
