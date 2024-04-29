#!/usr/bin/env Rscript
## ---------------------------
## Script name: SourceLoadGitHub.R
## Purpose: Function or importing scripts from Github for executing the Methylation Pipeline
## Date Created: March 18, 2023
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

gb <- globalenv(); assign("gb", gb)

LoadGitHubScripts <- function(ghRepo, scriptList){scripts = file.path(ghRepo, scriptList)
  return(lapply(scripts, function(i){message("Sourcing: ", i); devtools::source_url(i)}))}

mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R"

scriptList <- c("LoadInstallPackages.R", "SetRunParams.R", "MakeSampleSheet.R",
                "CopyInputs.R", "CopyOutput.R", "pipelineHelper.R", "CustomRuns.R")

rmdScripts <- c("ClassTables.R", "MLH1_Functions.R", "PipeLineU.R", "RedcapOutput.R", "TsneFunctions.R", "cnvggplotly.R")

LoadGitHubScripts(mainHub, scriptList); LoadGitHubScripts(file.path(mainHub,"Report-Scripts"), rmdScripts)

gb$reportMd <- reportMd <- file.path(fs::path_home(),"report.Rmd")

clear_sarc_dir <- function(outputFi) {
  if (file.exists(file.path(getwd(), outputFi))) {
      file.remove(file.path(getwd(), outputFi))
  }
  
  file.list <- dir(getwd(), pattern = ".html|.idat", full.names = T)
  
  if (length(file.list) > 0) {
      file.remove(file.list)
  }
}

if(!file.exists(reportMd)){
    if(!require("curl")){try(install.packages("curl", dependencies = T, ask = F), silent=T)}
    url = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/report.Rmd"
    curl::curl_download(url, reportMd, quiet = TRUE, mode = "wb", handle = curl::new_handle())
}
