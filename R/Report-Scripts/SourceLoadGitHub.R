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

rmdScripts <- c("ClassTables.R", "MLH1_Functions.R", "RedcapOutput.R")

LoadGitHubScripts(mainHub, scriptList)
LoadGitHubScripts(file.path(mainHub,"Report-Scripts"), rmdScripts)

gb$reportMd <- reportMd <- "/Volumes/CBioinformatics/Methylation/EPIC_V2_report_2.Rmd"

clear_sarc_dir <- function(outputFi) {
  if (file.exists(file.path(getwd(), outputFi))) {
      file.remove(file.path(getwd(), outputFi))
  }
  file.list <- dir(getwd(), pattern = ".html|.idat", full.names = T)
  if (length(file.list) > 0) {
    response <- readline(prompt = "Do you want to remove previous idat and html files in the current directory? (yes/no): ")
    if (tolower(response) == "yes") {
      message("Deleting existing .idat and .html files...")
      file.remove(file.list)
    }
  }
}


SetKnitProgress <- function() {
    library("knitr")
    library("progressr")
    library("cli")
    knitr::opts_knit$set(verbose = TRUE)
    options(knitr.package.verbose = FALSE)
    progressr::handlers(global = TRUE)
    progressr::handlers("cli")
    options(rmarkdown.render.message = FALSE)
    options(
        knitr.progress.fun = function(total, labels) {
            p <- progressr::progressor(
                steps = total,
                along = labels,
                trace = TRUE,
                auto_finish = FALSE,
                on_exit = FALSE
            )
            list(
                update = function(i) {
                    p(message = sprintf('Current Chunk: %s......', labels[i]),
                      class = "sticky")
                    newLab <-
                        ifelse(labels[i] == "", "In-line Code", labels[i])
                    message("Chunk #", i, ": ", newLab, " - ",
                            round((i / total) * 100), "% complete")
                },
                done = function() p(type = 'finish')
            )
        }
    )
}
