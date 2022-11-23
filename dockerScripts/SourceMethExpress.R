LoadGitHubScripts <- function(ghRepo, scriptList){
  scripts = file.path(ghRepo, scriptList)
  return(lapply(scripts, function(i){message("Sourcing: ", i); devtools::source_url(i)}))
}

# Source GitHub Scripts -------------------------------------------------------------------------
mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R"
scriptList <-
  c(
    "LoadInstallPackages.R",
    "SetRunParams.R",
    "MakeSampleSheet.R",
    "CopyInputs.R",
    "CopyOutput.R",
    "pipelineHelper.R",
    "CustomRuns.R"
  )

# Source Additional Scripts for Rmd Knit ---------------------------------------------------------
rmdScripts <- c("ClassTables.R", "MLH1_Functions.R", "PipeLineU.R",
                "RedcapOutput.R", "TsneFunctions.R", "cnvggplotly.R")

LoadGitHubScripts(mainHub, scriptList)
LoadGitHubScripts(file.path(mainHub,"Report-Scripts"), rmdScripts)
