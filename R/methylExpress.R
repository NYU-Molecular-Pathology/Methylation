#!/usr/bin/env Rscript
library("base")
gb <- globalenv(); assign("gb", gb)
args <- commandArgs(TRUE)

# Parameters Input trailing commandline -------------------------------------------------------
args[1] -> token      # <- NULL
args[2] -> runID      # <- 21-MGDM_TEST
args[3] -> selectRDs  # <- NULL
args[4] -> baseFolder # <- NULL
args[5] -> redcapUp   # <- TRUE
args[6] -> runLocal   # <- FALSE

if(!require("devtools")){install.packages("devtools", quiet=T)}

# Check Input Parameters -----------------------------------------------------------------------
CheckInputArg <- function(varValue, gb, defVal = NULL) {varStr <- deparse(substitute(varValue))
    if (length(varValue) == 0 | identical(varValue, NULL) | identical(varValue, "NULL")) {gb[[varStr]] <- varValue <- defVal} else{ varValue <- ifelse(is.na(varValue), NULL, varValue)}
    if(varStr=="token"){message("\n~~~~~~~~~~~~~~~~~~~~~Parameters input~~~~~~~~~~~~~~~~~~~~~")}
    message(varStr, ": " , ifelse(is.null(varValue), "NULL", varValue)); return(assign(varStr, varValue, envir = gb))
}

# Message and Check Input Args ------------------------------------------------------------------
CheckInputArg(token, gb); CheckInputArg(runID, gb); CheckInputArg(selectRDs, gb); CheckInputArg(baseFolder, gb)
CheckInputArg(redcapUp, gb, T); CheckInputArg(runLocal, gb, F)

stopifnot(!is.null(token)); stopifnot(!is.null(runID))
gb$ApiToken <- gb$token <- token; assign("token", token, envir = gb); assign("ApiToken", token, envir = gb)

LoadGitHubScripts <- function(ghRepo, scriptList){scripts <- file.path(ghRepo, scriptList)
    return(invisible(lapply(scripts, function(i){message("Sourcing: ", i); devtools::source_url(i)})))
}

# Source GitHub Scripts --------------------------------------------------------------------------
mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R"
scriptList <- c("LoadInstallPackages.R", "SetRunParams.R", "MakeSampleSheet.R", "CopyInputs.R", "CopyOutput.R", "pipelineHelper.R", "CustomRuns.R")
rmdScripts <- c("ClassTables.R", "MLH1_Functions.R", "PipeLineU.R", "RedcapOutput.R", "TsneFunctions.R", "cnvggplotly.R")
LoadGitHubScripts(mainHub, scriptList)
LoadGitHubScripts(file.path(mainHub,"Report-Scripts"), rmdScripts)

# Assign Parameters if Defined -------------------------------------------------------------------
baseFolder <- gb$CheckBaseFolderInput(baseFolder)
selectRDs <- gb$AssignArgs(runID, baseFolder, token, selectRDs, redcapUp, gb)
gb$reportMd <- reportMd <- "~/report.Rmd"

# Execute Pipeline Functions ----------------------------------------------------------------------
gb$PrepareRun(token, baseFolder, runID, runLocal=runLocal) # If running local set runLocal = TRUE
selectRDs <- gb$GetPriorityCases(selectRDs) # Prioritizes select RD-numbers and BN cases
gb$StartRun(selectRDs, emailNotify=T, redcapUp=redcapUp) # Can be changed to default false
# gb$MakeSarcomaReport()
