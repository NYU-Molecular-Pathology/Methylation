#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library("base")
gb <- globalenv(); assign("gb", gb)

if(!require("devtools")){install.packages("devtools", quiet=T)}

# Parameters Input trailing commandline ----------------------------------------------------
token <-      args[1]
runID <-      args[2]
selectRDs <-  args[3]
baseFolder <- args[4]
redcapUp <-   args[5]
runLocal <-   args[6]

osType <- Sys.info()[['sysname']]
if(osType=="Linux"){
  .libPaths("~/molecpathlab/production/Methylation/common_libs")
  if(!require("devtools")){install.packages("devtools", quiet=T)}
  baseFolder <- "~/molecpathlab/production/Methylation"
}

# Check Input Parameters ----------------------------------------------------
CheckInputArg <- function(varValue, gb, defVal = NULL) {
    varStr <- deparse(substitute(varValue))
    if (length(varValue) == 0 | identical(varValue, NULL) | identical(varValue, "NULL")) {
        gb[[varStr]] <- varValue <- defVal
    } else{ varValue <- ifelse(is.na(varValue), NULL, varValue)}
    message(varStr, ": " , ifelse(is.null(varValue), "NULL", varValue))
    return(assign(varStr, varValue, envir = gb))
}

message("\n~~~~~~~~~~~~~~~~~~~~~Parameters input~~~~~~~~~~~~~~~~~~~~~")
CheckInputArg(token, gb)
CheckInputArg(runID, gb)
CheckInputArg(selectRDs, gb)
CheckInputArg(baseFolder, gb)
CheckInputArg(redcapUp, gb, T)
CheckInputArg(runLocal, gb, F)
message("\n^^^^^^^^^^^^^^^^^^^^^Parameters input^^^^^^^^^^^^^^^^^^^^^\n")

# Cancel if no token or runID ----------------------------------------------------
stopifnot(!is.null(token)); stopifnot(!is.null(runID))
if(!is.null(baseFolder) & !identical(baseFolder, "NULL")) {
    message("Checking if custom run directory is valid: ", baseFolder, "\n")
    stopifnot("Directory input does not exist! Try creating it first with mkdir" = dir.exists(baseFolder) == T)
} else{baseFolder <- NULL}

# Source GitHub Scripts ----------------------------------------------------
mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/"
script.list = c("LoadInstallPackages.R", "SetRunParams.R", "MakeSampleSheet.R", "CopyInputs.R", "CopyOutput.R", "pipelineHelper.R")
scripts = paste0(mainHub, script.list)
lapply(scripts, function(i){message("Sourcing: ", i); devtools::source_url(i)})

# Define Parameters ----------------------------------------------------
selectRDs <- AssignArgs(runID, baseFolder, token, selectRDs, redcapUp, gb)

# Execute Functions ----------------------------------------------------
PrepareRun(token, baseFolder, runID, runLocal=runLocal) # If running local and  runLocal = TRUE
StartRun(selectRDs, emailNotify=T, redcapUp=redcapUp) # can change to default false
