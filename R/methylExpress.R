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

# Check BigPurple Enironment & Devtools -------------------------------------------------------
# if(Sys.info()[['sysname']]=="Linux"){
#   .libPaths("~/molecpathlab/production/Methylation/common_libs")
#   assign("baseFolder", "~/molecpathlab/production/Methylation", envir = gb)
# }
if(!require("devtools")){install.packages("devtools", quiet=T)}

# Check Input Parameters -----------------------------------------------------------------------
CheckInputArg <- function(varValue, gb, defVal = NULL) {
    varStr <- deparse(substitute(varValue))
    if (length(varValue) == 0 | identical(varValue, NULL) | identical(varValue, "NULL")) {
        gb[[varStr]] <- varValue <- defVal
    } else{ varValue <- ifelse(is.na(varValue), NULL, varValue)}
    message(varStr, ": " , ifelse(is.null(varValue), "NULL", varValue))
    return(assign(varStr, varValue, envir = gb))
}

# Message and Check Input Args ------------------------------------------------------------------
message("\n~~~~~~~~~~~~~~~~~~~~~Parameters input~~~~~~~~~~~~~~~~~~~~~")
CheckInputArg(token, gb); CheckInputArg(runID, gb)
CheckInputArg(selectRDs, gb); CheckInputArg(baseFolder, gb)
CheckInputArg(redcapUp, gb, T); CheckInputArg(runLocal, gb, F)

# Cancel if no token or runID -------------------------------------------------------------------
stopifnot(!is.null(token)); stopifnot(!is.null(runID))
if(!is.null(baseFolder) & !identical(baseFolder, "NULL")) {
    message("Checking if custom run directory is valid: ", baseFolder, "\n")
    stopifnot("Input directory does not exist! Create it with mkdir" = dir.exists(baseFolder) == T)
} else{baseFolder <- NULL}

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

# Unload Libraries to properly source v11 Package References -------------------------------------
unloadNamespace("mnp.v11b4"); unloadNamespace("mnp.v12b6"); unloadNamespace("sarc.v12b6")
loadNamespace("mnp.v11b6"); require("mnp.v11b6"); library("mnp.v11b6")

# Assign Parameters if Defined -------------------------------------------------------------------
selectRDs <- gb$AssignArgs(runID, baseFolder, token, selectRDs, redcapUp, gb)
gb$reportMd <- reportMd <- "~/report.Rmd"

# Execute Pipeline Functions ----------------------------------------------------------------------
gb$PrepareRun(token, baseFolder, runID, runLocal=runLocal) # If running local set runLocal = TRUE

GetPriorityCases <-
  function(selectRDs, samSheet = "samplesheet.csv", kwd="BN0") {
    csvFi <- read.csv(file.path(getwd(), samSheet))
    BN00 <- which(stringr::str_detect(csvFi$MP_num, kwd))
    if (length(BN00) > 0) {
      selectRDs <- c(selectRDs, csvFi$Sample_Name[BN00])
      message("Prioritizing these cases first:\n",
              paste(capture.output(selectRDs), collapse=" "))
    }
    return(selectRDs)
  }

selectRDs <- GetPriorityCases(selectRDs)

gb$StartRun(selectRDs, emailNotify=T, redcapUp=redcapUp) # Can be changed to default false
