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
message("\n~~~~~~~~~~~~~~~~~~~~~Parameters input~~~~~~~~~~~~~~~~~~~~~\n")
CheckInputArg <- function(varValue, gb, defVal = NULL) {
    varStr <- deparse(substitute(varValue))
    if (length(varValue) == 0 | identical(varValue, NULL) | identical(varValue, "NULL")) {
        gb[[varStr]] <- defVal
    } else{ varValue <- ifelse(is.na(varValue), NULL, varValue)}
    message(varStr, ": " , ifelse(is.null(varValue), "NULL", varValue))
    return(assign(varStr, varValue, envir = gb))
}

CheckInputArg(token, gb)
CheckInputArg(runID, gb)
CheckInputArg(selectRDs, gb)
CheckInputArg(baseFolder, gb)
CheckInputArg(redcapUp, gb)
CheckInputArg(runLocal, gb)

# Cancel if no token or runID
stopifnot(!is.null(token)); stopifnot(!is.null(runID))

if(!is.null(baseFolder)) {
    message("Checking if custom run directory is valid:\n", baseFolder, "\n")
    if (dir.exists(baseFolder) == F) {
        message("Directory does not exist, try creating it first:\n mkdir ", baseFolder,"\n")
        stopifnot(dir.exists(baseFolder) == TRUE)
    }
}

# Source GitHub Scripts ----------------------------------------------------
mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/"
script.list = c(
    "all_installer.R", "LoadInstallPackages.R",
    "SetRunParams.R", "MakeSampleSheet.R",
    "CopyInputs.R", "CopyOutput.R",
    "pipelineHelper.R", "CustomRuns.R"
)
scripts = paste0(mainHub, script.list)
lapply(scripts, function(i){message("Sourcing: ", i); devtools::source_url(i)})

# Define Parameters ----------------------------------------------------
gb$defineParams(methDir = baseFolder, baseDir = baseFolder, ApiToken = token)
gb$setVar("runID", runID)

# Executes the functions in order to setup a run
PrepareRun <- function(token, baseFolder=NULL, runID, runLocal=F){
    gb$msgFunName(paste0(mainHub,"methylExpress.R"),"PrepareRun")
    if(runLocal==F){
        gb$checkMounts()
        gb$checkValidRun(runID)
        baseFolder <- gb$SetBaseFolder(token, baseFolder, runID)
        setwd(file.path(baseFolder, runID))
        gb$copyWorksheetFile(runID = runID) # copies the xlsm file
        gb$readSheetWrite(runID = runID) # reads xlsm and generates input .csv samplesheet
        gb$get.idats() # Copy idat files to current folder from molecular and snuderlabspace to cwd
        gb$moveSampleSheet(baseFolder, runID) #copies outputs temp to desktop for QC.Rmd
    } else{
        baseFolder <- gb$SetBaseFolder(token, baseFolder, runID)
        setwd(file.path(baseFolder, runID))
        gb$RunLocalIdats(runID, token)
    }
}

# Starts run knitting html reports loop
StartRun <- function(selectRDs = NULL, emailNotify = T, redcapUp = T) {
    gb$msgFunName(paste0(mainHub,"methylExpress.R"),"startRun")
    # Re-order sample report generation for priority
    if (!is.null(selectRDs)) {runOrder <- gb$reOrderRun(selectRDs) }else{runOrder <- NULL}
    gb$makeReports.v11b6(
        skipQC = F,            # Don't skip QC generation
        email = emailNotify,   # to email after Run complete
        cpReport = T,          # Flag to copy files to network drive
        selectSams = runOrder, # Prioritize specific RD-numbers
        redcapUp = redcapUp    # Flag to import files to REDCap
    )
}

if(!is.null(selectRDs)){selectRDs <- stringr::str_split(selectRDs, ",")}
assign("redcapUp", redcapUp); assign("redcapUpload", redcapUp)
# Execute Functions ----------------------------------------------------
PrepareRun(token, baseFolder, runID, runLocal=runLocal) # If running local and  runLocal = TRUE
StartRun(selectRDs, emailNotify=T, redcapUp=redcapUp) # can change to default false
