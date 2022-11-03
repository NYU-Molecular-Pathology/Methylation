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
gb$PrepareRun(token, baseFolder, runID, runLocal=runLocal) # If running local and  runLocal = TRUE

require(compiler)
enableJIT(3)
mainPage = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/Report-Scripts/"
scripts <- c("ClassTables.R", "MLH1_Functions.R", "PipeLineU.R",
             "RedcapOutput.R", "TsneFunctions.R", "cnvggplotly.R")
scripts <- paste0(mainPage, scripts)
lapply(scripts, function(i) {message("Sourcing: ", i);devtools::source_url(i)})

unloadNamespace("mnp.v11b4")
unloadNamespace("mnp.v12b6")
unloadNamespace("sarc.v12b6")
loadNamespace("mnp.v11b6")
require("mnp.v11b6")
library("mnp.v11b6")

do_report <-  function(data = NULL, genCn=F) {
  msgFunName(pipeLnk,"do_report")
  msgParams("data")
  #data = data[1, ]
  if(!is.null(data)){
    dat <- getRunData(data)
    message(paste0(capture.output(dat), collapse="\n"))
    RGsetEpic <- getRGset(getwd(), dat$senLi)
    RGset <- RGsetEpic[,1]
    if(genCn==T){generateCNVpng(RGsetEpic,dat$sampleID)}
    msgRunUp(dat$sampleID, dat$run_id, dat$senLi)
    tryCatch(
      expr = {
        rmarkdown::render(
          reportMd,
          "html_document",
          dat$outFi,
          getwd(),
          quiet = FALSE,
          params = list(
            token = gb$ApiToken,
            rundata = dat,
            RGsetEpic = RGsetEpic,
            knitDir = getwd()
          )
        )
      },
      error=function(e){
        beepr::beep(1)
        message(bkRed("Report Generation Failed:"),"\n", dat$outFi)
        message("The following error returned:\n", e)
        stopifnot("Report Generation Failed"=FALSE)
      }, finally=message("\nRunning next sample\n")
    )
  } else {message(bkRed("Data is NULL, check your SampleSheet.csv"))}
}

gb$StartRun(selectRDs, emailNotify=T, redcapUp=redcapUp) # can change to default false
