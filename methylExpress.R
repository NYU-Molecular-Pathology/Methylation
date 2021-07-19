#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library("base"); gb <- globalenv(); assign("gb", gb)

# Main Parameters trailing commandline
token<-args[1]; runID<-args[2]; selectRDs<-args[3]
if(length(selectRDs)==0){selectRDs=NULL}else{if(is.na(selectRDs))selectRDs=NULL}

# Check Input Params
stopifnot(!is.null(token))
params <- c(token, runID, selectRDs)
invisible(lapply(params, function(i){message(i, " = ", get(i), " (",typeof(i),")\n")}))

# Paths to the GitHub Repo files
mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/"
script.list = c("LoadInstall_new.R","SetRunParams.R","CopyInputs.R","CopyOutput.R", "pipelineHelper.R")

# Source GitHub Scripts
scripts = paste0(mainHub, script.list)
invisible(lapply(scripts, function(i) {devtools::source_url(i)}))

# Execute Methylation Run
if(!is.null(runID)){gb$setVar("runID", runID)}
gb$prepareRun(token); gb$startRun(selectRDs)


