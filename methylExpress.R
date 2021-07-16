#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library("base")

# Main Parameters
token<-args[1]; runID<-args[2]; selectRDs<-args[3]
if(length(selectRDs)==0){selectRDs=NULL}else{if(is.na(selectRDs))selectRDs=NULL}

# Check Params
stopifnot(!is.null(token))
invisible(lapply(base::ls(), function(i){message(i, " = ", get(i), " (",typeof(i),")\n")}))
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

# Source Github Scripts
mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/"
scripts = paste0(mainHub, c("LoadInstall_new.R", "pipelineHelper.R"))
gb <- globalenv(); assign("gb", gb)
invisible(lapply(scripts, function(i) {devtools::source_url(i)}))

# Execute Methylation Run
if(!is.null(runID)){gb$setVar("runID", runID)}
gb$prepareRun(token) #print(gb$prepareRun)
gb$startRun(selectRDs)