#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library(base)
# Main Parameters
token <- paste0(args[1])
runID <- paste0(args[2]) # Set to NULL to get latest
selectRDs <- ifelse(length(args[3])>0, args[3], NULL)

# Source Github Scripts
mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/"
scripts = paste0(mainHub, c("LoadInstall_new.R", "pipelineHelper.R"))
gb <- globalenv(); assign("gb", gb)
invisible(lapply(scripts, function(i) {devtools::source_url(i)}))

# Execute Methylation Run
if(!is.null(runID)){gb$setVar("runID", runID)}
gb$prepareRun(token) #print(gb$prepareRun)
gb$startRun(selectRDs)

