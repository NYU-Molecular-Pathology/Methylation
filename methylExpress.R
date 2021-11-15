#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library("base"); gb <- globalenv(); assign("gb", gb)

# If you have issues with conda and r, try running the following:
# conda create -n r_env r-essentials r-base
# conda activate r_env

# Main Parameters trailing commandline
token<-args[1]; runID<-args[2]; selectRDs<-args[3]
if(length(selectRDs)==0){
    selectRDs=NULL
    }else{
        if(is.na(selectRDs)){selectRDs=NULL}
    }

# Check Input Params
cat("\n~~~~~~Parameters input:",token, runID, selectRDs,sep='\n')
message("~~~~~~~~~~~~~~~~~~~~~"); stopifnot(!is.null(token)); stopifnot(!is.null(runID))

# Paths to the GitHub Repo files
mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/"
script.list = c("LoadInstall_new.R","SetRunParams.R", "CopyInputs.R","CopyOutput.R", "pipelineHelper.R")

# Source GitHub Scripts
scripts = paste0(mainHub, script.list)
invisible(lapply(scripts, function(i) {devtools::source_url(i)}))


gb$defineParams() # Define any custom directories here

# Execute Methylation Run
if(!is.null(runID)){gb$setVar("runID", runID)}

# Executes the functions in order to setup a run
prepareRun <- function(token){
    runValid <- gb$checkValidRun(gb$runID)
    message("Is the runID valid? ", runValid)
    if(!runValid){
        message(crayon::bgRed$white$bold("runID",gb$runID,"is not valid"))
        message(crayon::bgBlue$white$bold(paste0(gb$runID,".xlsm"),"not found in worksheets folder"))
        stopifnot(runValid)
        }
    methylPath <- gb$setRunDir(gb$runID)
    message("Working directory set to:"); cat(crayon::bgGreen(methylPath)); setwd(methylPath)
    gb$setVar("ApiToken", token) # assign the ApiToken & print params
    gb$copyWorksheetFile(runID = gb$runID) # copies the xlsm file
    gb$readSheetWrite() # reads xlsm and generates input .csv samplesheet
    gb$get.idats() # Copy idat files to current folder from molecular and snuderlabspace to cwd
    gb$moveSampleSheet(gb$methDir) #copies outputs temp to desktop for QC.Rmd
    #gb$classifierInstall(instNew = F, rmpkg = F) # Loads pipeline or installs new
}

#gb$startRun <- function(selectRDs=NULL, runID=NULL, emailNotify=T){
#    if(!is.null(selectRDs)){
#        sampleOrder <- reOrderRun(selectRDs) # Re-order sample report generation for priority
#        makeReports.v11b6(skipQC=F, email=emailNotify, cpReport=T, selectSams=sampleOrder, redcapUp=T)
#    } else {makeReports.v11b6(skipQC=F, email=emailNotify, cpReport=F, selectSams=NULL, redcapUp=T)}
#}

prepareRun(token); gb$startRun(selectRDs)
