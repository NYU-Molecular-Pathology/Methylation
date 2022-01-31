#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library("base"); gb <- globalenv(); assign("gb", gb)

if(!require("devtools")){install.packages("devtools")}

# If you have issues with conda and r, try running the following:
# conda create -n r_env r-essentials r-base
# conda activate r_env

# Main Parameters trailing commandline
token<-args[1]; runID<-args[2]; selectRDs<-args[3]
baseFolder <- args[4] #NULL

# Check Parameters Input
if(length(selectRDs)==0){selectRDs=NULL}else {
    if(is.na(selectRDs)){selectRDs=NULL}
}

if(length(baseFolder)==0){
    baseFolder<-NULL
}else {
    if(is.na(baseFolder)){baseFolder<-NULL}
}

# Check Input Params
cat("\n~~~~~~Parameters input:")
message("token: ", token)
message("runID: ", runID)
message("selectRDs: ", selectRDs)
message("baseFolder: ", baseFolder)
message("~~~~~~~~~~~~~~~~~~~~~")

# Cancel if no token or runID
stopifnot(!is.null(token)); stopifnot(!is.null(runID))

if(!is.na(baseFolder) & !is.null(baseFolder)){
    message("Trying custom run directory from input:","\n", baseFolder,"\n")
    isValid <- dir.exists(baseFolder)
    message("Checking if directory exists: ", isValid)
    if(isValid==F){
        message("Directory does not exist, trying to create path:\n", baseFolder,"\n")

        tryCatch(
            expr={dir.create(baseFolder)},
            warning=function(er){
                message("\n","One warning during directory creation","\n")
                message("\n",er,"\n")
            },
            error=function(er){
                message("\n","An Error during directory creation:","\n")
                message("\n",er,"\n")
            },
            finally={
                stopifnot(dir.exists(baseFolder))
            }
        )
    }
}

# Paths to the GitHub Repo files
mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/"
script.list = c("LoadInstall_new.R","SetRunParams.R", "CopyInputs.R","CopyOutput.R", "pipelineHelper.R")

# Source GitHub Scripts
scripts = paste0(mainHub, script.list)
invisible(lapply(scripts, function(i) {devtools::source_url(i)}))

# Define any custom directories here, default is NULL
#mnp.pk.loc = NULL, ApiToken = NULL, methDir = NULL, clinDrv = NULL, rschOut = NULL, clinOut = NULL,
#rsch.idat = NULL, clin.idat = NULL, QC_file = NULL, isMC = T, baseDir = NULL, runID = NULL

gb$defineParams(
    methDir = baseFolder,
    baseDir = baseFolder,
    ApiToken = token
)

# Execute Methylation Run
if(!is.null(runID)){gb$setVar("runID", runID)}

# Executes the functions in order to setup a run
prepareRun <- function(token,baseFolder){
    runValid <- gb$checkValidRun(gb$runID)
    message("Is the runID valid? ", runValid)
    if(!runValid){
        message(crayon::bgRed$white$bold("runID",gb$runID,"is not valid"))
        message(crayon::bgBlue$white$bold(paste0(gb$runID,".xlsm"),"not found in worksheets folder"))
        stopifnot(runValid)
    }
    if(is.null(baseFolder)){
        baseFolder <- "/Volumes/CBioinformatics/Methylation/Clinical_Runs"
        }
    gb$methDir <- baseFolder
    gb$baseDir <- baseFolder
    methylPath <- gb$setRunDir(gb$runID, workFolder = baseFolder)
    if(str_detect(baseFolder, pattern="Desktop")==T){
        warning("Trying to run methylation from Desktop working directory is not allowed")
        message("Try setting baseFolder to '~/Documents/' instead")
        stopifnot(str_detect(baseFolder, pattern="Desktop")==F)
        }
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

prepareRun(token, baseFolder)
gb$startRun(selectRDs)
