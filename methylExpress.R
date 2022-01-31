#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library("base")
gb <- globalenv(); assign("gb", gb)

if(!require("devtools")){install.packages("devtools")}

# Main Parameters trailing commandline
token<-args[1]; runID<-args[2]; selectRDs<-args[3]
baseFolder <- args[4] #NULL

# Check Parameters Input
if(length(selectRDs)==0 | identical(selectRDs,NULL) | identical(selectRDs,"NULL")
){selectRDs=NULL}else {if(is.na(selectRDs)){selectRDs=NULL}}
if(length(baseFolder)==0 | identical(baseFolder,NULL) | identical(baseFolder,"NULL")){
    gb$baseFolder<-NULL
}else {if(is.na(baseFolder)){gb$baseFolder<-NULL}}

# Check Input Params
message("\n~~~~~~~~~~~~~~~~~~~~~Parameters input~~~~~~~~~~~~~~~~~~~~~\n")
message("token: ", token,"\nrunID: ", runID,"\nselectRDs: ", selectRDs, "\nbaseFolder: ", baseFolder)
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

# Cancel if no token or runID
stopifnot(!is.null(token)); stopifnot(!is.null(runID))

if(!is.na(baseFolder) & !is.null(baseFolder) & !identical(baseFolder,NULL)){
    message("Trying custom run directory from input:","\n", baseFolder,"\n")
    isValid <- dir.exists(baseFolder)
    message("Checking if directory exists: ", isValid)
    if(isValid==F){message("Directory does not exist, trying to create path:\n", baseFolder,"\n")
        tryCatch(
            expr={dir.create(baseFolder)},
            warning=function(er){message("\n","One warning during directory creation","\n",er)},
            error=function(er){message("\n","An Error during directory creation:","\n",er,"\n")},
            finally={stopifnot(dir.exists(baseFolder))}
        )
    }
}else {baseFolder=NULL}

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
gb$prepareRun <- function(token, baseFolder=NULL){
    if(is.null(baseFolder)){
        gb$baseFolder <- "/Volumes/CBioinformatics/Methylation/Clinical_Runs"
        gb$methDir <- "/Volumes/CBioinformatics/Methylation/Clinical_Runs"
        gb$baseDir <- "/Volumes/CBioinformatics/Methylation/Clinical_Runs"
    }

    if(length(baseFolder)>0){
        if(str_detect(baseFolder, pattern="Desktop")==T){
            warning("Trying to run methylation from Desktop working directory is not allowed")
            message("Try setting baseFolder to '~/Documents/' instead")
            stopifnot(str_detect(baseFolder, pattern="Desktop")==F)
        }
    }
    runValid <- gb$checkValidRun(gb$runID)
    message("Is the runID valid? ", runValid)
    if(!runValid){
        message(crayon::bgRed$white$bold("runID",gb$runID,"is not valid"))
        message(crayon::bgBlue$white$bold(paste0(gb$runID,".xlsm"),"not found in worksheets folder"))
        stopifnot(runValid)
    }
  
    methylPath <- gb$setRunDir(gb$runID, workFolder = baseFolder)
    message("Working directory set to:")
    cat(crayon::bgGreen(methylPath))
    message(" ")
    setwd(methylPath)
    gb$setVar("ApiToken", token) # assign the ApiToken & print params
    gb$copyWorksheetFile(runID = gb$runID) # copies the xlsm file
    gb$readSheetWrite() # reads xlsm and generates input .csv samplesheet
    gb$get.idats() # Copy idat files to current folder from molecular and snuderlabspace to cwd
    gb$moveSampleSheet(gb$methDir) #copies outputs temp to desktop for QC.Rmd
    #gb$classifierInstall(instNew = F, rmpkg = F) # Loads pipeline or installs new
}

gb$prepareRun(token, baseFolder)

if(!is.null(selectRDs)){selectRDs <- stringr::str_split(selectRDs, ",")}
gb$startRun(selectRDs)
