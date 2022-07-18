#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library("base")
gb <- globalenv(); assign("gb", gb)

if(!require("devtools")){install.packages("devtools")}

# Main Parameters trailing commandline
token <-      args[1]
runID <-      args[2]
selectRDs <-  args[3]
baseFolder <- args[4] #NULL
redcapUpload <- args[5]

# Check Parameters Input
if(length(selectRDs)==0 | identical(selectRDs,NULL) | identical(selectRDs,"NULL")
){selectRDs=NULL}else {if(is.na(selectRDs)){selectRDs=NULL}}

if(length(baseFolder)==0 | identical(baseFolder,NULL) | identical(baseFolder,"NULL")
  ){gb$baseFolder<-NULL}else {if(is.na(baseFolder)){gb$baseFolder<-NULL}}

if(length(redcapUpload)==0 | identical(redcapUpload,NULL) | identical(redcapUpload,"NULL")){
    gb$redcapUpload<-T
}else {
  if(is.na(redcapUpload)){
    redcapUpload<-T
  }else{
    redcapUpload <- as.logical(args[5])
  }
}

# Check Input Params
message("\n~~~~~~~~~~~~~~~~~~~~~Parameters input~~~~~~~~~~~~~~~~~~~~~\n")
message("token: ", token,"\nrunID: ", runID,"\nselectRDs: ", selectRDs, "\nbaseFolder: ", baseFolder, "\nredcapUpload: ", redcapUpload)
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

# Cancel if no token or runID
stopifnot(!is.null(token)); stopifnot(!is.null(runID))

if(!is.null(baseFolder) & !identical(baseFolder,NULL)){
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
mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/"
script.list = c("LoadInstall_new.R","SetRunParams.R","MakeSampleSheet.R", "CopyInputs.R","CopyOutput.R", "pipelineHelper.R")

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
        if(stringr::str_detect(baseFolder, pattern="Desktop")==T){
            warning("Trying to run methylation from Desktop working directory is not allowed")
            message("Try setting baseFolder to '~/Documents/' instead")
            stopifnot(stringr::str_detect(baseFolder, pattern="Desktop")==F)
        }
    }
    runValid <- gb$checkValidRun(gb$runID)

    methylPath <- gb$setRunDir(gb$runID, workFolder = baseFolder)
    message("Working directory set to:")
    cat(crayon::bgGreen(methylPath))
    gb$workFolder <- baseFolder
    gb$setVar("workFolder", baseFolder)
    message(" ")
    setwd(methylPath)
    gb$setVar("ApiToken", token) # assign the ApiToken & print params
    gb$copyWorksheetFile(runID = gb$runID) # copies the xlsm file
    gb$readSheetWrite() # reads xlsm and generates input .csv samplesheet
    gb$get.idats() # Copy idat files to current folder from molecular and snuderlabspace to cwd
    gb$moveSampleSheet(gb$methDir) #copies outputs temp to desktop for QC.Rmd
    message("workFolder:", gb$workFolder)
    #gb$classifierInstall(instNew = F, rmpkg = F) # Loads pipeline or installs new
}

gb$prepareRun(token, baseFolder)

if(!is.null(selectRDs)){selectRDs <- stringr::str_split(selectRDs, ",")}

gb$startRun <- function(selectRDs=NULL, emailNotify=T, redcapUp=T){
    gb$msgFunName("https://github.com/NYU-Molecular-Pathology/Methylation/edit/main/methylExpress.R","startRun")
    message("Default Params: selectRDs=NULL, emailNotify=T")
    nullVal = ifelse(is.null(selectRDs),"NULL", selectRDs)
    gb$msgParams(nullVal, emailNotify)

    if(!is.null(selectRDs)){
        sampleOrder <- gb$reOrderRun(selectRDs) # Re-order sample report generation for priority
        gb$makeReports.v11b6(skipQC=F, email=emailNotify, cpReport=T, selectSams=sampleOrder, redcapUp=T)
    } else {
        gb$makeReports.v11b6(skipQC=F, email=emailNotify, cpReport=T, selectSams=NULL, redcapUp=T)
    }
}

assign("redcapUpload", redcapUpload)
gb$startRun(selectRDs, emailNotify=T,redcapUp=T) # can change to default false
