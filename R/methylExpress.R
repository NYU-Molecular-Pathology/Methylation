#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library("base")
gb <- globalenv(); assign("gb", gb)

if(!require("devtools")){install.packages("devtools", quiet=T)}

# Parameters Input trailing commandline ----------------------------------------------------
token <-      args[1]
runID <-      args[2]
selectRDs <-  args[3]
baseFolder <- args[4] #NULL
redcapUpload <- args[5]

# Check Input Parameters ----------------------------------------------------
if(length(selectRDs)==0 | identical(selectRDs,NULL) | identical(selectRDs,"NULL")
){selectRDs=NULL}else {if(is.na(selectRDs)){selectRDs=NULL}}

if(length(baseFolder)==0 | identical(baseFolder,NULL) | identical(baseFolder,"NULL")){
    gb$baseFolder<-NULL
}else {
        if(is.na(baseFolder)){gb$baseFolder<-NULL}
    }
if(length(redcapUpload)==0 | identical(redcapUpload,NULL) | identical(redcapUpload,"NULL")){
    gb$redcapUpload<-T
}else {
    if(is.na(redcapUpload)){redcapUpload<-T
    }else{
        redcapUpload <- as.logical(args[5])
        }
}

# Message input on console
message("\n~~~~~~~~~~~~~~~~~~~~~Parameters input~~~~~~~~~~~~~~~~~~~~~\n")
message("token: ", token,"\nrunID: ", runID,"\nselectRDs: ", selectRDs, "\nbaseFolder: ", baseFolder, "\nredcapUpload: ", redcapUpload, "\n")

# Cancel if no token or runID
stopifnot(!is.null(token)); stopifnot(!is.null(runID))

if(!is.null(baseFolder) & !identical(baseFolder,NULL)){
    message("Trying custom run directory from input:","\n", baseFolder,"\n")
    isValid <- dir.exists(baseFolder)
    message("Checking if directory exists: ", isValid)
    if(isValid==F){
        message("Directory provided does not exist, try creating it first:\n",
                "mkdir ", baseFolder,"\n")
        stopifnot(isValid==TRUE)
    }
}else {baseFolder=NULL}


# Source GitHub Scripts ----------------------------------------------------
mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/"
script.list = c(
    "LoadInstallPackages.R",
    "SetRunParams.R",
    "MakeSampleSheet.R",
    "CopyInputs.R",
    "CopyOutput.R",
    "pipelineHelper.R"
)
scripts = paste0(mainHub, script.list)
lapply(scripts, function(i) {
    message("Sourcing: ", i)
    devtools::source_url(i)}
    )

# Define Parameters ----------------------------------------------------
gb$defineParams(
    methDir = baseFolder,
    baseDir = baseFolder,
    ApiToken = token
)
gb$setVar("runID", runID)

# Functions ----------------------------------------------------
checkBaseFolder <- function(baseFolder){
    if(is.null(baseFolder)){
        gb$baseFolder <- "/Volumes/CBioinformatics/Methylation/Clinical_Runs"
        gb$methDir <- "/Volumes/CBioinformatics/Methylation/Clinical_Runs"
        gb$baseDir <- "/Volumes/CBioinformatics/Methylation/Clinical_Runs"
    }else{
      gb$baseFolder <- baseFolder
    }
    if(length(baseFolder)>0){
        if(stringr::str_detect(baseFolder, pattern="Desktop")==T){
            warning("Trying to run methylation from Desktop working directory is not allowed")
            message("Try setting baseFolder to '~/Documents/' instead")
            stopifnot(stringr::str_detect(baseFolder, pattern="Desktop")==F)
        }
    }
    return(gb$baseFolder)
}

# Sets the working folder directory
SetBaseFolder <- function(token, baseFolder){
    methylPath <- gb$setRunDir(runID=gb$runID, workFolder = baseFolder)
    message("Working directory set to:\n", crayon::bgGreen(methylPath), "\n")
    gb$methDir <- gb$workFolder <- baseFolder
    gb$setVar("workFolder", baseFolder)
    gb$setVar("ApiToken", token) # assign the ApiToken & print params
}

GetLocalData <- function(rg){
    dat <- data.frame(
        sampleID = paste0(rg),
        bnumber = "NONE",
        senLi = paste0(rg),
        run_id = paste0(gb$runID),
        mp_number = "NONE",
        tech = "NONE",
        tech2 = "NONE",
        outFi = paste0(rg, ".html")
    )
    return(dat)
}

loop_local <- function(RGSet){
    for (rg in colnames(RGSet)) {
        thisSam <- RGSet[, rg]
        dat <- GetLocalData(rg)
        sentrix=dat$senLi
        RGsetEpic<-RGset<-thisSam
        rmarkdown::render(
            reportMd, "html_document", dat$outFi, getwd(), quiet = FALSE,
            params = list(token = gb$ApiToken, rundata = dat)
        )
        }
}

RunLocalIdats <- function(runID){
    if(!file.exists(file.path(getwd(), paste0(runID,".xlsm")))){
        idatFiles <- dir(path = getwd(), pattern = ".idat", full.names = T)
        idatBase <- unique(substring(idatFiles, 1, nchar(idatFiles) - 9))
        RGSet<-minfi::read.metharray(basenames =idatBase, force=TRUE)
        loop_local(RGSet)
    }else{
        gb$readSheetWrite(runID = runID)
        gb$moveSampleSheet(baseFolder, runID)
    }
}


# Executes the functions in order to setup a run
PrepareRun <- function(token, baseFolder=NULL, runID, runLocal=F){
    if(runLocal==F){
        gb$checkMounts()
        gb$checkValidRun(runID)
    }
    baseFolder <- checkBaseFolder(baseFolder)
    SetBaseFolder(token, baseFolder)
    setwd(file.path(baseFolder, runID))
    if(runLocal==F) {
        gb$copyWorksheetFile(runID = runID) # copies the xlsm file
        gb$readSheetWrite(runID = runID) # reads xlsm and generates input .csv samplesheet
        gb$get.idats() # Copy idat files to current folder from molecular and snuderlabspace to cwd
        gb$moveSampleSheet(baseFolder, runID) #copies outputs temp to desktop for QC.Rmd
    } else{
        RunLocalIdats(runID)
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
assign("redcapUpload", redcapUpload)

# Execute Functions ----------------------------------------------------
PrepareRun(token, baseFolder, runID) # If running local and  runLocal = TRUE
StartRun(selectRDs, emailNotify=T, redcapUp=T) # can change to default false
