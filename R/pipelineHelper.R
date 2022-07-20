#!/usr/bin/env Rscript
gb <- globalenv(); assign("gb", gb)
dsh = "-----------"
bky <- function(txtVar,...){crayon::black$bgYellow$bold(txtVar,...)}
bkRed <- function(txtVar,...){crayon::bgRed$bold$white(txtVar,...)}
bkgrn <- function(...){crayon::black$bgGreen$bold(paste0(...))}
bkBlu <- function(...){crayon::bgBlue$bold$white(paste0(...))}

# Global Variables ----------------------------------
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
reportMd <- "~/report.Rmd" # From curl github download
QC_file <- "~/Methyl_QC.Rmd" # From curl github download
pipeLnk <- "https://github.com/NYU-Molecular-Pathology/Methylation/edit/main/pipelineHelper.R"
predictionPath <- "/Volumes/CBioinformatics/Methylation/in_house/mnp.v116/mnp.v11b6/data/rfpred.v11b6.RData"
QC_file <- "~/Methyl_QC.Rmd"

# Message Functions and Variables for Debugging ----------------------------------
msgFunName <- function(pthLnk, funNam){
    message("\nExecuting function: ", bky(funNam), " from RScript in:\n", pthLnk)
}

msgParams <- function(...){
    message("\n", crayon::bgGreen("Params passed:"),
            "\n", paste(..., collapse = "_"), "\n")
}

# Helper function to return the index of priority selected samples first
reOrderRun <- function(selectRDs, sh=NULL){
    msgFunName(pipeLnk,"reOrderRun"); msgParams(selectRDs, sh)

    if(is.null(selectRDs)){return(NULL)}
    if(is.null(sh)){sh<-"samplesheet.csv"}
    allRd <- as.data.frame(read.csv(sh))
    runFirst <- which(allRd[,1] %in% selectRDs)
    runAfter <- which(!(allRd[,1] %in% selectRDs))
    return(c(runFirst,runAfter))
}

# Saves the methyl CNV as a png file in the cwd
generateCNVpng <- function(RGsetEpic, sampleName) {
    msgFunName(pipeLnk,"generateCNVpng")

    imgName <- paste(sampleName, "cnv.png", sep="_")
    MsetEpic <- minfi::preprocessRaw(RGsetEpic)
    png(filename=imgName,width=1820, height=1040)
    mnp.v11b6::MNPcnvplot(MsetEpic)
    dev.off(); message("File saved:\n",imgName,"\n")
}

getRGset <- function(runPath, sentrix){
    msgFunName(pipeLnk,"getRGset")

    barcode = stringr::str_split_fixed(sentrix, "_",2)[1]
    RGsetEpic <- minfi::read.metharray(file.path(runPath, sentrix), verbose = T, force = T)
    aEpic=c(array="IlluminaHumanMethylationEPIC", annotation="ilm10b4.hg19")
    a450k=c(array="IlluminaHumanMethylation450k", annotation="ilmn12.hg19")
    if (barcode >= as.numeric("204220033000")) {RGsetEpic@annotation=aEpic}
    if (RGsetEpic@annotation['array']=="IlluminaHumanMethylation450k") {RGsetEpic@annotation=a450k}
    return(RGsetEpic)
}

CopyQCrmd <- function(runID){
    msgFunName(pipeLnk, "CopyQCrmd")
    if (!file.exists(QC_file)) {message(crayon::bgRed("QC_file.rmd not found:"), "\n", QC_file)}
    #fs::file_copy(QC_file, getwd(), overwrite = T)
    qcFileName = paste0(runID, "_QC.html") # output file name
    if(file.exists(file.path(getwd(), qcFileName))) {
        message(qcFileName, "Already Exists! Detete file to output new QC")
        return(NULL)
    }
    file.copy(from = QC_file, to = file.path(getwd(),basename(QC_file)), copy.mode = F, overwrite = F)
    return(qcFileName)
    }


# QC REPORT maker: knits the QC RMD file
generateQCreport <- function(runID=NULL) {
    msgFunName(pipeLnk, "generateQCreport")

    if (is.null(runID)){runID<-paste0(basename(getwd()))}
    qcFileName <- CopyQCrmd(runID)
    if(!is.null(qcFileName)){
        outQCpath <- file.path(getwd(), qcFileName)
        rmdToKnit <- dir(getwd(), "*QC.Rmd", full.names = T)[1]
        rmarkdown::render(
            rmdToKnit,
            output_file = outQCpath,
            quiet = F,
            params = list(runID = runID)
        )
        qcCache <- stringr::str_replace_all(string = rmdToKnit, ".Rmd", "_cache")
        unlink(qcCache, recursive = T) #clear cache
        gb$uploadToRedcap(outQCpath, F)
    }
}

# Sends an email notification that the run is complete from redcap admin
launchEmailNotify <- function(runID){
    msgFunName(pipeLnk,"launchEmailNotify")
    isMC = sjmisc::str_contains(runID, "MGDM") | sjmisc::str_contains(runID, "MC")
    com <- ifelse(isMC==T, "sample_qc", "sample_research") # research or clinical notification
    record = data.frame(record_id=paste0(runID,"_QC"), comments=com)
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
    RCurl::postForm(
        gb$apiLink,
        token = gb$ApiToken,
        content = 'record',
        format = 'json',
        type = 'flat',
        data = datarecord
        )
    cat(crayon::white$bgBlue$bold("Check email to confirm Email Notification Created"), sep="\n")
}

# FUN: Creates the QC record for the current run on redcap if it does not exist
CreateRedcapRecord <- function(runID=NULL, recordWord="QC"){
    msgFunName(pipeLnk,"CreateRedcapRecord")
    if(is.null(runID)){runID<-paste0(basename(gb$workDir))}

    record = c(
        record_id = paste0(runID, "_", recordWord),
        run_number = runID
        )
    cntrl <- jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
    RCurl::postForm(
        gb$apiLink,
        token = gb$ApiToken,
        content = 'record',
        format = 'json',
        type = 'flat',
        data = cntrl
    )
    message(dsh, "Created Record: ", record[1], dsh, "\n")
}

# Check if the QC File will be read
checkRunOutput <- function(runID) {
    msgFunName(pipeLnk,"checkRunOutput")
    redcapFi <- paste0(runID,"_Redcap.csv")

    csvLocation <- file.path(fs::path_home(),"Desktop", runID, redcapFi)
    if (!file.exists(csvLocation)) {
        message(bkRed("File not found:")," ", csvLocation,"\n",
                bkRed("QC Summary Table will not Knit"),"\n")
    } else{
        outputPath <- getwd()
        if(file.exists(file.path(getwd(), redcapFi))){
            outputPath <- paste0(Sys.Date(), "_", redcapFi)
        }
        fs::file_copy(path=csvLocation, new_path=outputPath, overwrite = T)
    }
}

# gets rid of desktop files if run is successful
tidyUpFiles <- function(runID){
    msgFunName(pipeLnk,"tidyUpFiles")
    deskDir <- file.path(fs::path_home(),"Desktop",runID)
    backupD <- file.path(gb$methDir,"csvRedcap")
    if(!dir.exists(backupD)){dir.create(backupD)}
    file.copy(deskDir, backupD, overwrite=T, recursive = T, copy.mode = T)
    unlink(deskDir,T,T)
}

msgProgress <- function(msg, i, samList) {
    nOfTotal <- paste(i, "of", length(samList), dsh)
    startEnd <- ifelse(msg == 1, "Now Running", "Completed Report")
    if(msg==1){
        cat("\n", bky(dsh, startEnd , nOfTotal), sep = "\n")
        }else{
            cat("\n", bkgrn(dsh, startEnd , nOfTotal), sep = "\n")
        }
}

msgRunUp <- function(sampleID,run_id,senLi){
    cat(
        "\nCurrent Sample:\n",
        "samplename_data=", bky(sampleID),
        " run_id=", bky(run_id),
        "\npathEpic:", bky(file.path(getwd(), senLi)), "\n"
        )
}

msgSamSheet <- function(samSh) {
    if (length(samSh) > 1) {
        message(crayon::black$bgRed$bold("Multiple samplesheets found:"))
        cat(samSh, sep = "\n")
        samSh <- samSh[stringr::str_detect(samSh, pattern = "\\$", negate = T)]
        message(crayon::black$bgGreen("Reading the following .xlsm in current directory:"))
    } else{message(crayon::black$bgGreen("Reading the following .xlsm in current directory:"))}
    cat(samSh, sep = "\n")
    stopifnot(!length(samSh) > 1)
    return(samSh)
}

getRunData <- function(data) {
    runDt <- data.frame(
        sampleID = paste0(data[, 1]),
        bnumber = paste0(data[, 2]),
        senLi = paste0(data[, 5]),
        run_id = paste0(data[, 7]),
        mp_number = paste0(data[, 8]),
        tech = paste0(data[, 9]),
        tech2 = paste0(data[, 10]),
        outFi = paste0(data[, 1], ".html")
    )
    return(runDt)
}

NameControl <- function(data, runId) {
    library("data.table")
    if (any(data[, 1] %like% 'control')) {
        cntrl <- which(data[, 1] %like% 'control') #DNA_Number
        data[cntrl, 1] <- paste0(runId, "_control")
    } else{
        if (any(data[, 2] %like% 'control')) {
            cntrl <- which(data[, 2] %like% 'control') #DNA_Number
            data[cntrl, 1] <- paste0(runId, "_control")
        } else{
            warning('No word "control" in RD-number found in samplesheet')
        }
        return(data)
    }
}

ReadSamSheet <- function(samList){
    msgFunName(pipeLnk, "ReadSamSheet")
    msgParams("samList","=",samList)

    samSh <- gb$GrabSampleSheet()
    xlSheets <- readxl::excel_sheets(samSh)
    redSheet <- as.integer(which(grepl("REDCap",xlSheets)==T))
    message("Excel sheet names:")
    print(xlSheets)
    if(length(redSheet)==0){
        warning("REDCap sheet not found in Sheet names, setting sheet #3")
        redSheet <- 3
        }
    message(bkgrn("Sheet Index containing 'REDCap_Import':"), " ", redSheet,"\n")
    samSh <- readxl::read_excel(samSh, sheet=redSheet, range = "A1:M97", col_types = c("text"))
    message(bkgrn("SampleSheet:"))
    samplesSheet <- as.data.frame(samSh)[samList, 1:13]
    print(samplesSheet)
    return(samplesSheet)
}

# FUN: Parses the WetLab .xlsm sheet in the current directory
checkSamSh <- function(samList){
    msgFunName(pipeLnk, "checkSamSh")
    msgParams("samList","=",samList)
    require(rmarkdown)
    wksh <- ReadSamSheet(samList)
    #wksh <- NameControl(wksh, wksh$run_number[1])
    message(bkgrn(wksh),":"," ")
    print(wksh)
    stopifnot(!is.null(wksh))
    rownames(wksh)<- wksh[,1]
    return(wksh)
}

getRunList <- function(data, samList){
    msgFunName(pipeLnk, "getRunList")
    toRun <- unlist(lapply(samList, FUN=function(i){
        nOutDir = file.path(getwd(), paste0(data[i, 1], ".html"))
        if (file.exists(nOutDir)) {
            message(bky(basename(nOutDir), "already exists! Skipping..."),"\n")
        }else{return(i)}
    }))
    return(toRun)
}

# Helper function called by makeReports.v11b6 to generate the HTML report
do_report <- function(data = NULL, genCn=F) {
    msgFunName(pipeLnk,"do_report")
    msgParams("data")
    print(data)
    if(!is.null(data)){
        dat <- getRunData(data)
        RGsetEpic <- getRGset(getwd(), dat$senLi)
        if(genCn==T){generateCNVpng(RGsetEpic,dat$sampleID)}
        msgRunUp(dat$sampleID,dat$run_id,dat$senLi)
        tryCatch(
            expr={rmarkdown::render(
                reportMd, "html_document", dat$outFi, getwd(), quiet=FALSE,
                params = list(token=gb$ApiToken, rundata=dat))},
            error=function(e){
                message(bkRed("Report Generation Failed:"),"\n", dat$outFi,"\n")
                message(bkRed("The following error returned:"),"\n", e)
            }, finally=message("\nRunning next sample\n")
        )
    } else {message(bkRed("Data is NULL, check your SampleSheet.csv"))}
}

# FUN: Iterates over each sample in the csv file to generate a report
loopRender <- function(samList = NULL, data, redcapUp = T){
    msgFunName(pipeLnk, "loopRender")
    msgParams("samList = NULL, data, redcapUp = T")
# Debug: data <- read.csv("samplesheet.csv", strip.white=T)

    stopifnot(!is.null(data))
    #data <- NameControl(data, data$RunID[1])
    if (is.null(samList)) {
        samList <- 1:length(data$Sample_Name!=0)
    }
    wksh <- checkSamSh(samList)
    toRun <- getRunList(data, samList)
    for (i in toRun) {
        msgProgress(1, i, samList)
        do_report(data = data[i, ], gb$genCn)
        msgProgress(2, i, samList)
        if (redcapUp == T) {
            sh_Dat <- wksh[data[i, 1],]
            gb$importSingle(sh_Dat)
            }
    }
    message(crayon::black$bgGreen$bold(dsh, "RUN COMPLETE", dsh))
}

#' REPORT: Generates Html reports to cwd with samplesheet.csv
#' @param runPath The location of samplesheet.csv and idats
#' @param sheetName name of samplesheet if it is not "samplesheet.csv"
#' @param selectSams vector of integer index of specific samples to run from samplesheet
#' @param genCn also generate the CNV plot as a PNG file
#' @param skipQC default is false set as true to skip QC generation
#' @param email default is true, set to false to avoid email notification
#' @param cpReport default false, set TRUE to copy the reports to the Zdrive/research output Dir
#' @param redcapUp default is true, flag will upload output html files and dataframe to redcap
makeReports.v11b6 <- function(runPath = NULL,
                              sheetName = "samplesheet.csv",
                              selectSams = NULL,
                              genCn = F,
                              skipQC = F,
                              email = T,
                              cpReport = T,
                              redcapUp = T) {
    msgFunName(pipeLnk,"makeReports.v11b6")

    assign("genCn", genCn, envir = gb)
    data <- read.csv(sheetName, strip.white=T)
    runID <- paste0(data$RunID[1])

    if(file.exists(predictionPath)){
    message("\nLoading data...\n",predictionPath,"\n")
    load(predictionPath)
    }

    #CreateRedcapRecord(runID,"control")
    loopRender(selectSams, data, redcapUp)
    checkRunOutput(runID)
    if(skipQC == F){
        CreateRedcapRecord(runID)
        generateQCreport()
    }
    if(grepl("TEST",runID)){cpReport=F;redcapUp=F;email=F}
    if(cpReport==T){file.list <- gb$copy2outFolder(gb$clinDrv, runID)}
    if (redcapUp == T) {
        file.list <- dir(pattern = ".html", full.names = T)
        gb$uploadToRedcap(file.list, T)
    }
    if(email==T){launchEmailNotify(runID)}
    #beepr::beep(4)
    tidyUpFiles(runID)
}

# FUN: Checks if all the paths are accessible to the Rscript location
checkMounts <- function(){
    msgFunName(pipeLnk,"checkMounts")

    # List of three mount paths needed to run the pipleine
    critialMnts <- c("/Volumes/CBioinformatics/Methylation",
                     "/Volumes/molecular/",
                     "/Volumes/snudem01labspace/idats")
    researchMount = "smb://research-cifs.nyumc.org/Research/"
    molecularDrive = "smb://shares-cifs.nyumc.org/apps"
    failMount <-
        lapply(critialMnts, function(driveMount){ifelse(!dir.exists(driveMount), return(T), return(F))})
    if(any(failMount==T)){
        toFix <- paste(critialMnts[which(failMount==T)])
        cat("PATH does not exist, ensure network drive is mounted:", bkRed(toFix),"\n")
        cat(
            "You may need to re-mount one of the following paths:",
            bkBlu(researchMount,"CBioinformatics/"),
            bkgrn(researchMount,"snudem01lab/snudem01labspace"),
            bkBlu(molecularDrive,"/acc_pathology/molecular"), sep="\n"
            )
        stopifnot(!any(failMount==T))
    }
}

