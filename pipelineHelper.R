#!/usr/bin/env Rscript
gb <- globalenv(); assign("gb", gb)
dsh = "-----------"
bky <- function(txtVar,...){crayon::black$bgYellow$bold(txtVar,...)}
bkRed <- function(txtVar,...){crayon::bgRed$bold$white(txtVar,...)}

# Global Variables ----------------------------------
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
reportMd <- "~/report.Rmd" # From curl github download
QC_file <- "~/Methyl_QC.Rmd" # From curl github download
pipeLnk <- "https://github.com/NYU-Molecular-Pathology/Methylation/edit/main/pipelineHelper.R"
predictionPath <- "/Volumes/CBioinformatics/Methylation/in_house/mnp.v116/mnp.v11b6/data/rfpred.v11b6.RData"

msgFunName <- function(pthLnk, funNam){message("Executing function: ", bky(funNam), " from RScript in:\n", pthLnk)}
msgParams <- function(...){cat("\nParams passed: ", crayon::bgGreen(paste(..., sep = " , ")))}

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

# QC REPORT maker: knits the QC RMD file
generateQCreport <- function(runID=NULL) {
    msgFunName(pipeLnk, "generateQCreport")

    if (is.null(runID)){runID<-paste0(basename(gb$workDir))}
    if (!file.exists(QC_file)) {
        message(crayon::bgRed("QC_file.rmd not found:"), "\n", QC_file)
    }
    fs::file_copy(QC_file, getwd(), overwrite = T)
    currQc = dir(getwd(), "*QC.Rmd", full.names = T)
    qcFile = paste0(runID, "_QC.html") # output file name
    if(file.exists(file.path(getwd(), qcFile))) {
        message(qcFile, "Already Exists!  Skipping render...")
    } else{
        qcFile <- file.path(dirname(currQc), qcFile)
        rmarkdown::render(
            currQc,
            output_file = qcFile,
            quiet = T,
            params = list(runID = runID)
        )
        currQc <- stringr::str_replace_all(string = currQc, ".Rmd", "_cache")
        unlink(currQc, recursive = T) #clear cache
    }
    gb$uploadToRedcap(qcFile,F)
}

# Sends an email notification that the run is complete from redcap admin
launchEmailNotify <-function(runID){
    msgFunName(pipeLnk,"launchEmailNotify")
    isMC = sjmisc::str_contains(runID, "MGDM") | sjmisc::str_contains(runID, "MC")
    com <- ifelse(isMC==T, "sample_qc", "sample_research") # research or clinical notification
    record = data.frame(record_id=paste0(runID,"_QC"),comments=com)
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
    res<-RCurl::postForm(gb$apiLink,token=gb$ApiToken,content='record',format='json',type='flat',data=datarecord)
    cat(crayon::white$bgBlue$bold("Check email to confirm Email Notification Created"), sep="\n")
    cat(res)
}

# FUN: Creates the QC record for the current run on redcap if it does not exist
create.QC.record <- function(runID=NULL){
    msgFunName(pipeLnk,"create.QC.record")

    if(is.null(runID)){runID<-paste0(basename(gb$workDir))}
    record = c(record_id=paste0(runID,"_QC"), run_number=runID)
    qcdata <- jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
    RCurl::postForm(msgParams,token=gb$ApiToken,content='record',format='json', type='flat',data=qcdata)
    message(dsh,"Created QC Record",dsh)
    print(qcdata)
}

# Check if the QC File will be read
checkRunOutput <- function(runID) {
    msgFunName(pipeLnk,"checkRunOutput")
    csvLocation <- file.path(fs::path_home(),"Desktop",runID,paste0(runID,"_Redcap.csv"))
    if (!file.exists(csvLocation)) {
        message(bkRed("File not found:")," ", csvLocation,"\n",bkRed("QC Summary Table will not Knit"),"\n")
    } else{
        fs::file_copy(path=csvLocation, new_path=getwd(), overwrite = T)
    }
}

# gets rid of desktop files if run is successful
tidyUpFiles <- function(runID){
    msgFunName(pipeLnk,"tidyUpFiles")
    deskDir <- file.path("~/Desktop",runID)
    backupD <- file.path(gb$methDir,"csvRedcap")
    if(!dir.exists(backupD)){dir.create(backupD)}
    file.copy(deskDir, backupD,overwrite=T, recursive = T, copy.mode = T)
    unlink(deskDir,T,T)
}

msgProgress <- function(msg, i, samList) {
    nOfTotal <- paste(i, "of", length(samList), dsh)
    startEnd <- ifelse(msg == 1, "Now Running", "Completed Report")
    cat("\n", bky(dsh, startEnd , nOfTotal), sep = "\n")
}

msgRunUp <- function(sampleID,run_id,senLi){
    message(
        "Current Sample:\nsamplename_data=",
        bky(sampleID),"\nrun_id=",bky(run_id),
        "\npathEpic:\n", bky(file.path(getwd(), senLi)), "\n"
        )
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

# Helper function called by makeReports.v11b6 to generate the HTML report
do_report <-function(data = NULL, genCn=F) {
    msgFunName(pipeLnk,"do_report")
    msgParams("data\n", data)
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
                message(bkRed("Report Generation Failed:"),"\n", outFi,"\n")
                message(bkRed("The following error returned:"),"\n", e)
            }, finally=message("\nRunning next sample\n")
        )
    } else {message(bkRed("Data is NULL, check your SampleSheet.csv"))}
}


msgSamSheet <- function(samSh) {
    if (length(samSh) > 1) {
        message(crayon::black$bgRed$bold("Multiple samplesheets found:"))
        cat(samSh, sep = "\n")
        samSh <- samSh[stringr::str_detect(samSh, pattern = "\\$", negate = T)]
    } else{message("Reading the following .xlsm in current directory:")}
    cat(samSh, sep = "\n")
    stopifnot(!length(samSh) > 1)
    return(samSh)
}

# FUN: Parses the WetLab .xlsm sheet in the current directory
checkSamSh <- function(samList){
    msgFunName(pipeLnk, "checkSamSh")
    require(rmarkdown)
    samSh <- msgSamSheet(dir(path = getwd(), ".xlsm", full.names = T))
    shNames <- which(grepl("REDCap_Import",readxl::excel_sheets(samSh)))
    samSh <- readxl::read_excel(samSh, sheet=shNames, range = "A1:M97",col_types = c("text"))
    wksh <- as.data.frame(samSh)[1:length(samList), 1:13]
    rownames(wksh)<- wksh$Sample_Name
    stopifnot(!is.null(wksh))
    return(wksh)
}

getRunList <- function(data, samList){
    msgFunName(pipeLnk, "getRunList")
    toRun <- unlist(lapply(samList, FUN=function(i){
        nOutDir = file.path(gb$methDir, gb$runID, paste0(data[i, 1], ".html"))
        if (file.exists(nOutDir)) {
            message(bky(basename(nOutDir), "already exists! Skipping..."),"\n")
        }else{return(i)}
    }))
    return(toRun)
}


# FUN: Iterates over each sample in the csv file to generate a report
loopRender <- function(samList = NULL, data, redcapUp = T){
    msgFunName(pipeLnk, "loopRender")
    msgParams("samList = NULL, data, redcapUp = T")

    stopifnot(!is.null(data))
    if (is.null(samList)) {samList = 1:length(data$SentrixID_Pos)}
    wksh <- checkSamSh(samList)
    toRun <- getRunList(data, samList)
    for (i in toRun) {
        msgProgress(1, i, samList)
        do_report(data = data[i, ], gb$genCn)
        msgProgress(2, i, samList)
        if (redcapUp == T) {gb$importSingle(wksh[data[i, 1],])}
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
makeReports.v11b6<-function(runPath=NULL,sheetName="samplesheet.csv",selectSams=NULL,genCn=F,
                            skipQC=F,email=T,cpReport=T,redcapUp=T){
    msgFunName(pipeLnk,"makeReports.v11b6")

    assign("genCn",genCn, envir = gb)

    data <- read.csv(sheetName, strip.white=T)
    runID <- paste0(data$RunID[1])
    load(predictionPath)
    loopRender(selectSams, data, redcapUp)
    checkRunOutput(runID)
    if(skipQC == F){
        create.QC.record(runID)
        generateQCreport()
    }
    if(grepl("TEST",runID)){cpReport=F;redcapUp=F;email=F}
    if(cpReport==T){file.list <- gb$copy2outFolder(gb$clinDrv, runID)}
    if(redcapUp==T){file.list <- dir(pattern="*.html", full.names = T); gb$uploadToRedcap(file.list,T)}
    if(email==T){
        launchEmailNotify(runID)
        #beepr::beep(4)
    }
    tidyUpFiles(runID)
}

# Function to just run a default clinical run without changes, input selectRDs to prioritize samples running first
startRun <- function(selectRDs=NULL, emailNotify=T){
    msgFunName(pipeLnk,"startRun")
    msgParams("selectRDs=NULL, emailNotify=T")
    msgParams(selectRDs,emailNotify)

    if(!is.null(selectRDs)){
        sampleOrder <- reOrderRun(selectRDs) # Re-order sample report generation for priority
        makeReports.v11b6(skipQC=F, email=T, cpReport=T, selectSams=sampleOrder, redcapUp=T)
    } else {
        makeReports.v11b6(skipQC=F, email=T, cpReport=T, selectSams=NULL, redcapUp=T)
    }
}

# FUN: Checks if all the paths are accessible to the Rscript location
checkMounts <- function(){
    msgFunName(pipeLnk,"checkMounts")

    # List of three mount paths needed to run the pipleine
    critialMnts <- c("/Volumes/CBioinformatics/Methylation",
                     "/Volumes/molecular/MOLECULAR LAB ONLY", "/Volumes/snudem01labspace/idats")
    failMount <- lapply(critialMnts, function(driveMount){
        ifelse(!dir.exists(driveMount),return(T),return(F))})
    if(any(failMount==T)){
        toFix <- paste(critialMnts[which(failMount==T)])
        cat("PATH does not exist, ensure path is mounted:", bkRed(toFix),"\n")
        cat("You may need to mount one of the following paths:",
            crayon::white$bgBlue$bold(
                "smb://research-cifs.nyumc.org/Research/CBioinformatics/"),
            crayon::white$bgCyan$bold(
                "smb://research-cifs.nyumc.org/Research/snudem01lab/snudem01labspace"),
            crayon::white$bgGreen$bold(
                "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular"), sep="\n"
        )
        stopifnot(!any(failMount==T))
    }
}

# Executes the functions in order to setup a run
prepareRun <- function(token, baseFolder=NULL){
    msgFunName(pipeLnk,"prepareRun")

    runValid <- gb$checkValidRun(gb$runID)
    message("Is the runID valid? ", runValid)
    if(!runValid){
        message(crayon::bgRed$white$bold("runID",gb$runID,"is not valid"))
        message(crayon::bgBlue$white$bold(paste0(gb$runID,".xlsm"),"not found in worksheets folder"))
        stopifnot(runValid)
    }
    if(is.null(baseFolder)){baseFolder <- "/Volumes/CBioinformatics/Methylation/Clinical_Runs"}
    gb$methDir <- baseFolder
    methylPath <- gb$setRunDir(gb$runID, baseFolder)
    message("\n","Working directory set to:","\n")
    cat(crayon::bgGreen(methylPath))
    setwd(methylPath)
    gb$setVar("ApiToken", token) # assign the ApiToken & print params
    gb$copyWorksheetFile(runID = gb$runID) # copies the xlsm file
    gb$readSheetWrite() # reads xlsm and generates input .csv samplesheet
    gb$get.idats() # Copy idat files to current folder from molecular and snuderlabspace to cwd
    gb$moveSampleSheet(gb$methDir) #copies outputs temp to desktop for QC.Rmd
    #gb$classifierInstall(instNew = F, rmpkg = F) # Loads pipeline or installs new
}

checkMounts()
#gb$defineParams()
