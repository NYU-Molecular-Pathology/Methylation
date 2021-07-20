#!/usr/bin/env Rscript
gb <- globalenv(); assign("gb", gb)
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"

# Helper function to return the index of priority selected samples first
reOrderRun <- function(selectRDs, sh=NULL){
    if(is.null(selectRDs)){return(NULL)}; if(is.null(sh)){sh<-"samplesheet.csv"}
    allRd <- as.data.frame(read.csv(sh))
    runFirst <- which(allRd[,1] %in% selectRDs)
    runAfter <- which(!(allRd[,1] %in% selectRDs))
    return(c(runFirst,runAfter))
}

# Saves the methyl CNV as a png file in the cwd
generateCNVpng <- function(RGsetEpic, sampleName) {
    imgName <- paste(sampleName, "cnv.png", sep="_")
    MsetEpic <- minfi::preprocessRaw(RGsetEpic)
    png(filename=imgName,width=1820, height=1040)
    mnp.v11b6::MNPcnvplot(MsetEpic)
    dev.off(); message("File saved:\n",imgName,"\n")
}

getRGset <- function(runPath, sentrix){
    barcode = stringr::str_split_fixed(sentrix, "_",2)[1]
    RGsetEpic <- minfi::read.metharray(file.path(runPath, sentrix), verbose = T, force = T)
    aEpic=c(array="IlluminaHumanMethylationEPIC", annotation="ilm10b4.hg19")
    a450k=c(array="IlluminaHumanMethylation450k", annotation="ilmn12.hg19")
    if (barcode >= as.numeric("204220033000")) {RGsetEpic@annotation=aEpic}
    if (RGsetEpic@annotation['array']=="IlluminaHumanMethylation450k") {RGsetEpic@annotation=a450k}
    return(RGsetEpic)
}

# Helper function called by makeReports.v11b6 to generate the HTML report
do_report <-function(data = NULL, genCn=F) {
    require(rmarkdown)
    if(!is.null(data)){
        samplename_data = paste0(data[,1])
        sentrix_pos_list = (data[,5])
        tech = paste0(data[,9])
        tech2 = paste0(data[,10])
        mp_number = paste0(data[,8])
        run_id = paste0(data[,7])
        bnumber = paste0(data[,2])
        runPath = getwd()
        barcode = as.numeric(data[,3])
        pathEpic = file.path(runPath, sentrix_pos_list)
        RGsetEpic <- getRGset(runPath, sentrix_pos_list)
        RGset = RGsetEpic
        sampleID=paste0(samplename_data)
        FFPE = NULL
        output_dir = getwd()
        sample = 1
        output_file <- paste0(sampleID,".html")
        if(genCn==T){generateCNVpng(RGsetEpic,sampleID)}
        message(paste("Now running:", samplename_data, run_id, barcode, pathEpic,"\n"))
        rmarkdown::render(
            system.file("report.Rmd", package="mnp.v11b6"),
            output_dir = output_dir,
            output_format = "html_document",
            output_file = output_file,
            clean = T
        )
    } else {message("your data is null")}
}

# QC REPORT maker: knits the QC RMD file
generateQCreport <- function(runID=NULL, qc=NULL) {
    runID<-gb$ckNull(nullVar = runID, subVar=gb$runID, deparse(substitute(runID,env=.GlobalEnv)))
    QC_file <- system.file('Methyl_QC.Rmd', package = "mnp.v11b6")
    if (!file.exists(QC_file)){message("Check Working directory, QC_file.rmd not found")}
    fs::file_copy(QC_file, getwd(), overwrite = T)
    currentQC = dir(getwd(),"*QC.Rmd", full.names=T)
    qcFile = paste0(runID,"_QC.html") # output file name
    rmarkdown::render(currentQC, output_file=file.path(dirname(currentQC), qcFile), params=list(runID=runID))
    currentQC <- stringr::str_replace_all(string=currentQC, ".Rmd", "_cache")
    unlink(currentQC, recursive=T) #clear cache
}

# Sends an email notification that the run is complete from redcap admin
launchEmailNotify <-function(runID){
    isMC = sjmisc::str_contains(runID, "MGDM") | sjmisc::str_contains(runID, "MC")
    rcon <- redcapAPI::redcapConnection(apiLink, gb$ApiToken)
    ur=paste0(rcon$url);tk=rcon$token
    sam_id = paste0(runID,"_QC")
    com <- ifelse(isMC==T, "sample_qc", "sample_research") # research or clinical notification
    record = data.frame(record_id=sam_id,comments=com)
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T);print(datarecord)
    res<-RCurl::postForm(ur,token=tk,content='record',format='json',type='flat',data=datarecord)
    cat(crayon::black$bgGreen$bold("Email Notification Created\n"))
    cat(crayon::white$bgBlue$bold("Check email to confirm run notifcation\n"));cat(res)
}

# FUN: Creates the QC record for the current run on redcap if it does not exist
create.QC.record <- function(runID=NULL){
    rcon <- redcapAPI::redcapConnection(apiLink, gb$ApiToken)
    uri=paste0(rcon$url); tk=rcon$token
    if(is.null(runID)){runID<-paste0(basename(gb$workDir))}
    record = c(record_id=paste0(runID,"_QC"), run_number=runID)
    qcdata <- jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
    rr<-RCurl::postForm(uri,token=tk,content='record',format='json', type='flat',data=qcdata)
    rr;message("Created QC Record"); print(qcdata)
}

# Check if the QC File will be read
checkRunOutput <- function(runID) {
    location <- file.path("~/Desktop",runID,paste0(runID,"_Redcap.csv"))
    if (!file.exists(location)) {
        message("File not found: "); cat(location,sep="\n"); message("QC Summary Table will not Knit")
    } else {fs::file_copy(path=location, new_path=getwd(), overwrite = T)}
}

#' REPORT: Generates Html reports to cwd with samplesheet.csv
#' @param runPath The location of samplesheet.csv and idats
#' @param sheetName name of samplesheet if it is not "samplesheet.csv"
#' @param selectSams vector of integer index of specific samples to run from samplesheet
#' @param genCn also generate the CNV plot as a PNG file
#' @param skipQC default is false set as true to skip QC generation
#' @param email default is true, set to false to avoid email notification
#' @param cpReport default is false, set to true to copy the reports to the Zdrive or research output directories
#' @param redcapUp default is true, flag will upload output html files and dataframe to redcap
makeReports.v11b6<-function(runPath=NULL,sheetName=NULL,selectSams=NULL,genCn=F,skipQC=F,email=T,cpReport=F,redcapUp=T){
    dsh <- "-----------"
    bky <- function(txtVar,...){crayon::black$bgYellow$bold(txtVar,...)}
    if (is.null(runPath)) {runPath=gb$workDir}
    if (is.null(sheetName)) {sheetName="samplesheet.csv"}
    data <- read.csv(sheetName, strip.white=T)
    runID<-paste0(data$RunID[1])
    normList <- 1:length(as.character(data$SentrixID_Pos))
    if(is.null(selectSams)){samList <-normList} else {samList <-selectSams}
    for (i in samList) {
        outFileN = paste0(data[i,1],".html")
        outPathN = file.path(gb$workFolder,gb$runID,outFileN)
        if(file.exists(outPathN)){
            cat(bky(outFileN, "exists skipping sample","\n"))
            next
        } else {
            cat(bky("\n",dsh,"Now Running", i, "of", length(samList),dsh),sep="\n")
            do_report(data=data[i, ], genCn)
        }
    }
    cat(crayon::black$bgGreen$bold(dsh,"RUN COMPLETE",dsh),sep="\n")
    checkRunOutput(runID)
    if(skipQC==F){
        create.QC.record(runID)
        generateQCreport()
         # creates a redcap QC record and Knits the QC RMD file
    if(cpReport==T){file.list <- gb$copy2outFolder(gb$clinDrv, runID)}
    if(redcapUp==T){file.list <- dir(pattern="*.html", full.names = T); gb$uploadToRedcap(file.list)}
    if(email==T){launchEmailNotify(runID)}
    deskDir <- file.path("~/Desktop",runID)
    unlink(deskDir,T,T)
    }
}

# Function to just run a default clinical run without changes, input selectRDs to prioritize samples running first
startRun <- function(selectRDs=NULL, runID=NULL, emailNotify=T){
    if(!is.null(selectRDs)){
        sampleOrder <- reOrderRun(selectRDs) # Re-order sample report generation for priority samples first
        makeReports.v11b6(skipQC=F, email=emailNotify, cpReport=T, selectSams=sampleOrder, redcapUp=T)
    } else {makeReports.v11b6(skipQC=F, email=emailNotify, cpReport=F, selectSams=NULL, redcapUp=T)}
}

# FUN: Checks if all the paths are accessible to the Rscript location
checkMounts <- function(){
    # List of three mount paths needed to run the pipleine
    critialMnts <- c("/Volumes/CBioinformatics/jonathan",
                     "/Volumes/molecular/MOLECULAR LAB ONLY", "/Volumes/snudem01labspace/idats")
    failMount <- lapply(critialMnts, function(driveMount){
        ifelse(!dir.exists(driveMount),return(T),return(F))})
    if(any(failMount==T)){
        toFix <- paste(critialMnts[which(failMount==T)])
        cat("PATH does not exist, ensure path is mounted:", crayon::white$bgRed$bold(toFix),"\n")
        cat("You must mount each of the following paths:\n",
            crayon::white$bgGreen$bold(
                "smb://research-cifs.nyumc.org/Research/CBioinformatics/\n",
                "smb://research-cifs.nyumc.org/Research/snudem01lab/snudem01labspace\n",
                "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular\n")
        )
        stopifnot(!any(failMount==T))
    }
}

# Executes the functions in order to setup a run
prepareRun <- function(token){
    runValid <- gb$checkValidRun(gb$runID)
    message("Is the runID valid? ", runValid)
    if(!runValid){message("runID ",gb$runID," is not valid"); stopifnot(runValid)}
    methylPath <- gb$setRunDir(gb$runID)
    message("Working directory set to:"); cat(crayon::bgGreen(methylPath)); setwd(methylPath)
    gb$setVar("ApiToken", token) # assign the ApiToken & print params
    gb$copyWorksheetFile(runID = gb$runID) # copies the xlsm file
    gb$readSheetWrite() # reads xlsm and generates input .csv samplesheet
    gb$get.idats() # Copy idat files to current folder from molecular and snuderlabspace to cwd
    gb$moveSampleSheet(gb$methDir) #copies outputs temp to desktop for QC.Rmd
    gb$install.or.load(instNew = F, rmpkg = F) # Loads pipeline or installs new
}

checkMounts()
gb$defineParams()
