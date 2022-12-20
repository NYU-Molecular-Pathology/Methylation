#!/usr/bin/env Rscript
gb <- globalenv(); assign("gb", gb)
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
dsh = "-----------"
bky <- function(txtVar,...){crayon::black$bgYellow$bold(txtVar,...)}
bkRed <- function(txtVar,...){crayon::bgRed$bold$white(txtVar,...)}
reportMd <- "~/report.Rmd" #system.file("report.Rmd", package="mnp.v11b6")
#/Volumes/CBioinformatics/Methylation/in_house/mnp.v116/mnp.v11b6/inst
QC_file <- "~/Methyl_QC.Rmd" #system.file('Methyl_QC.Rmd', package = "mnp.v11b6")
#/Volumes/CBioinformatics/Methylation/in_house/mnp.v116/mnp.v11b6/inst/

pipeLnk <- "https://github.com/NYU-Molecular-Pathology/Methylation/edit/main/pipelineHelper.R"
predictionPath <- "/Volumes/CBioinformatics/Methylation/in_house/mnp.v116/mnp.v11b6/data/rfpred.v11b6.RData"

msgFunName <- function(pthLnk, funNam){message("Executing function: ", crayon::black$bgYellow(funNam), " from RScript in:\n", pthLnk)}
msgParams <- function(...){message("Params passed: ", crayon::bgGreen(paste(..., sep = ",")))}

# Helper function to return the index of priority selected samples first
reOrderRun <- function(selectRDs, sh=NULL){
    msgFunName(pipeLnk,"reOrderRun")

    if(is.null(selectRDs)){return(NULL)}; if(is.null(sh)){sh<-"samplesheet.csv"}
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
generateQCreport <- function(runID=NULL, qc=NULL) {
    msgFunName(pipeLnk, "generateQCreport")
  runID <- gb$ckNull(nullVar = runID,
                     subVar = gb$runID,
                     deparse(substitute(runID, env = .GlobalEnv)))
  if (!file.exists(QC_file)) {
    message(crayon::bgRed("Check Working directory, QC_file.rmd not found"))
  }
  fs::file_copy(QC_file, getwd(), overwrite = T)
  currentQC = dir(getwd(), "*QC.Rmd", full.names = T)
  qcFile = paste0(runID, "_QC.html") # output file name
  if(file.exists(file.path(getwd(), qcFile))){
      message(qcFile, "Already Exists!  Skipping render...")
  }else{
    qcFile <- file.path(dirname(currentQC), qcFile)
    rmarkdown::render(
        currentQC,
        output_file = qcFile,
        quiet=T,
        params = list(runID = runID)
    )
  currentQC <- stringr::str_replace_all(string = currentQC, ".Rmd", "_cache")
  unlink(currentQC, recursive = T) #clear cache
      }
  gb$uploadToRedcap(qcFile,F)
}

# Sends an email notification that the run is complete from redcap admin
launchEmailNotify <-function(runID){
    msgFunName(pipeLnk,"launchEmailNotify")

    isMC = sjmisc::str_contains(runID, "MGDM") | sjmisc::str_contains(runID, "MC")
    rcon <- redcapAPI::redcapConnection(apiLink, gb$ApiToken)
    ur=paste0(rcon$url);tk=rcon$token
    sam_id = paste0(runID,"_QC")
    com <- ifelse(isMC==T, "sample_qc", "sample_research") # research or clinical notification
    record = data.frame(record_id=sam_id,comments=com)
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
    res<-RCurl::postForm(ur,token=tk,content='record',format='json',type='flat',data=datarecord)
    cat(crayon::black$bgGreen$bold("Email Notification Created"), sep="\n")
    cat(crayon::white$bgBlue$bold("Check email to confirm run notifcation"), sep="\n")
    cat(res)
}

# FUN: Creates the QC record for the current run on redcap if it does not exist
create.QC.record <- function(runID=NULL){
    msgFunName(pipeLnk,"create.QC.record")

    rcon <- redcapAPI::redcapConnection(apiLink, gb$ApiToken)
    uri=paste0(rcon$url); tk=rcon$token
    if(is.null(runID)){runID<-paste0(basename(gb$workDir))}
    record = c(record_id=paste0(runID,"_QC"), run_number=runID)
    qcdata <- jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
    rr<-RCurl::postForm(uri,token=tk,content='record',format='json', type='flat',data=qcdata)
    rr
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

msgProgress <- function(msg,i,samList){
    if(msg==1){
        cat("\n",bky(dsh, "Now Running", i, "of", length(samList), dsh),sep = "\n")
    }else{
        cat(bky("\n",dsh,"Completed Report",i,"of",length(samList),dsh),sep = "\n")
    }
}

# Helper function called by makeReports.v11b6 to generate the HTML report
do_report <-function(data = NULL, genCn=F) {
    msgFunName(pipeLnk,"do_report")

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
        sampleID = paste0(samplename_data)
        FFPE = NULL
        outDir = getwd()
        sample = 1
        outFi <- paste0(sampleID,".html")
        if(genCn==T){generateCNVpng(RGsetEpic,sampleID)}
        msgUpdate <- paste0("Current Sample:\n", "samplename_data=", samplename_data, " ",
                           "run_id=", run_id," ", "barcode=", barcode," ", "pathEpic:\n", pathEpic)
        message(bky(msgUpdate), "\n")
        tryCatch(
            expr={rmarkdown::render(reportMd, "html_document", outFi, outDir, quiet=FALSE,
                                   params = list(token=gb$ApiToken))},
            error=function(e){
                message(bkRed("Report Generation Failed:"),"\n", outFi,"\n")
                message(bkRed("The following error returned:"),"\n", e)
                }, finally=message("\nRunning next sample\n")
        )
    } else {message(bkRed("Your data is null, check your SampleSheet.csv"))}
}


# FUN: Iterates over each sample in the csv file to generate a report
loopRender <- function(samList = NULL, data,redcapUp=T) {
        msgFunName(pipeLnk, "loopRender")
        msgParams("samList = NULL, data")
        if (is.null(samList)) {
            samList <- 1:length(as.character(data$SentrixID_Pos))
        }
        require(rmarkdown)
        samSh <- dir(path = getwd(), full.names = T, ".xlsm")
    if(length(samSh)>1){
    warning("Multiple samplesheets found:\n")
    print(samSh)
    removeTemp <- stringr::str_detect(samSh,pattern = "\\$",negate = T)
    samSh <- samSh[removeTemp]
  }
        message("Reading the following .xlsm in current directory:", samSh)
        sampleNumb <- gb$getTotalSamples()
        sampleNumb = as.integer(sampleNumb)
        sh_Dat <- as.data.frame(
            readxl::read_excel(samSh, sheet = 3,range = "A1:N97",col_types = c("text")))[, 1:13]
        stopifnot(!is.null(sh_Dat))
        for (i in samList) {
            outFileN = paste0(data[i, 1], ".html")
            outPathN = file.path(gb$methDir, gb$runID, outFileN)
            message("outPathN = file.path(gb$methDir,gb$runID,outFileN)", "\n", outPathN)
            if (file.exists(outPathN)) {
                cat(bky(outFileN, "already exists! Skipping sample"),"\n")
                next
            } else {
                msgProgress(1,i,samList)
                do_report(data = data[i,], gb$genCn)
                msgProgress(2,i,samList)
                if(redcapUp==T){
                sh_Dat = sh_Dat[1:sampleNumb,]
                currSam <- sh_Dat[, 1] == data[i, 1]
                gb$importSingle(sh_Dat = sh_Dat[currSam, ])
                    }
            }
        }
        message(crayon::black$bgGreen$bold(dsh, "RUN COMPLETE", dsh))
        #beepr::beep(3)
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
makeReports.v11b6<-function(runPath=NULL,sheetName=NULL,selectSams=NULL,genCn=F,
                            skipQC=F,email=T,cpReport=T,redcapUp=T){
    msgFunName(pipeLnk,"makeReports.v11b6")

    assign("genCn",genCn, envir = gb)
    if (is.null(sheetName)) {sheetName="samplesheet.csv"}
    data <- read.csv(sheetName, strip.white=T)
    runID <- paste0(data$RunID[1])
    load(predictionPath)
    loopRender(selectSams, data, redcapUp)
    checkRunOutput(runID)
    if (skipQC == F) {
        if(redcapUp==T){
        create.QC.record(runID)
        }
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

#checkMounts()
#gb$defineParams()
