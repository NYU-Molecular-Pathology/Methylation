#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library("base"); gb <- globalenv(); assign("gb", gb)

if(!require("devtools")){install.packages("devtools")}
if(!require("minfi")){BiocManager::install("minfi")}

# If you have issues with conda and r, try running the following:
# conda create -n r_env r-essentials r-base
# conda activate r_env

# Main Parameters trailing commandline
token<-args[1]; runID<-args[2]; runPath<-args[3]

apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
dsh = "-----------"
bky <- function(txtVar,...){crayon::black$bgYellow$bold(txtVar,...)}
bkRed <- function(txtVar,...){crayon::bgRed$bold$white(txtVar,...)}

grabYear<- function(yr) {
    rnum <- NULL
    if(nchar(yr)>2){rnum <- substring(yr, 3)}else{rnum <- yr}
    if(nchar(yr)>0){rnum <- paste0("20",rnum)}else{rnum}
    return(rnum)
}

setVar <- function(valueName,val){return(assign(valueName, val, envir=.GlobalEnv))}
assignVar <- function(varStr, assignedVal){
    return(
        tryCatch(
            expr = {if(!is.null(get(varStr))){cat("")}},
            error = {setVar(varStr,assignedVal)},
            finally = message(varStr," = ",assignedVal)
            )
        )
}

# Checks if variable is null and assigns value
ckNull <- function(nullVar, subVar, varName){
    if (is.null(nullVar)){
        setVar(as.character(varName),as.character(subVar))
        return(paste0(subVar))} else {return(paste0(nullVar))}
}


# FUN: copies the molecular or research lab Worksheet xlsm to cwd
copyWorksheetFile <- function(runID=NULL, runYear=NULL) {
    if (is.null(runID)){runID=paste0(basename(getwd()))} else {runID=runID}
    stopifnot(!is.null(runID))
    #if (is.null(runYear)){
    yr <- stringr::str_split_fixed(runID,"-",2)[,1]
    runYear=grabYear(yr)
    #   }
    isMC = sjmisc::str_contains(runID, "MGDM")|sjmisc::str_contains(runID, "MC")
    assign("isMC", isMC); message("\nIs methylation run Clinical? ",isMC)
    rschDrv <- "/Volumes/snudem01labspace/Methylation_Worksheets"
    mountLoc <- ifelse(isMC,file.path(gb$clinDrv,"WORKSHEETS"),rschDrv)
    message("\nCopying file from: "); cat(mountLoc,"\n")
    if (dir.exists(mountLoc)) {
        fileLoc <- file.path(mountLoc,runYear,paste0(runID,".xlsm"))
        if (file.exists(fileLoc)) {
            message("\n\nCopying worksheet from Worksheets Folder:")
            if(!file.exists(file.path(getwd(),basename(fileLoc)))){
                fs::file_copy(fileLoc, getwd(),overwrite=T)}
        } else {
            message("\nFile not found:\n",fileLoc)
            message("\nMake sure your path is correct, try print(gb$copyWorksheetFile)\n")
        }
    } else {
        message("\nFolder not found:\n",mountLoc)
        message("\nEnsure the path to this network Volume is correct\n")
    }
}

# Sets the methylation run directory named by the new run name
setRunDir <- function(runID=NULL, workFolder=NULL){
    runID <- gb$ckNull(runID, paste0(basename(getwd())), deparse(substitute(runID,env=gb)))
    workFolder <- gb$ckNull(workFolder, gb$methDir, deparse(substitute(workFolder,env=gb)))
    newRun <- file.path(workFolder, runID)
    assign("newRunPath", newRun)
    if(grepl("TEST",runID)){
        if(dir.exists(newRun)){unlink(newRun, T, T)}
        dir.create(newRun)
        try(unlink(file.path("~/Desktop",runID), T, T),silent = T)
    }
    if(!dir.exists(newRun)){
        base::dir.create(newRun, mode = "777", recursive=T)
        Sys.chmod(newRun, mode = "777")
        cat("Creating folder: ",newRun)
        gb$setDirectory(newRun)
    } else {gb$setDirectory(newRun)}
    return(newRun)
}

# Returns Total Sample Count in the run
getTotalSamples <- function(){
    templateDir = "Clinical_Methylation/methylation_run_TEMPLATE_new.xlsm"
    thisSh <- dir(getwd(), "*.xlsm")
    temp <- stringi::stri_detect_fixed(thisSh, "~$")
    thisSh <- thisSh[!temp]
    if(length(thisSh)==0){
        message("No .xlsm sheet, defaulting to NONE")
                          print(getwd())
                          return("NONE")
                         }
    worksheet <- readxl::read_excel(thisSh[1], col_names="Total", range="B4:B4")
    if (length(worksheet) == 0) {
        message("Samplesheet ", thisSh[1]," is invalid format, manually edit")
        message("Try copying the template:\n", templateDir)
    } else {message("\nTotal sample count found is: ", worksheet[1])}
    return(paste0(worksheet[1]))
}

# FUN: translates the xlsm excel file to the .csv samplesheet for redcap and minfi
writeSampleSheet <- function(df, samplesheet_ID, bn = NULL, sampleName, dnaNumber, Sentrix) {
    if (is.null(bn)) {bn = file.path(getwd(), df$Sentrix)}
    samplesheet_csv = data.frame(
        Sample_Name =  df[, sampleName],
        DNA_Number =df[, dnaNumber],
        Sentrix_ID = samplesheet_ID[, 1],
        Sentrix_Position = samplesheet_ID[, 2],
        SentrixID_Pos = df[, Sentrix],
        Basename = paste0(bn),
        RunID = df$Batch,
        MP_num = df$MP_number,
        tech = df$Tech,
        tech2 = df$Tech2,
        Date = df$Date
    )
    write.csv(samplesheet_csv, file = "samplesheet.csv", quote = F, row.names = F)
}

# FUN: Reads the csv samplesheet for minfi input
readSampleSheet <- function(runID=F, totalSam=F, wks=F) {
    file.list <- dir(path=getwd(), "*.xlsm")
    temps <- stringi::stri_detect_fixed(file.list, "~$")
    file.list <- file.list[!temps]
    sampleSheet <- paste0(file.list[1])
    message(paste0("Reading worksheet named: ", sampleSheet))
    worksheet <- readxl::read_excel(sampleSheet, sheet=2, col_names=T, col_types="text", trim_ws=T)
    wsDate <-  as.data.frame(readxl::read_excel(sampleSheet, sheet=1, col_names=F, range="F4:F4", trim_ws=T))[1]
    names(wsDate)="Date"
    colnames(worksheet)
    worksheet$Date <- paste0(wsDate$Date[1])
    if (runID == T) {return(worksheet$Project[1])}
    if (totalSam == T){
        runNum <- readxl::read_excel(sampleSheet, sheet=1, col_names=F, range="B4")
        runNum <- as.numeric(runNum)
        return(runNum)
    }
    if (wks == T) {return(worksheet)}
}

# FUN: reads the .xlsm worksheet and outputs the .csv methyl experiment for MINFI
readSheetWrite <- function(sampleNumb= NULL, runID = NULL) {
    if(is.null(sampleNumb)){sampleNumb<-getTotalSamples()}
    if(is.null(runID)){runID<-paste0(basename(getwd()))}
    current.run.Folder <- getwd()
    message("\nCurrent Run folder:", current.run.Folder,"\n")
    ww1=crayon::bgRed("No tech name found: check df$Tech in samplesheet")
    ww2=crayon::bgRed("No tech name found: check df$MP_number in samplesheet")
    if (!file.exists("samplesheet.csv")) {
        sampleNumb=as.integer(sampleNumb)
        worksheet=readSampleSheet(wks=T)
        hdrs = dimnames(worksheet)
        print(data.frame(HEADERS=hdrs[[2]]))
        sampleName=paste(hdrs[[2]][9])
        Sentrix=paste(hdrs[[2]][1])
        dnaNumber=paste(hdrs[[2]][10])
        mpnum=paste(hdrs[[2]][8])
        Var.names=c("~~~~~~~",sampleName,Sentrix, dnaNumber, mpnum);cat("\n\n")
        print(data.frame(Var.names))
        df=as.data.frame(worksheet)[1:sampleNumb, ]
        samplesheet_ID=as.data.frame(stringr::str_split_fixed(df[, Sentrix], "_", 2))
        bn <- file.path(current.run.Folder, df[, Sentrix])
        message("Basename layout:\n",bn[1])
        note=paste(df$Notes[1])
        df$Notes <- note
        if (is.null(df$Tech)){warning(ww1); df$Tech <- "NA"}
        if (is.null(df$MP_number)){warning(ww2); df$MP_number <- "none"}
        writeSampleSheet(df, samplesheet_ID=samplesheet_ID,bn=bn, sampleName, dnaNumber, Sentrix)
    } else {
        msgCvs = "samplesheet.csv already exists! To Create a new csv file, Delete the existing samplesheet.csv"
        message("\n", crayon::white$bgGreen(msgCvs),"\n")
    }
}

# FUN: Returns a list of idat files given an idat drive location -
getAllFiles <- function(idatDir, csvNam=NULL) {
    if(!is.null(csvNam)){ssheet=read.csv(csvNam,strip.white=T)
    barcode=as.vector(ssheet$Sentrix_ID)} else {
        ssheet=read.csv(csvNam,strip.white=T)
        barcode=as.vector(ssheet$Sentrix_ID)}
    allFi=NULL
    senP=ssheet$SentrixID_Pos
    for (i in 1:length(idatDir)){
        dirNames = file.path(idatDir[i], barcode)
        green.files <- file.path(dirNames, paste0(senP, "_Grn.idat"))
        red.files <- file.path(dirNames, paste0(senP, "_Red.idat"))
        allFi <- c(allFi, green.files,red.files)
    }
    return(allFi)
}

# FUN: Copies samplesheet to Desktop folder
moveSampleSheet <- function(methDir, runID=NULL) {
    if (is.null(runID)){runID=basename(getwd())}
    baseFolder=getwd()
    if(!dir.exists(baseFolder)){dir.create(baseFolder)}
    currDir=methDir
    endDir = file.path(baseFolder,paste0(runID,"_samplesheet.csv"))
    #fs::file_copy(path=paste0(currDir,"/samplesheet.csv"),new_path=baseFolder,overwrite=T)
    #file.rename(from=file.path(baseFolder,"samplesheet.csv"), to=endDir)
}

#  Copy idats and Worksheets creation
writeFromRedcap <- function(df, samplesheet_ID, bn = NULL) {
    if (is.null(bn)) {bn = file.path(getwd(), df$barcode_and_row_column)}
    message("~~~Writing from redcap samplesheet.csv:")
    names(df)
    samplesheet_csv = data.frame(
        Sample_Name = df[, "record_id"],
        DNA_Number = df[,"b_number"],
        Sentrix_ID = samplesheet_ID[, 1],
        Sentrix_Position = samplesheet_ID[, 2],
        SentrixID_Pos = df[, "barcode_and_row_column"],
        Basename = paste0(bn),
        RunID = df$run_number,
        MP_num = df$accession_number,
        tech = df$primary_tech,
        tech2 = df$second_tech,
        Date = df$arrived
    )
    print(samplesheet_csv)
    write.csv(samplesheet_csv,file = "samplesheet.csv",quote = F,row.names = F)
}

#' FUN: Returns dataframe of redcap search using default worksheet header and fields
search.redcap <- function(rd_numbers, token=NULL, flds=NULL) {
    if(!require("redcapAPI")){install.packages("redcapAPI", dependencies = T, type="both",ask=F);library("redcapAPI")}
    if(is.null(token)){message("You must provide an ApiToken!")};stopifnot(!is.null(token))
    rcon <- redcapAPI::redcapConnection(gb$apiLink, token)
    if (is.null(flds)){flds = c("record_id","b_number","primary_tech","second_tech","run_number",
                                "barcode_and_row_column","accession_number","arrived")}
    result <- redcapAPI::exportRecords(rcon,records = rd_numbers,fields = flds,dag = F,factors = F,
                                       labels = F,dates = F, form_complete_auto = F,format = 'csv')
    return(as.data.frame(result))
}

# FUN: Copies .idat files to your directory and saves samplesheet.csv
get.rd.info <- function(rd_numbers=NULL, token=NULL, sh_name=NULL){
    if (is.null(rd_numbers)){message("Input RD-numbers using get.rd.info(rd_numbers)")}
    if (is.null(sh_name)) {sh_name = "samplesheet.csv"}
    if(is.null(token)){message("You must provide an ApiToken!")};stopifnot(!is.null(token))
    result <- search.redcap(rd_numbers, token, NULL)
    samplesheet_ID = as.data.frame(stringr::str_split_fixed(result[,"barcode_and_row_column"],"_",2))
    writeFromRedcap(result, samplesheet_ID) # writes API export as minfi dataframe sheet
    get.idats(csvNam = sh_name)  # copies idat files from return to current directory
    return(result)
}

# Helper function to return the index of priority selected samples first
reOrderRun <- function(selectRDs, sh=NULL){
    if(is.null(selectRDs)){return(NULL)}; if(is.null(sh)){sh<-"samplesheet.csv"}
    allRd <- as.data.frame(read.csv(sh))
    runFirst <- which(allRd[,1] %in% selectRDs)
    runAfter <- which(!(allRd[,1] %in% selectRDs))
    return(c(runFirst,runAfter))
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
    reportMd <- "report.Rmd"
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
        msgUpdate <- paste0("Current Sample:\n", "samplename_data=", samplename_data, " ", "run_id=", run_id," ", "barcode=", barcode," ", "pathEpic:\n", pathEpic)
        message(bky(msgUpdate), "\n")
        tryCatch(
            expr={rmarkdown::render(reportMd, "html_document", outFi, outDir)},
            error=function(e){
                message(bkRed("Report Generation Failed:"),"\n", outFi,"\n")
                message(bkRed("The following error returned:"),"\n", e)
            }, finally=message("\nRunning next sample\n")
        )
    } else {message(bkRed("Your data is null, check your SampleSheet.csv"))}
}

# QC REPORT maker: knits the QC RMD file
generateQCreport <- function(runID=NULL, qc=NULL) {
    QC_file <- 'Methyl_QC.Rmd' #system.file('Methyl_QC.Rmd', package = "mnp.v11b6")
    runID<-gb$ckNull(nullVar = runID, subVar=gb$runID, deparse(substitute(runID,env=.GlobalEnv)))
    if (!file.exists(QC_file)){message(crayon::bgRed("Check Working directory, QC_file.rmd not found"))}
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
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
    res<-RCurl::postForm(ur,token=tk,content='record',format='json',type='flat',data=datarecord)
    cat(crayon::black$bgGreen$bold("Email Notification Created"), sep="\n")
    cat(crayon::white$bgBlue$bold("Check email to confirm run notifcation"), sep="\n")
    cat(res)
}

# FUN: Creates the QC record for the current run on redcap if it does not exist
create.QC.record <- function(runID=NULL){
    rcon <- redcapAPI::redcapConnection(apiLink, gb$ApiToken)
    uri=paste0(rcon$url); tk=rcon$token
    if(is.null(runID)){runID<-paste0(basename(gb$workDir))}
    record = c(record_id=paste0(runID,"_QC"), run_number=runID)
    qcdata <- jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
    rr<-RCurl::postForm(uri,token=tk,content='record',format='json', type='flat',data=qcdata)
    rr;message(dsh,"Created QC Record",dsh); print(qcdata)
}

# Check if the QC File will be read
checkRunOutput <- function(runID) {
    csvLocation <- file.path(fs::path_home(),"Desktop",runID,paste0(runID,"_Redcap.csv"))
    if (!file.exists(csvLocation)) {
        message(bkRed("File not found:")," ", csvLocation,"\n",bkRed("QC Summary Table will not Knit"),"\n")
    } else{
        fs::file_copy(path=csvLocation, new_path=getwd(), overwrite = T)
    }
}

# gets rid of desktop files if run is successful
tidyUpFiles <- function(runID){
    deskDir <- file.path(getwd())
    backupD <- file.path(gb$methDir,"csvRedcap")
    if(!dir.exists(backupD)){dir.create(backupD)}
    file.copy(deskDir, backupD,overwrite=T, recursive = T, copy.mode = T)
    #unlink(deskDir,T,T)
}

loopRender <- function(samList = NULL, data){
    if(is.null(samList)){
        samList<-1:length(as.character(data$SentrixID_Pos))
    }
    require(rmarkdown)
    require(mnp.v11b6)
    for (i in samList) {
        outFileN = paste0(data[i,1],".html")
        outPathN = file.path(getwd(),outFileN)
        message(outPathN)
        if(file.exists(outPathN)){
            cat(bky(outFileN, "already exists! Skipping sample"),"\n")
            next
        } else {
            cat("\n",bky(dsh,"Now Running", i, "of", length(samList),dsh),sep="\n")
            do_report(data=data[i, ], gb$genCn)
            cat(bky("\n",dsh,"Completed Report", i, "of", length(samList),dsh),sep="\n")
        }
    }
    cat(crayon::black$bgGreen$bold(dsh,"RUN COMPLETE",dsh),sep="\n")
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
                            skipQC=F,email=T,cpReport=F,redcapUp=T){
    assign("genCn",genCn, envir = gb)
    if (is.null(sheetName)) {sheetName="samplesheet.csv"}
    data <- read.csv(file.path(getwd(), sheetName), strip.white=T)
    runID <- paste0(data$RunID[1])
    loopRender(selectSams, data)
    checkRunOutput(runID)
    if (skipQC == F) {
        create.QC.record(runID)
        generateQCreport()
    }
    if(grepl("TEST",runID)){cpReport=F;redcapUp=F;email=F}
    if(cpReport==T){file.list <- gb$copy2outFolder(gb$clinDrv, runID)}
    if(redcapUp==T){file.list <- dir(pattern="*.html", full.names = T); gb$uploadToRedcap(file.list)}
    if(email==T){
        launchEmailNotify(runID)
        #beepr::beep(4)
    }
    tidyUpFiles(runID)
}

# Function to just run a default clinical run without changes, input selectRDs to prioritize samples running first
startRun <- function(selectRDs=NULL, runID=NULL, emailNotify=T){
    if(!is.null(selectRDs)){
        sampleOrder <- reOrderRun(selectRDs) # Re-order sample report generation for priority
        makeReports.v11b6(skipQC=F, email=emailNotify, cpReport=T, selectSams=sampleOrder, redcapUp=T)
    } else {
        makeReports.v11b6(skipQC=F, email=emailNotify, cpReport=F, selectSams=NULL, redcapUp=T)
    }
}

# Executes the functions in order to setup a run
gb$prepareRun <- function(token,runID){
    methylPath <- getwd()
    message("\n","Working directory set to:")
    cat(crayon::bgGreen(methylPath))
    setwd(methylPath)
    gb$setVar("ApiToken", token) # assign the ApiToken & print params
    gb$readSheetWrite() # reads xlsm and generates input .csv samplesheet
    gb$methDir <- methylPath
    gb$moveSampleSheet(methDir, runID=runID) #copies outputs temp to desktop for QC.Rmd
    #gb$classifierInstall(instNew = F, rmpkg = F) # Loads pipeline or installs new
}

gb$startRun <- function(selectRDs=NULL, runID=NULL, emailNotify=T){
    if(!is.null(selectRDs)){
        sampleOrder <- reOrderRun(selectRDs) # Re-order sample report generation for priority
        makeReports.v11b6(skipQC=F, email=emailNotify, cpReport=T, selectSams=sampleOrder, redcapUp=T)
    } else {makeReports.v11b6(skipQC=F, email=emailNotify, cpReport=F, selectSams=NULL, redcapUp=T)}
}

mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/"
script.list = c("LoadInstall_new.R","SetRunParams.R", "CopyInputs.R","CopyOutput.R")

# Source GitHub Scripts
scripts = paste0(mainHub, script.list)
invisible(lapply(scripts, function(i) {devtools::source_url(i)}))

#gb$prepareRun(token,runID)
#gb$startRun()
