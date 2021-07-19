#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE);gb <- globalenv(); assign("gb", gb)
# FUN: Loads all the packages used in the RMD Methylation QC file
checkQCpkg <- function(){
    methylQCpacks <- c(
        "kableExtra","magick","webshot","plyr","ggplot2","knitr","reshape2",
        "data.table","DT","plotly", "MethylAid","minfi","scales","IlluminaHumanMethylation450kmanifest",
        "IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
        "IlluminaHumanMethylationEPICmanifest",
        "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
        "Biobase", "RColorBrewer","limma","ggfortify","Rtsne",
        "qdapTools","gplots","readxl","stringr","ggrepel","Polychrome",
        "tinytex","gridExtra","rmarkdown", "BiocParallel", "grid"
    )
    easypackages::packages(methylQCpacks,prompt=F)
    reqPkg <- list("ggplot2","pals","stringr","scales","grid")
    invisible(lapply(reqPkg, FUN = function(X) {do.call("require", list(X))}))
}

# Helper functions to get and set global variables
setVar <- function(valueName,val){return(assign(valueName, val, envir=.GlobalEnv))}
assignVar <- function(varStr, assignedVal){
    return(tryCatch(expr = {if(!is.null(get(varStr))){message(varStr," = ",assignedVal)}},error = {setVar(varStr,assignedVal)}))}

ckNull <- function(nullVar, subVar, varName){if (is.null(nullVar)){setVar(as.character(varName),as.character(subVar))
    return(paste0(subVar))} else {return(paste0(nullVar))}
}

# FUN: Loads the main packages and dependencies checks if any are not installed
loadClassifierPacks <- function(){
    ldPkg <- tryCatch(
        expr={checkQCpkg()},
        error=function(cond){message("error in loading QC package dependency:\n");message(cond)},
        warning=function(cond){message("Warning in loading QC package dependency:\n");message(cond)}
    )
    return(ldPkg)
}
getDefaults <- function() {
    cbVol = "/Volumes/CBioinformatics"
    moVol = "/Volumes/molecular"
    rsVol = "/Volumes/snudem01labspace"
    defaultParams <- data.frame(
        mnp.pk.loc = paste0(file.path(cbVol, "Methylation/in_house/mnp.v116/mnp.v11b6")),
        ApiToken = "",
        methDir = paste0(file.path(cbVol, "jonathan/Rprojects/Clinical_Methylation")),
        clinDrv = paste0(file.path(moVol, "MOLECULAR LAB ONLY/NYU-METHYLATION")),
        rschOut = paste0(file.path(cbVol, "jonathan/Rprojects/Research_runs/")),
        clinOut = paste0(file.path(moVol, "MOLECULAR/MethylationClassifier")),
        rsch.idat = paste0(file.path(rsVol,"idats")),
        clin.idat = paste0(file.path(moVol, "MOLECULAR/iScan")),
        QC_file = paste0(system.file('Methyl_QC.Rmd', package = "mnp.v11b6")),
        baseDir = paste0(file.path(cbVol, "jonathan/Rprojects")),
        stringsAsFactors=F
    )
    return(defaultParams)
}
getSetvars <- function() {
    assignedVars <- data.frame(
        mnp.pk.loc = gb$mnp.pk.loc,
        ApiToken = gb$ApiToken,
        methDir = gb$methDir,
        clinDrv = gb$clinDrv,
        rschOut = gb$rschOut,
        clinOut = gb$clinOut,
        rsch.idat = gb$rsch.idat,
        clin.idat = gb$clin.idat,
        QC_file = gb$QC_file,
        baseDir = gb$baseDir,
        stringsAsFactors = F
    )
    return(assignedVars)
}

# Sets the default parameters for a methylation Run
defineParams <- function(mnp.pk.loc = NULL, ApiToken = NULL, methDir = NULL,
                         clinDrv = NULL, rschOut = NULL, clinOut = NULL,
                         rsch.idat = NULL, clin.idat = NULL, QC_file = NULL, isMC = T,
                         baseDir = NULL, runID = NULL){
    defVars <- getDefaults()
    inVars <- list(mnp.pk.loc, ApiToken, methDir, clinDrv, rschOut, clinOut,rsch.idat, clin.idat, QC_file, baseDir)
    i=1:length(inVars)
    invisible(lapply(i,function(x){if(!is.null(inVars[[x]])){setVar(names(defVars[x]), inVars[[x]])}}))
    loadClassifierPacks()
    i=1:length(defVars)
    invisible(
        sapply(i,FUN=function(i){ y <- defVars[i]
        assignVar(names((y)), paste0(defVars[,i]))})
    )
    if(!isMC) {methDir=rschOut;assign("workDir","/Volumes/CBioinformatics/jonathan/Rprojects/Research_runs/")}
}

# Changes the working directory using the system CD command
setDirectory <- function(foldr) {
    bsDir = paste("cd", foldr); mm2 = crayon::white$bgRed("Location Not Found:", foldr)
    if (dir.exists(foldr)) {system(bsDir);setwd(foldr);assign("workDir", foldr)} else{warning(mm2)}
}

# Returns a text string of the latest modified Run name, if isMC=False then research directory is returned
listMolecularSheets <- function(isMC=T, getAll=F) {
    researchWorksheets <- "/Volumes/snudem01labspace/Methylation_Worksheets"
    if(isMC){
        wsPath <- file.path(gb$clinDrv,"WORKSHEETS",format(Sys.time(),"%Y"))
        prevMC <- dir(path=wsPath, pattern="MGDM", full.names=T)
    }else{
        wsPath <- file.path(researchWorksheets,format(Sys.time(),"%Y"))
        prevMC <- dir(path=wsPath, pattern="MR", full.names=T)
    }
    if(getAll==T){return(prevMC)}
    newestFile <- which.max(file.info(prevMC)$mtime)
    prevMC<- sub(".xlsm","",basename(prevMC))
    newestRun = paste0(prevMC[newestFile])
    cat(crayon::bgCyan("List of Runs Found:\n"));cat(prevMC,sep="\n")
    cat(crayon::black$bgYellow("Newest Run Found:", crayon::red$bold(paste0(newestRun))))
    return(newestRun)
}

checkValidRun <- function(runID){
    isMC = sjmisc::str_contains(runID, "MGDM")|sjmisc::str_contains(runID, "MC")
    ws.list <- listMolecularSheets(isMC,getAll=T)
    paste0(runID,".xlsm") %in% basename(ws.list)
    found<- stringr::str_detect(ws.list,pattern = runID)
    ifelse(file.exists(ws.list[found]),T,F)
}

#prints the variables currently assigned
printParams <-function(){
    varsNames<- getSetvars()
    cat(crayon::white$bgGreen("\nThese are your variables:\n"))
    i=NULL; varCol=crayon::white$bgBlue$bold("Variable:")
    valCol = crayon::white$bgCyan$bold("Assigned value:")
    invisible(foreach::foreach(i=1:length(varsNames))%do%{cat(paste("\n",varCol, names(varsNames)[i], "\n",valCol, varsNames[1,i],"\n"))})
    funArgs <- crayon::white$bgBlue$bold('gb$setVar("varNameString", varValue)')
    cat("\nUse the function", funArgs ,"to change any default variables above if needed\n")
}

# FUN: copies the molecular or research lab Worksheet xlsm to cwd
copyWorksheetFile <- function(runID=NULL, runYear=NULL) {
    if (is.null(runYear)){runYear=paste0(format(Sys.Date(), "%Y"))}
    if (is.null(runID)){runID=paste0(basename(getwd()))}else{runID=runID}
    stopifnot(!is.null(runID))
    isMC = sjmisc::str_contains(runID, "MGDM")|sjmisc::str_contains(runID, "MC")
    assign("isMC", isMC); message("Is Clinical methylation:",isMC)
    rschDrv <- "/Volumes/snudem01labspace/Methylation_Worksheets"
    mountLoc <- ifelse(isMC,file.path(gb$clinDrv,"WORKSHEETS"),rschDrv)
    message("\nCopying file from: "); cat(mountLoc,"\n")
    if (dir.exists(mountLoc)) {
        fileLoc <- file.path(mountLoc,runYear,paste0(runID,".xlsm"))
        if (file.exists(fileLoc)) {
            message("\n\nCopying worksheet from Worksheets Folder:")
            if(!file.exists(file.path(getwd(),basename(fileLoc)))){
                fs::file_copy(fileLoc, getwd())}
        }else{
            message("File not found:");cat(fileLoc)
            message("\nMake sure your path is correct, try print(gb$copyWorksheetFile)")
        }
    } else{
        message("Folder not found:\n"); cat(mountLoc)
        message("\nEnsure your volume path is correct, check print(copyWorksheetFile)")
    }
}

# Sets the methylation run directory named by the new run name
setRunDir <- function(runID=NULL, workFolder=NULL){
    runID <- ckNull(runID, paste0(basename(getwd())), deparse(substitute(runID,env=.GlobalEnv)))
    workFolder <- ckNull(workFolder, gb$methDir, deparse(substitute(workFolder,env=.GlobalEnv)))
    newRun <- file.path(workFolder, runID)
    assign("newRunPath", newRun)
    if(runID=="21-MGDM_TEST" & dir.exists(newRun)){unlink(newRun, T, T);dir.create(newRun)
        try(unlink("~/Desktop/21-MGDM_TEST",recursive = T),silent = T)
        }
    if(!dir.exists(newRun)){
        dir.create(newRun);cat("creating folder: ",newRun)
        setDirectory(newRun)
        }else{setDirectory(newRun)}
    return(newRun)
}

# Returns Total Sample Count in the run
getTotalSamples <- function(){
    templateDir = "Clinical_Methylation/methylation_run_TEMPLATE_new.xlsm"
    thisSh <- dir(getwd(), "*.xlsm")
    temp <- stringi::stri_detect_fixed(thisSh, "~$")
    thisSh <- thisSh[!temp]
    if(length(thisSh)==0){print("No .xlsm sheet, defaulting to 16");return(16)}
    worksheet <- readxl::read_excel(thisSh[1], col_names="Total", range="B4:B4")
    if (length(worksheet) == 0) {
        message("Samplesheet ", thisSh[1]," is invalid format, manually edit")
        fs::file_copy(path=file.path(gb$baseDir, templateDir), new_path=getwd())
    } else {message("Total sample count found is: ", worksheet[1])}
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

# Check if the QC File will be read
checkRunOutput <- function(runID) {
    location <- file.path("~/Desktop",runID,paste0(runID,"_Redcap.csv"))
    if (!file.exists(location)) {
        message("File not found: ")
        cat(location,sep="\n")
        message("QC Summary Table will not Knit")
    } else{fs::file_copy(path=location, new_path=getwd(), overwrite = T)}
}

# FUN: Reads the csv samplesheet for minfi input
readSampleSheet <- function(runID=F, totalSam=F, wks=F) {
    file.list <- dir(path=getwd(), "*.xlsm")
    temps <- stringi::stri_detect_fixed(file.list, "~$")
    file.list <- file.list[!temps]
    sampleSheet <- paste0(file.list[1])
    message(paste0("Reading worksheet named: ", sampleSheet))
    worksheet <- readxl::read_excel(sampleSheet, sheet=2, col_names=T, col_types="text", trim_ws=T)
    wsDate <-  as.data.frame(
        readxl::read_excel(sampleSheet, sheet=1, col_names=F, range="F4:F4", trim_ws=T))[1]
    names(wsDate)="Date"
    colnames(worksheet)
    worksheet$Date <- paste0(wsDate$Date[1])
    if (runID  ==  T) {return(worksheet$Project[1])}
    if (totalSam  ==  T){
        runNum <- readxl::read_excel(sampleSheet, sheet=1, col_names=F, range="B4")
        runNum <- as.numeric(runNum)
        return(runNum)
    }
    if (wks  ==  T) {return(worksheet)}
}

# Sends an email notification that the run is complete from redcap admin
launchEmailNotify <-function(runID){
    isMC = sjmisc::str_contains(runID, "MGDM") | sjmisc::str_contains(runID, "MC")
    rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", gb$ApiToken)
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
    rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", gb$ApiToken)
    uri=paste0(rcon$url); tk=rcon$token
    if(is.null(runID)){runID<-paste0(basename(gb$workDir))}
    record = c(record_id=paste0(runID,"_QC"), run_number=runID)
    qcdata <- jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
    rr<-RCurl::postForm(uri,token=tk,content='record',format='json', type='flat',data=qcdata)
    rr;message("Created QC Record"); print(qcdata)
}

# Imports the xlsm sheet 3 data
importRedcapStart <- function(nfldr){
    rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", gb$ApiToken)
    uri=paste0(rcon$url); tk=rcon$token
    samSh <- dir(path=getwd(), full.names=T, ".xlsm")
    sampleNumb <- getTotalSamples()
    sh_Dat <-as.data.frame(readxl::read_excel(samSh,sheet=3,range="A1:N97", col_types=c("text")))[,1:13]
    sampleNumb=as.integer(sampleNumb);sh_Dat = sh_Dat[1:sampleNumb,]
    runName <- gb$runID
    filnm = paste0(sh_Dat$record_id, "_cnv.png")
    pathNam = file.path(nfldr,paste0(runName,"_CNVs"),filnm)
    rms=paste(c("control_","low_"),collapse ='|')
    pathNam <- stringr::str_replace(pathNam, rms, "")
    pathNam <- stringr::str_replace(pathNam, "//", "/")
    sh_Dat$cnv_file_path <- pathNam
    if (!is.null(sh_Dat)){
        for(n in 1:nrow(sh_Dat)){
            print(n);record=c(sh_Dat[n,])
            datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T);print(datarecord)
            res<-RCurl::postForm(uri,token=tk,content='record',format='json',type='flat',data=datarecord)
            message("Record Uploaded"); print(res)
        }
    }
}

# FUN: reads the .xlsm worksheet and outputs the .csv methyl experiment for MINFI
readSheetWrite <- function(sampleNumb= NULL, runID = NULL) {
    if(is.null(sampleNumb)){sampleNumb<-getTotalSamples()}
    if(is.null(runID)){runID<-paste0(basename(getwd()))}
    current.run.Folder <- file.path(gb$methDir,runID)
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
        head(df)
        samplesheet_ID=as.data.frame(stringr::str_split_fixed(df[, Sentrix], "_", 2))
        print(samplesheet_ID)
        bn <- file.path(current.run.Folder, df[, Sentrix])
        message("Basename layout: ");message(bn[1])
        note=paste(df$Notes[1])
        df$Notes <- note
        if (is.null(df$Tech)){warning(ww1); df$Tech <- "NA"}
        if (is.null(df$MP_number)){warning(ww2); df$MP_number <- "none"}
        writeSampleSheet(df, samplesheet_ID=samplesheet_ID,bn=bn, sampleName, dnaNumber, Sentrix)
    } else{cat(crayon::white$bgGreen(
        "samplesheet.csv already exists! To Create a new csv, Delete existing samplesheet.csv\n"
    ))}
}

# FUN: Returns a list of idat files given an idat drive location -
getAllFiles <- function(idatDir, csvNam=NULL) {
    if(!is.null(csvNam)){ssheet=read.csv(csvNam,strip.white=T)
    barcode=as.vector(ssheet$Sentrix_ID)} else{
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

# FUN: Copies .idat files to your current directory using sample sheet
copyBaseIdats <- function(allFi) {
    cat(crayon::white$bgCyan("Copying idats to current directory..."),"\n")
    fs::file_copy(allFi, file.path(getwd()))
    idcs = basename(allFi)
    idatsCopied <- idcs[idcs != ""]
    success = file.exists(idatsCopied)
    message(".idat files that failed to copy:")
    if ((length(success) == 0)) {cat("none","\n")} else{print(idatsCopied[!success])}
}

# Helper FUN called during copying idats to notify if a network mount is not found
warnMount <- function(idat.dir){
    cat(crayon::bgRed("\nDirectory not found, make sure idat folder location is mounted:\n", idat.dir))
}

# FUN: Returns a list of idat files that exist on Molecular and Snuderl lab drives -
get.idats <-function(csvNam = "samplesheet.csv"){
    rsch.idat = gb$rsch.idat;clin.idat = gb$clin.idat
    if(!dir.exists(rsch.idat)){warnMount(rsch.idat)}; if(!dir.exists(clin.idat)){warnMount(clin.idat)}
    stopifnot(dir.exists(rsch.idat)|dir.exists(clin.idat))
    if (file.exists(csvNam)) {
        allFi <- getAllFiles(idatDir = c(rsch.idat, clin.idat), csvNam = csvNam)
        allFi = allFi[file.exists(allFi)]
        if (length(allFi) > 0) {
            message("Files found: "); print(allFi)
            cur.idat <- dir(pattern = "*.idat")
            bcds <- paste0(basename(allFi))
            if (all(bcds %in% cur.idat)) {message(".idat files already copied")}
            if (length(cur.idat) < length(allFi)) {copyBaseIdats(allFi[!(bcds %in% cur.idat)])}
        } else{message("No .idat files found! Check worksheet and input folder path")}
    } else {message(paste("Cannot find your sheet named:", csvNam))}
}

# FUN: Copies samplesheet to Desktop folder
moveSampleSheet <- function(methDir, runID=NULL) {
    if (is.null(runID)){runID=paste0(basename(getwd()))}
    baseFolder=paste0("~/Desktop/",runID,"/")
    if(!dir.exists(baseFolder)){dir.create(baseFolder)}
    currDir=paste0(methDir, "/",runID); endDir = paste0(baseFolder, runID, "_samplesheet.csv")
    fs::file_copy(path=paste0(currDir,"/samplesheet.csv"),new_path=baseFolder,overwrite=T)
    file.rename(from=paste0(baseFolder,"samplesheet.csv"), to=endDir)
}

# FUN: Generate CNV image
saveCNVhtml <- function(data) {
    # get RGSet for sample
    RGsetEpic <- minfi::read.metharray(file.path(getwd(),'idats',unlist(data[4])), verbose = F, force = T)
    # Save the CNV as a png file
    MsetEpic <- minfi::preprocessRaw(RGsetEpic)
    sample_id <- unlist(data[1])
    xx <- mnp.v11b6::MNPcnv(MsetEpic,sex = unlist(data[6]),main = sample_id)
    hg <- mnp.v11b6::MNPcnvggplotly(xx,getTables = F)
    htmlwidgets::saveWidget(hg, paste(sample_id, "_cnv.html", sep = ""))
}

# Saves the methyl CNV as a png file in the cwd
generateCNVpng <- function(RGsetEpic, sampleName) {
    imgName <- paste(sampleName, "cnv.png", sep="_")
    MsetEpic <- minfi::preprocessRaw(RGsetEpic)
    png(filename=imgName,width=1820, height=1040)
    mnp.v11b6::MNPcnvplot(MsetEpic)
    dev.off()
}

# Check if your run uses the EPIC B5 or B4 Manifest
checkBarcode <- function(RGsetEpic, barcode) {
    aEpic=c(array="IlluminaHumanMethylationEPIC", annotation="ilm10b4.hg19")
    #bEpic=c(array="IlluminaHumanMethylationEPIC", annotation="ilm10b5.hg38")
    a450k=c(array="IlluminaHumanMethylation450k", annotation="ilmn12.hg19")
    if (barcode >= as.numeric("204220033000")) {
        RGsetEpic@annotation=aEpic;message("Barcode is new EPIC B5 manifest")
    }
    if (RGsetEpic@annotation['array']=="IlluminaHumanMethylation450k") {
        RGsetEpic@annotation=a450k; message("Barcode is 450k using B4 manifest")
    }else {message("Barcode is Default manifest")}
    return(RGsetEpic)
}

# QC REPORT maker: knits the QC RMD file
generateQCreport <- function(runID=NULL, qc=NULL) {
    runID<-ckNull(nullVar = runID, subVar=gb$runID, deparse(substitute(runID,env=.GlobalEnv)))
    QC_file <- system.file('Methyl_QC.Rmd', package = "mnp.v11b6")
    if (!file.exists(QC_file)){message("Check Working directory, QC_file.rmd not found")}
    fs::file_copy(QC_file, getwd(), overwrite = T)
    currentQC = dir(getwd(),"*QC.Rmd", full.names=T)
    qcFile = paste0(runID,"_QC.html") # output file name
    rmarkdown::render(currentQC, output_file=file.path(dirname(currentQC), qcFile), params=list(runID=runID))
    currentQC <- stringr::str_replace_all(string=currentQC, ".Rmd", "_cache")
    unlink(currentQC, recursive=T) #clear cache
}

# FUN: Clean wd Returns Directory Junk Files and Copies redcap file
removeExtras <- function(runID) {
    folder.list = list.dirs(path=getwd(), full.names=T, recursive=T)
    temp.folder.list <- folder.list[stringi::stri_detect_fixed(folder.list, "_files")]
    thisReport <- paste0(runID, "_Redcap.csv")
    file.list<-dir(path=paste0("~/Desktop/",runID,"/"), thisReport, full.names=T)
    fs::file_copy(file.list, getwd(), overwrite=T)
    print(as.data.frame(temp.folder.list))
    unlink(temp.folder.list, recursive=T)
    return(temp.folder.list)
}

# Helper archive function: creates a "previous" folder when reports already exists
save.prev.folder <- function(prevs,oldFi){
    cat(crayon::white$bgRed("Output folder already exists, moving existing reports to previous folder"))
    dir.create(prevs); fs::file_copy(path=oldFi, new_path=prevs)
}

# FUN: Copies Reports to Z drive
copy.to.clinical <- function(clinOut, runID, runYear) {
    newFolder <- file.path(clinOut, runYear,runID); message(newFolder)
    if (!dir.exists(newFolder)) {dir.create(newFolder)}
    if (dir.exists(newFolder)) {
        oldFi=dir(path=newFolder, full.names=T)
        prevs=file.path(newFolder,"previous")
        if (length(oldFi) > 0) {save.prev.folder(prevs,oldFi)}
        mmm1 <- paste("Copying Reports to output folder:\n", newFolder)
        cat(crayon::white$bgBlue$bold(mmm1),"\n")
        fi2copy <- dir(getwd(), pattern="*.html", full.names = T)
        cat(crayon::white$bgGreen$bold("Now copying html reports..."),"\n")
        fs::file_copy(path=fi2copy,newFolder)
    }
}

# REDCap: API call & Upload
# uploads the redcap classifier values must convert to JSON first
importDesktopCsv <- function(rcon,samsheet=NULL) {
    rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", gb$ApiToken)
    ur=paste0(rcon$url);tk=rcon$token
    if(is.null(samsheet)){
        samsheet=dir(path=getwd(),full.names=T,"_Redcap.csv",recursive=F)}else{samsheet=samsheet}
    if (length(samsheet) < 1) {message("Redcap headers csv file not found")}
    if (length(samsheet) == 1) {
        data<-read.csv(samsheet, stringsAsFactors=F)
        if(any(duplicated(data$record_id))){
            message("Removing duplicates in redcap dataframe")
            data = data[!duplicated(data$record_id),]
        }
        for (n in 1:nrow(data)) {
            datarecord = jsonlite::toJSON((as.list(data[n,])), auto_unbox=T);print(datarecord)
            message("\n~~Record Uploaded:")
            RCurl::postForm(ur, token = tk, content = 'record', format = 'csv', type = 'flat', data = datarecord)
            cat(redcapAPI::importRecords(rcon,data,"normal","ids",returnData = F))
        }
    }else{message("no redcap file found")}
}

# Copy Output cnv Files if generated
copy.cnv.files <- function(newFolder, runID, runYear=NULL) {
    if (is.null(runYear)){runYear=paste0(format(Sys.Date(), "%Y"))}
    cnv_folder <- file.path(newFolder,paste0(runID,"_CNVs/"))
    cnvNames <- dir(path=getwd(), full.names=T, "*_cnv.png")
    if(length(cnvNames)>2){
        message(paste0("Copying PNG to: ", cnv_folder));print(as.data.frame(cnvNames))
        if (!dir.exists(cnv_folder)){dir.create(cnv_folder)}
        if (dir.exists(cnv_folder)){fs::file_copy(cnvNames,cnv_folder)}
    }
}

# Uploads any created cnv png files to redcap database
uploadCnPng <- function() {
    rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", gb$ApiToken)
    samSh <- dir(path=getwd(), full.names=T, ".xlsm")
    sampleNumb <- getTotalSamples()
    sh_Dat <-as.data.frame(readxl::read_excel(samSh,sheet=3,range="A1:N97", col_types=c("text")))
    sampleNumb=as.integer(sampleNumb)
    sh_Dat = sh_Dat[1:sampleNumb,]
    #runName=basename(getwd())
    filnm =paste0(sh_Dat$tm_number,"_", sh_Dat$b_number, "_cnv.png")
    pathNam=file.path(getwd(),filnm)
    rms=paste(c("control_","low_"),collapse ='|')
    pathNam <- stringr::str_replace(pathNam, rms, "")
    pathNam <- stringr::str_replace(pathNam, "//", "/")
    records <- sh_Dat$record_id
    for (idx in 1:length(records)) {
        pth = pathNam[idx]
        recordName = paste0(records[idx])
        message("Importing Record:");cat(recordName, sep = "\n");cat(pth, sep = "\n")
        redcapAPI::importFiles(rcon, pth, recordName, field="methyl_cn", repeat_instance=1)
    }
}

# FUN: Copies Reports to Z drive
copy2outFolder <-function(clinDrv = NULL, runID, runYear = NULL) {
    if (is.null(runYear)) {runYear = paste0(format(Sys.Date(), "%Y"))}
    if (is.null(clinDrv)) {clinDrv <- gb$clinDrv}
    isMC = sjmisc::str_contains(runID, "MGDM") | sjmisc::str_contains(runID, "MC")
    
    researchOutDir = file.path("/Volumes/snudem01labspace/FINAL_PDF_Reports_Brain",runID)
    clinicalOutDir = file.path(clinDrv, "Results", runYear, runID)
    runYear = ifelse(isMC, paste0("20", stringr::str_split_fixed(runID, "-", 2)[1]), runYear)
    newFolder <- ifelse(isMC==T, clinicalOutDir, researchOutDir)
    
    cat(crayon::white$bgCyan("Output Folder is:\n", newFolder))
    if (!dir.exists(newFolder)) {dir.create(newFolder)}
    oldFi = dir(path = newFolder, full.names = T)
    prevs = file.path(newFolder, "previous")
    
    if (length(oldFi) > 0) {save.prev.folder(prevs, oldFi)} # saves any old files
    file.list = dir(path = getwd(), "*.html", full.names = T)
    newFolder <- base::rep(newFolder,length(file.list))
    message("\nCopying Existing Reports to Folder...\n")
    fs::file_copy(file.list, newFolder, overwrite = T)
    if (isMC) {
        cnList <- dir(getwd(), "_cnv.png",recursive = F)
        hasCn <- length(cnList) > 2
        if (hasCn) {
            copy.cnv.files(newFolder, runID, runYear)
            uploadCnPng()
            }
        clinOut = file.path(stringr::str_split_fixed(clinDrv, " ", 2)[1],"MethylationClassifier")
        importRedcapStart(clinOut)
        copy.to.clinical(clinOut, runID, runYear)
    }
    return(file.list)
}

# Creates QC record and uploads reports to redcap
uploadToRedcap <- function(file.list) {
    rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", gb$ApiToken)
    message("\nFiles to Import:\n"); print(file.list)
    importDesktopCsv(rcon)
    recordName <-stringr::str_replace_all(string = paste0(file.list), ".html", "")
    runIDs <- rep(gb$runID, length(recordName))
    records <- basename(recordName)
    for (idx in 1:length(records)) {
        pth = file.list[idx]
        recordName = paste0(records[idx])
        runID = paste0(runIDs[idx])
        message("\nImporting Record:")
        cat(recordName, sep = "\n")
        cat(pth, sep = "\n")
        data = data.frame(record_id = recordName, run_number = runID)
        redcapAPI::importRecords(
            rcon, data, overwriteBehavior = "normal",
            returnContent = "ids",
            returnData = T
        )
        redcapAPI::importFiles(
            rcon = rcon,
            file = pth,
            record = recordName,
            field = "classifier_pdf",
            repeat_instance = 1
        )
    }
}

# Helper function called by makeReports.v11b6 to generate the HTML report
do_report <-function(data = NULL) {
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
        RGsetEpic <- minfi::read.metharray(pathEpic, verbose = T, force = T)
        aEpic=c(array="IlluminaHumanMethylationEPIC", annotation="ilm10b4.hg19")
        a450k=c(array="IlluminaHumanMethylation450k", annotation="ilmn12.hg19")
        if (barcode >= as.numeric("204220033000")) {RGsetEpic@annotation=aEpic}
        if (RGsetEpic@annotation['array']=="IlluminaHumanMethylation450k") {RGsetEpic@annotation=a450k}
        RGset = RGsetEpic
        sampleID=paste0(samplename_data)
        FFPE = NULL
        output_dir = getwd()
        sample = 1
        output_file <- paste0(sampleID,".html")
        message(paste("Now running:", samplename_data, run_id, barcode, pathEpic,"\n"))
        rmarkdown::render(
            system.file("report.Rmd",package="mnp.v11b6"),
            output_dir = output_dir,
            output_format = "html_document",
            output_file = output_file,
            clean = T
        )

    }else{message("your data is null")}
}

# Helper function to return the index of priority selected samples first
reOrderRun <- function(selectRDs, sh=NULL){
    if(is.null(selectRDs)){return(NULL)}; if(is.null(sh)){sh<-"samplesheet.csv"}
    allRd <- as.data.frame(read.csv(sh))
    runFirst <- which(allRd[,1] %in% selectRDs)
    runAfter <- which(!(allRd[,1] %in% selectRDs))
    return(c(runFirst,runAfter))
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
    if(is.null(selectSams)){samList <-normList}else{samList <-selectSams}
    for (i in samList) {
        outFileN = paste0(data[i,1],".html")
        outPathN = file.path(gb$workFolder,gb$runID,outFileN)
        if(file.exists(outPathN)){
            cat(bky(outFileN, "exists skipping sample","\n"))
            next
        }else{
            cat(bky("\n",dsh,"Now Running", i, "of", length(samList),dsh),"\n")
            do_report(data=data[i, ])
        }
    }
    cat(crayon::black$bgGreen$bold(dsh,"RUN COMPLETE",dsh),sep="\n")
    checkRunOutput(runID)
    if(skipQC==F){
        create.QC.record(runID)
        generateQCreport()
         # creates a redcap QC record and Knits the QC RMD file
    if(cpReport==T){file.list <- copy2outFolder(gb$clinDrv, runID)}
    if(redcapUp==T){file.list <- dir(pattern="*.html", full.names = T); uploadToRedcap(file.list)}
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
    }else{makeReports.v11b6(skipQC=F, email=emailNotify, cpReport=F, selectSams=NULL, redcapUp=T)}
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
    }
}

prepareRun <- function(token){
    methylPath <- gb$setRunDir(gb$runID)
    message("Working directory set to:"); cat(crayon::bgGreen(methylPath)); setwd(methylPath)
    setVar("ApiToken", token) # assign the ApiToken & print params
    copyWorksheetFile(runID = gb$runID) # copies the xlsm file
    readSheetWrite() # reads xlsm and generates input .csv samplesheet
    get.idats() # Copy idat files to current folder from molecular and snuderlabspace to cwd
    moveSampleSheet(gb$methDir) #copies outputs temp to desktop for QC.Rmd
    gb$install.or.load(instNew = F, rmpkg = F) # Loads pipeline or installs new
}

checkMounts()
defineParams()
