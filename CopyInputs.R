#!/usr/bin/env Rscript
gb <- globalenv(); assign("gb", gb)
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
cpInLnk = "https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/CopyInputs.R"

msgFunName <- function(pthLnk, funNam){
message("\nExecuting function: ", funNam, " from RScript in:\n", pthLnk,"\n")
}

grabYear<- function(yr) {
    msgFunName(cpInLnk, "grabYear")
    rnum <- NULL
    if(nchar(yr)>2){
        rnum <- substring(yr, 3)
    }else{rnum <- yr}
    if(nchar(yr)>0){
        rnum <- paste0("20",rnum)
    }else{rnum}
    return(rnum)
}

# Returns a text string of the latest modified Run name, if isMC=False then research directory is returned
listMolecularSheets <- function(isMC=T, getAll=F,runID=gb$runID) {
        msgFunName(cpInLnk, "listMolecularSheets")
    researchWorksheets <- "/Volumes/snudem01labspace/Methylation_Worksheets"
    runYear=grabYear(stringr::str_split_fixed(runID,"-",2)[,1])
    if(isMC){
        wsPath <- file.path(gb$clinDrv,"WORKSHEETS",runYear)
        prevMC <- dir(path=wsPath, pattern="MGDM", full.names=T)
    } else {
        wsPath <- file.path(researchWorksheets,runYear)
        prevMC <- dir(path=wsPath, pattern="MR", full.names=T)
    }
    if(getAll==T){
        return(prevMC)
    }else{
        newestFile <- which.max(file.info(prevMC)$mtime)
        prevMC <- sub(".xlsm","",basename(prevMC))
        newestRun = paste0(prevMC[newestFile])
        cat(crayon::bgCyan("List of Runs Found:\n"));cat(prevMC,sep="\n")
        cat(crayon::black$bgYellow("Newest Run Found:", crayon::red$bold(paste0(newestRun))))
        return(newestRun)
    }
}

# Verifies if the runID inputted is valid
checkValidRun <- function(runID){
    msgFunName(cpInLnk, "checkValidRun")
    isMC = sjmisc::str_contains(runID, "MGDM")|sjmisc::str_contains(runID, "MC")
    ws.list <- listMolecularSheets(isMC,getAll=T,runID)
    found <- paste0(runID,".xlsm") %in% basename(ws.list)
    return(found)
}

# FUN: copies the molecular or research lab Worksheet xlsm to cwd
copyWorksheetFile <- function(runID=NULL, runYear=NULL) {
    msgFunName(cpInLnk, "copyWorksheetFile")
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
        msgFunName(cpInLnk, "setRunDir")
    runID <- gb$ckNull(runID, paste0(basename(getwd())), deparse(substitute(runID,env=gb)))
    if(is.null(gb$methDir)){gb$methDir<- "/Volumes/CBioinformatics/Methylation/Clinical_Runs"}
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
    msgFunName(cpInLnk, "getTotalSamples")
    templateDir = "Clinical_Methylation/methylation_run_TEMPLATE_new.xlsm"
    thisSh <- dir(getwd(), "*.xlsm")
    temp <- stringi::stri_detect_fixed(thisSh, "~$")
    thisSh <- thisSh[!temp]
    if(length(thisSh)==0){print("No .xlsm sheet, defaulting to 16");return(16)}
    worksheet <- readxl::read_excel(thisSh[1], col_names="Total", range="B4:B4")
    if (length(worksheet) == 0) {
        message("Samplesheet ", thisSh[1]," is invalid format, manually edit")
        message("Try copying the template:\n", templateDir)
    } else {message("Total sample count found is: ", worksheet[1])}
    return(paste0(worksheet[1]))
}

# FUN: translates the xlsm excel file to the .csv samplesheet for redcap and minfi
writeSampleSheet <- function(df, samplesheet_ID, bn = NULL, sampleName, dnaNumber, Sentrix) {
    msgFunName(cpInLnk, "writeSampleSheet")
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
    msgFunName(cpInLnk, "readSampleSheet")
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

checkSampleSheet <- function(df){
    msgFunName(cpInLnk, "checkSampleSheet")
    ww1=crayon::bgRed("No tech name found: check df$Tech in samplesheet.csv assinging NA")
    ww2=crayon::bgRed("No tech name found: check df$MP_number in samplesheet.csv assinging NA")
    ww3=crayon::bgRed("Duplicated sample name found: check df$Sample_Name in samplesheet.csv assigning as None")
    ww4=crayon::bgRed("Some Samples are missing RD-numbers Check df$Sample_Name in samplesheet.csv")
    ww5=crayon::bgRed("Some Samples are missing B-numbers Check df$b_number in samplesheet.csv")
    if (is.null(df$Tech)) {
        warning(ww1)
        df$Tech <- "NA"
    }
    if (is.null(df$MP_number)) {
        warning(ww2)
        df$MP_number <- "none"
        print(df$MP_number)
    }
    if (is.null(df$b_number)) {
        warning(ww5)
        df$b_number <- "blank"
        print(df$b_number)
    }
    if (any(duplicated(df$Sample_Name))) {
        warning(ww3)
        print(df[, c(1, 3, 8:11)])
    }
    if (any(!grepl("RD-", df$Sample_Name, fixed = T))) {
        warning(ww4)
        print(df[, c(1, 3, 8:11)])
              }
    
    return(df)
}

# FUN: reads the .xlsm worksheet and outputs the .csv methyl experiment for MINFI
readSheetWrite <- function(sampleNumb= NULL, runID = NULL) {
    msgFunName(cpInLnk, "readSheetWrite")
    if(is.null(sampleNumb)){sampleNumb<-getTotalSamples()}
    if(is.null(runID)){runID<-paste0(basename(getwd()))}
    current.run.Folder <- file.path(gb$methDir,runID)
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
        df <- checkSampleSheet(df)
        writeSampleSheet(df, samplesheet_ID=samplesheet_ID,bn=bn, sampleName, dnaNumber, Sentrix)
    } else {
        msgCvs = "samplesheet.csv already exists! To Create a new csv file, Delete the existing samplesheet.csv"
        message("\n", crayon::white$bgGreen(msgCvs),"\n")
    }
}

# FUN: Returns a list of idat files given an idat drive location -
getAllFiles <- function(idatDir, csvNam=NULL) {
    msgFunName(cpInLnk, "getAllFiles")
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

# FUN: Copies .idat files to your current directory using sample sheet
copyBaseIdats <- function(allFi) {
    msgFunName(cpInLnk, "copyBaseIdats")
    cat(crayon::white$bgCyan("Copying idats to current directory..."),"\n")
    fs::file_copy(allFi, file.path(getwd()),overwrite=T)
    idcs = basename(allFi)
    idatsCopied <- idcs[idcs != ""]
    success = file.exists(idatsCopied)
    all(success)
    message(".idat files that failed to copy:")
    if (all(success)) {cat("none","\n")} else {print(idatsCopied[!success])}
}

# Helper FUN called during copying idats to notify if a network mount is not found
warnMount <- function(idat.dir){
    msgFunName(cpInLnk, "warnMount")
    cat(crayon::bgRed("Directory not found, ensure the idat folder location is accessible:"),idat.dir,sep="\n")
}

# FUN: Returns a list of idat files that exist on Molecular and Snuderl lab drives -
get.idats <-function(csvNam = "samplesheet.csv"){
    msgFunName(cpInLnk, "get.idats")
    rsch.idat <- gb$rsch.idat;clin.idat <- gb$clin.idat
    if(!dir.exists(rsch.idat)){warnMount(rsch.idat)}; if(!dir.exists(clin.idat)){warnMount(clin.idat)}
    stopifnot(dir.exists(rsch.idat)|dir.exists(clin.idat))
    if (file.exists(csvNam)) {
        allFi <- getAllFiles(idatDir = c(rsch.idat, clin.idat), csvNam = csvNam)
        allFi = allFi[file.exists(allFi)]
        if (length(allFi) > 0) {
            message("Files found: "); print(allFi)
            cur.idat <- dir(pattern = "*.idat$")
            bcds <- paste0(basename(allFi))
            if (all(bcds %in% cur.idat)) {message(".idat files already copied")}
            if (!all(bcds %in% cur.idat)) {copyBaseIdats(allFi[!(bcds %in% cur.idat)])} #length(cur.idat) < length(allFi)
        } else {message("No .idat files found! Check worksheet and input folder path")}
    } else {message(paste("Cannot find your sheet named:", csvNam))}
}

# FUN: Copies samplesheet to Desktop folder
moveSampleSheet <- function(methDir, runID=NULL) {
    msgFunName(cpInLnk, "moveSampleSheet")
    if (is.null(runID)){runID=paste0(basename(getwd()))}
    baseFolder=paste0("~/Desktop/",runID,"/")
    if(!dir.exists(baseFolder)){dir.create(baseFolder)}
    currDir=paste0(methDir, "/",runID); endDir = paste0(baseFolder, runID, "_samplesheet.csv")
    fs::file_copy(path=paste0(currDir,"/samplesheet.csv"),new_path=baseFolder,overwrite=T)
    file.rename(from=paste0(baseFolder,"samplesheet.csv"), to=endDir)
}

#  Copy idats and Worksheets creation
writeFromRedcap <- function(df, samplesheet_ID, bn = NULL) {
     msgFunName(cpInLnk, "writeFromRedcap")
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
    write.csv(
        samplesheet_csv,
        file = "samplesheet.csv",
        quote = F,
        row.names = F
    )
}

#' FUN: Returns dataframe of redcap search using default worksheet header and fields
search.redcap <- function(rd_numbers, token=NULL, flds=NULL) {
    msgFunName(cpInLnk, "search.redcap")
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
    msgFunName(cpInLnk, "get.rd.info")
    if (is.null(rd_numbers)){message("Input RD-numbers using get.rd.info(rd_numbers)")}
    if (is.null(sh_name)) {sh_name = "samplesheet.csv"}
    if(is.null(token)){message("You must provide an ApiToken!")};stopifnot(!is.null(token))
    result <- search.redcap(rd_numbers, token, NULL)
    samplesheet_ID = as.data.frame(stringr::str_split_fixed(result[,"barcode_and_row_column"],"_",2))
    writeFromRedcap(result, samplesheet_ID) # writes API export as minfi dataframe sheet
    get.idats(csvNam = sh_name)  # copies idat files from return to current directory
    
    return(result)
}
