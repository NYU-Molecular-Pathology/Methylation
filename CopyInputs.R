#!/usr/bin/env Rscript
gb <- globalenv(); assign("gb", gb)
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"

# Returns a text string of the latest modified Run name, if isMC=False then research directory is returned
listMolecularSheets <- function(isMC=T, getAll=F) {
    researchWorksheets <- "/Volumes/snudem01labspace/Methylation_Worksheets"
    if(isMC){
        wsPath <- file.path(gb$clinDrv,"WORKSHEETS",format(Sys.time(),"%Y"))
        prevMC <- dir(path=wsPath, pattern="MGDM", full.names=T)
    } else {
        wsPath <- file.path(researchWorksheets,format(Sys.time(),"%Y"))
        prevMC <- dir(path=wsPath, pattern="MR", full.names=T)
    }
    if(getAll==T){return(prevMC)}
    newestFile <- which.max(file.info(prevMC)$mtime)
    prevMC <- sub(".xlsm","",basename(prevMC))
    newestRun = paste0(prevMC[newestFile])
    cat(crayon::bgCyan("List of Runs Found:\n"));cat(prevMC,sep="\n")
    cat(crayon::black$bgYellow("Newest Run Found:", crayon::red$bold(paste0(newestRun))))
    return(newestRun)
}

# Verifies if the runID inputted is valid
checkValidRun <- function(runID){
    isMC = sjmisc::str_contains(runID, "MGDM")|sjmisc::str_contains(runID, "MC")
    ws.list <- listMolecularSheets(isMC,getAll=T)
    found <- paste0(runID,".xlsm") %in% basename(ws.list)
    return(found)
}

# FUN: copies the molecular or research lab Worksheet xlsm to cwd
copyWorksheetFile <- function(runID=NULL, runYear=NULL) {
    if (is.null(runYear)){runYear=paste0(format(Sys.Date(), "%Y"))}
    if (is.null(runID)){runID=paste0(basename(getwd()))} else {runID=runID}
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
        } else {
            message("File not found:");cat(fileLoc)
            message("\nMake sure your path is correct, try print(gb$copyWorksheetFile)")
        }
    } else {
        message("Folder not found:\n"); cat(mountLoc)
        message("\nEnsure your volume path is correct, check print(copyWorksheetFile)")
    }
}

# Sets the methylation run directory named by the new run name
setRunDir <- function(runID=NULL, workFolder=NULL){
    runID <- gb$ckNull(runID, paste0(basename(getwd())), deparse(substitute(runID,env=.GlobalEnv)))
    workFolder <- gb$ckNull(workFolder, gb$methDir, deparse(substitute(workFolder,env=.GlobalEnv)))
    newRun <- file.path(workFolder, runID)
    assign("newRunPath", newRun)
    if(runID=="21-MGDM_TEST" & dir.exists(newRun)){unlink(newRun, T, T);dir.create(newRun)
        try(unlink("~/Desktop/21-MGDM_TEST",recursive = T),silent = T)
    }
    if(!dir.exists(newRun)){
        dir.create(newRun);cat("creating folder: ",newRun)
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
        samplesheet_ID=as.data.frame(stringr::str_split_fixed(df[, Sentrix], "_", 2))
        bn <- file.path(current.run.Folder, df[, Sentrix])
        message("Basename layout:\n",bn[1])
        note=paste(df$Notes[1])
        df$Notes <- note
        if (is.null(df$Tech)){warning(ww1); df$Tech <- "NA"}
        if (is.null(df$MP_number)){warning(ww2); df$MP_number <- "none"}
        writeSampleSheet(df, samplesheet_ID=samplesheet_ID,bn=bn, sampleName, dnaNumber, Sentrix)
    } else {cat(crayon::white$bgGreen(
        "samplesheet.csv already exists! To Create a new csv, Delete existing samplesheet.csv\n"
    ))}
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

# FUN: Copies .idat files to your current directory using sample sheet
copyBaseIdats <- function(allFi) {
    cat(crayon::white$bgCyan("Copying idats to current directory..."),"\n")
    fs::file_copy(allFi, file.path(getwd()))
    idcs = basename(allFi)
    idatsCopied <- idcs[idcs != ""]
    success = file.exists(idatsCopied)
    all(success)
    message(".idat files that failed to copy:")
    if (all(success)) {cat("none","\n")} else {print(idatsCopied[!success])}
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
            cur.idat <- dir(pattern = "*.idat$")
            bcds <- paste0(basename(allFi))
            if (all(bcds %in% cur.idat)) {message(".idat files already copied")}
            if (length(cur.idat) < length(allFi)) {copyBaseIdats(allFi[!(bcds %in% cur.idat)])}
        } else {message("No .idat files found! Check worksheet and input folder path")}
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
    write.csv(
        samplesheet_csv,
        file = "samplesheet.csv",
        quote = F,
        row.names = F
    )
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