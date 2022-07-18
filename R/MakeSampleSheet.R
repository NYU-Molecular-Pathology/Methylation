#!/usr/bin/env Rscript
gb <- globalenv(); assign("gb", gb)
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
cpInLnk2 = "https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/MakeSampleSheet.R"
rschSheets = "/Volumes/snudem01labspace/Methylation_Worksheets"

msgFunName <- function(pthLnk, funNam){
    message("\nExecuting function: ", crayon::black$bgYellow(funNam), " from RScript in:\n", pthLnk,"\n")
}

msgParams <- function(...){message("\nParams passed: ", crayon::bgGreen(paste(..., sep = ",")), "\n")}

grabYear<- function(yr) {
    msgFunName(cpInLnk2, "grabYear")
    msgParams("yr",yr)

    rnum <- NULL
    if(nchar(yr)>2){rnum <- substring(yr, 3)}else{rnum <- yr}
    if(nchar(yr)>0){rnum <- paste0("20",rnum)}else{rnum}
    return(rnum)
}

# Returns a text string of the latest modified Run name, if isMC=False then research directory is returned
listMolecularSheets <- function(isMC=T, getAll=F,runID=gb$runID) {
    msgFunName(cpInLnk2, "listMolecularSheets")
    msgParams("isMC=",isMC,"getAll=",getAll,"runID=",gb$runID)

    runYear=grabYear(stringr::str_split_fixed(runID,"-",2)[,1])
    if(isMC){
        wsPath <- file.path(gb$clinDrv,"WORKSHEETS",runYear)
        prevMC <- dir(path=wsPath, pattern="MGDM", full.names=T)
    } else {
        wsPath <- file.path(gb$rschSheets,runYear)
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
    msgFunName(cpInLnk2, "checkValidRun")
    isMC = sjmisc::str_contains(runID, "MGDM")|sjmisc::str_contains(runID, "MC")
    ws.list <- listMolecularSheets(isMC,getAll=T,runID)
    found <- paste0(runID,".xlsm") %in% basename(ws.list)
    message("Is the runID valid? ", found)
    if(!found){
        message(crayon::bgRed$white$bold("runID",gb$runID,"is not valid"))
        message(crayon::bgBlue$white$bold(paste0(gb$runID,".xlsm"),"not found in worksheets folders:"))
        message("/Volumes/molecular/MOLECULAR LAB ONLY/NYU-METHYLATION/WORKSHEETS/")
        message("/Volumes/snudem01labspace/Methylation_Worksheets/")
        message("Stopping pipeline. Check worksheet directory if the file is present and saved as .xlsm and not .xlsx")
    }
    stopifnot(found)
    return(found)
}

# FUN: copies the molecular or research lab Worksheet xlsm to cwd
copyWorksheetFile <- function(runID=NULL, runYear=NULL) {
    msgFunName(cpInLnk2, "copyWorksheetFile")
    msgParams(runID, runYear)

    if (is.null(runID)){runID=paste0(basename(getwd()))} else {runID=runID}
    stopifnot(!is.null(runID))
    yr <- stringr::str_split_fixed(runID,"-",2)[,1]
    runYear=grabYear(yr)
    isMC = sjmisc::str_contains(runID, "MGDM")|sjmisc::str_contains(runID, "MC")
    assign("isMC", isMC); message("\nIs methylation run Clinical? ",isMC)

    mountLoc <- ifelse(isMC,file.path(gb$clinDrv,"WORKSHEETS"),gb$rschSheets)
    message("\nCopying file from: "); cat(mountLoc,"\n")
    if (dir.exists(mountLoc)) {
        fileLoc <- file.path(mountLoc,runYear,paste0(runID,".xlsm"))
        if (file.exists(fileLoc)) {
            message("\n\nCopying worksheet from Worksheets Folder:")
            if(!file.exists(file.path(getwd(),basename(fileLoc)))){
                fs::file_copy(fileLoc, getwd(), overwrite = F)}
        } else {
            message("\nFile not found:\n",fileLoc)
            message("\nMake sure your path is correct, try print(gb$copyWorksheetFile)\n")
        }
    } else {
        message("\nFolder not found:\n",mountLoc)
        message("\nEnsure the path to this network Volume is correct\n")
    }
}

GrabSampleSheet <- function(){
    samSh <- dir(path=getwd(), full.names=T, ".xlsm")
    if(length(samSh)>1){
        message("Multiple samplesheets found:\n")
        print(samSh)
        samSh <- samSh[stringr::str_detect(samSh,pattern = "\\$",negate = T)]
    }
    message("Using following samplesheet:\n", samSh[1])
    stopifnot(length(samSh)>0)
    return(samSh[1])
}

# Returns Total Sample Count in the run
getTotalSamples <- function(thisSh=NULL){
    msgFunName(cpOutLnk, "getTotalSamples")

    templateDir = "Clinical_Methylation/methylation_run_TEMPLATE.xlsm"
    thisSh <-  ifelse(is.null(thisSh), GrabSampleSheet(), thisSh)
    thisSh <- thisSh[!stringi::stri_detect_fixed(thisSh, "~$")]
    if(length(thisSh)==0){print("No .xlsm sheet, defaulting to 16 total samples");return(16)}
    worksheet <- suppressMessages(readxl::read_excel(thisSh[1], col_names="Total", range="B4:B4"))
    if (length(worksheet) == 0) {
        warning("Samplesheet ", thisSh[1]," is invalid format or no integer in Cell B4 found.")
        message("Manually edit samplesheet to fix-Try copying data into the template file:\n", templateDir)
        stopifnot(length(worksheet)>0)
    } else {message("Total sample count found is: ", worksheet[1])}
    totNumb <- paste0(worksheet[1])
    return(as.integer(totNumb))
}

# FUN: Reads the csv samplesheet for minfi input
readSampleSheet <- function(runID=F, totalSam=F, wks=F) {
    msgFunName(cpInLnk2, "readSampleSheet")
    msgParams("runID","totalSam","wks")
    msgParams(runID, totalSam, wks)

    file.list <- dir(path=getwd(), "*.xlsm")
    temps <- stringi::stri_detect_fixed(file.list, "~$")
    file.list <- file.list[!temps]
    sampleSheet <- paste0(file.list[1])
    message(paste0("Reading worksheet named: ", sampleSheet))

    worksheet <-
        suppressMessages(readxl::read_excel(
            sampleSheet,
            sheet = 2,
            col_names = T,
            col_types = "text",
            trim_ws = T
        ))
    wsDate <-
        suppressMessages(as.data.frame(
            readxl::read_excel(
                sampleSheet,
                sheet = 1,
                col_names = F,
                range = "F4:F4",
                trim_ws = T
            )
        )[1])
    names(wsDate)="Date"
    message("DATE: ", wsDate$Date)
    colnames(worksheet)
    worksheet$Date <- paste0(wsDate$Date[1])
    if (runID == T) {return(worksheet$Project[1])}
    if (totalSam == T){
        message("Reading cell B4 of ", sampleSheet, "...")
        runNum <- suppressMessages(readxl::read_excel(sampleSheet, sheet=1, col_names=F, range="B4"))
        runNum <- as.numeric(runNum)
        message("Total Samples: ", runNum)
        if(is.null(runNum)|runNum==0){
            warning("Total Number of samples is 0 or Null in worksheet, check cell B4 of .xlsm")
        }
        return(runNum)
    }
    if (wks == T) {return(worksheet)}
}

checkSampleSheet <- function(df){
    msgFunName(cpInLnk2, "checkSampleSheet")

    ww1=crayon::bgRed("No tech name found: check df$Tech in samplesheet.csv assinging NA")
    ww2=crayon::bgRed("No tech name found: check df$MP_number in samplesheet.csv assinging NA")
    ww3=crayon::bgRed("Duplicated sample name found: check df$Sample_Name in samplesheet.csv assigning as None")
    ww4=crayon::bgRed("Some samples are missing RD-numbers Check df$Sample_Name in samplesheet.csv")
    ww5=crayon::bgRed("Some samples are missing B-numbers Check df$b_number in samplesheet.csv")
    ww6=crayon::bgRed("Some samples are missing RD-Numbers or are 0! Check samplesheet.csv")
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
    isMissing <- df$Sample_Name=="0"|is.na(df$Sample_Name)
    if (any(isMissing)) {
        warning(ww6)
        print(df[, c(1, 3, 8:11)])
        stopifnot(any(isMissing)==FALSE)
    }
    dupes <- duplicated(df$Sample_Name)
    if (any(dupes)) {
        warning(ww3)
        message("Duplicated Sample_Name(s): ", paste(df$Sample_Name[dupes]))
        print(df[, c(1, 3, 8:11)])
        stopifnot(any(dupes)==FALSE)
    }
    if (any(stringr::str_count(df$Sample_Name, "control|RD-")==0)) {
        warning(ww4)
        print(df[, c(1, 3, 8:11)])
    }
    if (nrow(df)<8){
        warning("samplesheet.csv contains less that 8 samples!  Check .xslm sheet tab 'raw_labels' for any errors")
        print(df[,1:11])
    }

    return(df)
}

# FUN: translates the xlsm excel file to the .csv samplesheet for redcap and minfi
writeSampleSheet <- function(df, samplesheet_ID, bn = NULL, sampleName, dnaNumber, Sentrix) {
    msgFunName(cpInLnk2, "writeSampleSheet")
    msgParams("df", "samplesheet_ID", samplesheet_ID)
    cat("")
    msgParams("bn",bn)
    msgParams('sampleName',sampleName)
    msgParams('dnaNumber',dnaNumber)
    msgParams("Sentrix", Sentrix)
    if (is.null(bn)) {
        #bn = file.path(getwd(), df$Sentrix)
        bn = file.path(gb$methDir, df$Batch, df$Sentrix)
    }
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

# FUN: Checks header names of .xlsm sheet for mismatch
checkHeaders <- function(Var.names){
    cat("\n\n")
    print(data.frame(Var.names))
    var.default= c("Sample_Name","Sentrix_ID","b_number","MP_number")
    if(any(Var.names!=var.default)){
        stop("Headers in .xlsm 'raw_input' tab are missing or mis-matched:\n",
             paste(var.default[1:4], collapse=" "))
    }
}

# FUN: reads the .xlsm worksheet and outputs the .csv methyl experiment for MINFI
readSheetWrite <- function(sampleNumb= NULL, runID = NULL) {
    msgFunName(cpInLnk2, "readSheetWrite")
    msgParams("sampleNumb","runID")
    msgParams(sampleNumb,runID)

    if(is.null(sampleNumb)){sampleNumb<-getTotalSamples()}
    if(is.null(runID)){runID<-paste0(basename(getwd()))}
    current.run.Folder <- file.path(gb$methDir,runID)
    if (!file.exists("samplesheet.csv")) {
        sampleNumb=as.integer(sampleNumb)
        if(!is.integer(sampleNumb) | sampleNumb<8 ){
        stop("Check samplesheet .xlsm cell B4 for valid integer total arrays")
        }
        worksheet=readSampleSheet(wks=T)
        hdrs = dimnames(worksheet)[[2]]

        print(data.frame(HEADERS=hdrs))
        sampleName=paste(hdrs[9]) # "Sample_Name"
        Sentrix=paste(hdrs[1]) # "Sentrix_ID"
        dnaNumber=paste(hdrs[10]) # "b_number"
        mpnum=paste(hdrs[8]) # "MP_number"

        Var.names=c(sampleName,Sentrix, dnaNumber, mpnum)
        checkHeaders(Var.names)

        df = as.data.frame(worksheet)[1:sampleNumb, ]
        samplesheet_ID=as.data.frame(stringr::str_split_fixed(df[, Sentrix], "_", 2))
        bn <- file.path(current.run.Folder, df[, Sentrix])
        message("Basename layout:\n",bn[1])
        df$Notes <- paste(df$Notes[1])
        df <- checkSampleSheet(df)
        writeSampleSheet(df, samplesheet_ID=samplesheet_ID,bn=bn, sampleName, dnaNumber, Sentrix)
        } else {
            msgCvs = "samplesheet.csv already exists! To Create a new csv file, Delete the existing samplesheet.csv"
            message("\n", crayon::white$bgGreen(msgCvs),"\n")
            }
}
