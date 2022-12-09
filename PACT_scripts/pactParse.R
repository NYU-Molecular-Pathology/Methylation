#!/usr/bin/env Rscript
library("base")
args <- commandArgs(TRUE)
gb <- globalenv(); assign("gb", gb)

# Main arguments input in comandline
token <- args[1]
inputSheet <- args[2]
runID <- args[3]
#token=NULL;inputSheet=NULL
dsh="\n================"
dsh2="\n==========================\n"

# Displays the Input args -----
message(dsh,"\nParameters input",dsh)
message("token: ",token)
message("PACT Run: ", inputSheet,"\n")
message("Run ID: ", runID,"\n")

# Check Input Params -----
stopifnot(exists("token")|!is.null(token))
stopifnot(exists("inputSheet")|!is.null(inputSheet))

# FUN: Checks if z-drive is accessible to the Rscript
checkMounts <- function(){
    molecDrive = "/Volumes/molecular/MOLECULAR LAB ONLY"
    zDrive = "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular"
    failMount <- ifelse(dir.exists(molecDrive),T,F)
    if(failMount!=T){
        cat("\nPATH does not exist, ensure path is mounted:\n")
        cat(crayon::white$bgRed$bold(molecDrive))
        cat("\nYou must mount the network Z-drive path:\n")
        cat(crayon::white$bgRed$bold(zDrive),"\n")
        stopifnot(!any(failMount==T))
    } else {message("Z-drive path is accessible")}
}

GetPhilipsColumns <- function() {
    filterColumns <- c(
        "MRN",
        "Test Name",
        "Tumor Specimen ID",
        "Normal Specimen ID",
        "Test Number",
        "Epic Order Number",
        "Diagnosis for interpretation",
        "Tumor DNA/RNA Number",
        "Normal DNA/RNA Number",
        "Tumor Percentage"
    )
    return(filterColumns)
}

# Function to silently load library without conflict warnings
libLoad <- function(libName) {
    lib.opts <- list(package = libName, character.only = T, verbose = F, warn.conflicts = F)
    suppressWarnings(suppressPackageStartupMessages(do.call(library, c(lib.opts))))
    message(libName, " ...loaded successful")
}

# Functions to load or install missing required packages -----
loadPacks <- function(){
    if(suppressWarnings(!require("redcapAPI"))){
        params=list('nutterb/redcapAPI', dependencies=T, upgrade="always", type="source")
        do.call(devtools::install_github,c(params))
    } else{libLoad("redcapAPI")}
    pkgs = c("foreach","jsonlite","RCurl","readxl","stringr")
    rlis = getOption("repos");rlis["CRAN"] = "http://cran.us.r-project.org"
    options(repos = rlis)
    invisible(lapply(pkgs, function(pk){
        if(suppressWarnings(!require(pk, character.only=T))){
            install.packages(pk,dependencies=T, verbose=T, repos="http://cran.us.r-project.org", type="both")
            libLoad(pk)
        }else{libLoad(pk)}}))
}

# Returns Path to xlsx file -----
getExcelPath <- function(inputSheet, pathType=1){
    if(stringr::str_detect(inputSheet,"/")==T){return(inputSheet)}
    drive = file.path("", "Volumes", "molecular", "MOLECULAR LAB ONLY")
    folder = file.path("NYU PACT Patient Data", "Workbook")
    runyr <- stringr::str_split_fixed(inputSheet, "-", 3)[, 2]
    ending <-ifelse(pathType==1,".xlsm","_FinalExportedList.xlsx")
    inputFiPath <- file.path(drive, folder, paste0("20", runyr),
                             inputSheet,
                             paste0(inputSheet, ending))
    return(inputFiPath)
}

# Removes and fixes newlines, commas, and blanks from samplesheet
sanitizeSheet <- function(sheetVals){
    mainSheet <- sheetVals[!is.na(sheetVals[,1]),]
    mainSheet$Tumor_Type <- gsub(" ", "-", mainSheet$Tumor_Type)
    for(i in 1:ncol(mainSheet)){
        mainSheet[,i] <- sapply(mainSheet[,i], function(x) { gsub("[\r\n]", "", x) })
        mainSheet[,i] <- sapply(mainSheet[,i], function(x) { gsub(",", "", x) })
        mainSheet[,i] <- sapply(mainSheet[,i], function(x) { gsub(" ", "", x) })
    }
    mainSheet$Paired_Normal[mainSheet$Paired_Normal==0|is.na(mainSheet$Paired_Normal)] <-""
    mainSheet$Tumor_Type[mainSheet$Tumor_Type==0] <- "NA"
    mainSheet[,1:16] <- sapply(mainSheet[,1:16], function(x) { gsub("\\\\", "-", x) })
    controlNames <- "NTC_H20|SC_SERACARE|NC_HAPMAP"
    controlSamples <- grepl(pattern=controlNames, mainSheet$Sample_Name)
    if(table(controlSamples)[['TRUE']]!=3){
        warning("There are not 3 control samples, either NTC_H20, SC_SERACARE, or NC_HAPMAP is missing or added")
    }else{
        controlIndexes <- which(controlSamples==T)
        mainSheet[controlIndexes,'Paired_Normal'] <- ""
    }
    if(any(duplicated(mainSheet$Sample_ID))){
        warning("There are duplicated Sample_ID in the SampleSheet")
        print(mainSheet$Sample_ID[duplicated(mainSheet$Sample_ID)])
    }
    return(mainSheet)
}

# Reads the genders and outputs a .tsv file on desktop
WritePhilipsGender <- function(mainSheet, inputFi, shNames){
    sh2 <- which(grepl("Philips", shNames, ignore.case = T))[1]
    cnvSheet <- mainSheet[,1:15]
    runID <- cnvSheet[1,"Sample_Project"]
    philipVals <- as.data.frame(readxl::read_excel(inputFi, sheet = sh2, skip = 3, col_types = "text"))
    cnvSheet$Gender <- philipVals$Gender[match(cnvSheet$Test_Number, philipVals$`Test Number`)]
    cnvPath <- "/Volumes/molecular/Molecular/MethylationClassifier/CNV_PNG"
    cnvPath <- file.path(cnvPath,paste0(runID,".tsv"))
    message("Writing table to: ", cnvPath)
    write.table(cnvSheet,quote=F, sep='\t', file=cnvPath,row.names=F)
}

GetExcelData <- function(inputFi, sheetNum, shRange, toSkip=0, cm=F){
    sheetData <- suppressMessages(as.data.frame(
        readxl::read_excel(
            inputFi,
            sheet = sheetNum,
            na = "",
            range = shRange,
            col_types = "text",
            col_names = cm,
            skip=toSkip
        )
    ))
    sheetData[is.na(sheetData)] <- ""
    return(sheetData)
}

WriteFileHeader <- function(inputFi){
    pactName <- stringr::str_split_fixed(basename(inputFi), ".xl", 2)[1,1]
    td <- as.character(paste0(format(Sys.time(), "%m/%d/%y")))
    header1 <- c(
        "[Header]",
        "IEMFileVersion",
        "Investigator",
        "Project Name",
        "Experiment Name",
        "Date",
        "Workflow",
        "Application",
        "Assay",
        "Description",
        "Chemistry",
        "[Reads]",
        "151",
        "151",
        "[Settings]",
        "AdapterSequenceRead1",
        "AdapterSequenceRead2"
    )

    header2 <- c(
        "",
        "4",
        "Name",
        pactName,
        pactName,
        td,
        "GenerateFASTQ",
        "FASTQ Only",
        "KAPA-IDT",
        "Description",
        "Capture Enrichment",
        "",
        "",
        "",
        "",
        "AGATCGGAAGAGCACACGTC",
        "AGATCGGAAGAGCGTCGTGT"
    )
    theHeader <- cbind(header1,header2)
    return(theHeader)

}


GetTypeIndex <- function(samNumber, rawSheetData){
    return(
        unlist(lapply(samNumber, function(samName){return(which(samName == rawSheetData$`DNA #`))}))
        )
}

CheckControlRows <- function(rawSheetData, allBnumber){
    controlIndex <- which(rawSheetData$`DNA #` %in% allBnumber == F)
    message(crayon::bgRed("The following samples are Controls or not in Philips:"),"\n",
            paste(rawSheetData$`DNA #`[controlIndex], collapse="\n"),"\n")
    rawSheetData$`Type & Tissue`[controlIndex] <- "Filler"
    isControl <- unlist(lapply(rawSheetData$`Accession#`, function(sam){ sam=="NTC" | sam=="SC" | sam=="NC"}))
    rawSheetData$`Type & Tissue`[isControl] <- "Control"
    message(crayon::bgGreen("Assigned Values:"))
    message(paste(rawSheetData$`DNA #`[controlIndex], rawSheetData$`Type & Tissue`[controlIndex], sep = " = ", collapse="\n"))
    return(rawSheetData)
}

FixPairedList <- function(philipsExport, rawSheetData){
    tumorIndex <- GetTypeIndex(philipsExport$`Tumor DNA/RNA Number`, rawSheetData)
    normalIndex <- GetTypeIndex(philipsExport$`Normal DNA/RNA Number`, rawSheetData)
    rawSheetData$`Type & Tissue` <- ""
    rawSheetData$`Type & Tissue`[tumorIndex] <- "Tumor"
    rawSheetData$`Type & Tissue`[normalIndex] <- "Normal"
    allBnumber <- c(philipsExport$`Tumor DNA/RNA Number`, philipsExport$`Normal DNA/RNA Number`)
    if(!all(allBnumber %in% rawSheetData$`DNA #`)==T){
        warning("Not all B-numbers from PhilipsExport are found in the SampleSheet!")
    }
    rawSheetData <- CheckControlRows(rawSheetData, allBnumber)
    return(rawSheetData$`Type & Tissue`)
}

AddSampleIndexes <- function(pairedList, rawSheetData, philipsExport){
    tissueType <- FixPairedList(philipsExport, rawSheetData)
    sheetTumors <- which(tissueType == "Tumor")
    sheetNormals <- which(tissueType == "Normal" & tissueType != "Control")
    mainSheet <- data.frame(matrix("", nrow = nrow(pairedList), ncol = 0))
    mainSheet$Sample_Name <- mainSheet$Sample_ID <- paste(pairedList[, 1])
    mainSheet$Paired_Normal <- ""
    mainSheet$Paired_Normal[sheetTumors] <- paste(mainSheet$Sample_ID[sheetNormals])
    mainSheet$I7_Index_ID <- paste(rawSheetData$I7_Index_ID)
    mainSheet$index <- paste(rawSheetData$index)
    mainSheet$Specimen_ID <- paste(rawSheetData$`Accession#`)
    mainSheet$Tumor_Content <- mainSheet$Test_Number <- mainSheet$EPIC_ID <- "0"
    mainSheet$Description <- mainSheet$Tumor_Type <- ""
    return(mainSheet)
}


BindUnpairedRows <- function(rawSheetData, pairedList, runID){
    accessions <- rawSheetData$`Accession#`
    controls <- which(rawSheetData$`Type & Tissue` == "Control")
    rowsMissing <- !unlist(lapply(accessions[-controls], function(x){any(grepl(x , pairedList))}))
    doubleBars <- unlist(lapply(pairedList, function(X){stringr::str_detect(X, pattern = "__")}))
    if(any(doubleBars)){
        message(crayon::bgRed("The following rows are missing data:"),"\n",
                paste0(capture.output(pairedList[doubleBars]), collapse = "\n"))
    }
    accToBind <- accessions[rowsMissing]
    if(length(accToBind)==0) {
        message(crayon::bgGreen("No additional paired sample rows to bind to sample sheet"))
        newRows <- NULL
    } else{
        message("Binding rows:","\n", paste0(capture.output(accToBind), collapse = "\n"))
        newRows <- unlist(lapply(1:length(accToBind), function(acc){
            return(paste(0, runID, accToBind[acc], rawSheetData$`DNA #`[acc], sep = "_"))
        }))
    }
    controlRows <- unlist(lapply(1:length(accessions[controls]), function(acc) {
        dnaNum <- which(accessions == accessions[controls][acc])
        return(paste(0, runID, accessions[controls][acc], rawSheetData$`DNA #`[dnaNum], sep = "_"))
    }))
    message(crayon::bgBlue("Binding controls:"),"\n", paste0(capture.output(controlRows), collapse = "\n"))
    pairedList <- data.frame("Sample_ID"=unlist(c(pairedList, newRows, controlRows)))
    rownames(pairedList) <- 1:nrow(pairedList)
    return(pairedList)
}

GetPairedList <- function(philipsExport, runID){
    pairedList <- paste(unlist(lapply(X = 1:length(philipsExport$`Tumor Specimen ID`), function(acc) {
        tumorSam <- paste(
            philipsExport[,"Epic Order Number"][acc],
            runID,
            philipsExport$`Tumor Specimen ID`[acc],
            philipsExport$`Tumor DNA/RNA Number`[acc],
            sep = "_"
        )
        normalSam <- paste(
            philipsExport[,"Epic Order Number"][acc],
            runID,
            philipsExport[,"Normal Specimen ID"][acc],
            philipsExport$`Normal DNA/RNA Number`[acc],
            sep = "_"
        )
        return(rbind(tumorSam, normalSam))
    })))
    return(pairedList)
}


MatchIndex <- function(list1, list2){
    return(unlist(lapply(list1, function(x) {which(x==list2)})))
}

GetIndexMatch <- function(rawSheetData, philipsExport){
    accessions <- rawSheetData$`Accession#`
    philipsN <- c(philipsExport$`Normal Specimen ID`)
    philipsT <- c(philipsExport$`Tumor Specimen ID`)

    idx <- data.frame(
        philipsT = MatchIndex(accessions, philipsT),
        philipsN = MatchIndex(accessions, philipsN),
        wetLabT = MatchIndex(philipsT, accessions),
        wetLabN = MatchIndex(philipsN, accessions)
    )
    return(idx)
}

AddRunChipColumns <- function(mainSheet, runID, pact_run){
    mainSheet$Run_Number <- runID
    splitRun <- stringr::str_split_fixed(runID, "_", 4)
    mainSheet$Sequencer_ID <- splitRun[1, 2]
    mainSheet$Chip_ID <- splitRun[1, 4]
    mainSheet$Sample_Project <- paste0(pact_run)
    return(mainSheet)
}

AddPairedColumn <- function(mainSheet, sheetColumn, philipsColumn, idx){
    mainSheet[idx$wetLabN, sheetColumn] <- philipsColumn[idx$philipsN]
    mainSheet[idx$wetLabT, sheetColumn] <- philipsColumn[idx$philipsT]
    return(mainSheet)
}


FixLastColumns <- function(mainSheet, rawSheetData){
    mainSheet$Tumor_Type[rawSheetData$`Type & Tissue`=="Normal"] <- ""
    mainSheet$Description <- paste0(rawSheetData[,'Description'])
    mainSheet$GenomeFolder <- as.character("PhiX\\Illumina\\RTA\\Sequence\\WholeGenomeFASTA")
    dupes <- base::anyDuplicated(mainSheet$I7_Index_ID)
    if(length(dupes)>0 & dupes!=0){
        message(crayon::bgRed("The following rows are duplicated and will be removed:"), "\n",
                paste0(capture.output(mainSheet[dupes,]), collapse = "\n"))
        mainSheet <- mainSheet[-dupes,]
        row.names(mainSheet) <- 1:nrow(mainSheet)
    }
    mainSheet$Tumor_Content[mainSheet$Paired_Normal==""] <- 0
    return(mainSheet)
}


BuildMainSheet <- function(philipsExport, rawSheetData, runID, pact_run) {
    epicOrder <- philipsExport[,'Epic Order Number']
    testNumber <-  philipsExport[,'Test Number']
    tumorPercent <- philipsExport[,'Tumor Percentage']
    diagColumn <- philipsExport[,'Diagnosis for interpretation']

    pairedList <- GetPairedList(philipsExport, runID)
    pairedList <- BindUnpairedRows(rawSheetData, pairedList, runID)
    mainSheet <- AddSampleIndexes(pairedList, rawSheetData, philipsExport)
    idx <- GetIndexMatch(rawSheetData, philipsExport)
    mainSheet <- AddPairedColumn(mainSheet, "EPIC_ID", epicOrder, idx)
    mainSheet <- AddPairedColumn(mainSheet, "Test_Number", testNumber, idx)
    mainSheet <- AddRunChipColumns(mainSheet, runID, pact_run)
    mainSheet <- AddPairedColumn(mainSheet, "Tumor_Content", tumorPercent, idx)
    mainSheet <- AddPairedColumn(mainSheet, "Tumor_Type", diagColumn, idx)
    mainSheet <- FixLastColumns(mainSheet, rawSheetData)
    return(mainSheet)
}

GetSheetHeading <- function(inputFi){
    shNames <- readxl::excel_sheets(inputFi)
    sheetHead <- WriteFileHeader(inputFi)
    sheetHead <- rbind(sheetHead,c("",""),c("[Data]",""))
    return(sheetHead)
}

GetPhilipsData <- function(inputFi){
    shNames <- readxl::excel_sheets(inputFi)
    sh2 <- which(grepl("Philips", shNames, ignore.case = T))[1]
    philipsExport <- GetExcelData(inputFi, sh2, NULL, 3, cm=T)
    filterColumns <- GetPhilipsColumns()
    philipsExport <- philipsExport[,filterColumns]
    blankOrder <- philipsExport[,"Epic Order Number"] == ""
    if(any(blankOrder)) {
        warning("Some Philips Samples are missing MRNs, defaulting to 0")
        philipsExport[,"Epic Order Number"][blankOrder] <- 0
    }
    blankSpecimen <- philipsExport[,"Normal Specimen ID"] == ""
    if(any(blankSpecimen)) {
        missingSam <- philipsExport$`Tumor Specimen ID`[blankSpecimen]
        message(crayon::bgRed("The following Philips Samples are missing a Normal Specimen ID:\n"),
                paste0(capture.output(missingSam), collapse = "\n"))
    }
    philipsExport[philipsExport==""] <- 0
    return(philipsExport)
}

GetExcelData <- function(inputFi, sheetNum, shRange, toSkip=0, cm=F){
    sheetData <- suppressMessages(as.data.frame(
        readxl::read_excel(
            inputFi,
            sheet = sheetNum,
            na = "",
            range = shRange,
            col_types = "text",
            col_names = cm,
            skip=toSkip
        )
    ))
    sheetData[is.na(sheetData)] <- ""
    return(sheetData)
}

GrabRunNumber <- function(inputFi){
    shNames <- readxl::excel_sheets(inputFi)
    sh <- which(grepl("PACT-", shNames, ignore.case = T))[1]
    rawSheetData <- GetExcelData(inputFi, sh, shRange="A6:X200", cm=T)
    run_number <- try(rawSheetData[which(rawSheetData$`DNA #`=="Run ID:"), 8], silent=T)
    return(run_number)
}

GetRawSamplesheet <- function(inputFi){
    shNames <- readxl::excel_sheets(inputFi)
    sh <- which(grepl("PACT-", shNames, ignore.case = T))[1]
    msgRd <- paste0('Reading Excel Sheet named \"', shNames[sh],'\" from file:')
    message(crayon::bgGreen(msgRd),'\n',inputFi)
    rawSheetData <- GetExcelData(inputFi, sh, shRange="A6:X200", cm=T)
    toDrop <- which(rawSheetData[, "DNA #"]=="HAPMAP")[1]
    rawSheetData <- rawSheetData[1:toDrop,]
    return(rawSheetData)
}

WriteMainSheet <- function(mainSheet, sheetHead){
    outFile <- file.path("~","Desktop",paste(mainSheet[1,"Run_Number"],"SampleSheet.csv",sep="-"))
    write.table(sheetHead, sep=",", file=outFile, row.names=F, col.names=F, quote=F)
    message("Writing file output: ", outFile)
    suppressWarnings(write.table(mainSheet,sep=",", file=outFile, row.names=F, col.names=T, append=T,quote=F))
    return(outFile)
}


AltParseFormat <- function(inputFi, runID){
    rawSheetData <- GetRawSamplesheet(inputFi)
    philipsExport <- GetPhilipsData(inputFi)
    pact_run <- stringr::str_split_fixed(base::basename(inputFi), ".xls", 2)[1,1]
    mainSheet <-  BuildMainSheet(philipsExport, rawSheetData, runID, pact_run)
    return(mainSheet)
}

# Parses xlsx file and writes as csv file -----
parseExcelFile <- function(inputFi, runID = NULL){
    shNames <- readxl::excel_sheets(inputFi)
    message(paste0(capture.output(data.frame(`Sheet names in Workbook` = shNames)), collapse = "\n"))
    sh <- which(grepl("SampleSheet", shNames, ignore.case = T))[1]
    if(!is.na(sh)){
        sheetHead <- GetExcelData(inputFi, sheetNum=shNames[sh], shRange="A1:B17")
        sheetHead <- rbind(sheetHead,c("",""),c("[Data]",""))
        sheetVals <- GetExcelData(inputFi, shNames[sh], NULL, 19, T)
        } else{
        sheetHead <- GetSheetHeading(inputFi)
        sheetVals <- AltParseFormat(inputFi, runID)
        }
    mainSheet <- sanitizeSheet(sheetVals)
    try(WritePhilipsGender(mainSheet,inputFi, shNames), silent=T)
    outFile <- WriteMainSheet(mainSheet, sheetHead)
    return(c(runID = mainSheet[1, "Sample_Project"], outFile = outFile))
}

# Filters list of possible files in the directory for worksheet
filterFiles <- function(potentialFi) {
    wbFiles <- stringr::str_which(basename(potentialFi), pattern = "xlsm")
    if (length(wbFiles) == 0) {
        message("No .xlsm worksheet found. Checking if .xlsx file exists:\n")
        wbFiles <- stringr::str_which(basename(potentialFi), pattern = "book")
    }
    potentialFi <- potentialFi[wbFiles]
    potentialFi <- potentialFi[!stringr::str_detect(potentialFi, "\\$")]
    return(potentialFi)
}

CheckOtherFiles <- function(inputFi, runID){
    message(crayon::bgRed("The PACT run worksheet was not found:"),"\n", inputFi, "\n")
    inputFi <- gsub(paste0("/", basename(inputFi)), "", inputFi)
    potentialFi <- list.files(path = inputFi, full.names = T)
    if(length(potentialFi)>1){
        message(crayon::bgRed("Checking other existing files:"), "\n")
        print(potentialFi)
        potentialFi <- filterFiles(potentialFi)
    }
    if (file.exists(potentialFi[1])) {
        message(crayon::bgGreen("Now trying to read:"),"\n",potentialFi[1],"\n")
        pfile <- potentialFi[1]
    }else{
        message("\n",crayon::bgRed("The PACT run worksheet was not found:"),"\n", potentialFi[1])
        message("\nTrying:", potentialFi[2], "\n")
        stopifnot(file.exists(potentialFi[2]))
        pfile <- potentialFi[2]
    }
    return(suppressMessages(parseExcelFile(inputFi=pfile, runID)))
}

PostRedcapCurl <- function(rcon, datarecord, retcon = 'ids') {
    message(crayon::bgMagenta("REDCap Output:"))
    tryCatch(
        expr = {RCurl::postForm(
            rcon$url, token = rcon$token, content = 'record', format = 'json',
            type = 'flat', data = datarecord, returnContent = retcon, returnFormat = 'csv')
            },
        error = function(e) {message("REDCap push failed!\n", e)},
        finally = message(datarecord)
    )
}

# Generate Email notification and attach csv file
emailNotify <- function(record, rcon){
    record$pact_csv_email <- "pact_csv_email"
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
    PostRedcapCurl(rcon, datarecord)
    message("\n",dsh2,"Email Notification Created",dsh2)
}

# Connect to REDCap and send email attachments of csv file ----
pushToRedcap <- function(outVals, token=NULL) {
    stopifnot(!is.null(token) & !is.null(outVals))
    runID <- outVals[[1]]
    outFile <- outVals[[2]]
    rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", token)
    record = data.frame(record_id = runID, pact_run_number = runID)
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox = T)
    PostRedcapCurl(rcon, datarecord, retcon='nothing')
    tryCatch(
        redcapAPI::importFiles(rcon = rcon, file = outFile, record = runID, field = "pact_csv_sheet", repeat_instance = 1),
        error=function(e){message("REDCap file upload error failed:\n", e)}
        )
    emailNotify(record, rcon)
}

# Gets dataframe and saves as CSV file -----
writeSampleSheet <- function(inputSheet, token, runID = NULL) {
    inputFi <- getExcelPath(inputSheet)
    if (file.exists(inputFi)) {
        outVals <- parseExcelFile(inputFi, runID)
    } else {
        outVals <-  CheckOtherFiles(inputFi, runID)
    }
    if(exists("outVals")){pushToRedcap(outVals, token)}
}

loadPacks()
checkMounts()
writeSampleSheet(inputSheet, token, runID)
