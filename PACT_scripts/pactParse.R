#!/usr/bin/env Rscript
## ---------------------------
## Script name: pactParse.R
## Purpose: source of global scripts and generate PACT -SampleSheet.csv file
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------

library("base"); args <- commandArgs(TRUE); gb <- globalenv(); assign("gb", gb)
dsh<-"\n================"
dsh2<-"\n==========================\n"

# Main arguments input in comandline (Uncomment to Debug or run Locally) -----------------------
args[1] -> token        #<- '8XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' = APITOKEN
args[2] -> inputSheet   #<- 'PACT-YY-##'
args[3] -> runID        #<- 'YYMMDD_NB######_0###_AH7ABCDEFG'

# Display & Verify the Input args --------------------------------------------------------------
message(dsh,"\nParameters input",dsh)
message("token:    ", token, "\n", "PACT Run: ", inputSheet,"\n","Run ID:   ", runID,"\n")
stopifnot(exists("token") & !is.null(token) & exists("inputSheet") & !is.null(inputSheet))

MsgDF <- function(data){
    return(message(paste0(capture.output(as.data.frame(data)), collapse = "\n")))
}

# FUN: Checks if z-drive is accessible to the Rscript ------------------------------------------
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

# FUN: Returns hardcoded Philips Export Tab Column Names----------------------------------------
GetPhilipsColumns <- function() {
    filterColumns <- c(
        "MRN", "Test Name", "Tumor Specimen ID", "Normal Specimen ID", "Test Number",
        "Epic Order Number", "Diagnosis for interpretation", "Tumor DNA/RNA Number",
        "Normal DNA/RNA Number","Tumor Percentage"
    )
    return(filterColumns)
}

# Function to silently load library without conflict warnings ----------------------------------
libLoad <- function(libName) {
    lib.opts <- list(package = libName, character.only = T, verbose = F, warn.conflicts = F)
    suppressWarnings(suppressPackageStartupMessages(do.call(library, c(lib.opts))))
    message(libName, " ...loaded successful")
}

# Functions to load or install missing required packages ---------------------------------------
loadPacks <- function(){
    if(!require("devtools")){install.packages("devtools", quiet=T, dependencies=T)}
    if(suppressWarnings(!require("redcapAPI"))){
        params=list('nutterb/redcapAPI', dependencies=T, upgrade="always", type="source")
        do.call(devtools::install_github, c(params))
    } else{libLoad("redcapAPI")}
    pkgs = c("foreach","jsonlite","RCurl","readxl","stringr")
    rlis = getOption("repos");rlis["CRAN"] = "http://cran.us.r-project.org"
    options(repos = rlis)
    invisible(lapply(pkgs, function(pk){
        if(suppressWarnings(!require(pk, character.only=T))){
            install.packages(pk, dependencies=T, verbose=T, repos="http://cran.us.r-project.org", type="both")
            libLoad(pk)
        }else{libLoad(pk)}}))
    if(paste(utils::packageVersion("redcapAPI")) != "2.7.4"){
        install.packages("redcapAPI", ask=F, update=T, dependencies=T)
    }

}


# FUN: Returns Path to xlsx file ---------------------------------------------------------------
getExcelPath <- function(inputSheet, pathType = 1, isNxtSeq=F){
    if(stringr::str_detect(inputSheet,"/")==T){return(inputSheet)}
    drive = file.path("", "Volumes", "molecular", "MOLECULAR LAB ONLY")
    folder = file.path("NYU PACT Patient Data", "Workbook")
    runyr <- stringr::str_split_fixed(inputSheet, "-", 3)[, 2]
    yearDir <- paste0("20", runyr)
    ending <- ifelse(pathType==1,".xlsm","_FinalExportedList.xlsx")
    xlFi <- paste0(inputSheet, ending)
    inputFiPath <- file.path(drive, folder, yearDir, inputSheet, xlFi)
    if(isNxtSeq == T){
        folder <- "Validations/PACT new i7-NextSeq2000/Wet Lab/Workbook"
        inputFiPath <- file.path(drive, folder, yearDir, xlFi)
    }
    return(inputFiPath)
}


MsgChangesMade <- function(mainSheet, ptrn=","){
    for(i in 1:ncol(mainSheet)){
        theRows <- mainSheet[,i]
        if(any(is.na(theRows))){
            theRows[which(is.na(theRows))] <- ""
        }
        theMatch <- stringr::str_detect(theRows, pattern=ptrn)
        if(any(theMatch)){
            theCommas <- which(theMatch==T)
            theSym <- paste0('"', ptrn, '"')
            message(paste('The following', theSym, 'will be removed from rows',
                          paste(theCommas, collapse=", "), 'in the column:', colnames(mainSheet)[i]))
            MsgDF(theRows[theCommas])
        }
    }
}


# Removes and fixes newlines, commas, and blanks from samplesheet ------------------------------
sanitizeSheet <- function(mainSheet){
    mainSheet <- mainSheet[!is.na(mainSheet[,1]),]
    spacedTxt <- stringr::str_detect(mainSheet$Tumor_Type, " ")
    if(any(spacedTxt)){
        message("The following Tumor_Type have spaces:")
        MsgDF(mainSheet$Tumor_Type[spacedTxt])
    }
    mainSheet$Tumor_Type <- gsub(" ", "-", mainSheet$Tumor_Type)

    MsgChangesMade(mainSheet)
    MsgChangesMade(mainSheet, " ")
    MsgChangesMade(mainSheet, "[\r\n]")

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
        warning("There are not 3 control samples, either NTC_H20, SC_SERACARE, or NC_HAPMAP is missing OR there are extra controls added in this run")
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

# Reads the genders and outputs a tsv file to /Molecular/MethylationClassifier/CNV_PNG ---------
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
            skip = toSkip
        )))
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
    mislabelled <- philipsExport$`Tumor DNA/RNA Number` %in% philipsExport$`Normal DNA/RNA Number`
    if(any(mislabelled)){
        message(crayon::bgRed("There are samples potentially mislabelled! in Philips:"),"\n")
        MsgDF(philipsExport$`Tumor DNA/RNA Number`[mislabelled])
    }
    tumorIndex <- GetTypeIndex(philipsExport$`Tumor DNA/RNA Number`, rawSheetData)
    normalIndex <- GetTypeIndex(philipsExport$`Normal DNA/RNA Number`, rawSheetData)
    rawSheetData$`Type & Tissue` <- ""
    rawSheetData$`Type & Tissue`[tumorIndex] <- "Tumor"
    rawSheetData$`Type & Tissue`[normalIndex] <- "Normal"
    allBnumber <- c(philipsExport$`Tumor DNA/RNA Number`, philipsExport$`Normal DNA/RNA Number`)
    if(!all(allBnumber %in% rawSheetData$`DNA #`)==T){
        message(crayon::bgRed("Not all B-numbers from PhilipsExport tab are in the SampleSheet DNA # Column!"))
    }
    rawSheetData <- CheckControlRows(rawSheetData, allBnumber)
    return(rawSheetData$`Type & Tissue`)
}

AddSampleIndexes <- function(pairedList, rawSheetData, philipsExport){
    tissueType <- FixPairedList(philipsExport, rawSheetData)
    sheetTumors <- which(tissueType == "Tumor")
    sheetControl <- which(tissueType != "Control")
    sheetNormals <- which(tissueType == "Normal" & tissueType != "Control")
    mainSheet <- data.frame(matrix("", nrow = nrow(pairedList), ncol = 0))
    mainSheet$Sample_Name <- mainSheet$Sample_ID <- paste(pairedList[, 1])
    mainSheet$Paired_Normal <- ""

    if(length(mainSheet$Paired_Normal[sheetTumors])!=length(paste(mainSheet$Sample_ID[sheetNormals]))){
        suppressWarnings(mainSheet$Paired_Normal[sheetTumors] <- paste(mainSheet$Sample_ID[sheetNormals]))
        theSampleIdEx <- paste(mainSheet$Sample_ID[sheetNormals])
        extraIndex <- which(theSampleIdEx %in% mainSheet$Paired_Normal[sheetTumors] == F)
        message(crayon::bgRed(paste(length(extraIndex),
                                    "Samples are extra normals and will need to be added manually:")))
        thenNext <- theSampleIdEx[extraIndex]
        MsgDF(thenNext)
        newNormList <- paste(mainSheet$Sample_ID[sheetNormals])[-extraIndex]
        mainSheet$Paired_Normal[sheetTumors] <- newNormList
        #missingTumorNorm <- stringr::str_ends(mainSheet$Sample_ID, "0_0")
        missingTumorNorm <- which(stringr::str_ends(mainSheet$Sample_ID,  paste(thenNext,collapse="|")))
        if(length(missingTumorNorm)>0){
            message(crayon::bgRed("There are samples missing tumor normal pairs and will be dropped:"))
            MsgDF(mainSheet[missingTumorNorm,1])
            mainSheet <- mainSheet[-missingTumorNorm,]
            rownames(mainSheet) <- 1:nrow(mainSheet)
        }
    }else{
        mainSheet$Paired_Normal[sheetTumors] <- paste(mainSheet$Sample_ID[sheetNormals])
    }
    if (length(rawSheetData$I7_Index_ID) != nrow(mainSheet)){

    }
    mainSheet$I7_Index_ID <- paste(rawSheetData$I7_Index_ID)
    mainSheet$index <- paste(rawSheetData$index)
    mainSheet$Specimen_ID <- paste(rawSheetData$`Accession#`)
    mainSheet$Tumor_Content <- mainSheet$Test_Number <- mainSheet$EPIC_ID <- "0"
    mainSheet$Description <- mainSheet$Tumor_Type <- ""
    return(mainSheet)
}


FixDuplicateControls <- function(controlRows) {
    message("Fixing duplicated controls:\n", paste(controlRows[duplicated(controlRows)], collapse = "\n"))
    
    for (unique_string in unique(controlRows)) {
        indices <- which(controlRows == unique_string)
        n <- length(indices)
        if (n > 1) {
            new_values <- sapply(seq_len(n) - 1, function(x) {
                substr(unique_string, 1, 1) <- as.character(x)
                return(unique_string)
            })
            controlRows[indices] <- rev(new_values)
        }
    }
    
    return(controlRows)
}


BindUnpairedRows <- function(rawSheetData, pairedList, runID) {
    accessions <- rawSheetData$`Accession#`
    dnaNumbers <- rawSheetData$`DNA #`
    controls <- which(str_detect(rawSheetData$`Type & Tissue`, "Control"))
    missRows <- vapply(accessions[-controls], function(x) !any(str_detect(pairedList, x)), logical(1))
    doubleBars <- vapply(pairedList, function(X) {str_detect(X, pattern = "__")}, logical(1))
    
    if(any(doubleBars)){
        message(crayon::bgRed("The following rows are missing data:"),"\n",
                paste(pairedList[doubleBars], collapse = "\n"))
    }
    
    newRows <- NULL
    
    if(any(missRows)) {
        message(crayon::bgGreen("Binding additional rows/filler:"), "\n", 
                paste(accessions[missRows], collapse = "\n"))
        extraRow <- accessions == accessions[missRows]
        newRows <-  paste0(accessions[missRows], "_", runID, "_", accessions[extraRow], "_", dnaNumbers[extraRow]))
    } else {
        message(crayon::bgGreen("No additional fillers or paired sample rows to bind to sample sheet"))
    }
    
    cntrlT <- rawSheetData$Test_Number == rawSheetData$Test_Number[controls]
    controlRows <- paste0(rawSheetData$Test_Number[controls], "_", runID, "_", accessions[cntrlT], "_", dnaNumbers[cntrlT])
    message(crayon::bgBlue("Binding controls:"),"\n", paste(controlRows, collapse = "\n"))
    
    if(any(duplicated(controlRows))) {
        controlRows <- FixDuplicateControls(controlRows)
        message(crayon::bgBlue("New control names:"), paste(controlRows, collapse = "\n"))
    }
    
    pairedList <- data.frame("Sample_ID" = c(pairedList, newRows, controlRows))
    rownames(pairedList) <- seq_len(nrow(pairedList))
    
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


DisplayMissing <- function(ngsMissing, valMissing, type) {
    message(crayon::bgRed("Some samples are missing tumor/normal pairs:"))
    message(paste(ngsMissing, "Philips", type, "=" , valMissing, collapse = "\n")) 
}

CheckMissingPairs <- function(ngsNumbers, philipsT, philipsN) {
    listCheck <- list("Normal" = philipsN, "Tumor" = philipsT)
    
    for (type in names(listCheck)) {
        val <- listCheck[[type]]
        miss <- val == 0 | is.na(val)
        if (any(miss)) {DisplayMissing(ngsNumbers[miss], val[miss], type)}
    }
    
    dupes <- philipsN %in% philipsT
    if (any(dupes)) {
        message(crayon::bgRed("Philips samples have the same tumor/normal accession#:"))
        MsgDF(ngsNumbers[dupes])
        MsgDF(philipsN[dupes])
    }
}


CheckDupesMatch <- function(nameList, secondCol) {
    matchIdx <- MatchIndex(secondCol, nameList)
    if (any(duplicated(matchIdx))) {
        message(crayon::bgRed("The sample sheet may contain duplicate accession Numbers!"))
        dupesIdx <- matchIdx[duplicated(matchIdx)]
        message("Duplicate at row Index: ", dupesIdx)
        MsgDF(nameList[dupesIdx])
        MsgDF(secondCol[dupesIdx])
        #return(unique(matchIdx))
        return(matchIdx)
    } else{
        return(matchIdx)
    }
}


CheckTotalIndexes <- function(tumorSam, normalSam, ngsNumbers, sheetType = "Philips"){
    if(length(tumorSam) != length(normalSam)){
        message("Total Tumors: ", length(tumorSam))
        message("Total Normals: ", length(normalSam))
        tumorLonger <- length(tumorSam) > length(normalSam)
        totalAdd <- abs(length(tumorSam) - length(normalSam))
        message(crayon::bgRed(paste("Not all",sheetType,"samples are paired with a tumor or normal!")))
        if(tumorLonger==F){
            totalDrop <- normalSam %in% tumorSam
            message(crayon::bgRed(paste(totalAdd, "more Normals than Tumors listed in worksheet!")))
            MsgDF(data.frame(
                NGS=ngsNumbers[!totalDrop], Potential.Extra=philipsN[!totalDrop]))
            tumorSam <- c(tumorSam, rep(0,totalAdd))
        }else{
            totalDrop <- tumorSam %in% normalSam
            message(crayon::bgRed(paste(totalAdd, "more Tumors than Normals listed in worksheet!")))
            MsgDF(data.frame(
                NGS=ngsNumbers[!totalDrop], Potential.Extra=philipsT[!totalDrop]))
            normalSam <- c(normalSam, rep(0,totalAdd))
        }
        message("Extra ", length(which(totalDrop==F)), " cases...")
    }

}

GetIndexMatch <- function(rawSheetData, philipsExport){
    message("Running GetIndexMatch function in pactParse.R")
    accessions <- rawSheetData$`Accession#`
    philipsN <- c(philipsExport$`Normal Specimen ID`)
    philipsT <- c(philipsExport$`Tumor Specimen ID`)
    ngsNumbers <- philipsExport$`Test Number`

    CheckMissingPairs(ngsNumbers, philipsT, philipsN)

    philipsIdxT <- CheckDupesMatch(philipsT, accessions)
    philipsIdxN <- CheckDupesMatch(philipsN, accessions)

    wetLabIdxT <- CheckDupesMatch(accessions, philipsT)
    wetLabIdxN <- CheckDupesMatch(accessions, philipsN)

    if(length(philipsIdxT) != length(philipsIdxN)){
        message("Total Tumors: ", length(philipsIdxT))
        message("Total Normals: ", length(philipsIdxN))
        tumorLonger <- length(philipsIdxT) > length(philipsIdxN)
        totalAdd <- abs(length(philipsIdxT) - length(philipsIdxN))
        message(crayon::bgRed("Not all Philips samples are paired with a tumor or normal!"))
        if(tumorLonger==F){
            totalDrop <- philipsIdxN %in% philipsIdxT
            message(crayon::bgRed(paste(totalAdd, "more Normals than Tumors listed in worksheet!")))
            MsgDF(data.frame(
                NGS=ngsNumbers[!totalDrop], Potential.Extra=philipsN[!totalDrop]))
            philipsIdxT <- c(philipsIdxT, rep(0,totalAdd))
        }else{
            totalDrop <- philipsIdxT %in% philipsIdxN
            message(crayon::bgRed(paste(totalAdd, "more Tumors than Normals listed in worksheet!")))
            MsgDF(data.frame(
                NGS=ngsNumbers[!totalDrop], Potential.Extra=philipsT[!totalDrop]))
            philipsIdxN <- c(philipsIdxN, rep(0,totalAdd))
        }
        message("Extra ", length(which(totalDrop==F)), " cases...")
    }

    if(length(wetLabIdxT) != length(wetLabIdxN)){
        message("Total Tumors: ", length(wetLabIdxT))
        message("Total Normals: ", length(wetLabIdxN))
        message("Wetlab Tumor/Normal samples are not paired correctly!")
        tumorLonger <- length(wetLabIdxT) > length(wetLabIdxN)
        totalAdd <- abs(length(wetLabIdxT) - length(wetLabIdxN))
        if(tumorLonger==F){
            totalDrop <- wetLabIdxN %in% wetLabIdxT
            #totalDrop <- 1:length(philipsIdxT)
            message(crayon::bgRed("More Normals than Tumors listed in worksheet!"))
            MsgDF(data.frame(Potential.Extra=philipsN[totalDrop]))
            wetLabIdxT <- c(wetLabIdxT, rep(0,totalAdd))
            message("Extra ", length(which(totalDrop==T)), " cases...")
        }else{
            totalDrop <- wetLabIdxT %in% wetLabIdxN
            message(crayon::bgRed("More Tumors than Normals listed in worksheet!"))
            MsgDF(data.frame(Potential.Extra=philipsT[!totalDrop]))
            wetLabIdxN <- c(wetLabIdxN, rep(0,totalAdd))
            message("Extra ", length(which(totalDrop==F)), " cases...")
        }

    }

    idx <- data.frame(
        philipsT = philipsIdxT,
        philipsN = philipsIdxN,
        wetLabT = wetLabIdxT,
        wetLabN = wetLabIdxN
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
        message(crayon::bgRed("The following rows are duplicated and will be removed:"), "\n")
        MsgDF(mainSheet[dupes,])
        mainSheet <- mainSheet[-dupes,]
        row.names(mainSheet) <- 1:nrow(mainSheet)
    }
    mainSheet$Tumor_Content[mainSheet$Paired_Normal==""] <- 0
    return(mainSheet)
}


BuildNoPhilips <- function(rawSheetData, runID, pact_run) {
    message("Generating SampleSheet without Philips Data: ALL cases are validation!")
    
    seqId <- stringr::str_split_fixed(runID, "_", 4)
    whichNormal <- rawSheetData$`Type & Tissue` %in% c("Normal", "Control")
    onlyNormals <- rawSheetData$`Type & Tissue` == "Normal"
    concat_id <- paste(0, runID, rawSheetData$`Accession#`, rawSheetData$`DNA #`, sep="_")
    
    mainSheet <- data.frame(
        Sample_ID = concat_id,
        Sample_Name = concat_id,
        Paired_Normal = ifelse(!whichNormal, concat_id[onlyNormals], ""),
        I7_Index_ID = rawSheetData$I7_Index_ID,
        index = rawSheetData$index,
        Specimen_ID = rawSheetData$`Accession#`,
        EPIC_ID = 0,
        Test_Number = rawSheetData$Test_Number,
        Tumor_Content = 0,
        Tumor_Type = ifelse(rawSheetData$`Tumor Type` == "NA", "", rawSheetData$`Tumor Type`),
        Description = rawSheetData$Description,
        Run_Number = runID,
        Sequencer_ID = seqId[1, 2],
        Chip_ID = seqId[1, 4],
        Sample_Project = pact_run,
        GenomeFolder = "PhiX-Illumina-RTA-Sequence-WholeGenomeFASTA",
        stringsAsFactors = FALSE
    )
    
    return(mainSheet)
}


BuildMainSheet <- function(philipsExport, rawSheetData, runID, pact_run) {
    epicOrder <- philipsExport[,'Epic Order Number']
    testNumber <- philipsExport[,'Test Number']
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

GetPhilipsData <- function(inputFi){
    shNames <- readxl::excel_sheets(inputFi)
    sh2 <- which(grepl("Philips", shNames, ignore.case = T))[1]
    philipsExport <- GetExcelData(inputFi, sh2, NULL, 3, cm=T)
    if(nrow(philipsExport)==0){
        warning("No PhilipsExport tab data found!")
        return(NULL)
    }
    philipsExport <- philipsExport[philipsExport$`Test Name`!="",]
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


GrabRunNumber <- function(inputFi){
    shNames <- readxl::excel_sheets(inputFi)
    sh <- which(grepl("PACT-", shNames, ignore.case = T))[1]
    if(is.na(sh)){
        sh = 1
    }
    rawSheetData <- GetExcelData(inputFi, sh, shRange="A6:X200", cm=T)
    runID_row <- which(rawSheetData$`DNA #`=="Run ID:")
    runID_col <- which(stringr::str_detect(paste(rawSheetData[runID_row, ]),"NB551709|NB501073|VH01471"))
    run_number <- try(paste(rawSheetData[runID_row, runID_col]), silent=T)
    if(length(run_number) == 0){
        message(crayon::bgRed('Keyword "Run ID:" not found in SampleSheet, defaulting to script input RUNID'))
        run_number <- ""
    }
    return(run_number[1])
}

GetRawSamplesheet <- function(inputFi){
    shNames <- readxl::excel_sheets(inputFi)
    sh <- which(grepl("PACT-", shNames, ignore.case = T))[1]
    if(is.na(sh)){
        message(crayon::bgRed('Did not detect "PACT" in Excel sheetnames, defaulting to reading sheet 1'))
        sh = 1
    }
    msgRd <- paste0('Reading Excel Sheet named \"', shNames[sh],'\" from file:')
    message(crayon::bgGreen(msgRd),'\n',inputFi)
    rawSheetData <- GetExcelData(inputFi, sh, shRange="A6:X200", cm=T)
    #toDrop <- which(rawSheetData[, "DNA #"]=="HAPMAP")[1]
    toDrop <- which(rawSheetData[,15] == "")[1] - 1
    rawSheetData <- rawSheetData[1:toDrop,]
    return(rawSheetData)
}

WriteMainSheet <- function(mainSheet, sheetHead){
    xlRunId <- mainSheet[1, "Run_Number"]
    runIdFi <- paste(xlRunId,  "SampleSheet.csv", sep = "-")
    outFile <- file.path(fs::path_home(), "Desktop", runIdFi)
    write.table(sheetHead, sep=",", file=outFile, row.names=F, col.names=F, quote=F)
    message("Writing file output, check your samplesheet here: ", outFile)
    if(any(mainSheet$Tumor_Type == "NA")){
        mainSheet[mainSheet$Tumor_Type == "NA", "Tumor_Type"] <- ""
    }
    suppressWarnings(write.table(mainSheet, sep=",", file=outFile, row.names=F, col.names=T, append=T, quote=F))
    return(outFile)
}

AltParseFormat <- function(inputFi, runID){
    rawSheetData <- GetRawSamplesheet(inputFi)
    sheetRunID <- GrabRunNumber(inputFi)
    if(sheetRunID != runID){
        message(crayon::bgRed("SampleSheet and Input RUNID do not match!"))
        message("Samplesheet RunID: ",'"', crayon::bgRed(sheetRunID),'"')
        message("Your input RunID: ", crayon::bgGreen(runID))
    }
    philipsExport <- GetPhilipsData(inputFi)
    pact_run <- stringr::str_split_fixed(base::basename(inputFi), ".xls", 2)[1,1]
    if(is.null(philipsExport)){
        mainSheet <- BuildNoPhilips(rawSheetData, sheetRunID, pact_run)
    }else{
        mainSheet <- BuildMainSheet(philipsExport, rawSheetData, sheetRunID, pact_run)
    }
    return(mainSheet)
}

# Parses xlsx file and writes as csv file -----
parseExcelFile <- function(inputFi, runID = NULL){
    shNames <- readxl::excel_sheets(inputFi)
    MsgDF(data.frame(`Sheet names in Workbook` = shNames))
    sh <- which(grepl("SampleSheet", shNames, ignore.case = T))[1]
    if(!is.na(sh)){
        sheetHead <- GetExcelData(inputFi, sheetNum=shNames[sh], shRange="A1:B17")
        sheetHead <- rbind(sheetHead,c("",""),c("[Data]",""))
        mainSheet <- GetExcelData(inputFi, shNames[sh], NULL, 19, T)
    } else{
        sheetHead <- GetSheetHeading(inputFi)
        mainSheet <- AltParseFormat(inputFi, runID)
    }
    mainSheet <- sanitizeSheet(mainSheet)
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
    message("\n", dsh2, "Email Notification Created", dsh2)
}


callApiFileCsv <- function(rcon, recordName, fiPath, ovwr = T) {

    if (ovwr == F) {
        writeLogFi(recordName)
    } else{
        fld <- "pact_csv_sheet"
        body <- list(
            token = rcon$token,
            content = 'file',
            action = 'import',
            record = recordName,
            field = fld,
            file = httr::upload_file(fiPath),
            returnFormat = 'csv'
        )
        res <-
            tryCatch(
                httr::POST(url = rcon$url, body = body, config = rcon$config),
                error = function(cond){
                    list(status_code = "200")
                }
            )
        if(res$status_code=="200"){
            message("REDCap file upload successful: ", fiPath)
        }else{
            message("REDCap file upload failed: ", fiPath)
        }
    }
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
    callApiFileCsv(rcon, runID, outFile, ovwr = T)
    tryCatch(
        redcapAPI::importFiles(rcon = rcon, file = outFile, record = runID, field = "pact_csv_sheet", repeat_instance = 1),
        error=function(e){message("REDCap file upload error failed:\n", e)}
    )
    emailNotify(record, rcon)
}


# Gets dataframe and saves as CSV file -----
writeSampleSheet <- function(inputSheet, token, runID = NULL) {
    isPath <- stringr::str_detect(inputSheet, .Platform$file.sep) == T
    isNxtSeq <- ifelse(stringr::str_detect(inputSheet, "2000"), T, F)
    if(isPath==F){
        inputFi <- getExcelPath(inputSheet, 1, isNxtSeq)
    }else{
        inputFi <- inputSheet
    }
    if (file.exists(inputFi)) {
        outVals <- parseExcelFile(inputFi, runID)
    } else {
        outVals <- CheckOtherFiles(inputFi, runID)
    }
    if(exists("outVals")){
        pushToRedcap(outVals, token)
    }
}

loadPacks()
checkMounts()

writeSampleSheet(inputSheet, token, runID)
