#!/usr/bin/env Rscript
## ---------------------------
## Script name: pactParse.R
## Purpose: source of global scripts and generate PACT -SampleSheet.csv file
## Author: Jonathan Serrano
## Version: 1.0.1
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------

library("base"); args <- commandArgs(TRUE); gb <- globalenv(); assign("gb", gb)

# Main arguments input in comandline (Uncomment to Debug or run Locally) -----------------------
args[1] -> token        #<- '8XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' = APITOKEN
args[2] -> inputSheet   #<- 'PACT-YY-##'
args[3] -> runID        #<- 'YYMMDD_NB######_0###_AH7ABCDEFG'

# Display & Verify the Input args --------------------------------------------------------------
message("\n================\nParameters input\n================")
message("token:    ", token)
message("PACT Run: ", inputSheet)
message("Run ID:   ", runID)
stopifnot(exists("token") & !is.null(token))
stopifnot(exists("inputSheet") & !is.null(inputSheet))

MsgDF <- function(data) {
    return(message(paste0(capture.output(as.data.frame(data)), collapse = "\n")))
}

boldRed <- function(txt) {return(crayon::white$bgRed$bold(txt))}


# FUN: Checks if z-drive is accessible to the Rscript ------------------------------------------
checkMounts <- function() {
    molecDrive = "/Volumes/molecular/MOLECULAR LAB ONLY"
    zDrive = "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular"
    failMount <- ifelse(dir.exists(molecDrive), T, F)
    if (failMount != T) {
        mountMsg <- paste(
            "PATH does not exist, ensure path is mounted:", boldRed(molecDrive),
            "You must mount the network Z-drive path:", boldRed(zDrive), sep = "\n"
        )
        message("\n", mountMsg, "\n")
        stopifnot(!any(failMount == T))
    } else {
      message("Z-drive path is accessible")
    }
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
    message(libName, "...loaded successful")
}


SilentInstall <- function(pkgs){
  rlis = getOption("repos")
  rlis["CRAN"] = "http://cran.us.r-project.org"
  options(repos = rlis)

  invisible(lapply(pkgs, function(pk) {
    if (suppressWarnings(!requireNamespace(pk, quietly = T))) {
      install.packages(pk, dependencies = T, verbose = T,
                       repos = "http://cran.us.r-project.org", type = "both")
    }
    libLoad(pk)
  }))
}


# Functions to load or install missing required packages ---------------------------------------
loadPacks <- function() {
    SilentInstall("devtools")

    if (!requireNamespace("redcapAPI", quietly = T)) {
        params = list('nutterb/redcapAPI', dependencies = T, upgrade = "always", type = "source")
        do.call(devtools::install_github, c(params))
    }
    libLoad("redcapAPI")

    pkgs <- c("foreach", "jsonlite", "RCurl", "readxl", "stringr")
    SilentInstall(pkgs)

    if (!(utils::packageVersion("redcapAPI") >= "2.7.4")) {
        install.packages("redcapAPI", ask = F, update = T, dependencies = T)
    }
}


MsgChangesMade <- function(mainSheet, badChars = c("[\r\n]", ",", " ")) {
    for (ptrn in badChars) {
        for (i in seq_along(mainSheet)) {
            col <- ifelse(is.na(mainSheet[[i]]), "", mainSheet[[i]])
            idxMatch <- which(stringr::str_detect(col, pattern = ptrn))
            if (length(idxMatch) > 0) {
                message(sprintf('The following "%s" will be removed from the column: %s',
                                ptrn, names(mainSheet)[i]))
                MsgDF(mainSheet[idxMatch, i])
            }
        }
    }
}


CountControls <- function(mainSheet, controlNames="NTC_H20|SC_SERACARE|NC_HAPMAP") {
    controlSamples <- grepl(pattern = controlNames, mainSheet$Sample_Name)
    if (table(controlSamples)[['TRUE']] != 3) {
        message(crayon::bgRed("There are not 3 Control samples"))
        message("Either NTC_H20, SC_SERACARE, or NC_HAPMAP is missing")
        message("OR\nThere are extra controls added in this run:")
        MsgDF(mainSheet$Sample_Name[controlSamples])
    }else{
        controlIndexes <- which(controlSamples == T)
        mainSheet[controlIndexes,'Paired_Normal'] <- ""
    }
    return(mainSheet)
}


# Removes and fixes newlines, commas, and blanks from samplesheet ------------------------------
sanitizeSheet <- function(mainSheet) {
    MsgChangesMade(data.frame(Tumor_Type = mainSheet$Tumor_Type), " ")
    mainSheet$Tumor_Type <- gsub(" ", "-", mainSheet$Tumor_Type)
    MsgChangesMade(mainSheet)
    mainSheet <- data.frame(lapply(mainSheet, function(col) gsub("[\r\n, ]", "", col)))

    pn_blank <- mainSheet$Paired_Normal == 0 | is.na(mainSheet$Paired_Normal)
    mainSheet$Paired_Normal[pn_blank] <- ""
    mainSheet$Tumor_Type[mainSheet$Tumor_Type == 0] <- "NA"
    mainSheet[, 1:16] <- sapply(mainSheet[, 1:16], function(x) {gsub("\\\\", "-", x)})
    mainSheet <- CountControls(mainSheet)

    if (any(duplicated(mainSheet$Sample_ID))) {
        message(crayon::bgRed("There are duplicated Sample_ID in the SampleSheet"))
        MsgDF(mainSheet$Sample_ID[duplicated(mainSheet$Sample_ID)])
    }

    return(mainSheet)
}


# Reads the genders and outputs a tsv file to /Molecular/MethylationClassifier/CNV_PNG ---------
WritePhilipsGender <- function(mainSheet, inputFi, shNames) {
    sh2 <- which(grepl("Philips", shNames, ignore.case = T))[1]
    cnvSheet <- mainSheet[,1:15]
    runID <- cnvSheet[1,"Sample_Project"]
    philipVals <- as.data.frame(readxl::read_excel(inputFi, sheet = sh2, skip = 3, col_types = "text"))
    cnvSheet$Gender <- philipVals$Gender[match(cnvSheet$Test_Number, philipVals$`Test Number`)]
    cnvPath <- "/Volumes/molecular/Molecular/MethylationClassifier/CNV_PNG"
    cnvPath <- file.path(cnvPath, paste0(runID,".tsv"))
    message("Writing table to: ", cnvPath)
    write.table(cnvSheet, quote = F, sep = '\t', file = cnvPath, row.names = F)
}


GetExcelData <- function(inputFi, sheetNum, shRange, toSkip = 0, cm = F) {
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


WriteFileHeader <- function(inputFi) {
    pactName <- stringr::str_split_fixed(basename(inputFi), ".xl", 2)[1, 1]
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


GetTypeIndex <- function(samNumber, rawData) {
    return(
        unlist(lapply(samNumber, function(samName) {return(which(samName == rawData$`DNA #`))}))
    )
}


CheckControlRows2 <- function(rawData, allBnumber) {
    controlIndex <- which(rawData$`DNA #` %in% allBnumber == F)
    message(crayon::bgRed("The following samples are Controls or not in Philips:"), "\n",
            paste(rawData$`DNA #`[controlIndex], collapse = "\n"), "\n")

    rawData$`Type & Tissue`[controlIndex] <- paste(rawData$original_tissue[controlIndex], "Filler")
    isControl <- unlist(lapply(rawData$`Accession#`, function(sam) {
        sam == "NTC" | sam == "SC" | sam == "NC"
        }))
    rawData$`Type & Tissue`[isControl] <- "Control"
    message(crayon::bgGreen("Assigned Values:"))
    message(paste(
        rawData$`DNA #`[controlIndex],
        rawData$`Type & Tissue`[controlIndex],
        sep = " = ",
        collapse = "\n"
    ))
    return(rawData)
}


FixPairedList <- function(philipsExport, rawData) {
    mislabelled <- philipsExport$`Tumor DNA/RNA Number` %in% philipsExport$`Normal DNA/RNA Number`
    if (any(mislabelled)) {
        message(crayon::bgRed("There are samples potentially mislabelled! in Philips:"),"\n")
        MsgDF(philipsExport$`Tumor DNA/RNA Number`[mislabelled])
    }
    tumorIndex <- GetTypeIndex(philipsExport$`Tumor DNA/RNA Number`, rawData)
    normalIndex <- GetTypeIndex(philipsExport$`Normal DNA/RNA Number`, rawData)
    rawData$original_tissue <- rawData$`Type & Tissue`
    rawData$`Type & Tissue` <- ""
    rawData$`Type & Tissue`[tumorIndex] <- "Tumor"
    rawData$`Type & Tissue`[normalIndex] <- "Normal"
    allBnumber <- c(philipsExport$`Tumor DNA/RNA Number`, philipsExport$`Normal DNA/RNA Number`)
    if (!all(allBnumber %in% rawData$`DNA #`) == T) {
        b_txt <- "Not all B-numbers from PhilipsExport tab are in the SampleSheet DNA # Column!"
        message(crayon::bgRed(b_txt))
    }
    rawData <- CheckControlRows2(rawData, allBnumber)
    return(rawData$`Type & Tissue`)
}


FindMissingPairs <- function(extraNormals, philipsExport, sheetPairTumor, sheetPairNorm, mainSheet) {
    unmatchedSamples <- if (extraNormals) {
        message(crayon::bgRed("More NORMALS than Tumors on this run!"))
        which(!philipsExport$`Tumor DNA/RNA Number` %in% sheetPairTumor)
    } else {
        message(crayon::bgRed("More TUMORS than Normals on this run!"))
        which(!philipsExport$`Normal DNA/RNA Number` %in% sheetPairNorm)
    }
    message("Manually edit sample sheet!")

    missingSams <-
        paste(if (extraNormals) {
            sheetPairNorm[unmatchedSamples]
        } else {
            sheetPairTumor[unmatchedSamples]
        }, collapse = "|")

    dnaMissingPair <- stringr::str_detect(mainSheet$Sample_ID, pattern = missingSams)

    return(dnaMissingPair)
}


MessageMismatched <- function(extraIndices, extraNormals, mainSheet, dnaMissingPair) {
    blu_miss <- paste(
        "The following case(s) are missing a paired",
        if (extraNormals) {"Tumor:"} else {"Normal:"}
    )
    message(crayon::bgBlue(blu_miss))
    MsgDF(mainSheet$Sample_ID[dnaMissingPair])
    msg_tum_norm <-
        paste(length(extraIndices), "extra", if (extraNormals) {
            "Normal(s) may need to be paired manually to missing Tumor(s)"
        } else {
            "Tumor(s) may need to be paired manually to missing Normal(s)"
        })
    message(crayon::bgRed(msg_tum_norm))
}


ProcessExtraNormals <- function(normSamPair, missingCases, tumorSamPair, mainSheet, sheetTumors) {
    pairedNormals <- which(!normSamPair %in% missingCases)
    extraNormals <- normSamPair[-pairedNormals]

    message("This normal sample will not be paired with a tumor:\n", extraNormals)

    indexToAdd <- which(normSamPair %in% missingCases)
    indexToAdd <- sort(indexToAdd, decreasing = FALSE)

    # Insert normal as an unpaired tumor at this index
    for (idx in indexToAdd) {
        tumorSamPair <- append(tumorSamPair, normSamPair[idx], after = idx - 1)
    }

    tumorRowAdd <- which(mainSheet$Sample_ID %in% missingCases)
    indexToAdd <- sort(tumorRowAdd, decreasing = FALSE)

    for (idx in indexToAdd) {
        sheetTumors <- append(sheetTumors, idx)
    }

    for (xSam in extraNormals) {
        toDrop <- normSamPair == xSam
        normSamPair[toDrop] <- ""
    }

    mainSheet$Paired_Normal[sheetTumors] <- normSamPair
    mainSheet$Sample_ID[sheetTumors] <- mainSheet$Sample_ID[sheetTumors] <- tumorSamPair
    return(mainSheet)
}


# Insert "missing" at this index to manually pair
AppendMissingNotes <- function(indexToAdd, normSamPair, mainSheet, sheetTumors) {
    for (idx in indexToAdd) {
        normSamPair <-
            append(normSamPair, "MISSING_NORMAL_PAIR", after = idx - 1)
    }
    mainSheet$Paired_Normal[sheetTumors] <- normSamPair
    return(mainSheet)
}


AddUnmatchNormals <- function(sheetPairNorm, sheetPairTumor, mainSheet,
                              sheetTumors, sheetNormals, philipsExport) {
    extraNormals <- length(sheetPairNorm) > length(sheetPairTumor)
    dnaMissingPair <- FindMissingPairs(extraNormals, philipsExport, sheetPairTumor,
                                       sheetPairNorm, mainSheet)

    if (any(dnaMissingPair)) {
        extraIndices <- which(dnaMissingPair)
        MessageMismatched(extraIndices, extraNormals, mainSheet, dnaMissingPair)

        normSamPair <- mainSheet$Sample_ID[sheetNormals]
        tumorSamPair <- mainSheet$Sample_ID[sheetTumors]
        missingCases <- mainSheet$Sample_ID[dnaMissingPair]

        if (extraNormals == T) {
            mainSheet <- ProcessExtraNormals(normSamPair, missingCases, tumorSamPair,
                                             mainSheet, sheetTumors)
            return(mainSheet)
        } else{
            indexToAdd <- which(tumorSamPair %in% missingCases)
            mainSheet <- AppendMissingNotes(indexToAdd, normSamPair, mainSheet, sheetTumors)
            return(mainSheet)
        }
    }
    return(mainSheet)
}


AddSampleIndexes <- function(pairedList, rawData, philipsExport) {
        tissueType <- FixPairedList(philipsExport, rawData)
        sheetTumors <-
            which(stringr::str_detect(tissueType, pattern = "(?i)tumor"))
        sheetNormals <-
            which(stringr::str_detect(tissueType, pattern = "(?i)norm"))
        sheetControl <-
            which(stringr::str_detect(tissueType, pattern = "(?i)cont"))

        mainSheet <-
            data.frame(matrix("", nrow = nrow(pairedList), ncol = 0))
        mainSheet$Sample_Name <-
            mainSheet$Sample_ID <- paste(pairedList[, 1])
        mainSheet$Paired_Normal <- ""

        sheetPairTumor <- rawData$`DNA #`[sheetTumors]
        sheetPairNorm <- rawData$`DNA #`[sheetNormals]

        # Check if length mismatch between tumor normals and normal samples
        if (length(sheetPairTumor) != length(sheetPairNorm)) {
            mainSheet <-
                AddUnmatchNormals(
                    sheetPairNorm,
                    sheetPairTumor,
                    mainSheet,
                    sheetTumors,
                    sheetNormals,
                    philipsExport
                )
        } else {
            mainSheet$Paired_Normal[sheetTumors] <- mainSheet$Sample_ID[sheetNormals]
        }

        if (length(rawData$I7_Index_ID) != nrow(mainSheet)) {
            stop("length(rawData$I7_Index_ID) != nrow(mainSheet)")
        }
        mainSheet$I7_Index_ID <- paste(rawData$I7_Index_ID)
        mainSheet$index <- paste(rawData$index)
        mainSheet$Specimen_ID <- paste(rawData$`Accession#`)
        mainSheet$Tumor_Content <-
            mainSheet$Test_Number <- mainSheet$EPIC_ID <- "0"
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


CheckControlRows <- function(rawData, runID, newRows) {
    accessions <- rawData$`Accession#`
    dnaNumbers <- rawData$`DNA #`
    controls <- which(stringr::str_detect(rawData$`Type & Tissue`, "Cont|cont"))
    controlRows <-
        paste0(rawData$Test_Number[controls], "_",
               runID, "_", accessions[controls], "_", dnaNumbers[controls])
    hasControls <- controlRows %in% newRows

    if (any(duplicated(controlRows))) {
      controlRows <- FixDuplicateControls(controlRows)
      message(crayon::bgBlue("New control names:"), paste(controlRows, collapse = "\n"))
    }

    if (!any(hasControls)) {
      message(crayon::bgBlue("Binding controls:"),"\n", paste(controlRows, collapse = "\n"))
    } else{
      controlRows <- NULL
    }

    return(controlRows)
}


CheckMissingRows <- function(pairedList, rawData, runID) {
    dnaNumbers <- rawData$`DNA #`
    no_control <- stringr::str_detect(dnaNumbers, "H20|SERACARE|HAPMAP", negate = T)
    accessions <- rawData$`Accession#`[no_control]
    tst_number <- rawData$Test_Number[no_control]
    dnaNumbers <- dnaNumbers[no_control]
    missRows <- vapply(dnaNumbers, function(x) {!any(stringr::str_detect(pairedList, x))}, logical(1))
    if (any(missRows)) {
        message(crayon::bgGreen("Binding additional rows/filler:"), "\n",
                paste(accessions[missRows], collapse = "\n"))
        extraRow <- dnaNumbers %in% dnaNumbers[missRows]
        tst_number[missRows] <- 0
        newRows <- paste(tst_number[missRows],
                         runID,
                         accessions[extraRow],
                         dnaNumbers[extraRow], sep = "_")
    } else {
        message(crayon::bgGreen("No additional fillers or paired sample rows to bind to sample sheet"))
        newRows <- NULL
    }
    return(newRows)
}


CheckMissingData <- function(pairedList) {
    doubleBars <- vapply(pairedList, function(X) {
        stringr::str_detect(X, pattern = "__")}, logical(1))
    if (any(doubleBars)) {
        message(crayon::bgRed("The following rows are missing data:"),"\n",
                paste(pairedList[doubleBars], collapse = "\n"))
    }
}


BindUnpairedRows <- function(rawData, pairedList, runID) {
    CheckMissingData(pairedList)
    newRows <- CheckMissingRows(pairedList, rawData, runID)
    controlRows <- CheckControlRows(rawData, runID, newRows)
    new_paired_list <- data.frame("Sample_ID" = c(1:nrow(rawData)))
    b_numbers <- rawData$`DNA #`
    all_samples <- c(pairedList, newRows, controlRows)
    for (sam in 1:length(b_numbers)) {
        sam_idx <- which(stringr::str_detect(all_samples, b_numbers[sam]))
        if (length(sam_idx) > 1) {
            for(i in 1:length(sam_idx)) {
                currSam <- sam_idx[i]
                matchedIdx <- which(b_numbers == b_numbers[sam])
                newIdx <- matchedIdx[i]
                new_paired_list[newIdx,] <- all_samples[currSam]
            }
        }else{
            new_paired_list[sam,] <- all_samples[sam_idx]
        }
    }
    rownames(new_paired_list) <- seq_len(nrow(new_paired_list))
    return(new_paired_list)
}


GetPairedList <- function(philipsExport, runID) {
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


MatchIndex <- function(list1, list2) {
    return(unlist(lapply(list1, function(x) {which(x == list2)})))
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
    }

    return(matchIdx)
}


CheckTotalIndexes <- function(tumorSam, normalSam, ngsNumbers, sheetType = "Philips") {
    if (length(tumorSam) != length(normalSam)) {
        message("Total Tumors: ", length(tumorSam))
        message("Total Normals: ", length(normalSam))
        tumorLonger <- length(tumorSam) > length(normalSam)
        totalAdd <- abs(length(tumorSam) - length(normalSam))
        message(crayon::bgRed(paste("Not all",sheetType,"samples are paired with a tumor or normal!")))
        if (tumorLonger == F) {
            totalDrop <- normalSam %in% tumorSam
            message(crayon::bgRed(paste(totalAdd, "more Normals than Tumors listed in worksheet!")))
            MsgDF(data.frame(
                NGS = ngsNumbers[!totalDrop], Potential.Extra = philipsN[!totalDrop]))
            tumorSam <- c(tumorSam, rep(0,totalAdd))
        }else{
            totalDrop <- tumorSam %in% normalSam
            message(crayon::bgRed(paste(totalAdd, "more Tumors than Normals listed in worksheet!")))
            MsgDF(data.frame(
                NGS = ngsNumbers[!totalDrop], Potential.Extra = philipsT[!totalDrop]))
            normalSam <- c(normalSam, rep(0,totalAdd))
        }
        message("Extra ", length(which(totalDrop == F)),  " cases...")
    }
}


MsgTumNormErr <- function(tums, norms, sheetSide = "Philips"){
    message("Total Tumors: ", length(tums))
    message("Total Normals: ", length(norms))
    if (sheetSide == "Philips") {
        message(crayon::bgRed("Not all Philips samples are paired with a Tumor/Normal!"))
    }else{
        message(crayon::bgRed("Wetlab Tumor/Normal samples are not paired correctly!"))
    }
}


MsgTotalAdded <- function(totalAdd, ngsNumbers, totalDrop, samExtra, tumorLonger) {
    if (tumorLonger) {
        sheetMessage <- "more Tumors than Normals listed in worksheet!"
    } else{
        sheetMessage <- "more Normals than Tumors listed in worksheet!"
    }
    message(crayon::bgRed(paste(totalAdd, sheetMessage)))
    MsgDF(data.frame(NGS = ngsNumbers[!totalDrop], Potential.Extra = samExtra[!totalDrop]))
    message("Extra ", length(which(totalDrop == F)), " cases...")
}


GetIndexMatch <- function(rawData, philipsExport) {
    message("Running GetIndexMatch function in pactParse.R")
    accessions <- rawData$`Accession#`
    philipsN <- c(philipsExport$`Normal Specimen ID`)
    philipsT <- c(philipsExport$`Tumor Specimen ID`)
    ngsNumbers <- philipsExport$`Test Number`

    CheckMissingPairs(ngsNumbers, philipsT, philipsN)

    philipsIdxT <- CheckDupesMatch(philipsT, accessions)
    philipsIdxN <- CheckDupesMatch(philipsN, accessions)

    wetLabIdxT <- CheckDupesMatch(accessions, philipsT)
    wetLabIdxN <- CheckDupesMatch(accessions, philipsN)

    if (length(philipsIdxT) != length(philipsIdxN)) {
        MsgTumNormErr(philipsIdxT, philipsIdxN)
        tumorLonger <- length(philipsIdxT) > length(philipsIdxN)
        totalAdd <- abs(length(philipsIdxT) - length(philipsIdxN))
        if (tumorLonger == F) {
            totalDrop <- philipsIdxN %in% philipsIdxT
            philipsIdxT <- c(philipsIdxT, rep(0,totalAdd))
            MsgTotalAdded(totalAdd, ngsNumbers, totalDrop, philipsN, tumorLonger)
        }else{
            totalDrop <- philipsIdxT %in% philipsIdxN
            philipsIdxN <- c(philipsIdxN, rep(0, totalAdd))
            MsgTotalAdded(totalAdd, ngsNumbers, totalDrop, philipsT, tumorLonger)
        }
    }

    if (length(wetLabIdxT) != length(wetLabIdxN)) {
        MsgTumNormErr(wetLabIdxT, wetLabIdxN, "WetLab")
        tumorLonger <- length(wetLabIdxT) > length(wetLabIdxN)
        totalAdd <- abs(length(wetLabIdxT) - length(wetLabIdxN))
        if (tumorLonger == F) {
            totalDrop <- wetLabIdxN %in% wetLabIdxT
            wetLabIdxT <- c(wetLabIdxT, rep(0,totalAdd))
            MsgTotalAdded(totalAdd, ngsNumbers, totalDrop, wetLabIdxN, tumorLonger)
        } else{
            totalDrop <- wetLabIdxT %in% wetLabIdxN
            wetLabIdxN <- c(wetLabIdxN, rep(0, totalAdd))
            MsgTotalAdded(totalAdd, ngsNumbers, totalDrop, wetLabIdxT, tumorLonger)
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


AddRunChipColumns <- function(mainSheet, runID, pact_run) {
    mainSheet$Run_Number <- runID
    splitRun <- stringr::str_split_fixed(runID, "_", 4)
    mainSheet$Sequencer_ID <- splitRun[1, 2]
    mainSheet$Chip_ID <- splitRun[1, 4]
    mainSheet$Sample_Project <- paste0(pact_run)
    return(mainSheet)
}


AddPairedColumn <- function(mainSheet, sheetColumn, philipsColumn, idx) {
    mainSheet[idx$wetLabN, sheetColumn] <- philipsColumn[idx$philipsN]
    mainSheet[idx$wetLabT, sheetColumn] <- philipsColumn[idx$philipsT]
    return(mainSheet)
}


FixLastColumns <- function(mainSheet, rawData) {
    mainSheet$Tumor_Type[rawData$`Type & Tissue` == "Normal"] <- ""
    mainSheet$Description <- paste0(rawData[,'Description'])
    mainSheet$GenomeFolder <- as.character("PhiX\\Illumina\\RTA\\Sequence\\WholeGenomeFASTA")
    dupes <- base::anyDuplicated(mainSheet$I7_Index_ID)
    if (length(dupes) > 0 & dupes != 0) {
        message(crayon::bgRed("The following rows are duplicated and will be removed:"), "\n")
        MsgDF(mainSheet[dupes,])
        mainSheet <- mainSheet[-dupes,]
        row.names(mainSheet) <- 1:nrow(mainSheet)
    }
    mainSheet$Tumor_Content[mainSheet$Paired_Normal == ""] <- 0
    return(mainSheet)
}


BuildNoPhilips <- function(rawData, runID, pact_run) {
    message("Generating SampleSheet without Philips Data: ALL cases are validation!")
    seqId <- stringr::str_split_fixed(runID, "_", 4)
    whichNormal <- stringr::str_detect(rawData$`Type & Tissue`, "Norm|norm|cont|Cont")
    onlyNormals <- stringr::str_detect(rawData$`Type & Tissue`, "Norm|norm")
    concat_id <-
        data.frame(sid = paste(0, runID, rawData$`Accession#`, rawData$`DNA #`,
                               sep = "_"))
    pairedNorm <- rep("", nrow(concat_id))
    pairedNorm[!whichNormal] <- concat_id[onlyNormals, 1]

    rawData$`Tumor Type` <- ""
    mainSheet <- data.frame(
        Sample_ID = concat_id[,1],
        Sample_Name = concat_id[,1],
        Paired_Normal = pairedNorm,
        I7_Index_ID = rawData$I7_Index_ID,
        index = rawData$index,
        Specimen_ID = rawData$`Accession#`,
        EPIC_ID = 0,
        Test_Number = rawData$Test_Number,
        Tumor_Content = 0,
        Tumor_Type = rawData$`Tumor Type`,
        Description = rawData$Description,
        Run_Number = runID,
        Sequencer_ID = seqId[1, 2],
        Chip_ID = seqId[1, 4],
        Sample_Project = pact_run,
        GenomeFolder = "PhiX-Illumina-RTA-Sequence-WholeGenomeFASTA",
        stringsAsFactors = FALSE
    )
    mainSheet[!whichNormal, "Tumor_Content"] <- 70
    return(mainSheet)
}


BuildMainSheet <- function(philipsExport, rawData, runID, pact_run) {
    epicOrder <- philipsExport[,'Epic Order Number']
    testNumber <- philipsExport[,'Test Number']
    tumorPercent <- philipsExport[,'Tumor Percentage']
    diagColumn <- philipsExport[,'Diagnosis for interpretation']
    pairedList <- GetPairedList(philipsExport, runID)
    pairedList <- BindUnpairedRows(rawData, pairedList, runID)
    mainSheet <- AddSampleIndexes(pairedList, rawData, philipsExport)
    idx <- GetIndexMatch(rawData, philipsExport)
    mainSheet <- AddPairedColumn(mainSheet, "EPIC_ID", epicOrder, idx)
    mainSheet <- AddPairedColumn(mainSheet, "Test_Number", testNumber, idx)
    mainSheet <- AddRunChipColumns(mainSheet, runID, pact_run)
    mainSheet <- AddPairedColumn(mainSheet, "Tumor_Content", tumorPercent, idx)
    mainSheet <- AddPairedColumn(mainSheet, "Tumor_Type", diagColumn, idx)
    mainSheet <- FixLastColumns(mainSheet, rawData)
    return(mainSheet)
}


GetSheetHeading <- function(inputFi) {
    sheetHead <- WriteFileHeader(inputFi)
    sheetHead <- rbind(sheetHead, c("", ""), c("[Data]", ""))
    return(sheetHead)
}


GetPhilipsData <- function(inputFi) {
    shNames <- readxl::excel_sheets(inputFi)
    sh2 <- which(grepl("Philips", shNames, ignore.case = T))[1]
    philipsExport <- GetExcelData(inputFi, sh2, NULL, 3, cm = T)
    if (nrow(philipsExport) == 0) {
        warning("No PhilipsExport tab data found!")
        return(NULL)
    }
    if (length(colnames(philipsExport)) == 0) {
      stop("PhilipsExport tab is missing headers!")
    }
    philipsExport <- philipsExport[philipsExport$`Test Name` != "", ]
    filterColumns <- GetPhilipsColumns()
    philipsExport <- philipsExport[,filterColumns]
    blankOrder <- philipsExport[,"Epic Order Number"] == ""
    if (any(blankOrder)) {
        warning("Some Philips Samples are missing MRNs, defaulting to 0")
        philipsExport[, "Epic Order Number"][blankOrder] <- 0
    }
    blankSpecimen <- philipsExport[, "Normal Specimen ID"] == ""
    if (any(blankSpecimen)) {
        missingSam <- philipsExport$`Tumor Specimen ID`[blankSpecimen]
        message(crayon::bgRed("The following Philips Samples are missing a Normal Specimen ID:\n"),
                paste0(capture.output(missingSam), collapse = "\n"))
    }
    philipsExport[philipsExport == ""] <- 0
    return(philipsExport)
}


Err_runID <- function(runID) {
    err_msg <- crayon::bgRed('Keyword "Run ID:" not found in SampleSheet "DNA #" column')
    new_id <- paste('\nDefaulting to script input RUNID:', runID)
    message(err_msg, new_id)
}


GrabRunNumber <- function(inputFi, runID) {
    shNames <- readxl::excel_sheets(inputFi)
    sh <- which(grepl("PACT-", shNames, ignore.case = TRUE))[1]
    sh <- ifelse(is.na(sh), 1, sh)

    rawData <- GetExcelData(inputFi, sh, shRange = "A6:X200", cm = T)
    if (!"Test_Number" %in% colnames(rawData)) {
        stop("Input Sheet formatting is incorrect, column names are:\n",
             paste(colnames(rawData), collapse = "\n"))
    }
    lastRow <- which(rawData$Test_Number == 0)[1]
    sheetBottom <- rawData[c(lastRow:nrow(rawData)), ]

    # Locate 'Run ID:' and extract relevant run number
    runID_row <- which(stringr::str_detect(sheetBottom$`DNA #`, "Run ID:"))

    # If 'Run ID:' keyword not found, use provided runID as default
    if (length(runID_row) != 1) {
        Err_runID(runID)
        return(runID)
    }

    runID_find <- paste(sheetBottom[runID_row, ])
    runID_col <- which(stringr::str_detect(runID_find, "NB551709|NB501073|VH01471"))

    if (length(runID_col) == 0) {
        Err_runID(runID)
        return(runID)
    }
    run_number <- paste(sheetBottom[runID_row, runID_col])
    run_number <- stringr::str_replace_all(run_number, c("Run ID:" = "", " " = ""))
    return(run_number)
}


GetRawSamplesheet <- function(inputFi) {
    shNames <- readxl::excel_sheets(inputFi)
    sh <- which(grepl("PACT-", shNames, ignore.case = T))[1]
    if (is.na(sh)) {
        message(crayon::bgRed('Did not detect "PACT" in Excel sheetnames, defaulting to reading sheet 1'))
        sh = 1
    }
    msgRd <- paste0('Reading Excel Sheet named \"', shNames[sh], '\" from file:')
    message(crayon::bgGreen(msgRd),'\n',inputFi)
    rawData <- GetExcelData(inputFi, sh, shRange = "A6:X200",  cm = T)
    #toDrop <- which(rawData[, "DNA #"]=="HAPMAP")[1]
    toDrop <- which(rawData[, 15] == "")[1] - 1
    rawData <- rawData[1:toDrop, ]
    rawData$`DNA #` <- stringr::str_replace_all(rawData$`DNA #`, "_", "-")
    return(rawData)
}

WriteMainSheet <- function(mainSheet, sheetHead) {
    xlRunId <- mainSheet[1, "Run_Number"]
    runIdFi <- paste(xlRunId,  "SampleSheet.csv", sep = "-")
    outFile <- file.path(fs::path_home(), "Desktop", runIdFi)
    write.table(sheetHead, sep = ",", file = outFile, row.names = F, col.names = F, quote = F)
    message("Writing file output, check your samplesheet here: ", outFile)
    if (any(mainSheet$Tumor_Type == "NA")) {
        mainSheet[mainSheet$Tumor_Type == "NA", "Tumor_Type"] <- ""
    }
    suppressWarnings(
      write.table(
        mainSheet,
        sep = ",",
        file = outFile,
        row.names = F,
        col.names = T,
        append = T,
        quote = F
      )
    )
    return(outFile)
}


AltParseFormat <- function(inputFi, runID) {
    rawData <- GetRawSamplesheet(inputFi)
    sheetRunID <- GrabRunNumber(inputFi, runID)
    if (sheetRunID != runID) {
        message(crayon::bgRed("SampleSheet and Input RUNID do not match!"))
        message("Samplesheet RunID: ",'"', crayon::bgRed(sheetRunID),'"')
        message("Your input RunID: ", crayon::bgGreen(runID))
    }
    philipsExport <- GetPhilipsData(inputFi)
    pact_run <- stringr::str_split_fixed(base::basename(inputFi), ".xls", 2)[1,1]
    message(crayon::bgBlue(paste("PACT ID is:", pact_run)))
    if (is.null(philipsExport)) {
        mainSheet <- BuildNoPhilips(rawData, sheetRunID, pact_run)
    }else{
        mainSheet <- BuildMainSheet(philipsExport, rawData, sheetRunID, pact_run)
    }
    return(mainSheet)
}

# Parses xlsx file and writes as csv file -----
parseExcelFile <- function(inputFi, runID = NULL) {
    shNames <- readxl::excel_sheets(inputFi)
    MsgDF(data.frame(`Sheet names in Workbook` = shNames))
    sh <- which(grepl("SampleSheet", shNames, ignore.case = T))[1]
    if (!is.na(sh)) {
        sheetHead <- GetExcelData(inputFi, sheetNum = shNames[sh], shRange = "A1:B17")
        sheetHead <- rbind(sheetHead,c("",""),c("[Data]",""))
        mainSheet <- GetExcelData(inputFi, shNames[sh], NULL, 19, T)
    } else{
        sheetHead <- GetSheetHeading(inputFi)
        mainSheet <- AltParseFormat(inputFi, runID)
    }
    mainSheet <- sanitizeSheet(mainSheet)
    try(WritePhilipsGender(mainSheet,inputFi, shNames), silent = T)
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


CheckOtherFiles <- function(inputFi, runID) {
    notFoundMsg <- crayon::bgRed("The PACT run worksheet was not found:")
    message(notFoundMsg,"\n", inputFi, "\n")
    parentFolder <- dirname(inputFi)
    potentialFi <- list.files(path = parentFolder, full.names = T)
    if (length(potentialFi) > 1) {
        message(crayon::bgRed("Checking other existing files:"), "\n")
        message(paste(potentialFi, sep = "\n"))
        potentialFi <- filterFiles(potentialFi)
    }else{
        stop(paste("No PACT worksheet was found in the directory:", parentFolder))
    }
    if (file.exists(potentialFi[1])) {
        message(crayon::bgGreen("Now trying to read:"), "\n", potentialFi[1], "\n")
        pfile <- potentialFi[1]
    }else{
        nextMsg <- paste(notFoundMsg, potentialFi[1], "Trying:", potentialFi[2], sep = "\n")
        message(nextMsg)
        pfile <- potentialFi[2]
    }
    # TODO: Write check for mismatched year
    # runyr <- stringr::str_split_fixed(inputFi, "-", 3)[, 2]
    stopifnot(file.exists(pfile))
    outVals <- parseExcelFile(inputFi = pfile, runID)
    return(outVals)
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
emailNotify <- function(record, rcon) {
    record$pact_csv_email <- "pact_csv_email"
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox = T)
    PostRedcapCurl(rcon, datarecord)
    newDash <- "\n==========================\n"
    message(sprintf("\n%sEmail Notification Created%s", newDash, newDash))
}


callApiFileCsv <- function(rcon, recordName, fiPath, ovwr = T) {

    if (ovwr == F) {
        gb$writeLogFi(recordName)
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
        res <- tryCatch(
          httr::POST(url = rcon$url, body = body, config = rcon$config),
          error = function(e) {
            message(e)
            return(list(status_code = "400"))
            }
          )
        if (res$status_code == "200") {
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
    PostRedcapCurl(rcon, datarecord, retcon = 'nothing')
    callApiFileCsv(rcon, runID, outFile, ovwr = T)
    tryCatch(
      redcapAPI::importFiles(
        rcon = rcon,
        file = outFile,
        record = runID,
        field = "pact_csv_sheet",
        repeat_instance = 1
      ),
        error = function(e) {message("REDCap file upload error failed:\n", e)}
    )
    emailNotify(record, rcon)
}


determineRunType <- function(inputSheet) {
    if (stringr::str_detect(inputSheet, "2000|NextSeq550")) {
        runType <- "NextSeq2000"
    } else if (stringr::str_detect(inputSheet, "TEST")) {
        runType <- "test"
    } else if (stringr::str_detect(inputSheet, "ICL")) {
        runType <- "Illumina"
    } else {
        runType <- "regular"
    }

    return(runType)
}


# FUN: Returns Path to xlsx file ---------------------------------------------------------------
getExcelPath <- function(inputSheet, runType) {
    if (stringr::str_detect(inputSheet, .Platform$file.sep)) {
        return(inputSheet)
    }
    drive = file.path("", "Volumes", "molecular", "MOLECULAR LAB ONLY")
    folder = file.path("NYU PACT Patient Data", "Workbook")
    runyr <- stringr::str_split_fixed(inputSheet, "-", 3)[, 2]

    yearDir <- paste0("20", runyr)

    xlFi <- paste0(inputSheet, ".xlsm")
    worksheetPath <- file.path(drive, folder, yearDir, inputSheet, xlFi)

    if (runType == "NextSeq2000") {
        folder <- "Validations/PACT new i7-NextSeq2000/Wet Lab/Workbook"
        worksheetPath <- file.path(drive, folder, yearDir, xlFi)
    }

    if (runType == "test") {
        drive = "/Volumes/molecular/Molecular/Validation/PACT/Test_Sheets"
        worksheetPath <- file.path(drive, xlFi)
    }

    return(worksheetPath)
}


# Gets dataframe and saves as CSV file -----
writeSampleSheet <- function(inputSheet, token, runID = NULL) {
    outVals <- NULL
    runType <- determineRunType(inputSheet)
    worksheetPath <- getExcelPath(inputSheet, runType)

    if (file.exists(worksheetPath)) {
        outVals <- parseExcelFile(worksheetPath, runID)
    } else {
        outVals <- CheckOtherFiles(worksheetPath, runID)
    }

    if (!is.null(outVals)) {
        pushToRedcap(outVals, token)
    }
}

loadPacks()
checkMounts()

writeSampleSheet(inputSheet, token, runID)
