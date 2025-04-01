#!/usr/bin/env Rscript
## Script name: pactParse.R
## Purpose: source of global scripts and generate PACT -SampleSheet.csv file
## Author: Jonathan Serrano
## Date Created: August 16, 2021
## Version: 1.0.1
## Copyright (c) NYULH Jonathan Serrano, 2024

library("base"); args <- commandArgs(TRUE); gb <- globalenv(); assign("gb", gb)

# Function to assign parameter values based on args ---------------------------
assign_param <- function(args, idx) {
    if (idx == 4) {
        if (length(args) >= idx) {
            return(args[idx])
        } else {
            return("NONE")
        }
    }
    if (length(args) >= idx) args[idx] else stop("arg ", idx, " is NULL")
}

# Functions to print parameters  ) --------------------------------------------
brick_message <- function(phrase) {
    dsh <- paste(rep("-", nchar(phrase)), collapse = "")
    message(paste(dsh, phrase, dsh, sep = "\n"))
}

print_parameters <- function(token, input, runID, keywd) {
    brick_message("Parameters input")
    message("token:   ", token)
    message("PACT ID: ", input)
    message("Run ID:  ", runID, "\n")
    message("ILMN-#:  ", keywd, "\n")
}

# Assign values to variables & print ------------------------------------------
token <- assign_param(args, 1)
input <- assign_param(args, 2)
runID <- assign_param(args, 3)
keywd <- assign_param(args, 4)
print_parameters(token, input, runID, keywd)

# Helper functions to message as data frame and red background ----------------
MsgDF <- function(...) {
    return(message(paste0(
        capture.output(as.data.frame(...)), collapse = "\n"
    )))
}
mkRed <- function(...) {return(crayon::white$bgRed$bold(...))}

# FUN: Checks if z-drive is accessible to the RScript -------------------------
checkMounts <- function() {
    molecDrive = "/Volumes/molecular/MOLECULAR LAB ONLY"
    zDrive = "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular"
    failMount <- ifelse(dir.exists(molecDrive), T, F)
    if (failMount != T) {
        mountMsg <- paste(
            "PATH does not exist, ensure path is mounted:", mkRed(molecDrive),
            "You must mount the network Z-drive path:", mkRed(zDrive),
            sep = "\n"
        )
        message("\n", mountMsg, "\n")
        stopifnot(!any(failMount == T))
    } else {
        message("Z-drive path is accessible")
    }
}

# FUN: Returns hardcoded Philips Export Tab Column Names-----------------------
GetPhilipsColumns <- function() {
    filterColumns <- c(
        "MRN", "Test Name", "Tumor Specimen ID", "Normal Specimen ID",
        "Test Number","Epic Order Number", "Diagnosis for interpretation",
        "Tumor DNA/RNA Number", "Normal DNA/RNA Number","Tumor Percentage"
    )
    return(filterColumns)
}


# Function to silently load library without conflict warnings -----------------
libLoad <- function(libName) {
    lib.opts <- list(package = libName, character.only = T,
                     verbose = F, warn.conflicts = F)
    suppressWarnings(suppressPackageStartupMessages(do.call(library, c(lib.opts))))
    message(libName, "...loaded successful")
}

# Installs package(s) if not already installed then loads the library ---------
SilentInstall <- function(pkgs) {
    rlis = getOption("repos")
    rlis["CRAN"] = "http://cran.us.r-project.org"
    options(repos = rlis)

    invisible(lapply(pkgs, function(pk) {
        if (suppressWarnings(!requireNamespace(pk, quietly = T))) {
            install.packages(pk, dependencies = T, verbose = T,
                             repos = "http://cran.us.r-project.org",
                             type = "both")
        }
        libLoad(pk)
    }))
}


append_log <- function(txt, logFile = NULL) {
    if (is.null(logFile)) {
        logFile = paste0(Sys.Date(), "_run_log.tsv")
    }
    if (stringr::str_detect(logFile, .Platform$file.sep, T)) {
        logFile <- file.path(fs::path_home(), logFile)
    }
    message(paste("Check", logFile))
    write.table(
        txt,
        logFile,
        append = T,
        quote = F,
        sep = '\t',
        row.names = F,
        col.names = F
    )
}


write_to_log <- function(recordName, logFile = NULL) {
    if (is.null(logFile)) {
        logFile = paste0(Sys.Date(), "_pact_log.tsv")
    }
    txt = paste(recordName[1], "already has an data in REDCap!\n")
    append_log(txt, logFile)
}


# Functions to load or install missing required packages ----------------------
loadPacks <- function() {
    SilentInstall("devtools")

    if (!requireNamespace("redcapAPI", quietly = T)) {
        params = list('nutterb/redcapAPI', dependencies = T,
                      upgrade = "never", type = "source")
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


CountControls <- function(mainSheet,
                          controlNames = "NTC_H20|SC_SERACARE|NC_HAPMAP") {
    controlSamples <- grepl(pattern = controlNames, mainSheet$Sample_Name)
    if (table(controlSamples)[['TRUE']] != 3) {
        message(crayon::bgRed("There are not 3 Control samples"))
        message("Either NTC_H20, SC_SERACARE, or NC_HAPMAP is missing")
        message("OR\nThere are extra controls added in this run:")
        MsgDF(mainSheet$Sample_Name[controlSamples])
    } else{
        controlIndexes <- which(controlSamples == T)
        mainSheet[controlIndexes, 'Paired_Normal'] <- ""
    }
    return(mainSheet)
}


# Removes and fixes newlines, commas, and blanks from samplesheet -------------
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


# Outputs gender to a tsv file in /Molecular/MethylationClassifier/CNV_PNG ----
WritePhilipsGender <- function(mainSheet, inputFi, shNames) {
    sh2 <- which(grepl("Philips", shNames, ignore.case = T))[1]
    cnvSheet <- mainSheet[,1:15]
    runID <- cnvSheet[1,"Sample_Project"]
    philipVals <- as.data.frame(readxl::read_excel(
        inputFi, sheet = sh2, skip = 3, col_types = "text"))
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
    return(unlist(lapply(samNumber, function(samName) {
        return(which(samName == rawData$`DNA #`))
    })))
}


CheckControlRows2 <- function(rawData, allBnumber) {
    idx <- which(rawData$`DNA #` %in% allBnumber == F)
    dna <- rawData$`DNA #`[idx]

    message(mkRed("The following are Controls or not in Philips:"))
    message(paste(dna, collapse = "\n"))

    isControl <- rawData$`Accession#` %in% c("NTC", "SC", "NC") |
        grepl("Cont", rawData$`Type & Tissue`, ignore.case = TRUE)

    rawData$`Type & Tissue`[idx] <-
        paste(rawData$original_tissue[idx], "Filler")

    rawData$`Type & Tissue`[isControl] <- "Control"

    message(crayon::bgGreen("Assigned Values:"))
    message(paste(
        dna,
        rawData$`Type & Tissue`[idx],
        sep = " = ",
        collapse = "\n"
    ))
    return(rawData)
}


CheckBnumbers <- function(rawData, allBnumber) {
    if (all(allBnumber %in% rawData$`DNA #`) == F) {
        b_txt <-
            "Not all B-numbers from PhilipsExport are in the 'DNA #' Column!"
        message(crayon::bgRed(b_txt))
        missing_Bnumber <- setdiff(allBnumber, rawData$`DNA #`)
        MsgDF(missing_Bnumber)
    }
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
    rawData <- CheckControlRows2(rawData, allBnumber)
    CheckBnumbers(rawData, allBnumber)
    return(rawData$`Type & Tissue`)
}


FindMissingPairs <- function(extraNormals,
                             philipsExport,
                             sheetPairTumor,
                             sheetPairNorm,
                             mainSheet, rawData) {

    if (extraNormals) {
        message(crayon::bgRed("More NORMALS than Tumors on this run!"))
        unmatchedSamples <- which(!philipsExport$`Tumor DNA/RNA Number` %in% sheetPairTumor)
    } else {
        message(crayon::bgRed("More TUMORS than Normals on this run!"))
        unmatchedSamples <-  which(!philipsExport$`Normal DNA/RNA Number` %in% sheetPairNorm)
    }

    if (length(unmatchedSamples) == 0) {
        unique_values <- rawData$Test_Number[!duplicated(rawData$Test_Number) &
                                                 !duplicated(rawData$Test_Number, fromLast = TRUE)]

        if (length(unique_values) > 0) {
            unmatchedSamples <- which(rawData$Test_Number %in% unique_values)
            missingSams <- mainSheet$Sample_ID[unmatchedSamples]
            dnaMissingPair <- stringr::str_detect(mainSheet$Sample_ID, pattern = missingSams)
        } else {
            start_with_0 <- grepl("^0_", mainSheet$Sample_ID)
            missingSams <- mainSheet$Sample_ID[start_with_0]
            contains_any_b <- function(x, substrings) {
                if (any(sapply(substrings, function(sub) grepl(sub, x)))) {
                    return(x)
                }
            }
            extraSams <- unlist(lapply(missingSams, contains_any_b, substrings = sheetPairTumor))
            extraSams <- paste(extraSams, collapse = "|")
            dnaMissingPair <- stringr::str_detect(mainSheet$Sample_ID, pattern = extraSams)
        }

        return(dnaMissingPair)
    }

    missingSams <-
        paste(if (extraNormals) {
            sheetPairNorm[unmatchedSamples]
        } else {
            sheetPairTumor[unmatchedSamples]
        }, collapse = "|")

    dnaMissingPair <- stringr::str_detect(mainSheet$Sample_ID, pattern = missingSams)

    return(dnaMissingPair)
}


MessageMismatched <- function(extraIndices,
                              extraNormals,
                              mainSheet,
                              dnaMissingPair) {

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


ProcessExtraNormals <- function(normSamPair,
                                missingCases,
                                tumorSamPair,
                                mainSheet,
                                sheetTumors) {

    pairedNormals <- which(!normSamPair %in% missingCases)
    extraNormals <- normSamPair[-pairedNormals]

    message("This normal sample will not be paired with a tumor:\n", extraNormals)

    indexToAdd <- which(normSamPair %in% missingCases)
    indexToAdd <- sort(indexToAdd, decreasing = FALSE)

    # Insert normal as an unpaired tumor at this index
    for (idx in indexToAdd) {
        tumorSamPair <- append(tumorSamPair, normSamPair[idx], after = idx - 1)
    }

    tumorRowAdd <- which(tumorSamPair %in% missingCases)
    indexToAdd <- sort(tumorRowAdd, decreasing = FALSE)

    for (sam in missingCases) {
        new_idx <- which(mainSheet$Sample_ID == sam)
        sheetTumors <- append(sheetTumors, new_idx)
    }
    sheetTumors <- sort(sheetTumors)

    for (xSam in extraNormals) {
        toDrop <- normSamPair == xSam
        normSamPair[toDrop] <- "NORMAL_ONLY_NO_TUMOR_PAIR"
    }

    mainSheet$Paired_Normal[sheetTumors] <- normSamPair
    mainSheet$Sample_ID[sheetTumors] <- mainSheet$Sample_ID[sheetTumors] <- tumorSamPair
    return(mainSheet)
}


# Insert "missing" at this index to manually pair -----------------------------
AppendMissingNotes <- function(indexToAdd,
                               normSamPair,
                               mainSheet,
                               sheetTumors) {
    for (idx in indexToAdd) {
        normSamPair <-
            append(normSamPair, "MISSING_NORMAL_PAIR", after = idx - 1)
    }
    if (length(sheetTumors) != length(normSamPair)) {
        message("length(sheetTumors) != length(normSamPair)")
        MsgDF(mainSheet$Paired_Normal[sheetTumors])
    } else{
        mainSheet$Paired_Normal[sheetTumors] <- normSamPair
        message("Missing Indexes Added: ", paste(indexToAdd, collapse = " "))
        MsgDF(normSamPair)
    }

    return(mainSheet)
}


check_if_unpaired <- function(filler_cases, rawData, indexToAdd, not_philips) {
    for (idx in 1:nrow(filler_cases)) {
        curr_test <- filler_cases[idx, "Test_Number"]
        curr_bnum <- filler_cases[idx, "DNA #"]
        is_control <- stringr::str_detect(rawData$`Type & Tissue`, pattern = "(?i)cont")
        non_controls <- rawData[!is_control, ]
        match_sam <- which(non_controls$Test_Number == curr_test)
        is_paired <- length(match_sam) > 1
        if (is_paired) {
            msg1 <- paste(curr_bnum, "is not on Philips, but seems to be paired:")
            message(crayon::bgGreen(msg1))
            MsgDF(non_controls[match_sam, ])
        } else{
            msg2 <- paste(curr_bnum, "is not on Philips, and is not paired.")
            message(crayon::bgGreen(msg2))
            indexToAdd <- c(indexToAdd, not_philips[idx])
            message("This case must be paired manually")
        }
    }
    return(indexToAdd)
}


process_extra_tumors <- function(tumorSamPair,
                                 missingCases,
                                 philipsExport,
                                 sheetPairTumor,
                                 rawData,
                                 sheetTumors,
                                 normSamPair,
                                 mainSheet) {
    indexToAdd <- which(tumorSamPair %in% missingCases)
    all_tumors_philips <- philipsExport$`Tumor DNA/RNA Number`
    if (length(sheetPairTumor) > length(all_tumors_philips)) {
        message("More tumors in XLSheet: ", length(all_tumors_philips))
        message("Less tumors in Philips: ", length(sheetPairTumor))
        not_philips <- which(sheetPairTumor %in% all_tumors_philips == F)
        raw_tumors <- rawData[sheetTumors, ]
        filler_cases <- raw_tumors[not_philips,]
        indexToAdd <- check_if_unpaired(filler_cases, rawData, indexToAdd, not_philips)
    }
    if (length(sheetPairTumor) < length(all_tumors_philips)) {
        message("More tumors in Philips: ", length(all_tumors_philips))
        message("Less tumors in XLSheet: ", length(sheetPairTumor))
    }
    mainSheet <- AppendMissingNotes(indexToAdd, normSamPair, mainSheet, sheetTumors)
    return(mainSheet)
}


AddUnmatchNormals <- function(sheetPairNorm, sheetPairTumor, mainSheet,
                              sheetTumors, sheetNormals, philipsExport,
                              rawData) {
    extraNormals <- length(sheetPairNorm) > length(sheetPairTumor)
    dnaMissingPair <- FindMissingPairs(extraNormals, philipsExport, sheetPairTumor,
                                       sheetPairNorm, mainSheet, rawData)

    if (any(dnaMissingPair)) {
        extraIndices <- which(dnaMissingPair)
        MessageMismatched(extraIndices, extraNormals, mainSheet, dnaMissingPair)

        normSamPair <- mainSheet$Sample_ID[sheetNormals]
        tumorSamPair <- mainSheet$Sample_ID[sheetTumors]
        missingCases <- mainSheet$Sample_ID[dnaMissingPair]

        if (extraNormals == T) {
            mainSheet <- ProcessExtraNormals(normSamPair, missingCases, tumorSamPair,
                                             mainSheet, sheetTumors)
        } else{
            mainSheet <- process_extra_tumors(
                tumorSamPair,
                missingCases,
                philipsExport,
                sheetPairTumor,
                rawData,
                sheetTumors,
                normSamPair,
                mainSheet
            )
        }
    }
    return(mainSheet)
}


check_old_runs <- function(mainSheet, philipsExport) {
    has_missing <- which(mainSheet$Paired_Normal == "MISSING_NORMAL_PAIR")
    if (length(has_missing) > 0) {
        needed_sams <- mainSheet$Specimen_ID[has_missing]
        on_philips <- philipsExport$`Tumor Specimen ID` == needed_sams
        if (any(on_philips)) {
            norm_msg <- "Missing normal(s) may be paired with previous run(s):"
            message(crayon::bgBlue(norm_msg))
            old_ids <- paste(philipsExport$`Epic Order Number`[on_philips],
                             philipsExport$`Normal Run Number`[on_philips],
                             philipsExport$`Normal Specimen ID`[on_philips],
                             philipsExport$`Normal DNA/RNA Number`[on_philips],
                             sep = "_")
            message(paste(old_ids, collapse = "\n"))
        }
    }
}


AddSampleIndexes <- function(pairedList, rawData, philipsExport) {
    tissueType <- FixPairedList(philipsExport, rawData)
    sheetTumors <-
        which(stringr::str_detect(tissueType, pattern = "(?i)tumor"))
    sheetNormals <-
        which(stringr::str_detect(tissueType, pattern = "(?i)norm"))
    sheetControl <-
        which(stringr::str_detect(tissueType, pattern = "(?i)cont"))

    mainSheet <- data.frame(matrix("", nrow = nrow(pairedList), ncol = 0))
    mainSheet$Sample_Name <- mainSheet$Sample_ID <- paste(pairedList[, 1])
    mainSheet$Paired_Normal <- ""

    sheetPairTumor <- rawData$`DNA #`[sheetTumors]
    sheetPairNorm <- rawData$`DNA #`[sheetNormals]

    # Check if length is not equal between tumor and normal samples
    if (length(sheetPairTumor) != length(sheetPairNorm)) {
        message("Issue in AddSampleIndexes():")
        message(crayon::bgRed("length(sheetPairTumor) != length(sheetPairNorm)"))
        mainSheet <- AddUnmatchNormals(sheetPairNorm, sheetPairTumor,
                                       mainSheet, sheetTumors, sheetNormals,
                                       philipsExport, rawData)
    } else {
        mainSheet$Paired_Normal[sheetTumors] <-
            mainSheet$Sample_ID[sheetNormals]
    }

    if (length(rawData$I7_Index_ID) != nrow(mainSheet)) {
        message(crayon::bgRed("Check sample sheet I7_Index_ID column!"))
        stop("length(rawData$I7_Index_ID) != nrow(mainSheet)")
    }
    mainSheet$I7_Index_ID <- paste(rawData$I7_Index_ID)
    mainSheet$index <- paste(rawData$index)
    mainSheet$Specimen_ID <- paste(rawData$`Accession#`)
    mainSheet$Tumor_Content <- mainSheet$Test_Number <- mainSheet$EPIC_ID <- "0"
    mainSheet$Description <- mainSheet$Tumor_Type <- ""
    check_old_runs(mainSheet, philipsExport)
    return(mainSheet)
}

FixDuplicateControls <- function(controlRows) {
    message("Fixing duplicated controls:\n",
            paste(controlRows[duplicated(controlRows)], collapse = "\n"))
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


GetControlInfo <- function(rawData) {
    idx <- stringr::str_detect(rawData$`Type & Tissue`, "Cont|cont")
    numbs <- rawData[idx, c("DNA #", "Accession#", "Test_Number")]
    colnames(numbs) <- c("dna", "acc", "tst")
    return(numbs)
}


CheckControlRows <- function(rawData, runID, newRows) {
    numbs <- GetControlInfo(rawData)
    controlRows <- with(numbs, paste(tst, runID, acc, dna, sep = "_"))

    if (any(duplicated(controlRows))) {
        controlRows <- FixDuplicateControls(controlRows)
    }
    if (any(controlRows %in% newRows)) {
        message("The following control rows are already in newRows to add:")
        extraRows <- controlRows %in% newRows
        MsgDF(controlRows[extraRows])
        controlRows <- controlRows[!extraRows]
        message("Only the following will be considered controls:")
        MsgDF(controlRows)
    } else{
        message(crayon::bgBlue("New control names:"))
        MsgDF(controlRows)
    }

    return(controlRows)
}


GetAccessionInfo <- function(rawData) {
    no_control <- stringr::str_detect(rawData$`DNA #`,
                                      "H20|SERACARE|HAPMAP", negate = T)
    numbs <- rawData[no_control, c("DNA #", "Accession#", "Test_Number")]
    colnames(numbs) <- c("dna", "acc", "tst")
    return(numbs)
}


CheckMissingRows <- function(pairedList, rawData, runID) {
    numbs <- GetAccessionInfo(rawData)
    missRows <- sapply(numbs$dna, function(x) {!any(grepl(x, pairedList))})
    newRows <- NULL
    if (any(missRows)) {
        message(crayon::bgGreen("Binding additional rows/filler:"))
        MsgDF(data.frame(accession_numbers = numbs$acc[missRows]))
        xtraRow <- numbs$dna %in% numbs$dna[missRows]
        newRows <- paste(0, runID, numbs$acc[xtraRow], numbs$dna[xtraRow], sep = "_")
    } else {
        message(crayon::bgGreen("No additional fillers or sample rows to bind to output."))
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

    used_indices <- integer()
    unique_b <- unique(b_numbers)
    b_seen <- setNames(integer(length(unique_b)), unique_b)

    for (sam in seq_along(b_numbers)) {
        current_b <- b_numbers[sam]
        sam_idx <- which(
            str_detect(all_samples, fixed(current_b)) &
                !(seq_along(all_samples) %in% used_indices)
        )
        if (length(sam_idx) == 0) {
            next
        }
        for (idx in sam_idx) {
            b_seen[current_b] <- b_seen[current_b] + 1
            row_idx <- which(b_numbers == current_b)[b_seen[current_b]]
            if (length(row_idx) > 0) {
                new_paired_list[row_idx, ] <- all_samples[idx]
                used_indices <- c(used_indices, idx)
            }
            if (b_seen[current_b] >= sum(b_numbers == current_b)) {
                break
            }
        }
    }

    rownames(new_paired_list) <- NULL
    return(new_paired_list)
}


GetPairedList <- function(philipsExport, runID) {
    pairedList <- paste(unlist(lapply(
        X = 1:length(philipsExport$`Tumor Specimen ID`), function(acc) {
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
        message(crayon::bgRed(
            "The sample sheet may contain duplicate accession Numbers!"))
        dupesIdx <- matchIdx[duplicated(matchIdx)]
        message("Duplicated in row(s): ", paste(dupesIdx, collapse = " "))
        df_x <- setNames(as.data.frame(nameList[dupesIdx]), deparse(substitute(nameList)))
        dupes <- secondCol %in% nameList[dupesIdx]
        df_y <- setNames(as.data.frame(secondCol[dupes]), deparse(substitute(secondCol)))
        MsgDF(df_x)
        MsgDF(df_y)
    }
    return(matchIdx)
}


MsgTumNormErr <- function(tums, norms, sheetSide = "Philips") {
    message("Total Tumors: ", length(tums))
    message("Total Normals: ", length(norms))
    if (sheetSide == "Philips") {
        message(crayon::bgRed(
            "Not all Philips samples are paired with a Tumor/Normal!"))
    } else{
        message(crayon::bgRed(
            "Wetlab Tumor/Normal samples are not paired correctly!"))
    }
}


CheckTotalIndexes <- function(tumorSam, normalSam, ngsNumbers, sheetType = "Philips") {
    if (length(tumorSam) != length(normalSam)) {
        message("Total Tumors: ", length(tumorSam))
        message("Total Normals: ", length(normalSam))
        tumorLonger <- length(tumorSam) > length(normalSam)
        totalAdd <- abs(length(tumorSam) - length(normalSam))
        message(crayon::bgRed(
            paste("Not all", sheetType,
                  "samples are paired with a tumor or normal!")
        ))
        if (tumorLonger == F) {
            totalDrop <- normalSam %in% tumorSam
            message(crayon::bgRed(paste(
                totalAdd, "more Normals than Tumors listed in worksheet!"
            )))
            MsgDF(data.frame(NGS = ngsNumbers[!totalDrop], Potential.Extra = philipsN[!totalDrop]))
            tumorSam <- c(tumorSam, rep(0, totalAdd))
        } else{
            totalDrop <- tumorSam %in% normalSam
            message(crayon::bgRed(paste(
                totalAdd, "more Tumors than Normals listed in worksheet!"
            )))
            MsgDF(data.frame(
                NGS = ngsNumbers[!totalDrop], Potential.Extra = philipsT[!totalDrop]))
            normalSam <- c(normalSam, rep(0,totalAdd))
        }
        message("Extra ", length(which(totalDrop == F)),  " cases...")
    }
}


MsgTotalAdded <- function(totalAdd, ngsNumbers, totalDrop, samExtra, tumorLonger) {
    if (tumorLonger) {
        sheetMessage <- "more Tumors than Normals listed in worksheet!"
    } else{
        sheetMessage <- "more Normals than Tumors listed in worksheet!"
    }
    message(crayon::bgRed(paste(totalAdd, sheetMessage)))

    if (all(totalDrop == F)) {
        ngsNeedsNormal <- ngsNumbers[which(samExtra == 0)]
        return(MsgDF(data.frame(
            NGS = ngsNeedsNormal, Missing_Normal = rep(1, length(ngsNeedsNormal))
        )))
    }

    if (tumorLonger) {
        MsgDF(data.frame(NGS = ngsNumbers[!totalDrop], Potential.Extra = samExtra[!totalDrop]))
        message("Extra ", length(which(totalDrop == F)), " case(s)...")
    } else{
        MsgDF(data.frame(NGS = ngsNumbers[totalDrop], Potential.Extra = samExtra[totalDrop]))
        message("Extra ", length(which(totalDrop == T)), " case(s)...")
    }
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
        } else{
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
            MsgTotalAdded(totalAdd, ngsNumbers, totalDrop, wetLabIdxN, tumorLonger)
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
    mainSheet$GenomeFolder <-
        as.character("PhiX\\Illumina\\RTA\\Sequence\\WholeGenomeFASTA")
    dupes <- base::anyDuplicated(mainSheet$I7_Index_ID)
    if (length(dupes) > 0 & dupes != 0) {
        message(crayon::bgRed(
            "The following rows are duplicated and will be removed:"), "\n"
            )
        MsgDF(mainSheet[dupes,])
        mainSheet <- mainSheet[-dupes,]
        row.names(mainSheet) <- 1:nrow(mainSheet)
    }
    mainSheet$Tumor_Content[mainSheet$Paired_Normal == ""] <- 0
    return(mainSheet)
}


GetExtraNormalRows <- function(rawData){
    non_paired_sam <- !duplicated(rawData$Test_Number) &
        !duplicated(rawData$Test_Number, fromLast = TRUE)
    unique_values <- rawData$Test_Number[non_paired_sam]
    mkRed("The following samples are unpaired:")
    message(paste(unique_values, collapse = "\n"))
    unmatched_idx <- which(rawData$Test_Number %in% unique_values)
    return(unmatched_idx)
}


BuildNoPhilips <- function(rawData, runID, pact_run) {
    message(
        "Generating SampleSheet without Philips Data: ALL are validation!"
        )
    seqId <- stringr::str_split_fixed(runID, "_", 4)
    whichNormal <- stringr::str_detect(rawData$`Type & Tissue`, "Norm|norm|cont|Cont")
    onlyNormals <- stringr::str_detect(rawData$`Type & Tissue`, "Norm|norm")
    concat_id <-
        data.frame(sid = paste(0, runID, rawData$`Accession#`, rawData$`DNA #`, sep = "_"))
    pairedNorm <- rep("", nrow(concat_id))
    whichTumor <- stringr::str_detect(rawData$`Type & Tissue`, "Tum|tum")

    if (sum(whichTumor) != sum(onlyNormals)) {
        message("Not all tumors and normals are paired")
        if (sum(whichTumor) > sum(onlyNormals)) {
            true_indices <- which(whichTumor)[1:sum(onlyNormals)]
            pairedNorm[true_indices] <- concat_id[onlyNormals, 1]
        }else {
            message("Run has multiple normal-only samples")
            true_indices <- which(whichTumor)
            extra_idx <- GetExtraNormalRows(rawData)

            pairsToFind <- rawData$`Accession#`[extra_idx]
            newPairs <- concat_id[onlyNormals, 1]

            unmatched <- unlist(lapply(pairsToFind, function(sam) {
                which(grepl(sam , newPairs))
            }))
            newPairs[unmatched] <- ""  # Un-pair Normal-only samples
            true_indices <- sort(c(true_indices, extra_idx))
            pairedNorm[true_indices] <- newPairs
        }
    } else{
        pairedNorm[whichTumor] <- concat_id[onlyNormals, 1]
    }

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
    mainSheet[!whichNormal, "Tumor_Type"] <- "Not Available"
    return(mainSheet)
}



fix_filler_zeros <- function(rawData, mainSheet) {
    is_filler <- stringr::str_detect(rawData$`Type & Tissue`,
                                     pattern = "(?i)filler")
    if (any(is_filler)) {
        message("Run has filler cases to be output with 0_ below:")
        MsgDF(rawData[is_filler, 1:5])
        to_be_zero <- rawData$Test_Number[is_filler]
        toReplace0 <- which(mainSheet$Test_Number %in% to_be_zero)
        for (idx in toReplace0) {
            curr_vals <- mainSheet[idx, c("Sample_ID", "Sample_Name")]
            new_end <- stringr::str_split_fixed(curr_vals, "_", 2)[, 2]
            new_zeros <- paste0("0_", new_end)
            mainSheet[idx, c("Sample_ID", "Sample_Name")] <- new_zeros
            normals2fix <- grepl(curr_vals[[1]], mainSheet[, "Paired_Normal"] )
            if (any(normals2fix)) {
                curr_norm <- mainSheet[normals2fix, "Paired_Normal"]
                new_end <- stringr::str_split_fixed(curr_norm, "_", 2)[, 2]
                new_norm <- paste0("0_", new_end)
                mainSheet[normals2fix, "Paired_Normal"] <- new_norm
            }
        }
    }
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
    mainSheet <- fix_filler_zeros(rawData, mainSheet)
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

    if ( is.na(sh2)) {
        warning("No PhilipsExport tab data found!")
        return(NULL)
    }

    philipsExport <- GetExcelData(inputFi, sh2, NULL, 3, cm = T)

    if (nrow(philipsExport) == 0) {
        warning("No PhilipsExport tab data found!")
        return(NULL)
    }

    filterColumns <- GetPhilipsColumns()

    if (!all(filterColumns %in% colnames(philipsExport))) {
        stop("PhilipsExport tab is missing headers!")
    }
    if ("Normal Run Number" %in% colnames(philipsExport)) {
        filterColumns <- c(filterColumns, "Normal Run Number")
    }
    philipsExport <- philipsExport[philipsExport$`Test Name` != "", ]
    philipsExport <- philipsExport[, filterColumns]
    blankOrder <- philipsExport[,"Epic Order Number"] == ""

    if (any(blankOrder)) {
        warning("Some Philips Samples are missing MRNs, defaulting to 0")
        philipsExport[, "Epic Order Number"][blankOrder] <- 0
    }

    blankSpecimen <- philipsExport[, "Normal Specimen ID"] == ""

    if (any(blankSpecimen)) {
        missingSam <- philipsExport$`Tumor Specimen ID`[blankSpecimen]
        message(crayon::bgRed(
            "The following Philips Samples are missing a Normal Specimen ID:\n"
        ))
        MsgDF(missingSam)
    }

    philipsExport[philipsExport == ""] <- 0
    return(philipsExport)
}


Err_runID <- function(runID) {
    err_msg <- crayon::bgRed(
        'Keyword "Run ID:" not found in SampleSheet "DNA #" column'
        )
    new_id <- paste('\nDefaulting to script input RUNID:', runID)
    message(err_msg, new_id)
}


GrabSheetBottom <- function(inputFi) {
    shNames <- readxl::excel_sheets(inputFi)
    sh <- which(grepl("PACT-", shNames, ignore.case = TRUE))[1]
    sh <- ifelse(is.na(sh), 1, sh)

    rawData <- GetExcelData(inputFi, sh, shRange = "A6:X200", cm = T)
    if (!"Test_Number" %in% colnames(rawData)) {
        stop("Input Sheet formatting is incorrect, column names are:\n",
             paste(colnames(rawData), collapse = "\n"))
    }
    lastRow <- which(rawData$Test_Number == 0)[1]
    if (is.na(lastRow)) {
        stop(mkRed("Cannot find ending of samplesheet!"))
    }
    sheetBottom <- rawData[c(lastRow:nrow(rawData)), ]
    return(sheetBottom)
}


GrabRunNumber <- function(inputFi, runID) {
    key_run_ids <- "NB551709|NB501073|VH01471"
    sheetBottom <- GrabSheetBottom(inputFi)
    runID_row <- which(stringr::str_detect(sheetBottom$`DNA #`, "Run ID:"))

    # If 'Run ID:' keyword not found, use provided runID as default
    if (length(runID_row) != 1) {
        message(mkRed("Run ID: not found in column 'DNA #, trying column 8..."))
        runID_row <- which(stringr::str_detect(sheetBottom[, 8], "Run ID:"))
        if (length(runID_row) != 1) {
            Err_runID(runID)
            return(runID)
        } else{
            message("runID_row found in column 8! row: ", runID_row)
        }
    }

    runID_find <- paste(sheetBottom[runID_row, ])
    runID_col <- which(stringr::str_detect(runID_find, key_run_ids))

    if (length(runID_col) == 0) {
        Err_runID(runID)
        return(runID)
    }

    run_number <- paste(sheetBottom[runID_row, runID_col])
    run_number <-
        stringr::str_replace_all(run_number, c("Run ID:" = "", " " = ""))

    return(run_number)
}


GetRawSamplesheet <- function(inputFi) {
    shNames <- readxl::excel_sheets(inputFi)
    sh <- which(grepl("PACT-", shNames, ignore.case = T))[1]
    if (is.na(sh)) {
        miss_msg <-
            '"PACT" not detected in Excel sheetnames, defaulting to sheet 1'
        message(crayon::bgRed(miss_msg))
        sh = 1
    }
    msgRd <-
        paste0('Reading Excel Sheet named \"', shNames[sh], '\" from file:')
    message(crayon::bgGreen(msgRd),'\n',inputFi)
    rawData <- GetExcelData(inputFi, sh, shRange = "A6:X200",  cm = T)
    toDrop <- which(rawData[, 15] == "")[1] - 1
    rawData <- rawData[1:toDrop, ]
    rawData$`DNA #` <- stringr::str_replace_all(rawData$`DNA #`, "_", "-")
    return(rawData)
}


WriteMainSheet <- function(mainSheet, sheetHead) {
    csv_file <- paste(mainSheet[1, "Run_Number"], "SampleSheet.csv", sep = "-")
    out_path <- file.path(fs::path_home(), "Desktop", csv_file)

    message("Writing output, check your samplesheet here:\n", out_path)

    if (any(mainSheet$Tumor_Type == "NA")) {
        mainSheet[mainSheet$Tumor_Type == "NA", "Tumor_Type"] <- ""
    }

    write.table(
        sheetHead,
        sep = ",",
        file = out_path,
        row.names = F,
        col.names = F,
        quote = F
    )

    suppressWarnings(
        write.table(
            mainSheet,
            sep = ",",
            file = out_path,
            row.names = F,
            col.names = T,
            append = T,
            quote = F
        )
    )
    return(out_path)
}


# Checks if the user input RunID matches worksheet RunID ----------------------
CheckRunIDMatch <- function(sheetRunID, runID) {
    if (sheetRunID != runID) {
        message(mkRed("SampleSheet and Input RUNID do not match!"))
        message(sprintf("Samplesheet RunID: %s\nUser input RunID:  %s",
                        crayon::bgBlue(sheetRunID), crayon::bgGreen(runID)))
        message("Using runID from inside worksheet: ", sheetRunID)
    }
}


# Parses File when no samplesheet.csv tab is detected -------------------------
AltParseFormat <- function(worksheetPath, runID) {
    rawData <- GetRawSamplesheet(worksheetPath)
    sheetRunID <- GrabRunNumber(worksheetPath, runID)
    CheckRunIDMatch(sheetRunID, runID)
    philipsExport <- GetPhilipsData(worksheetPath)
    input_base <- base::basename(worksheetPath)
    pact_run <- stringr::str_split_fixed(input_base, ".xls", 2)[1,1]
    message(crayon::bgBlue(paste("PACT ID is:", pact_run)))

    if (is.null(philipsExport) | nrow(philipsExport) == 0) {
        mainSheet <- BuildNoPhilips(rawData, sheetRunID, pact_run)
    } else{
        mainSheet <-
            BuildMainSheet(philipsExport, rawData, sheetRunID, pact_run)
    }
    return(mainSheet)
}


PostRedcapCurl <- function(rcon, datarecord, retcon = 'ids') {
    message(crayon::bgMagenta("REDCap Output:"))
    tryCatch(
        expr = {RCurl::postForm(
            rcon$url, token = rcon$token, content = 'record', format = 'json',
            type = 'flat', data = datarecord,
            returnContent = retcon, returnFormat = 'csv')
        },
        error = function(e) {message("REDCap push failed!\n", e)},
        finally = message(datarecord)
    )
}


# Generate Email notification and attach csv file -----------------------------
emailNotify <- function(record, rcon) {
    record$pact_csv_email <- "pact_csv_email"
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox = T)
    PostRedcapCurl(rcon, datarecord)
    brick_message("Email Notification Created")
}

# Calls API for CSV file ------------------------------------------------------
validateInputs <- function(rcon, recordName, fiPath) {
    if (is.null(rcon) || !is.list(rcon)) {
        stop("rcon should be a list")
    }
    if (!is.character(recordName) || length(recordName) != 1) {
        stop("recordName should be a single character string.")
    }
    if (!file.exists(fiPath)) {
        stop("The file specified as fiPath does not exist.")
    }
}

# Uploads CSV file to REDCap --------------------------------------------------
callApiFileCsv <- function(rcon, recordName, fiPath) {
    validateInputs(rcon, recordName, fiPath)
    
    upload_success <- try(redcapAPI::importFileToRecord(
        rcon,
        file = fiPath,
        record = recordName,
        event = NULL,
        field = "pact_csv_sheet"
    ), silent = T)
    
    if (upload_success == T) {
        message("File uploaded to REDCap successfully:\n", fiPath)
    } else {
        message("File failed to upload to REDCap in callApiFileCsv:\n", fiPath)
    }
}


# Connect to REDCap and send email attachments of csv file --------------------
pushToRedcap <- function(outVals, token = NULL) {

    stopifnot(!is.null(token) & !is.null(outVals))
    runID <- outVals[[1]]
    outFile <- outVals[[2]]
    api_url <- "https://redcap.nyumc.org/apps/redcap/api/"

    rcon <- redcapAPI::redcapConnection(api_url, token)
    record = data.frame(record_id = runID, pact_run_number = runID)
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox = T)

    PostRedcapCurl(rcon, datarecord, retcon = 'nothing')
    callApiFileCsv(rcon, runID, outFile)

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


replace_strings <- function(sheet_df, old_string, new_string) {
    sheet_df %>%
        dplyr::mutate_all(~str_replace_all(., fixed(old_string), new_string))
}


check_keyword <- function(keywd) {
    ilm_msg <- "ILMN samples were detected on this run, but no number provided!"
    pact_sh <- "/Volumes/CBioinformatics/PACT/parsepact.sh"
    ilm_num <- "ILMNVAL-#"
    cmd_msg <- paste(pact_sh, input, runID, ilm_num)
    if (stringr::str_detect(keywd, "-") == F || nchar(keywd) != 9) {
        if (nchar(keywd) != 9) {message("Keyword is not 9 characters: ", keywd)}
        message(mkRed(ilm_msg), " ", keywd)
        message("For example, execute with third arg:\n", cmd_msg)
        stop(mkRed("Third argument needed: ILMNVAL-#"))
    }
}


MakeValidationSheet <- function(sheetHead, mainSheet, has_validation) {
    if (!is.null(keywd)) {
        if (substr(keywd, 1, 1) == "-") {
            keywd <- substring(keywd, 2)
        }
    }
    check_keyword(keywd)
    mainSheet_val <- mainSheet[has_validation, ]
    pact_id <- mainSheet_val$Sample_Project[1]
    run_num <- mainSheet_val$Run_Number[1]
    new_pact_id <- paste0(pact_id, "-" , keywd)
    new_run_num <- paste0(run_num, "-", keywd)
    mainSheet_val <- replace_strings(mainSheet_val, pact_id, new_pact_id)
    mainSheet_val <- replace_strings(mainSheet_val, run_num, new_run_num)
    sheetHead_val <- sheetHead
    sheetHead_val[4, 2] <- new_pact_id
    sheetHead_val[5, 2] <- new_pact_id
    outFile <- WriteMainSheet(mainSheet_val, sheetHead_val)
    outValsVal <- c(runID = mainSheet_val[1, "Sample_Project"],
                    outFile = outFile)
    pushToRedcap(outValsVal, token)
    mainSheet <- mainSheet[!has_validation, ]
    row.names(mainSheet) <- 1:nrow(mainSheet)
    return(mainSheet)
}


# Parses xlsx file and writes as csv file -------------------------------------
parseExcelFile <- function(worksheetPath, runID = NULL) {
    shNames <- readxl::excel_sheets(worksheetPath)
    MsgDF(data.frame(`Sheet names in Workbook` = shNames))
    sheetHead <- GetSheetHeading(worksheetPath)
    mainSheet <- AltParseFormat(worksheetPath, runID)
    mainSheet <- sanitizeSheet(mainSheet)
    try(WritePhilipsGender(mainSheet,worksheetPath, shNames), silent = T)
    has_validation <- stringr::str_detect(mainSheet$Sample_Name, "-ILC-VAL_")

    if (any(has_validation)) {
        mainSheet <- MakeValidationSheet(sheetHead, mainSheet, has_validation)
    }

    outFile <- WriteMainSheet(mainSheet, sheetHead)
    return(c(runID = mainSheet[1, "Sample_Project"], outFile = outFile))
}


# Filters list of possible files in the directory for worksheet ---------------
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
    } else{
        stop(paste("No PACT worksheet was found in the directory:", parentFolder))
    }
    if (file.exists(potentialFi[1])) {
        message(crayon::bgGreen("Now trying to read:"), "\n", potentialFi[1], "\n")
        pfile <- potentialFi[1]
    } else{
        nextMsg <- paste(notFoundMsg, potentialFi[1], "Trying:", potentialFi[2], sep = "\n")
        message(nextMsg)
        pfile <- potentialFi[2]
    }
    # TODO: Write check for mismatched year
    # runyr <- stringr::str_split_fixed(inputFi, "-", 3)[, 2]
    stopifnot(file.exists(pfile))
    outVals <- parseExcelFile(pfile, runID)
    return(outVals)
}


determineRunType <- function(input) {
    if (stringr::str_detect(input, "2000|NextSeq550")) {
        runType <- "NextSeq2000"
    } else if (stringr::str_detect(input, "TEST")) {
        runType <- "test"
    } else if (stringr::str_detect(input, "ILC")) {
        runType <- "Illumina"
    } else if (stringr::str_detect(input, "TMB")) {
        runType <- "TMB"
    } else {
        runType <- "regular"
    }
    return(runType)
}


# FUN: Returns Path to xlsx file ----------------------------------------------
getExcelPath <- function(input, runType) {
    if (stringr::str_detect(input, .Platform$file.sep)) {
        return(input)
    }
    drive = file.path("", "Volumes", "molecular", "MOLECULAR LAB ONLY")
    folder = file.path("NYU PACT Patient Data", "Workbook")
    runyr <- stringr::str_split_fixed(input, "-", 3)[, 2]
    yearDir <- paste0("20", runyr)
    xlFi <- paste0(input, ".xlsm")
    worksheetPath <- file.path(drive, folder, yearDir, input, xlFi)

    if (runType == "NextSeq2000") {
        folder <- "Validations/PACT new i7-NextSeq2000/Wet Lab/Workbook"
        worksheetPath <- file.path(drive, folder, yearDir, xlFi)
    }
    if (runType == "test") {
        drive = "/Volumes/molecular/Molecular/Validation/PACT/Test_Sheets"
        worksheetPath <- file.path(drive, xlFi)
    }
    if (runType == "Illumina") {
        drive = "/Volumes/molecular/Molecular/Validation/PACT/Test_Sheets"
        worksheetPath <- file.path(drive, xlFi)
    }
    if (runType == "TMB") {
        drive = "/Volumes/molecular/MOLECULAR LAB ONLY/Validations/TMB"
        worksheetPath <- file.path(drive, xlFi)
    }
    return(worksheetPath)
}


# Gets dataframe and saves as CSV file ----------------------------------------
writeSampleSheet <- function(input, token, runID = NULL) {
    outVals <- NULL
    runType <- determineRunType(input)
    worksheetPath <- getExcelPath(input, runType)
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

writeSampleSheet(input, token, runID)
