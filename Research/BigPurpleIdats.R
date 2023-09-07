#!/usr/bin/env R
## ---------------------------
## Script name: idats_to_BigPurple.R
## Purpose: source of global scripts, copy idat files to HPC, & create minfi samplesheet.csv
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------

# Set package installation options
library("base"); gb <- globalenv(); assign("gb", gb)
options(repos = c(CRAN = "http://cran.us.r-project.org"))

# Input Arguments
args <- commandArgs(TRUE)
try(args[1] ->> token)
try(args[2] ->> inputSheet)
try(args[3] ->> copyToFolder)

# Helper function to return shared drive names
NameDrive <- function(volumeName, driveSpace){
    return(file.path("","Volumes", volumeName, driveSpace))
}

# Default Paths
apiLink <<- "https://redcap.nyumc.org/apps/redcap/api/"
rsch.idat <<- NameDrive("snudem01labspace","idats")
clin.idat <<- NameDrive("molecular", "MOLECULAR/iScan")

# REDcap Heading Fields to pull for SampleSheet
flds <<- c("record_id", "b_number", "tm_number", "accession_number", "block", "diagnosis", "organ", "tissue_comments", "run_number", "nyu_mrn")

# Helper function to suppress messages and warnings
supM <- function(sobj) suppressMessages(suppressWarnings(sobj))

loadOrInstallPackages <- function(pkg){
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, dependencies = TRUE, verbose = TRUE, ask = FALSE)
    }
    library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
}

loadOrInstallPackages("devtools")

# Package loading function
loadPacks <- function() {
    pkgs <- c("data.table", "foreach", "openxlsx", "jsonlite", "RCurl", "readxl", "stringr", "tidyverse", "crayon", "redcapAPI", "remotes")
    lapply(pkgs, loadOrInstallPackages)
}

WarnMounts <- function(idat.dir) {
    if (!dir.exists(idat.dir)) {
        stop(crayon::bgRed("Directory not found, ensure the idat path is accessible:"), "\n", idat.dir)
    }
}

checkMounts <- function() {
    molecDrive = "/Volumes/molecular/MOLECULAR LAB ONLY"
    zDrive = "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular"

    if (dir.exists(molecDrive)) {
        message("\n", crayon::bgGreen("Z-drive path is accessible"), "\n")
    } else {
        cat("\nPATH does not exist, ensure path is mounted:\n")
        cat(crayon::white$bgRed$bold(molecDrive))
        cat("\nYou must mount the network Z-drive path:\n")
        cat(crayon::white$bgRed$bold(zDrive), "\n")
        stop("Z-drive path is inaccessible.")
    }
}

grabAllRecords <- function(flds, rcon) {
    dbCols <- redcapAPI::exportRecordsTyped(
        rcon,
        fields = flds,
        labels = FALSE,
        dates = FALSE,
        survey = FALSE,
        dag = FALSE,
        factors = FALSE,
        form_complete_auto = FALSE
    )
    return(as.data.frame(dbCols))
}

search.redcap <- function(rd_numbers, token=NULL, flds=NULL) {
    if (is.null(token)){stop("You must provide an ApiToken!")}
    rcon <- redcapAPI::redcapConnection(gb$apiLink, token)
    flds <-
        if (is.null(flds)) {
            c(
                "record_id",
                "b_number",
                "primary_tech",
                "second_tech",
                "run_number",
                "barcode_and_row_column",
                "accession_number",
                "arrived"
            )
        } else {
            flds
        }
    result <-
        redcapAPI::exportRecords(
            rcon,
            records = rd_numbers,
            fields = flds,
            dag = F,
            factors = F,
            labels = F,
            dates = F,
            form_complete_auto = F,
            format = 'csv'
        )
    return(as.data.frame(result))
}

getDefaults <- function() {
    cbVol <- switch(Sys.info()[['sysname']], "Darwin"="/Volumes/CBioinformatics/Methylation", "Linux"="~/molecpathlab/production/Methylation")
    moVol = "/Volumes/molecular"
    rsVol = "/Volumes/snudem01labspace"
    defaultParams <- data.frame(
        mnp.pk.loc = file.path(cbVol, "classifiers/mnp.v11b6"),
        # other parameters as in the original code
        stringsAsFactors=F
    )
    return(defaultParams)
}

setVar <- function(valueName,val){
    return(assign(valueName, val, envir=.GlobalEnv))
}

assignVar <- function(varStr, assignedVal){
    if (exists(varStr, envir = .GlobalEnv)) {
        message(varStr," = ", assignedVal)
    } else {
        setVar(varStr, assignedVal)
    }
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

defineParams <- function(mnp.pk.loc = NULL,
                         ApiToken = NULL,
                         methDir = NULL,
                         clinDrv = NULL,
                         rschOut = NULL,
                         clinOut = NULL,
                         rsch.idat = NULL,
                         clin.idat = NULL,
                         QC_file = NULL,
                         isMC = T,
                         baseDir = NULL,
                         runID = NULL,
                         loadClassifier = T) {
    defVars <- getDefaults()
    inVars <- list(mnp.pk.loc, ApiToken, methDir, clinDrv, rschOut, clinOut, rsch.idat, clin.idat, QC_file, baseDir)

    non_null_indices <- !sapply(inVars, is.null)
    defVars[non_null_indices] <- inVars[non_null_indices]

    lapply(names(defVars), function(x) {
        if (!is.null(defVars[[x]])) setVar(x, defVars[[x]])
    })

    sapply(names(defVars), function(i) assignVar(i, defVars[[i]]))

    cbVol <- switch(Sys.info()[['sysname']], "Darwin"="/Volumes/CBioinformatics/", "Linux"="~/molecpathlab/production/")

    if (!isMC) {
        methDir = rschOut
        assign("workDir", cbVol)
    }
}


setDirectory <- function(foldr) {
    bsDir = paste("cd", foldr)
    if (dir.exists(foldr)) {
        system(bsDir)
        setwd(foldr)
        assign("workDir", foldr)
    } else{
        warning(paste("Location Not Found:", foldr))
    }
}

SourceFunctions <- function(workingPath = NULL) {
    git_url <- "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/CopyInputs.R"
    workingPath <- if (is.null(workingPath)) getwd() else workingPath
    invisible(devtools::source_url(url = git_url))
    setDirectory(foldr = workingPath)
    return(defineParams(loadClassifier = F))
}

readInfo <- function(inputSheet) {
    rds <- if (endsWith(inputSheet, ".csv")) {
        read.delim(inputSheet, header = F, sep = ",", colClasses = "character")[, 1]
    } else {
        readxl::read_excel(inputSheet, col_names = F)[, 1]
    }
    rds <- rds[!is.na(rds)]
    return(rds)
}


WriteSamSheet200 <- function(samplesheet_csv, outputFi, rowBatch = 200){
    nRows <- nrow(samplesheet_csv)
    if (nRows > rowBatch) {
        num_iterations <- ceiling(nRows / rowBatch)
        for (i in seq_len(num_iterations)) {
            start_row <- (i - 1) * rowBatch + 1
            end_row <- min(i * rowBatch, nRows)
            outputFi_sub <- gsub(".csv$", paste0("_set_", i, ".csv"), outputFi)
            write.csv(samplesheet_csv[start_row:end_row,],
                      outputFi_sub,
                      quote = FALSE,
                      row.names = FALSE)
        }
    }
}


makeSampleSheet <- function(df, samplesheet_ID, bn = NULL, outputFi = "samplesheet_og.csv") {
    bn <- if (is.null(bn)) file.path(getwd(), df$barcode_and_row_column) else bn
    df <- df[!is.na(df$barcode_and_row_column), ]

    samplesheet_csv <- data.frame(
        Sample_Name = df[, "record_id"],
        DNA_Number = df[, "b_number"],
        Sentrix_ID = samplesheet_ID[, 1],
        Sentrix_Position = samplesheet_ID[, 2],
        SentrixID_Pos = df[, "barcode_and_row_column"],
        Basename = paste0(bn),
        RunID = df$run_number,
        MP_num = df$accession_number,
        Date = df$arrived
    )

    toDrop <- stringr::str_detect(samplesheet_csv$SentrixID_Pos, "DUPLICATE")

    if (any(toDrop)) {
        message("Dropping duplicated samples!")
        otherCsv <- samplesheet_csv[toDrop, ]
        write.csv(otherCsv, file = "duplicated_samples.csv", quote = FALSE, row.names = FALSE)
    }

    samplesheet_csv <- samplesheet_csv[!toDrop, ]

    write.csv(samplesheet_csv, outputFi, quote = FALSE, row.names = FALSE)
    WriteSamSheet200(samplesheet_csv, outputFi)
}


Find_copy_idats <- function(rd_numbers, token, copyIdats = T, outputFi = "samplesheet_og.csv", idatPath = NULL) {
    idatPath <- if (is.null(idatPath)) file.path(getwd(), "idats") else idatPath
    result_raw <- search.redcap(rd_numbers, token)
    result <- result_raw[!is.na(result_raw$barcode_and_row_column), ]
    res_barcodes <- result[, "barcode_and_row_column"]
    samplesheet_ID <- as.data.frame(stringr::str_split_fixed(res_barcodes, "_", 2))

    if (nrow(samplesheet_ID) == 0) {
        message("The RD-numbers entered do not have sentrix IDs!",
                "Either samples have not been run or do not have SentrixID in REDCap:",
                paste(capture.output(result_raw), collapse = "\n"))
    }

    stopifnot(nrow(samplesheet_ID) > 0)

    makeSampleSheet(result, samplesheet_ID, bn = NULL, outputFi = outputFi)

    if (copyIdats) {
        GetIdatsFromDrives(csvNam = outputFi, runDir = idatPath)
    }
}


RsyncBigPurple <- function(allFi, idatPath = NULL) {
    idatPath <- if (is.null(idatPath)){gb$copyToFolder}else{idatPath}
    msgFunName(cpInLnk, "RsyncBigPurple")
    cat(sprintf("%s\n", crayon::white$bgCyan("Copying idats to current directory...")))

    userNam <- Sys.info()[["user"]]
    outputCmd <- sprintf("%s@bigpurple.nyumc.org:%s", userNam, copyToFolder)

    logFile <- "rsync_logs.txt"
    if (!file.exists(logFile)) {
        file.create(logFile)
    }

    # Read existing logs to determine files that have already been copied
    existingFiles <- readLines(logFile, warn = FALSE)

    for (iFile in allFi) {
        if (!(iFile %in% existingFiles)) {
            copyCmd <- sprintf("rsync -e ssh %s %s", iFile, outputCmd)
            message(sprintf("Copying file:\n%s", copyCmd))
            systemStatus <- try(system(copyCmd), silent = TRUE)

            # If rsync command is successful, append the filename to the log
            if (inherits(systemStatus, "try-error") == FALSE) {
                write(iFile, file = logFile, append = TRUE)
            }
        } else {
            message(sprintf("File %s is already copied according to logs. Skipping.", iFile))
        }
    }
}

           
LoopIdatFiles <- function(idatsToCopy, tRows=200) {
    totalFiles <- length(idatsToCopy)
    if (totalFiles > tRows) {
        message("More than ", tRows, " idats being copied!")
        message(sprintf("Looping through %d files at a time from total: %d", tRows, totalFiles))
        num_iterations <- ceiling(totalFiles / tRows)
        for (i in 1:num_iterations) {
            current_files <- idatsToCopy[((i - 1) * tRows + 1) : min(i * tRows, totalFiles)]
            red_files <- stringr::str_replace_all(current_files, "_Grn", "_Red")
            current_files <- unique(c(red_files, current_files))
            message(sprintf("Copying %d out of %d sets of %d files", i, num_iterations, tRows))
            RsyncBigPurple(current_files)
        }
    } else {
        RsyncBigPurple(idatsToCopy)
    }
}

           
GetIdatsFromDrives <- function(csvNam = "samplesheet.csv", runDir = NULL) {
    runDir <- if (is.null(runDir)) getwd() else runDir
    extr.idat <- file.path(gb$rsch.idat, "External")

    lapply(c(rsch.idat, clin.idat), WarnMounts)

    stopifnot(file.exists(csvNam))

    ssheet <- read.csv(csvNam, strip.white = TRUE)
    allIdats <- getAllFiles(idatDir = c(rsch.idat, clin.idat), csvNam = csvNam)
    allFi <- allIdats[file.exists(allIdats)]
    if (length(allFi) == 0) {
        allFi <- GetExternalIdats(allFi, ssheet, extr.idat)
    }
    stopifnot(length(allFi) > 0)
    message("Files found: ")
    DataFrameMessage(allFi)

    allFi <- GetExternalIdats(allFi, ssheet, extr.idat)
    cur.idat <- basename(dir(path = runDir, pattern = "*.idat$", recursive = FALSE))

    idatsToCopy <- allFi[!(basename(allFi) %in% cur.idat)]

    if (length(idatsToCopy) > 0) {
        LoopIdatFiles(idatsToCopy)
    } else {
        message(".idat files already copied to run directory")
    }
}

if(Sys.info()[['sysname']] == "Darwin") {checkMounts()}
SourceFunctions()

if (length(inputSheet) > 0 && length(token) > 0 && !is.na(inputSheet) && !is.na(token)) {
    rd_numbers <- readInfo(inputSheet)
    Find_copy_idats(rd_numbers = rd_numbers, token = token, idatPath = copyToFolder)
}
