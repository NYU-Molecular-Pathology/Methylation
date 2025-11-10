#!/usr/bin/env Rscript
## Script name: CopyInputs.R
## Purpose: Functions to copy input files for a methylation run directory
## Date Created: June 13, 2022
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025

gb <- globalenv(); assign("gb", gb)
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
cpInLnk = "https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/CopyInputs.R"
defaultDir = "/Volumes/CBioinformatics/Methylation/Clinical_Runs"

# Helper function to debug: prints the name of the function on GitHub source
msgFunName <- function(pthLnk, funNam) {
    message("Executing function: ", crayon::black$bgYellow(funNam), " from RScript in:\n", pthLnk,"\n")
}

# Helper function to debug prints parameters passed to a function
msgParams <- function(...) {
    message("\nParam passed: ", crayon::bgGreen(paste(deparse(substitute(...)), "=", ...)), "\n")
}

# DEFAULT DIRECTORY PATHS -----------------------------------------------------
cbVol <- "/Volumes/CBioinformatics/Methylation"
moVol <- "/Volumes/molecular"
rsVol <- "/Volumes/snudem01labspace"

mnp.pk.loc = file.path(cbVol, "classifiers/mnp.v12epicv2")
methDir = file.path(cbVol, "Clinical_Runs")
clinDrv = file.path(moVol, "MOLECULAR LAB ONLY/NYU-METHYLATION")
rschOut = file.path(rsVol, "FINAL_PDF_Reports_Brain")
clinOut = file.path(moVol, "MOLECULAR/MethylationClassifier")
rsch.idat = file.path(rsVol, "idats")
clin.idat = file.path(moVol, "MOLECULAR/iScan")

# Changes the working directory using the system CD command
setDirectory <- function(foldr) {
    msgFunName(cpInLnk,"setDirectory")
    bsDir = paste("cd", foldr)
    mm2 = crayon::white$bgRed("Location Not Found:", foldr)
    
    if (dir.exists(foldr)) {
        system(bsDir)
        setwd(foldr)
        assign("workDir", foldr)
    }else{warning(mm2)}
}


# FUN: Creates a new run directory if it does not exist and chmods permission
CreateRunDir <- function(newRun) {
    if (endsWith(newRun, "/")) {
        newRun <- substr(newRun, 1, nchar(newRun) - 1)
    }
    
    message(crayon::bgGreen("New Run Path:"), "\n", newRun)
    
    if (!dir.exists(newRun)) {
        dir.create(newRun, recursive = T)
        Sys.chmod(newRun, "0777", use_umask = FALSE)
    }
    
    setDirectory(newRun)
    return(newRun)
}


# Sets the methylation run directory named by the new run name
setRunDir <- function(runID = NULL, workFolder = NULL) {
    msgFunName(cpInLnk, "setRunDir")
    msgParams(runID)
    msgParams(workFolder)
    
    workFolder <- ifelse(is.null(workFolder), defaultDir, workFolder)
    
    if (is.null(runID)) {
        runID <- paste0(basename(getwd()))
    }
    newRun <- file.path(workFolder, runID)
    if (grepl("TEST", runID)) {
        if (dir.exists(newRun)) {
            unlink(newRun, T, T)
        }
        dir.create(newRun)
        desk_path <- file.path(fs::path_home(), "Desktop", runID)
        try(unlink(desk_path, T, T), silent = T)
    }
    methylPath <- CreateRunDir(newRun)
    return(methylPath)
}


# Lists all idat files in the run directory and checks uniform sizes in Mb
check_idat_sizes <- function(runFolder) {
    msgFunName(cpInLnk, "check_idat_sizes")
    idat_files <- dir(path = runFolder,
                      pattern = "\\.idat$",
                      full.names = TRUE)
    idat_sizes <- round(file.size(idat_files) / 1e6, 1)
    mb_unique <- unique(idat_sizes)
    
    if (length(mb_unique) > 1) {
        mb_sizes <- c(14.4, 13.7, 10.3, 8.1)
        miss_copied <- !idat_sizes %in% mb_sizes
        if (any(miss_copied)) {
            file_names <- paste(basename(idat_files)[miss_copied])
            file_megaB <- paste(" ", idat_sizes[miss_copied], "MB")
            files_sizes <- c(rbind(file_names, file_megaB, "\n"))
            message(
                paste(
                    crayon::bgRed("Not all idat files are the same file size!"),
                    "Ensure idat files copied correctly & are uniform in size:\n",
                    sep = "\n"
                ),
                files_sizes
            )
        }
    }
}


# FUN: Returns a list of idat files given an idat drive location -
getAllFiles <- function(idatDir, csvNam = NULL) {
    msgFunName(cpInLnk, "getAllFiles")
    msgParams(idatDir)
    msgParams(csvNam)
    if (!is.null(csvNam)) {
        ssheet = read.csv(csvNam, strip.white = T)
        barcode = as.vector(ssheet$Sentrix_ID)
    } else {
        ssheet = read.csv("samplesheet.csv", strip.white = T)
        barcode = as.vector(ssheet$Sentrix_ID)
    }
    allFi = NULL
    senP = ssheet$SentrixID_Pos
    for (i in 1:length(idatDir)) {
        dirNames = file.path(idatDir[i], barcode)
        green.files <- file.path(dirNames, paste0(senP, "_Grn.idat"))
        red.files <- file.path(dirNames, paste0(senP, "_Red.idat"))
        allFi <- c(allFi, green.files, red.files)
    }
    return(allFi)
}


check_success_copy <- function(allFi) {
    idcs = basename(allFi)
    idatsCopied <- idcs[idcs != ""]
    success = file.exists(idatsCopied)
    all(success)
    message(".idat files that failed to copy:")
    if (all(success)) {
        cat("none", "\n")
    } else {
        print(idatsCopied[!success])
    }
}


MakeLogFile <- function(infoData, logFile) {
    message(crayon::bgBlue("~~~Message logged~~~"), "\n", infoData)
    message(crayon::bgGreen("To file:"), " ", logFile)
    write.table(infoData, file = logFile, append = T, quote = F,
                sep = '\t', row.names = F, col.names = F)
}


copyBaseIdats <- function(allFi, idatPath = NULL) {
    msgFunName(cpInLnk, "copyBaseIdats")
    suppressPackageStartupMessages(library("cli"))
    
    if (is.null(idatPath)) {
        idatPath <- getwd()
    }
    
    # Check read permission (necessary for copying)
    readable <- fs::file_access(allFi, mode = "read")
    if (any(readable == F)) {
        failed_reads <- allFi[!readable]
        infoData <- paste("Cannot read idat file:", failed_reads, collapse = "\n")
        MakeLogFile(infoData, "read_error_idat.txt")
        allFi <- allFi[readable]
    }
    
    cli::cli_progress_bar("Copying files", total = length(allFi))
    
    for (f in allFi) {
        tryCatch(
            fs::file_copy(f, file.path(getwd(), basename(f)), overwrite = TRUE),
            error = function(e) {
                infoData <- paste("Failed to copy:", f)
                cli::cli_alert_danger(infoData)
                MakeLogFile(infoData, "missing_idat_files.txt")
            }
        )
        cli::cli_progress_update()  # Updates progress
    }
    
    cli::cli_progress_done()  # Mark progress as complete
    check_success_copy(allFi)
}


# Helper function to notify if a network mount is not found
WarnMounts <- function(idat.dir){
    msgFunName(cpInLnk, "WarnMounts")
    if (!dir.exists(idat.dir)) {
        message(
            crayon::bgRed("Share drive not found, ensure path is accessible:"),
            "\n", idat.dir
        )
        stopifnot(dir.exists(idat.dir))
    }
}

# FUN: Helper function returns input as dataframe message
DataFrameMessage <- function(dat){
    return(message(paste0(capture.output(as.data.frame(dat)), collapse = "\n")))
}


# Function to ensure all idat files exist and are not missing -----------------
CheckIdatsReal <- function(ssheet, allFi) {
    basesNeeded = as.vector(ssheet$SentrixID_Pos)
    if (length(unique(basesNeeded)) * 2 != length(unique(basename(allFi)))) {
        message(crayon::bgRed("Still missing idat files not in External folder:"))
        themissed <- stringr::str_split_fixed(basename(allFi), "_", 3)[, 1:2]
        themissed <- paste(themissed[, 1], themissed[, 2], sep = "_")
        message("The following samples are missing:")
        missing_samples <- basesNeeded[!(basesNeeded %in% themissed)]
        DataFrameMessage(missing_samples)
        missing_samples_df <- data.frame(Missing_Samples = missing_samples)
        write.csv(
            missing_samples_df,
            "missing_idats_log.csv",
            row.names = F,
            quote = F
        )
        message(
            crayon::bgRed("Check the log file to see which idats were not found:"),
            " ", "missing_idats_log.csv"
        )
    }
}


VerifyIdatFound <- function(foundIdat, otherIdat, toBeFound, extr.idat){
    if (any(foundIdat) == F) {
        
        message(crayon::bgRed("Still missing idat files not in External folder:"))
        DataFrameMessage(toBeFound)
        return(NULL)
    }
    message(crayon::bgGreen("Found extra idats in External folder:"),
            " ",
            extr.idat)
    idatsToAdd <- otherIdat[foundIdat]
    DataFrameMessage(idatsToAdd)
    return(idatsToAdd)
}


ListMissedIdats <- function(allFi, basesNeeded){
    basesSplit <- stringr::str_split_fixed(basename(allFi), "_", 3)[, c(1, 2)]
    basesFound <- unique(paste0(basesSplit[, 1], "_", basesSplit[, 2]))
    stillMissing <- !(basesNeeded %in% basesFound)
    return(stillMissing)
}


GetIdats2Add <- function(toBeFound, extr.idat){
    message(crayon::bgRed("The following idats are missing:"))
    DataFrameMessage(toBeFound)
    message(crayon::bgGreen("Searching the External folder for more idats..."))
    redGreenFi <- paste0(rep(toBeFound, each = 2), c("_Grn.idat", "_Red.idat"))
    idatsToAdd <- file.path(extr.idat, redGreenFi)
    if (all(file.exists(idatsToAdd)) == F) {
        otherIdat <- dir(extr.idat, pattern = ".idat", full.names = T, recursive = T)
        toBeSearch <- paste(toBeFound, collapse = "|")
        foundIdat <- stringr::str_detect(otherIdat, pattern = toBeSearch)
        idatsToAdd <- VerifyIdatFound(foundIdat, otherIdat, toBeFound, extr.idat)
    }
    return(idatsToAdd)
}


# FUN: Checks for additional idat files if not in iScan or research folders
GetExternalIdats <- function(allFi, ssheet, extr.idat) {
    basesNeeded = as.vector(ssheet$SentrixID_Pos)
    if (length(unique(basesNeeded)) * 2 == length(unique(basename(allFi)))) {
        message("All idats detected in folders!")
        return(allFi)
    }
    message(
        crayon::bgRed("Still missing some idats! Checking External Folder:"),
        " ",
        extr.idat
    )
    
    if (length(allFi) > 0) {
        stillMissing <- ListMissedIdats(allFi, basesNeeded)
    } else{
        stillMissing <- basesNeeded %in% basesNeeded
        allFi <- NULL
    }
    
    if (any(stillMissing) == T) {
        message("Missing idats:\n", DataFrameMessage(ssheet[stillMissing, ]))
        idatsToAdd <- GetIdats2Add(basesNeeded[stillMissing], extr.idat)
        if (!is.null(allFi)) {
            allFi <- unique(c(allFi, setdiff(idatsToAdd, allFi)))
            CheckIdatsReal(ssheet, allFi)
        } else{
            allFi <- idatsToAdd
        }
    }
    return(allFi)
}


# FUN: Returns list of idat files that exist on Molecular and SnuderlLab drives -
get.idats <- function(csvNam = "samplesheet.csv", runDir = NULL) {
    msgFunName(cpInLnk, "get.idats")
    extr.idat <- file.path(rsch.idat, "External")
    WarnMounts(rsch.idat)
    WarnMounts(clin.idat)
    
    if (is.null(runDir)) {
        runDir <- getwd()
    }
    
    if (!file.exists(csvNam)) {
        message("Cannot find your sheet named:", csvNam)
        stopifnot(file.exists(csvNam))
    }
    
    ssheet = read.csv(csvNam, strip.white = T)
    allFi <- getAllFiles(idatDir = c(rsch.idat, clin.idat),
                         csvNam = csvNam)
    
    # Check existence
    exists <- fs::file_exists(allFi)
    allFi <- allFi[exists]
    
    if (length(allFi) == 0) {
        allFi <- GetExternalIdats(allFi, ssheet, extr.idat)
    }
    if (!length(allFi) > 0) {
        warning(crayon::bgRed("No .idat files found!"))
        message("Check worksheet for typos and if the barcode folder exists in the search path(s):")
        message(rsch.idat, "\nor\n", clin.idat)
        stop(
            crayon::bgRed(
                "No .idat files found for these sample(s)!  The case(s) may have not been run yet."
            )
        )
    }
    message("Files found: ")
    DataFrameMessage(allFi)
    allFi <- GetExternalIdats(allFi, ssheet, extr.idat)
    bcds <- paste0(basename(allFi))
    message("Checking if idats exist in run directory...")
    cur.idat <- basename(dir(
        path = runDir,
        pattern = "*.idat$",
        recursive = F
    ))
    idats_found <- bcds %in% cur.idat
    if (all(idats_found)) {
        message(".idat files already copied to run directory")
    } else{
        idats_needed <- allFi[!idats_found]
        copyBaseIdats(idats_needed)
    }
}


# FUN: Copies samplesheet to Desktop folder
moveSampleSheet <- function(methDir = NULL, runID = NULL, deskDir = NULL) {
    msgFunName(cpInLnk, "moveSampleSheet")
    if (is.null(runID)) {
        runID = paste0(basename(getwd()))
        message("Setting runID=", runID)
    }
    if (is.null(deskDir)) {
        deskDir = file.path(fs::path_home(), "Desktop", runID)
        message("Setting deskDir=", deskDir)
    }
    if (!dir.exists(deskDir)) {
        message("Creating Directory: ", deskDir)
        dir.create(deskDir)
    }
    if (is.null(methDir)) {
        methDir <- "/Volumes/CBioinformatics/Methylation/Clinical_Runs"
    }
    thisDir = file.path(methDir, runID)
    outFile = paste0(runID, "_samplesheet.csv")
    outputDir = file.path(deskDir, "samplesheet.csv")
    message("Copying from\nthisDir: ",
            thisDir,
            "\n",
            "to",
            "\noutputDir: ",
            outputDir)
    fs::file_copy(
        path = file.path(thisDir, "samplesheet.csv"),
        new_path = outputDir,
        overwrite = T
    )
    message("Renaming samplesheet in new folder as: ", outFile)
    file.rename(from = outputDir, to = file.path(deskDir, outFile))
}


#  Copy idats and Worksheets creation
writeFromRedcap <- function(df, samplesheet_ID, bn = NULL) {
    msgFunName(cpInLnk, "writeFromRedcap")
    if (is.null(bn)) {bn = file.path(getwd(), df$barcode_and_row_column)}
    message("~~~Writing from REDCap data samplesheet.csv:")
    names(df)
    samplesheet_csv = data.frame(
        Sample_Name = df[, "record_id"],
        DNA_Number = df[,"b_number"],
        Sentrix_ID = samplesheet_ID[, 1],
        Sentrix_Position = samplesheet_ID[, 2],
        SentrixID_Pos = df[, "barcode_and_row_column"],
        Basename = paste0(bn),
        RunID = df$run_number,
        MP_num = df$tm_number,
        tech = df$primary_tech,
        tech2 = df$second_tech,
        Date = df$arrived
    )
    DataFrameMessage(samplesheet_csv)
    write.csv(samplesheet_csv, file = "samplesheet.csv", quote = F, row.names = F)
}

# Ensures the correct version of redcapAPI is installed
check_REDCap_vers <- function(min_version = "2.7.4") {
    if (!"redcapAPI" %in% rownames(installed.packages())) {
        devtools::install_github('nutterb/redcapAPI', dependencies = TRUE,
                                 upgrade = "always", ask = F, type = "source")
    }
    current_vers <- as.character(utils::packageVersion("redcapAPI"))
    is_current <- utils::compareVersion(current_vers, min_version) >= 0
    if (!is_current) {
        if ("redcapAPI" %in% loadedNamespaces()) {
            try(unloadNamespace("redcapAPI"), TRUE)
        }
        devtools::install_github('nutterb/redcapAPI', dependencies = TRUE,
                                 upgrade = "always", ask = F, type = "source")
    }
    suppressPackageStartupMessages(library("redcapAPI", logical.return = TRUE))
}


#' FUN: Returns dataframe of REDCap search using a default header and fields
search.redcap <- function(rd_numbers,
                          token = NULL,
                          flds = NULL) {
    msgFunName(cpInLnk, "search.redcap")
    if (!require("redcapAPI")) {
        install.packages("redcapAPI", dependencies = T, ask = F)
    }
    check_REDCap_vers()
    if (is.null(token)) {
        message("You must provide an ApiToken!")
    }
    stopifnot(!is.null(token))
    rcon <- redcapAPI::redcapConnection(apiLink, token)
    if (is.null(flds)) {
        flds = c(
            "record_id",
            "b_number",
            "primary_tech",
            "second_tech",
            "run_number",
            "barcode_and_row_column",
            "accession_number",
            "tm_number",
            "arrived"
        )
    }
    result <- redcapAPI::exportRecordsTyped(
        rcon,
        records = rd_numbers,
        fields = flds,
        dag = F,
        factors = F,
        form_complete_auto = F,
        format = 'csv'
    )
    return(as.data.frame(result))
}


# FUN: Copies .idat files to your directory and saves samplesheet.csv
get.rd.info <- function(rd_numbers = NULL,
                        token = NULL,
                        sh_name = NULL) {
    msgFunName(cpInLnk, "get.rd.info")
    if (is.null(rd_numbers)) {
        message("Input RD-numbers using get.rd.info(rd_numbers)")
    }
    if (is.null(sh_name)) {
        sh_name = "samplesheet.csv"
    }
    if (is.null(token)) {
        message("You must provide an ApiToken!")
    }
    stopifnot(!is.null(token))
    result <- search.redcap(rd_numbers, token, NULL)
    samplesheet_ID = as.data.frame(stringr::str_split_fixed(result[, "barcode_and_row_column"], "_", 2))
    writeFromRedcap(result, samplesheet_ID) # writes API export as minfi dataframe sheet
    get.idats(csvNam = sh_name)  # copies idat files from return to current directory
    return(result)
}
