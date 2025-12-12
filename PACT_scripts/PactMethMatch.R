#!/usr/bin/env Rscript
## Script name: PactMethMatch.R
## Purpose: search REDCap for PACT samples with methylation & generate cnv PNG
## Date Created: September 2, 2021
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025

gb <- globalenv(); assign("gb", gb)

# Input Arguments -------------------------------------------------------------
args <- commandArgs(TRUE)
args[1] -> token       # REDCap Methylation API Token
args[2] -> inputSheet  # pactID

# Validate arguments
stopifnot(!is.na(token), !is.na(inputSheet))
readFlag <- grepl("\\.csv$", inputSheet)

message(paste("Now Running $HOME/PactMethMatch.R", token, inputSheet))

# REDCap Fields  --------------------------------------------------------------
meth_repo <- "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation"
cnv_outFolder = "/Volumes/molecular/Molecular/MethylationClassifier/CNV_PNG"

flds <- c("record_id", "b_number", "tm_number", "accession_number", "block",
          "diagnosis", "organ", "tissue_comments", "run_number", "nyu_mrn",
          "qc_passed", "arrived")

main_pkgs <- c("data.table", "openxlsx", "jsonlite", "RCurl", "readxl",
               "stringr", "tidyverse", "crayon", "tinytex", "systemfonts",
               "remotes", "dplyr", "fs")

# Message Inputs --------------------------------------------------------------
message("\n================ Parameters input ================\n")
message("token: ", token, "\ninputSheet: ", inputSheet, "\n")

options(repos = c(CRAN = "https://cran.r-project.org"))

if (!"devtools" %in% rownames(installed.packages())) {
    install.packages("devtools", ask = FALSE, dependencies = TRUE)
}

library("devtools")

# FUN: Install Homebrew if not installed -------------------------------------
install_brew <- function() {
    message("Installing Homebrew...")
    system(
        '/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"'
        , wait = TRUE
    )
}

# FUN: Check if system command exists
command_exists <- function(cmd) nzchar(Sys.which(cmd))

# Ensures that brew is found in the system R PATH
fix_brew_path <- function() {
    # pick the correct brew binary in one line
    arm_brew <- "/opt/homebrew/bin/brew"
    x64_brew <- "/usr/local/bin/brew"

    if (!file.exists(arm_brew) && !file.exists(x64_brew)) install_brew()

    brew_path <- ifelse(file.exists(arm_brew), arm_brew, x64_brew)
    brew_dir <- dirname(brew_path)

    # update session PATH
    old_path <- strsplit(Sys.getenv("PATH", ""), ":", fixed = TRUE)[[1]]
    new_path <- unique(c(brew_dir, old_path))
    Sys.setenv(PATH = paste(new_path, collapse = ":"))

    # prepare .Renviron update
    renv_file <- file.path(Sys.getenv("HOME"), ".Renviron")
    entry <- paste0('PATH="', paste(new_path, collapse = ":"), '"')
    lines <- if (file.exists(renv_file)) {
        readLines(renv_file, warn = FALSE)
    } else {
        character()
    }

    # only write if the exact entry is missing
    if (!any(grepl(entry, lines, fixed = TRUE))) {
        idx   <- grep("^PATH=", lines)
        lines <- if (length(idx) > 0) {
            lines[idx] <- entry
            lines
        } else {c(lines, entry)}
        writeLines(lines, renv_file)
    }
}


# Install Homebrew and packages if necessary
ensure_homebrew <- function() {
    pkgs <- c("gcc", "llvm", "lld", "open-mpi", "pkgconf", "gdal", "proj",
              "apache-arrow")
    fix_brew_path()
    if (!command_exists("brew")) install_brew(); fix_brew_path()
    installed_pkgs <- system2("brew", c("list", "--formula"),
                              stdout = T, stderr = F)

    for (pkg in pkgs) {
        if (!(pkg %in% installed_pkgs)) {
            system2("brew", c("install", pkg), wait = TRUE)
        }
    }
}

# Returns the path to the brew module
get_prefix <- function(pkg = "") {
    system2("brew", c("--prefix", pkg), stdout = TRUE, stderr = FALSE)
}

# Set environment variables dynamically
set_env_vars <- function() {
    brew_prefix = get_prefix()
    llvm_path = get_prefix("llvm")
    arrow_path = get_prefix("apache-arrow")
    Sys.setenv(
        CC = file.path(llvm_path, "bin/clang"),
        CXX = file.path(llvm_path, "bin/clang++"),
        OBJC = file.path(llvm_path, "bin/clang"),
        LDFLAGS = paste(
            paste0("-L", file.path(llvm_path, "lib")),
            paste0("-L", file.path(llvm_path, "lib", "c++")),
            paste0("-L", file.path(llvm_path, "lib", "unwind")),
            paste0("-Wl,-rpath,", file.path(llvm_path, "lib", "c++")),
            paste0("-Wl,-rpath,", file.path(llvm_path, "lib", "unwind")),
            "-lunwind"),
        CPPFLAGS = paste0("-I", file.path(llvm_path, "include")),
        PKG_CFLAGS = paste(
            paste0("-I", file.path(brew_prefix, "include")),
            paste0("-I", file.path(arrow_path, "include"))),
        PKG_LIBS = paste(
            paste0("-L", file.path(brew_prefix, "lib")),
            paste0("-L", file.path(llvm_path, "lib")),
            paste0("-L", file.path(arrow_path, "lib"), " -larrow")),
        LD_LIBRARY_PATH = file.path(brew_prefix, "lib"),
        R_LD_LIBRARY_PATH = paste(
            file.path(brew_prefix, "lib"), file.path(llvm_path, "lib/c++"), sep = ":"),
        DYLD_LIBRARY_PATH = file.path(arrow_path, "lib")
    )
    if (command_exists("gfortran")) Sys.setenv(FC = Sys.which("gfortran"))
}

# Installs the pak package
install_pak <- function() {
    tryCatch(
        install.packages("pak", repos = sprintf(
            "https://r-lib.github.io/p/pak/stable/%s/%s/%s",
            .Platform$pkgType, R.Version()$os, R.Version()$arch)),
        error = function(e) {
            install.packages(
                "pak", ask = FALSE, dependencies = TRUE,
                repos = "https://packagemanager.rstudio.com/all/latest")
        }
    )
}


#' Ensures Required Packages are Loaded the `pak` package installs any missing
#' @param pkgs A character vector of package names
ensure_packages <- function(pkgs) {
    installed_pk <- rownames(installed.packages())
    missing_pkgs <- setdiff(pkgs, installed_pk)
    if (length(missing_pkgs) > 0) {
        message(missing_pkgs)
        if (!"pak" %in% installed_pk) install_pak()
        library("pak")
        for (pkg in missing_pkgs) {
            tryCatch(
                pak::pkg_install(pkg, ask = FALSE),
                error = function(e) {
                    install.packages(pkg, ask = FALSE, dependencies = TRUE)
                }
            )
        }
    }
    stopifnot(all(sapply(pkgs, function(pkg) {
        suppressWarnings(suppressPackageStartupMessages(library(
            pkg, mask.ok = TRUE, character.only = TRUE, logical.return = TRUE
        )))
    })))
}


# Function to setup compilers, load and install necessary packages ------------
check_pkg_install <- function() {
    ensure_homebrew()
    Sys.unsetenv(c("CC", "CXX", "OBJC", "LDFLAGS", "CPPFLAGS", "PKG_CFLAGS",
                   "PKG_LIBS", "LD_LIBRARY_PATH", "R_LD_LIBRARY_PATH"))
    set_env_vars()
    ensure_packages(main_pkgs)
    if (!"mnp.v12epicv2" %in% rownames(installed.packages())) {
        devtools::source_url(file.path(meth_repo, "refs/heads/main/R/all_installer.R"))
    }
    suppressWarnings(suppressPackageStartupMessages(library("mnp.v12epicv2")))
}


# Check REDCap Version --------------------------------------------------------
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


# Function to check if a network drive is mounted
checkMounts <- function() {
    molecDrive <- "/Volumes/molecular/MOLECULAR LAB ONLY"
    zDrive <- "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular"
    if (!dir.exists(molecDrive)) {
        message("Network share is not mounted:\n", crayon::bgRed(zDrive))
        stop("Molecular shared drive is not mounted")
    }
}


# API Call functions -----
grabAllRecords <- function(flds, rcon) {
    message("Pulling REDCap data...")
    params = list(rcon, fields = flds, survey = F, dag = F, factors = F, form_complete_auto = F)
    dbCols <- do.call(redcapAPI::exportRecordsTyped, c(params))

    db <- as.data.frame(dbCols)

    if (nrow(db) == 0) {
        message("REDCap API connection failed!\n",
                "Check the Database for non-ASCII characters and verify API Token: ", token)
        stopifnot(nrow(db) > 0)
    }

    return(db)
}


# Messages the item matched in the database with the NGS number
message_matched <- function(item, dbInfo, ngsNum, i) {
    match_log <- file.path(fs::path_home(), "Desktop", paste0(inputSheet, "_match_log.tsv"))
    match_line <- sprintf("Match found for '%s' (%s) for %s in: \"%s\" column",
                          item, dbInfo, ngsNum, i)
    message(match_line)
    cat(match_line, file = match_log, append = TRUE, sep = "\n")
}


# Searches the REDCap db against the queryList of items -----------------------
searchDb <- function(queryList, db) {
    res <- data.frame()
    for (idx in 1:length(queryList)) {
        item <- queryList[idx]
        ngsNum = paste(names(item)[1])
        for (i in colnames(db)) {
            ngsMatch <- which(grepl(item, db[, i]))
            if (length(ngsMatch) > 0) {
                message_matched(item, db[ngsMatch, i], ngsNum, i)
                dbMatch <- db[ngsMatch,]
                dbMatch$Test_Number <- ngsNum
                res <- rbind(res, dbMatch)
            }
        }
    }
    return(res)
}

# Filters out workbooks from xlsx and other files that may be in the directory
filterFiles <- function(potentialFi) {
    wbFiles <- grep("\\.xlsm$|book", basename(potentialFi), value = TRUE)
    if (!any(grepl("\\.xlsm$", wbFiles))) {
        message("\nNo .xlsm worksheet found. Checking .xlsx files and others...\n")
    }
    filteredFiles <- wbFiles[!grepl("\\$", wbFiles)]
    return(filteredFiles)
}

# Checks alternative directories if the file with expected name is not found
getAltPath <- function(inputFi) {
    message(crayon::bgRed("PACT run worksheet not found:"), "\n", inputFi)
    message("Checking other files in PACT folder: ", basename(dirname(inputFi)))
    potentialFi <- list.files(dirname(inputFi), full.names = TRUE)
    altFi <- filterFiles(potentialFi)
    if (length(altFi) > 0) {
        chosenFi <- altFi[1]
        message(crayon::bgGreen("Using this workbook instead:"), basename(chosenFi))
        return(chosenFi)
    } else {
        stop("No alternative file found.")
    }
}

# Returns the folder path to the PACT run where inputSheet is the PACT ID
getPactFolder <- function(inputSheet) {
    drive <- file.path("", "Volumes", "molecular", "MOLECULAR LAB ONLY")
    folder <- file.path(drive, "NYU PACT Patient Data", "Workbook")
    runyr <- stringr::str_split_fixed(inputSheet,"-",3)[,2]
    runFolder <- file.path(folder, paste0("20", runyr), inputSheet)
    stopifnot(dir.exists(runFolder))
    return(runFolder)
}

# Returns the file path to the PACT run worksheet
# Returns the file path to the PACT run worksheet
getFilePath <- function(inputSheet) {
    runFolder <- getPactFolder(inputSheet)
    inputFi <- file.path(runFolder, paste0(inputSheet, ".xlsm"))
    if (!file.exists(inputFi)) {
        inputFi <- getAltPath(inputFi)
    }
    message("Using PACT run worksheet file:\n", inputFi)
    return(inputFi)
}

# Parses the input file for the "PhilipsExport" tab
parseWorksheet <- function(inputFi) {
    sheet2Read <- "PhilipsExport"
    message("Reading the file:\n", inputFi)
    shNames <- readxl::excel_sheets(inputFi)
    message("Excel sheet names:\n", paste(shNames, collapse = "\n"))
    stopifnot(!is.null(shNames) & length(shNames) > 2)
    sh <- which(grepl(sheet2Read, shNames, ignore.case = T))[1]
    sheet_to_read <- shNames[sh]
    pact_columns <-
        c("Tumor Specimen ID",
          "Normal Specimen ID",
          "Tumor DNA/RNA Number",
          "MRN",
          "Test Number")
    vals2find <-  suppressMessages(as.data.frame(
        readxl::read_excel(
            inputFi,
            sheet = sheet_to_read,
            skip = 3,
            col_types = "text"
        )
    ))
    vals2find <- vals2find[, pact_columns]
    vals2find <- vals2find[!is.na(vals2find[,1]),]
    return(vals2find)
}

# Parses the input file depending on if the input is a csv file or a xlsx file path
getCaseValues <- function(inputSheet, readFlag) {
    isSamSheet <- stringr::str_detect(inputSheet, "-SampleSheet")
    isFilePath <- stringr::str_detect(inputSheet, .Platform$file.sep)

    if (readFlag && isSamSheet) {
        message("Parsing Data from Demux SampleSheet.csv file...")
        vals2find <- utils::read.csv(inputSheet, skip = 19)[, c(6, 7, 9)]
        vals2find <- as.data.frame(vals2find[!grepl("H20|SERACARE|HAPMAP", vals2find[, 2]),])
        return(vals2find)
    }

    if (readFlag && !isSamSheet) {
        message("Parsing Data from .csv file that is not a Demux SampleSheet...")
        vals2find <- read.csv(inputSheet)
        if (ncol(vals2find) > 1) {
            vals2find <- unlist(lapply(vals2find, identity))
        }
        return(as.data.frame(unique(vals2find[vals2find != ""])))
    }

    if (!readFlag && isFilePath) {
        message("Parsing Data from .xlsm/.xlsx file path...")
        vals2find <- parseWorksheet(inputSheet)
        return(vals2find)
    }
    # Default case: inputSheet is a PACT ID, get the file path
    message("Parsing Data from PACT RUN ID: ", inputSheet,
            " finding run worksheet...")
    pact_xlsm <- getFilePath(inputSheet)
    vals2find <- parseWorksheet(pact_xlsm)
    return(vals2find)
}

# Generates a query list of items to search for in the REDCap database
genQuery <- function(dbCol, vals2find) {
    currCol <- vals2find[, dbCol]
    toKeep <- which(currCol != 0 & !is.na(currCol) & currCol != "")
    q1 <- currCol[toKeep]
    stopifnot("Test Number" %in% colnames(vals2find))
    names(q1) <- vals2find$`Test Number`[toKeep]
    return(q1)
}

# Queries the REDCap database for the items in the query list
queryCases <- function(vals2find, db) {
    # Label each item with the NGS-number
    queryList <- unlist(lapply(1:ncol(vals2find), function(i) {
        genQuery(i, vals2find)
    }), use.names = TRUE)
    theTScases <- queryList[stringr::str_detect(queryList, "TS|TB|TC")]
    theTScases <- sapply(theTScases, function(x) {
        paste(stringr::str_split_fixed(x, "-", 3)[1, 1:2], collapse = "-")
    })
    queryList <- c(queryList, theTScases)
    queryList <- queryList[!duplicated(queryList)]
    methQuery <- searchDb(queryList, db)
    return(unique(methQuery))
}

# Generates the year path for the report link in the xlsx output file
get_year_paths <- function(output) {
    output_final <- output[!is.na(output$run_number),]
    yearSplit <- stringr::str_split_fixed(output_final$run_number, "-", 2)[, 1]
    yearSplit <- gsub("MC", "", yearSplit)
    yearPath <- lapply(yearSplit, function(yr) {
        if (nchar(yr) > 2) {
            yr <- substring(yr, 3)
        }
        paste0("20", yr)
    })
    return(yearPath)
}

# Appends the smb file paths to the report path columns in the output excel sheet
addOutputLinks <- function(output) {
    if (ncol(output) == 0) return(output)

    winpath <- "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular/Molecular/MethylationClassifier"
    yearPath <- get_year_paths(output)
    output$report_complete <- ifelse(!is.na(output$run_number), "YES", "NOT_YET_RUN")
    out_htmls <- paste0(output$record_id, ".html")
    output$'Report Link' <- file.path(winpath, yearPath, output$run_number, out_htmls)
    output$'Report Link'[is.na(output$run_number)] <- ""
    output$'Report Path' <- output$'Report Link'

    # Handle the case where the run number does not contain "MGDM"
    oldRun <- !grepl("MGDM", output$run_number)
    if (any(oldRun)) {
        output$report_complete[oldRun] <- "NOT_YET_RUN"
        output$`Report Link`[oldRun] <- ""
        output$`Report Path`[oldRun] <- ""
    }
    return(output)
}

# Finds and fills in missing NGS numbers in the output data frame for matching fields
FillMissingNGS <- function(output, vals2find) {
    rows2fill <- which(is.na(output$Test_Number))
    if (length(rows2fill) > 0) {
        for (row in rows2fill) {
            accessionNum <- output$accession_number[row]
            matchedIdx <- which(vals2find$`Tumor Specimen ID` == accessionNum)

            if (length(matchedIdx) > 0) {
                output$Test_Number[row] <- vals2find$`Test Number`[matchedIdx[1]]
            }
        }
    }
    rowsStillMissing <- which(is.na(output$Test_Number))
    if (length(rowsStillMissing) > 0) {
        warning("Some samples still missing NGS Numbers:\n",
                paste(output$record_id[rowsStillMissing], collapse = "\n"))
    }
    return(output)
}

# Corrects any NA or missing NGS numbers, and appends output links
modifyOutput <- function(output, vals2find) {
    if (!("Test Number" %in% colnames(vals2find))) {
        return(output)
    }
    if (length(vals2find$`Test Number`) == 0) {
        vals2find$`Test Number` <- ""
    }
    # Fixes NA's if any found matches did not have NGS numbers
    if (all(output$Test_Number %in% vals2find$`Test Number`)) {
        NGSmissing <- F
        message("All NGS Test Numbers Found in Methylation Database")
    } else{
        NGSmissing <- T
        message("Not all NGS do not have methylation")
    }

    if (NGSmissing == T) {
        output <- FillMissingNGS(output, vals2find)
    }

    output <- addOutputLinks(output)
    return(output)
}

# Returns the volume paths to the methylation reports in the dataframe
GetVolumePaths <- function(methData) {
    smb_path <- "smb://shares-cifs.nyumc.org/apps/acc_pathology"
    checkPaths <- stringr::str_replace_all(methData$`Report Path`, smb_path, "/Volumes")
    checkPaths <- checkPaths[checkPaths != "" & !is.na(checkPaths)]
    checkPaths <- checkPaths[stringr::str_detect(checkPaths, "MGDM")]
    return(checkPaths)
}

# Checks if paths to the methylation reports are valid and fixes any broken
CheckMethPaths <- function(methData) {
    for (i in 1:length(methData$`Report Path`)) {
        currPath <- methData$`Report Path`[i]
        currSplit <- stringr::str_split_fixed(currPath, "/", 11)[1, ]
        if (stringr::str_detect(currSplit[10], "MGDM") == FALSE) {
            next
        }
        runYear <- stringr::str_split_fixed(currSplit[10], "-", 2)[1, 1]
        runYear <- paste0("20", runYear)
        currSplit[9] <- runYear
        newPath <- paste(currSplit, collapse = "/")
        methData[i, "Report Path"] <- newPath
    }
    checkPaths <- GetVolumePaths(methData)
    anyPathsFalse <- file.exists(checkPaths) == FALSE
    if (any(anyPathsFalse)) {
        message("Fixing broken file paths...")
        toReplace <- basename(checkPaths[anyPathsFalse])
        mainDirs <- dirname(checkPaths[anyPathsFalse])
        mainDirs <- unique(mainDirs)
        for (x in 1:length(mainDirs)) {
            if (!dir.exists(mainDirs[x])) {
                correct_dir <- dir(
                    path = dirname(mainDirs[x]),
                    pattern = basename(mainDirs[x]),
                    full.names = TRUE
                )
                mainDirs[x] <- correct_dir
            }
        }
        for (missing in toReplace) {
            patt <- stringr::str_split_fixed(missing, ".html", 2)[1, 1]
            file_found <- dir(path = mainDirs, pattern = patt, full.names = TRUE)
            if (length(file_found) > 0) {
                toSwap <- which(grepl(missing, methData$`Report Path`))
                newPath <- stringr::str_replace(
                    methData$`Report Path`[toSwap],
                    missing,
                    basename(file_found)
                )
                message("Updating file path:\n", newPath)
                methData$`Report Path`[toSwap] <- newPath
            }
        }
        checkPaths <- GetVolumePaths(methData)
        anyPathMissed <- file.exists(checkPaths) == FALSE
        if (any(anyPathMissed)) {
            msg1 <- "Some paths to html reports need editing in MethylMatch.xlsx sheet!"
            msg2 <- "Fix the following paths in worksheet 'Report Path' column that do not exist:"
            message(crayon::bgRed(msg1))
            message(crayon::bgRed(msg2), "\n")
            message(paste(checkPaths[anyPathsFalse], collapse = "\n"), "\n")
        }
    }

    return(methData)
}

# Adds hyperlinks to the report links in the output excel file
addExcelLink <- function(output, fiLn, wb, runId) {
    x <- c(output$'Report Link'[fiLn])
    names(x) <- paste0(output$record_id[fiLn], ".html")
    class(x) <- "hyperlink"
    colNum <- which(colnames(output) == "Report Link")
    openxlsx::writeData(wb, sheet = runId, x = x, startCol = colNum, startRow = fiLn + 1)
}

# Writes the openxlsx workbook file to the Desktop
createXlFile <- function(runId, output) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, runId)
    output <- CheckMethPaths(methData = output)
    openxlsx::writeData(wb, sheet = runId, x = output)
    for (fiLn in 1:length(output$'Report Link')) {
        if (output$'Report Link'[fiLn] != '') {
            addExcelLink(output, fiLn, wb, runId)
        }
    }
    outFi <- file.path(fs::path_home(),"Desktop",
                       paste0(runId,"_MethylMatch.xlsx"))
    openxlsx::saveWorkbook(wb, outFi, overwrite = T)
    return(outFi)
}

# Uses Rcurl to post Json form to REDCap
postData <- function(rcon, record) {
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox = T)
    res <-
        RCurl::postForm(
            rcon$url,
            token = rcon$token,
            content = 'record',
            format = 'json',
            type = 'flat',
            data = datarecord,
            returnContent = 'nothing',
            returnFormat = 'csv'
        )
    cat(res)
}

# Uploads the xlsx file to REDCap and sends an email notification
emailFile <- function(runId, outFi, rcon) {
    record = data.frame(record_id = runId, run_number = runId)
    postData(rcon, record)
    isDone <-
        redcapAPI::exportRecordsTyped(
            rcon,
            factors = F,
            records = record$record_id,
            fields = c("record_id", "other_file")
        )
    if (length(isDone$other_file) == 0) {
        body <- list(
            token = rcon$token,
            content = 'file',
            action = 'import',
            record = runId,
            field = "other_file",
            file = httr::upload_file(outFi),
            returnFormat = 'csv'
        )
        res <-
            tryCatch(
                httr::POST(url = rcon$url, body = body, config = rcon$config),
                error = function(cond) {
                    list(status_code = "200")
                }
            )
        if (res$status_code == "200") {
            message("REDCap file upload failed:\n", outFi)
        }else{
            message("REDCap file upload successful:\n", outFi)
        }
    }
    record$comments <- "pact_sample_list_email"
    postData(rcon, record)
    message("------", "Email Notification Created", "------")
}

# Grabs the run ID from the input xlsx sheet name
grab_run_id <- function(readFlag, inputSheet) {
    runId <-
        ifelse(readFlag == T, substr(inputSheet, 1, nchar(inputSheet) - 4), inputSheet)
    if (stringr::str_detect(inputSheet, ".xls")) {
        runId <-
            substr(basename(inputSheet), 1, nchar(basename(inputSheet)) - 5)
    }
    return(runId)
}


process_values <- function(vals2find, db) {
    output <- data.frame()

    if (nrow(vals2find) != 0) {
        output <- queryCases(vals2find, db)
    }
    if (nrow(output) == 0) {
        warning("No Methylation Cases on this PACT run, generating blank file")
        output[1, ] <- "NONE"
    } else {
        rownames(output) <- 1:nrow(output)
    }
    if (ncol(output) == 0) {
        message(crayon::bgRed("No Methylation on this PACT run!"))
        return(output)
    }
    if (readFlag == T) {
        output <- modifyOutput(output, vals2find)
        write.csv(output, file = "meth_sample_data.csv", quote = F, row.names = F)
        return(output)
    }
    return(output)
}


# Grabs REDCap data and finds matches to inputSheet columns to fields
getOuputData <- function(token, flds, inputSheet, readFlag) {
    apiUrl = "https://redcap.nyumc.org/apps/redcap/api/"
    rcon <- redcapAPI::redcapConnection(apiUrl, token)

    vals2find <- getCaseValues(inputSheet, readFlag)

    if (class(vals2find) != "data.frame") {
        vals2find <- as.data.frame(vals2find)
    }

    # Get entire REDCap Database matrix
    db <- grabAllRecords(flds, rcon)
    output <- process_values(vals2find, db)

    runId <- grab_run_id(readFlag, inputSheet)
    output <- modifyOutput(output, vals2find)
    toDrop <- is.na(output$Test_Number)

    if (any(toDrop)) {
        hasNGS <- which(grepl("NGS", output$tm_number))
        message("NGS Number found in 'tm_number' field of REDCap!")
        if (length(hasNGS) > 0) {
            output$Test_Number[hasNGS] <- output$tm_number[hasNGS]
        }
    }

    outFi <- createXlFile(runId, output)
    emailFile(runId, outFi, rcon)
    return(output)
}

# Installs correct version of EPICv2 manifest and minfi from GitHub
minfi_install <- function() {
    Sys.setenv(R_COMPILE_AND_INSTALL_PACKAGES = "always")
    devtools::install_github(
        "mwsill/minfi",
        upgrade = "always",
        force = T, dependencies = T, type = "source", auth_token = NULL
    )
    devtools::install_github(
        "mwsill/IlluminaHumanMethylationEPICv2manifest",
        upgrade = "always",
        force = T, dependencies = T, type = "source", auth_token = NULL
    )
}

# Checks the correct package versions are installed and loads them
source_pkg_vers <- function() {
    minfiVers <- as.character(utils::packageVersion("minfi"))
    if (minfiVers != "1.43.1") minfi_install()

    v2_manifest <- "IlluminaHumanMethylationEPICv2manifest"
    epicVers <- as.character(utils::packageVersion(v2_manifest))
    if (epicVers != "0.1.0") minfi_install()

    v2Pkg_needed <- !"mnp.v12epicv2" %in% rownames(installed.packages())
    v2Con_needed <- !"conumee2.0" %in% rownames(installed.packages())
    if (v2Pkg_needed | v2Con_needed) {
        source(
            "/Volumes/CBioinformatics/Methylation/Rscripts/install_epic_v2_classifier.R"
        )
    }

    stopifnot(library("conumee2.0", logical.return = T))
    stopifnot(library("minfi", logical.return = T))
    stopifnot(library("IlluminaHumanMethylationEPICv2manifest", logical.return = T))
    stopifnot(library("mnp.v12epicv2", logical.return = T))
}

# Loads additional packages and functions for generating CNV PNGs
sourceFuns2 <- function() {
    mainHub <- file.path(meth_repo, "main", "R")
    script.list <- c("SetRunParams.R", "CopyInputs.R")
    scripts <- file.path(mainHub, script.list)
    invisible(lapply(scripts, function(i) {
        suppressPackageStartupMessages(devtools::source_url(i))
    }))
    gb$setDirectory(getwd())
    source_pkg_vers()
    return(gb$defineParams())
}

# Messages the RD-numbers with idats and sets the API token as global variable
msg_rd_num <- function(rds, token) {
    message("\nRD-numbers with idats:\n", paste(rds, collapse = "\n"))
    assign("rds", rds)
    message("--------", crayon::bgMagenta("Starting CNV PNG Creation"), "--------")
    ApiToken <- token
    assign("ApiToken", ApiToken)
}

# Generates a minfi samplesheet from REDCap RD numbers and copies idat files
get_pact_rds <- function(rd_numbers, token) {
    result <- gb$search.redcap(rd_numbers, token)
    result <- result[!is.na(result$barcode_and_row_column),]
    samplesheet_ID = as.data.frame(stringr::str_split_fixed(result[,"barcode_and_row_column"],"_",2))
    gb$writeFromRedcap(result, samplesheet_ID)  # writes API export as minfi dataframe sheet
    gb$get.idats()  # copies idat files from return to current directory
}

# Messages final output of png CNV file output success on the user Desktop
msgCreated <- function(mySentrix) {
    pngFiles <- paste0(file.path(fs::path_home(), "Desktop", mySentrix[, 1]),"_cnv.png")
    cnvMade <- file.exists(pngFiles)
    if (any(cnvMade == F)) {
        message("The following failed to be created:")
        print(pngFiles[!cnvMade])
        message("Try running again or check GitHub troubleshooting")
    }
    if (any(cnvMade == T)) {
        message("The following were created successfully:")
        print(pngFiles[cnvMade])
    }
}

# Function to check if CNV PNGs already exist on Z-drive to skip them
CheckIfPngExists <- function(rds, outFolder = NULL) {
    if (is.null(outFolder)) {
        outFolder <- cnv_outFolder
    }
    outFiles <- file.path(outFolder, paste0(rds, "_cnv.png"))
    finished <- file.exists(outFiles)
    if (any(finished)) {
        message(crayon::bgGreen("Cases below already exist and will be skipped:"))
        message(paste(capture.output(outFiles[finished]), collapse = '\n'))

        rds <- rds[!finished]

        if (length(rds) > 0) {
            message(crayon::bgGreen("The following CNV PNG will be generated:"))
            message(paste(capture.output(rds), collapse = '\n'))
        }
    }
    return(rds)
}

# Main function that generates the CNV PNG for a single idat file path
get_pact_cnv <- function(samName, samEpic, idatPath = getwd()) {
    sam_out_png <- file.path(fs::path_home(), "Desktop",
                             paste0(samName, "_cnv.png"))
    pathEpic <- file.path(idatPath, samEpic)
    RGsetEpic <- minfi::read.metharray(pathEpic, verbose = TRUE, force = TRUE)
    MsetEpic <- minfi::preprocessIllumina(RGsetEpic, bg.correct = TRUE,
                                          normalize = "controls")
    cnv_obj <- mnp.v12epicv2::MNPcnv(MsetEpic, sex = NULL, main = samName)
    chrAll <- paste0("chr", 1:22)
    message("Saving file to:\n", sam_out_png)
    png(filename = sam_out_png, width = 1820, height = 1040, res = 150)
    conumee2.0::CNV.genomeplot(
        cnv_obj,
        chr = chrAll,
        main = samName,
        bins_cex = "sample_level",
        cols = c("darkred", "salmon", "lightgrey", "lightgreen", "darkgreen")
    )
    invisible(dev.off())
}

# Function to generate CNV PNGs for all samples with sentrix IDs
generate_new_cnv <- function(targets) {
    has_idat <- targets[, "SentrixID_Pos"] %like% "_R0"
    targets <- targets[has_idat, ]
    all_installer <- file.path(meth_repo, "refs/heads/main/R/all_installer.R")

    if ("conumee2.0" %in% rownames(installed.packages()) == F) {
        devtools::source_url(all_installer)
    }
    if ("mnp.v12epicv2" %in% rownames(installed.packages()) == F) {
        devtools::source_url(all_installer)
    }
    if (!requireNamespace("mnp.v12epicv2", quietly = TRUE)) {
        stop(crayon::bgRed("The classifier is not installed no CNV will generate!"))
    }
    if (nrow(targets) > 0) {
        sam_names <- as.character(targets[, 1])
        sentrix.ids <- as.character(targets$SentrixID_Pos)
        mapply(get_pact_cnv, sam_names, sentrix.ids)
    } else{
        message("The RD-number(s) do not have idat files in REDCap:/n")
        print(targets)
    }
    msgCreated(targets)
    while (!is.null(dev.list())) {dev.off()}
}

# Reads the sample sheet and returns a dataframe of samples to be processed
GetSampleList <- function(rds, sampleSheet = "samplesheet.csv") {
    targets <- read.csv(sampleSheet)
    toDrop <- targets$Sample_Name %in% rds
    targets <- targets[toDrop, ]
    rownames(targets) <- 1:nrow(targets)
    return(targets)
}

# Loads additional functions from GitHub and a TryCatch to generate CNV PNGs
try_cnv_make <- function(rds, token) {
    msg_rd_num(rds, token)
    sourceFuns2()
    get_pact_rds(rds, token)
    targets <- GetSampleList(rds)

    tryCatch(
        expr = {
            generate_new_cnv(targets)
        },
        error = function(e) {
            message("The following error occured:\n", e)
            message("\nTry checking the troubleshooting section on GitHub:\n")
            git_url <- "https://github.com/NYU-Molecular-Pathology/Methylation/"
            message(git_url, "blob/main/PACT_scripts/README.md\n")
            stop(crayon::bgRed("CNV PNG generation failed"))
        }
    )
}

# Copies any CNV PNGs from the Desktop to the Molecular Z-drive folder
copy_output_png <- function(outFolder = NULL) {
    if (!"dplyr" %in% loadedNamespaces()) library("dplyr")

    if (is.null(outFolder)) {
        outFolder <- cnv_outFolder
    }
    desk <- file.path(fs::path_home(), "Desktop")
    the.cnvs <- dir(desk, "_cnv.png", full.names = T) %>% file.info() %>%
        tibble::rownames_to_column() %>% filter(as.Date(ctime) == Sys.Date()) %>%
        pull(rowname)
    if (length(the.cnvs) > 0) {
        savePath <- file.path(outFolder, basename(the.cnvs))
        message("\nCopying png files to Molecular folder:\n", outFolder, "\n")
        message(paste(capture.output(the.cnvs), collapse = '\n'))

        fs::file_copy(path = the.cnvs, new_path = savePath)

        if (any(!file.exists(savePath))) {
            message("The following failed to copy from the desktop:\n")
            print(basename(savePath[!file.exists(savePath)]))
            message(crayon::bgRed(
                "Manually copy any methylation CNV PNGs that failed to copy"))
        }
    } else{
        message("No CNV files found on Desktop to copy")
    }
}

# Checks the output generated to see if any CNV PNGs need to be created
queue_cnv_maker <- function(output, token) {
    rds <- output$record_id[output$report_complete == "YES"]
    rds <- rds[grep("^RD-", rds)]

    if (all(!is.null(rds)) == F | all(!is.na(rds)) == F | length(rds) == 0) {
        return(message(crayon::bgGreen(
            "The PACT run has no cases with methylation.")))
    }

    rds <- CheckIfPngExists(rds)

    if (length(rds) > 0) {
        try_cnv_make(rds, token)
        copy_output_png()
    } else{
        message(crayon::bgGreen(
            "No CNV png images to generate. Check the output directory."
            ))
    }
}

# MAIN Execution start -----
check_pkg_install()
check_REDCap_vers() # Check REDCap API version
checkMounts()
output <- getOuputData(token, flds, inputSheet, readFlag)

# CNV PNG Creation -------------------------------------
if (ncol(output) > 0) {
    queue_cnv_maker(output, token)
}
