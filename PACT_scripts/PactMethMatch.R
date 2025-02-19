#!/usr/bin/env Rscript
## Script name: PactMethMatch.R
## Purpose: search REDCap for PACT samples with methylation & generate cnv PNG
## Date Created: September 2, 2021
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024

gb <- globalenv(); assign("gb", gb)

# Input Arguments -------------------------------------------------------------
args <- commandArgs(TRUE)
token <- args[1]
inputSheet <- args[2]

# Validate arguments
stopifnot(!is.na(token), !is.na(inputSheet))
readFlag <- grepl("\\.csv$", inputSheet)

message(paste("Now Running $HOME/PactMethMatch.R", token, inputSheet))

# REDCap Fields  --------------------------------------------------------------
flds <- c("record_id", "b_number", "tm_number", "accession_number", "block",
          "diagnosis", "organ", "tissue_comments", "run_number", "nyu_mrn",
          "qc_passed")

main_pkgs <- c("data.table", "openxlsx", "jsonlite", "RCurl", "readxl",
               "stringr", "tidyverse", "crayon", "tinytex", "systemfonts",
               "remotes")

# Message Inputs --------------------------------------------------------------
message("\n================ Parameters input ================\n")
message("token: ", token, "\ninputSheet: ", inputSheet, "\n")

if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", ask = F, dependencies = T)
}

# FUN: Check if system command exists
command_exists <- function(cmd) nzchar(Sys.which(cmd))

# Ensures that brew is found in the system R PATH
fix_brew_path <- function() {
    # Define the target directories to add to the PATH
    target_dirs <- c("/opt/homebrew/bin", "/usr/local/bin")
    target_dirs <- target_dirs[dir.exists(target_dirs)]
    
    # Get the current PATH in the R session
    current_path <- strsplit(Sys.getenv("PATH"), ":")[[1]]
    
    # Add the target directories to the PATH if they are not already present
    new_path <- unique(c(target_dirs, current_path))
    Sys.setenv(PATH = paste(new_path, collapse = ":"))
    message("PATH updated for the current session.")
    
    # Persist the change in ~/.Renviron for future sessions
    renviron_file <- file.path(Sys.getenv("HOME"), ".Renviron")
    target_entry <- paste0('PATH="', paste(new_path, collapse = ":"), '"')
    
    if (file.exists(renviron_file)) {
        renviron_content <- readLines(renviron_file, warn = FALSE)
    } else {
        renviron_content <- character(0)
        writeLines(c(renviron_content, target_entry), renviron_file)
        message("PATH added to ~/.Renviron.")
    }
    
    if (!any(grepl(target_dirs, renviron_content))) {
        # Append the new PATH if no PATH is defined
        writeLines(c(renviron_content, target_entry), renviron_file)
        message("PATH added to ~/.Renviron.")
    } else {
        # Update the PATH entry if it exists
        #renviron_content <- sub("^PATH=.*", target_entry, renviron_content)
        #writeLines(renviron_content, renviron_file)
        message("PATH already updated in ~/.Renviron.")
    }
}

# Install Homebrew and packages if necessary
ensure_homebrew <- function() {
    pkgs <- c("gcc", "llvm", "lld", "open-mpi", "pkgconf", "gdal", "proj",
              "apache-arrow")
    fix_brew_path()
    
    if (!command_exists("brew")) {
        message("Installing Homebrew...")
        system('/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"', wait = TRUE)
    }
          
    installed_pkgs <- system2("brew", c("list", "--formula"),
                              stdout = T, stderr = F)

    for (pkg in pkgs) {
        if (!(pkg %in% installed_pkgs)) {
            message("Installing ", pkg, " via Homebrew...")
            system2("brew", c("install", pkg), wait = TRUE)
        }
    }
}

get_prefix <- function(pkg = "") {
    system2("brew", c("--prefix", pkg), stdout = TRUE, stderr = FALSE)
}

# Set environment variables dynamically
set_env_vars <- function() {
    brew_prefix = get_prefix()
    llvm_path = get_prefix("llvm")
    mpi_path = get_prefix("open-mpi")
    arrow_path = get_prefix("apache-arrow")
    Sys.setenv(
        CC = file.path(llvm_path, "bin/clang"),
        CXX = file.path(llvm_path, "bin/clang++"),
        OBJC = file.path(llvm_path, "bin/clang"),
        LDFLAGS = paste(
            paste0("-L", file.path(llvm_path, "lib")),
            paste0("-L", file.path(llvm_path, "lib/c++")),
            paste0("-Wl,-rpath,", file.path(llvm_path, "lib/c++")),
            paste0("-L", file.path(llvm_path, "lib/unwind"), "-lunwind")
        ),
        CPPFLAGS = paste0("-I", file.path(llvm_path, "include")),
        PKG_CFLAGS = paste(
            paste0("-I", file.path(brew_prefix, "include")),
            if (nzchar(arrow_path))
                paste0("-I", file.path(arrow_path, "include"))
            else
                NULL
        ),
        PKG_LIBS = paste(
            paste0("-L", file.path(brew_prefix, "lib")),
            paste0("-L", file.path(llvm_path, "lib")),
            if (nzchar(arrow_path)) {
                paste0("-L", file.path(arrow_path, "lib"), " -larrow")
            }else { NULL }
        ),
        LD_LIBRARY_PATH = file.path(brew_prefix, "lib"),
        R_LD_LIBRARY_PATH = paste(
            file.path(brew_prefix, "lib"),
            file.path(llvm_path, "lib/c++"),
            sep = ":"
        ),
        DYLD_LIBRARY_PATH = file.path(arrow_path, "lib"),
        PATH = paste(file.path(mpi_path, "bin"), Sys.getenv("PATH"), sep = ":")
    )
    if (command_exists("gfortran")) Sys.setenv(FC = Sys.which("gfortran"))
}


# FUN: Returns all packages that are not installed ----------------------------
check_needed <- function(pkgs) {
    neededPkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
    return(neededPkgs)
}


# Function to setup compilers, load and install necessary packages ------------
check_pkg_install <- function() {
    params <- list(
        dependencies = TRUE,
        ask = FALSE,
        update = "never",
        quiet = TRUE,
        repos = 'http://cran.us.r-project.org',
        Ncpus = 4
    )

    ensure_homebrew()
    Sys.unsetenv(c("CC", "CXX", "OBJC", "LDFLAGS", "CPPFLAGS", "PKG_CFLAGS",
                   "PKG_LIBS", "LD_LIBRARY_PATH", "R_LD_LIBRARY_PATH"))
    set_env_vars()

    needed_pkgs <- check_needed(main_pkgs)

    if (!requireNamespace("pak", quietly = T)) {
        install.packages("pak", dependencies = T, ask = F, type = "binary")
    }

    if (length(needed_pkgs) > 0) {
        for (n_pkg in needed_pkgs) {
            tryCatch({
                pak::pkg_install(n_pkg, ask = F)
            }, error = function(e) {
                message("Pak install failed, attempting binary installation.")
                params$type <- "binary"
                do.call(install.packages, c(list(pkgs = n_pkg), params))
            })
        }
    }

    sapply(main_pkgs, library, character.only = TRUE, logical.return = TRUE)

    if (!requireNamespace("mnp.v12epicv2", quietly = TRUE)) {
        devtools::source_url(
            "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/refs/heads/main/R/all_installer.R")
        source("/Volumes/CBioinformatics/Methylation/Rscripts/install_epic_v2_classifier.R")
    }
}


# Check REDCap Version --------------------------------------------------------
check_REDCap_vers <- function(min_version = "2.7.4") {
    if (!requireNamespace("redcapAPI", quietly = TRUE)) {
        devtools::install_github('nutterb/redcapAPI', dependencies = TRUE,
                                 upgrade = "always", ask = F, type = "source")
    }
    current_vers <- as.character(utils::packageVersion("redcapAPI"))
    is_current <- utils::compareVersion(current_vers, min_version) >= 0
    if (!is_current) {
        install.packages("redcapAPI", dependencies = T, ask = F, update = T)
    }
    library("redcapAPI", logical.return = TRUE)
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


message_matched <- function(item, dbInfo, ngsNum, i) {
    message(sprintf("Match found for '%s' (%s) for %s in: \"%s\" column", item, dbInfo, ngsNum, i))
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


filterFiles <- function(potentialFi) {
    wbFiles <- grep("\\.xlsm$|book", basename(potentialFi), value = TRUE)
    if (!any(grepl("\\.xlsm$", wbFiles))) {
        message("\nNo .xlsm worksheet found. Checking .xlsx files and others...\n")
    }
    filteredFiles <- wbFiles[!grepl("\\$", wbFiles)]
    return(filteredFiles)
}


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


getPactFolder <- function(inputSheet) {
    drive <- file.path("", "Volumes", "molecular", "MOLECULAR LAB ONLY")
    folder <- file.path(drive, "NYU PACT Patient Data", "Workbook")
    runyr <- stringr::str_split_fixed(inputSheet,"-",3)[,2]
    runFolder <- file.path(folder, paste0("20", runyr), inputSheet)
    stopifnot(dir.exists(runFolder))
    return(runFolder)
}


getFilePath <- function(inputSheet) {
    runFolder <- getPactFolder(inputSheet)
    inputFi <- file.path(runFolder, paste0(inputSheet, ".xlsm"))
    if (!file.exists(inputFi)) {
        inputFi <- getAltPath(inputFi)
    }
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


getCaseValues <- function(inputSheet, readFlag) {
    isSamSheet <- stringr::str_detect(inputSheet, "-SampleSheet")
    hasFileSep <- stringr::str_detect(inputSheet, .Platform$file.sep)

    if (readFlag && isSamSheet) {
        vals2find <- utils::read.csv(inputSheet, skip = 19)[, c(6, 7, 9)]
        return(as.data.frame(vals2find[!grepl("H20|SERACARE|HAPMAP", vals2find[, 2]),]))
    }

    if (readFlag && !isSamSheet) {
        vals2find <- read.csv(inputSheet)
        if (ncol(vals2find) > 1) {
            vals2find <- unlist(lapply(vals2find, identity))
        }
        return(as.data.frame(unique(vals2find[vals2find != ""])))
    }

    if (!readFlag && hasFileSep) {
        return(parseWorksheet(inputSheet))
    }

    return(parseWorksheet(getFilePath(inputSheet)))
}


genQuery <- function(dbCol,vals2find) {
    currCol <- vals2find[, dbCol]
    toKeep <- which(currCol != 0 & !is.na(currCol) & currCol != "")
    q1 <- currCol[toKeep]
    stopifnot("Test Number" %in% colnames(vals2find))
    names(q1) <- vals2find$`Test Number`[toKeep]
    return(q1)
}


queryCases <- function(vals2find, db) {
    queryList <- unlist(lapply(1:ncol(vals2find), function(i) {genQuery(i, vals2find)}), use.names = TRUE)
    theTScases <- queryList[stringr::str_detect(queryList, "TS|TB|TC")]
    theTScases <- sapply(theTScases, function(x) paste(stringr::str_split_fixed(x, "-", 3)[1, 1:2], collapse = "-"))
    queryList <- c(queryList, theTScases)
    queryList <- queryList[!duplicated(queryList)]
    methQuery <- searchDb(queryList, db)
    return(unique(methQuery))
}


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


addOutputLinks <- function(output) {
    if (ncol(output) == 0) {
        return(output)
    }

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


GetVolumePaths <- function(methData) {
    smb_path <- "smb://shares-cifs.nyumc.org/apps/acc_pathology"
    checkPaths <- stringr::str_replace_all(methData$`Report Path`, smb_path, "/Volumes")
    checkPaths <- checkPaths[checkPaths != "" & !is.na(checkPaths)]
    checkPaths <- checkPaths[stringr::str_detect(checkPaths, "MGDM")]
    return(checkPaths)
}

CheckMethPaths <- function(methData) {
    for (i in 1:length(methData$`Report Path`)) {
        currPath <- methData$`Report Path`[i]
        currSplit <- stringr::str_split_fixed(currPath, "/", 11)[1, ]
        if (stringr::str_detect(currSplit[10], "MGDM") == F) {
            next
        }
        runYear <- stringr::str_split_fixed(currSplit[10],"-", 2)[1,1]
        runYear <- paste0("20", runYear)
        currSplit[9] <- runYear
        newPath <- paste(currSplit, collapse = "/")
        methData[i, "Report Path"] <- newPath
    }
    checkPaths <- GetVolumePaths(methData)
    anyPathsFalse <- file.exists(checkPaths) == F
    if (any(anyPathsFalse)) {
        message("Fixing broken file paths...")
        toReplace <- basename(checkPaths[anyPathsFalse])
        mainDirs <- dirname(checkPaths[anyPathsFalse])
        if (!dir.exists(mainDirs)) {
            correct_dir <- dir(
                path = dirname(mainDirs),
                pattern = basename(mainDirs),
                full.names = T
            )
            mainDirs <- correct_dir
        }
        for (missing in toReplace) {
            patt <- stringr::str_split_fixed(missing, ".html", 2)[1,1]
            file_found <- dir(path = mainDirs, pattern = patt, full.names = T)
            if (length(file_found) > 0) {
                toSwap <- which(grepl(missing, methData$`Report Path`))
                newPath <- stringr::str_replace(methData$`Report Path`[toSwap], missing, basename(file_found))
                message("Updating file path:\n", newPath)
                methData$`Report Path`[toSwap] <- newPath
            }
        }
        checkPaths <- GetVolumePaths(methData)
        anyPathMissed <- file.exists(checkPaths) == F
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


addExcelLink <- function(output, fiLn, wb, runId) {
    x <- c(output$'Report Link'[fiLn])
    names(x) <- paste0(output$record_id[fiLn], ".html")
    class(x) <- "hyperlink"
    colNum <- which(colnames(output) == "Report Link")
    openxlsx::writeData(wb, sheet = runId, x = x, startCol = colNum, startRow = fiLn + 1)
}


createXlFile <- function(runId, output) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, runId)
    output <- CheckMethPaths(methData = output)
    openxlsx::writeData(wb, sheet = runId, x = output)
    for (fiLn in 1:length(output$'Report Link')) {
        if (output$'Report Link'[fiLn]!='') {
            addExcelLink(output, fiLn, wb, runId)
        }
    }
    outFi <- file.path(fs::path_home(),"Desktop",
                       paste0(runId,"_MethylMatch.xlsx"))
    openxlsx::saveWorkbook(wb, outFi, overwrite = T)
    return(outFi)
}


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


grab_run_id <- function(readFlag, inputSheet) {
    runId <-
        ifelse(readFlag == T, substr(inputSheet, 1, nchar(inputSheet) - 4), inputSheet)
    if (stringr::str_detect(inputSheet, ".xls")) {
        runId <-
            substr(basename(inputSheet), 1, nchar(basename(inputSheet)) - 5)
    }
    return(runId)
}


getOuputData <- function(token, flds, inputSheet, readFlag) {
    apiUrl = "https://redcap.nyumc.org/apps/redcap/api/"
    rcon <- redcapAPI::redcapConnection(apiUrl, token)

    vals2find <- getCaseValues(inputSheet, readFlag)
    if (class(vals2find) != "data.frame") {
        vals2find <- as.data.frame(vals2find)
    }

    db <- grabAllRecords(flds, rcon)

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


source_pkg_vers <- function() {
    minfiVers <- as.character(utils::packageVersion("minfi"))
    if (minfiVers != "1.43.1") {
        minfi_install()
    }

    v2_manifest <- "IlluminaHumanMethylationEPICv2manifest"
    epicVers <- as.character(utils::packageVersion(v2_manifest))
    if (epicVers != "0.1.0") {
        minfi_install()
    }

    v2Pkg_needed <- !requireNamespace("mnp.v12epicv2", quietly = T)
    v2Con_needed <- !requireNamespace("conumee2.0", quietly = T)

    if (v2Pkg_needed | v2Con_needed) {
        script_path <- "/Volumes/CBioinformatics/Methylation/Rscripts"
        v2_script <- "install_epic_v2_classifier.R"
        epicV2script <- file.path(script_path, v2_script)
        source(epicV2script)
    }

    stopifnot(library("conumee2.0", logical.return = T))
    stopifnot(library("minfi", logical.return = T))
    stopifnot(library("IlluminaHumanMethylationEPICv2manifest", logical.return = T))
    stopifnot(library("mnp.v12epicv2", logical.return = T))
}


sourceFuns2 <- function(workingPath = getwd()) {
    mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main"
    script.list <- c("R/SetRunParams.R", "R/CopyInputs.R")
    scripts <- file.path(mainHub, script.list)
    invisible(lapply(scripts, function(i) {
        suppressPackageStartupMessages(devtools::source_url(i))
    }))
    gb$setDirectory(workingPath)
    source_pkg_vers()
    return(gb$defineParams())
}


msg_rd_num <- function(rds, token) {
    message("\nRD-numbers with idats:\n", paste(rds, collapse = "\n"))
    assign("rds", rds)
    message("--------", crayon::bgMagenta("Starting CNV PNG Creation"), "--------")
    ApiToken <- token
    assign("ApiToken", ApiToken)
}


get_pact_rds <- function(rd_numbers, token) {
    result <- gb$search.redcap(rd_numbers, token)
    result <- result[!is.na(result$barcode_and_row_column),]
    samplesheet_ID = as.data.frame(stringr::str_split_fixed(result[,"barcode_and_row_column"],"_",2))
    gb$writeFromRedcap(result, samplesheet_ID)  # writes API export as minfi dataframe sheet
    gb$get.idats()  # copies idat files from return to current directory
}


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


TryCatchCNV <- function(mySentrix, sam, asPNG) {
    sampleName <- mySentrix[sam, 1]
    sentrix <- mySentrix[sam, "SentrixID_Pos"]
    message("\nGetting RGset for ", sentrix, "\n")
    RGsetEpic <- gb$grabRGset(getwd(), sentrix)
    tryCatch(
        expr = {gb$gen.cnv.png2(RGsetEpic, sampleName, asPNG)},
        error = function(e) {
            stop(crayon::bgRed(sprintf(
                "An error occured with %s png creation:", sampleName)),
                "\n", e)
        }
    )
}


loopCNV <- function(mySentrix, asPNG) {
    for (sam in rownames(mySentrix)) {
        sampleName <- mySentrix[sam, 1]
        fn = file.path(fs::path_home(), "Desktop", paste0(sampleName, "_cnv.png"))
        if (file.exists(fn)) {
            message("\nFile already exists, skipping:", fn, "\n")
        } else{
            TryCatchCNV(mySentrix, sam, asPNG)
        }
    }
}


CheckIfPngExists <- function(rds, outFolder = NULL) {
    if (is.null(outFolder)) {
        outFolder = "/Volumes/molecular/Molecular/MethylationClassifier/CNV_PNG"
    }
    outFiles <- file.path(outFolder, paste0(rds, "_cnv.png"))
    finished <- file.exists(outFiles)
    if (any(finished)) {
        message(crayon::bgGreen("The following are completed and will be skipped:"))
        message(paste(capture.output(outFiles[finished]), collapse = '\n'))

        rds <- rds[!finished]

        if (length(rds) > 0) {
            message(crayon::bgGreen("The following CNV PNG will be generated:"))
            message(paste(capture.output(rds), collapse = '\n'))
        }
    }
    return(rds)
}


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


generate_new_cnv <- function(targets) {
    has_idat <- targets[, "SentrixID_Pos"] %like% "_R0"
    targets <- targets[has_idat, ]
    if (!requireNamespace("conumee2.0", quietly = TRUE)) {
        devtools::source_url("https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/refs/heads/main/R/all_installer.R")
        source("/Volumes/CBioinformatics/Methylation/Rscripts/install_epic_v2_classifier.R")
    }
    if (!requireNamespace("mnp.v12epicv2", quietly = TRUE)) {
        devtools::source_url("https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/refs/heads/main/R/all_installer.R")
        source("/Volumes/CBioinformatics/Methylation/Rscripts/install_epic_v2_classifier.R")
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


GetSampleList <- function(rds, sampleSheet="samplesheet.csv") {
    targets <- read.csv(sampleSheet)
    toDrop <- targets$Sample_Name %in% rds
    targets <- targets[toDrop,]
    rownames(targets) <- 1:nrow(targets)
    return(targets)
}


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


copy_output_png <- function(outFolder = NULL) {
    library("dplyr")
    if (is.null(outFolder)) {
        outFolder = "/Volumes/molecular/Molecular/MethylationClassifier/CNV_PNG"
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
            "No CNV png images to generate. Check the output directory."))
    }
}

check_pkg_install()
check_REDCap_vers() # Check REDCap API version
checkMounts()
output <- getOuputData(token, flds, inputSheet, readFlag)

# CNV PNG Creation -------------------------------------
if (ncol(output) > 0) {
    queue_cnv_maker(output, token)
}
