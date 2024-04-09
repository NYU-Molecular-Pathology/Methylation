#!/usr/bin/env Rscript
## ---------------------------
## Script name: PactMethMatch.R
## Purpose: search REDCap for PACT samples with methylation & generate cnv PNG
## Date Created: September 2, 2021
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

library("base"); gb <- globalenv(); assign("gb", gb)

# Input Arguments -------------------------------------------------------------
args <- commandArgs(TRUE)
token <- args[1]
inputSheet <- args[2]

# Validate arguments
stopifnot(!is.na(token), !is.na(inputSheet))
readFlag <- grepl("\\.csv$", inputSheet)

# REDCap Fields  --------------------------------------------------------------
flds <- c("record_id", "b_number", "tm_number", "accession_number", "block",
          "diagnosis", "organ", "tissue_comments", "run_number", "nyu_mrn",
          "qc_passed")

# Message Inputs --------------------------------------------------------------
message("\n================ Parameters input ================\n")
message("token: ", token, "\ninputSheet: ", inputSheet, "\n")


if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", ask = F, dependencies = T)
}


fix_compiler_flags <- function(){
    # Check if brew installed
    if (is.na(Sys.which("brew")["brew"][[1]])) {
        message("Homebrew is not installed. Installing Homebrew...")
        cmd <-
            '/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"'
        system(cmd, wait = TRUE)
    }
    # Check if LLVM installed
    llvm_installed <- file.exists("/usr/local/opt/llvm/bin/clang")
    if (!llvm_installed) {
        message("LLVM is not installed. Installing LLVM via Homebrew...")
        system("brew install llvm", intern = TRUE, wait = TRUE)
    }
    Sys.setenv(CC = "")
    Sys.setenv(CFLAGS = "")
    Sys.setenv(CXX = "")
    Sys.setenv(CXXFLAGS = "")
    Sys.setenv(FC = "")
    Sys.setenv(FFLAGS = "")
    Sys.setenv(LDFLAGS = "")
    Sys.setenv(CPPFLAGS = "")
    Sys.setenv(SHLIB_CXXLD = "")
    Sys.setenv(SHLIB_LDFLAGS = "")
    Sys.setenv(OBJC = "")
    Sys.setenv(PATH = paste("/usr/local/opt/llvm/bin", Sys.getenv("PATH"), sep = ":"))
    Sys.setenv(CC = "/usr/local/opt/llvm/bin/clang")
    Sys.setenv(CXX = "/usr/local/opt/llvm/bin/clang++")
    Sys.setenv(CXX11 = "/usr/local/opt/llvm/bin/clang++")
    Sys.setenv(CXX14 = "/usr/local/opt/llvm/bin/clang++")
    Sys.setenv(CXX17 = "/usr/local/opt/llvm/bin/clang++")
    Sys.setenv(CXX1X = "/usr/local/opt/llvm/bin/clang++")
    Sys.setenv(OBJC = "/usr/local/opt/llvm/bin/clang")
    Sys.setenv(LDFLAGS = "-L/usr/local/opt/llvm/lib")
    Sys.setenv(CPPFLAGS = "-I/usr/local/opt/llvm/include")
}


# Function to load and install necessary packages -----------------------------
check_pkg_install <- function() {
    params <- list(
        dependencies = TRUE,
        ask = FALSE,
        update = "never",
        quiet = TRUE,
        repos = 'http://cran.us.r-project.org',
        Ncpus = 4
    )

    fix_compiler_flags()

    pkgs <- c("data.table", "openxlsx", "jsonlite", "RCurl",
              "readxl", "stringr", "tidyverse", "crayon", "tinytex",
              "systemfonts", "remotes")

    needed_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]

    if (length(needed_pkgs) > 0) {
        do.call(install.packages, c(list(pkgs = needed_pkgs), params))
    }

    sapply(pkgs, library, character.only = T, logical.return = T, quietly = T)
}


# Check REDCap Version --------------------------------------------------------
check_REDCap_vers <- function(min_version = "2.7.4") {
    if (!requireNamespace("redcapAPI", quietly = TRUE)) {
        devtools::install_github('nutterb/redcapAPI', dependencies = TRUE,
                                 upgrade = "always", ask = F, type = "source")
    }
    current_vers <- as.character(utils::packageVersion("redcapAPI"))
    if (utils::compareVersion(current_vers, min_version) < 0) {
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


message_matched <- function(item, dbInfo, ngsNum, i){
    message(sprintf("Match found for '%s' (%s) for %s in: \"%s\" column", item, dbInfo, ngsNum, i))
}

searchDb <- function(queryList, db){
    res <- data.frame()
    for (idx in 1:length(queryList)) {
        item <- queryList[idx]
        ngsNum = paste(names(item)[1])
        for (i in colnames(db)) {
            ngsMatch <- which(grepl(item, db[, i]))
            if(length(ngsMatch) > 0) {
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
    names(q1) <- vals2find$`Test Number`[toKeep]
    return(q1)
}


queryCases <- function(vals2find, db) {
    queryList <- unlist(lapply(1:ncol(vals2find), function(i) {genQuery(i, vals2find)}), use.names = TRUE)
    theTScases <- queryList[stringr::str_detect(queryList, "TS|TB|TC")]
    theTScases <- sapply(theTScases, function(x) paste(stringr::str_split_fixed(x, "-", 3)[1, 1:2], collapse = "-"))
    queryList <- unique(c(queryList, theTScases))
    methQuery <- searchDb(queryList, db)
    return(unique(methQuery))
}


get_year_paths <- function(output) {
    yearSplit <- stringr::str_split_fixed(output$run_number, "-", 2)[, 1]
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
    smb_path <- "smb://shares-cifs.nyumc.org/apps/acc_pathology"
    checkPaths <- stringr::str_replace_all(methData$`Report Path`, smb_path, "/Volumes")
    checkPaths <- checkPaths[checkPaths != "" & !is.na(checkPaths)]
    checkPaths <- checkPaths[stringr::str_detect(checkPaths, "MGDM")]
    anyPathsFalse <- file.exists(checkPaths) == F
    if (any(anyPathsFalse)) {
        message(crayon::bgRed("Some paths to html reports need editing in MethylMatch.xlsx sheet!"))
        message("Fix the following paths in worksheet 'Report Path' column that do not exist:\n\n")
        message(paste(checkPaths[anyPathsFalse], collapse = "\n"))
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
    openxlsx::writeData(wb, sheet = runId, x = output)
    output <- CheckMethPaths(methData = output)
    for (fiLn in 1:length(output$'Report Link')) {
        if (output$'Report Link'[fiLn]!='') {
            addExcelLink(output, fiLn, wb, runId)
        }
    }
    outFi <-file.path(fs::path_home(),"Desktop", paste0(runId,"_MethylMatch.xlsx"))
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


sourceFuns2 <- function(workingPath = NULL) {
    mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/"
    script.list <- c(
        "R/SetRunParams.R",
        "R/CopyInputs.R"
    )
    if (is.null(workingPath)) {workingPath = getwd()}
    scripts <- paste0(mainHub, script.list)
    invisible(lapply(scripts, function(i) {suppressPackageStartupMessages(devtools::source_url(i))}))
    gb$setDirectory(workingPath)
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
    gb$writeFromRedcap(result, samplesheet_ID) # writes API export as minfi dataframe sheet
    gb$get.idats()  # copies idat files from return to current directory
}


msgCreated <- function(mySentrix) {
    pngFiles <- paste0(file.path(fs::path_home(),"Desktop",mySentrix[, 1]),"_cnv.png")
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
            message(crayon::bgRed(sprintf("An error occured with %s png creation:", sampleName)), "\n", e)
            message(crayon::bgGreen("Trying next sample"))
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
    if(is.null(outFolder)){
        outFolder = "/Volumes/molecular/Molecular/MethylationClassifier/CNV_PNG"
    }
    outFiles <- file.path(outFolder, paste0(rds, "_cnv.png"))
    finished <- file.exists(outFiles)
    if (any(finished)) {
        message(crayon::bgGreen("The following samples are completed and will be skipped:"),
                "\n",
                paste(capture.output(outFiles[finished]), collapse = '\n')
        )
        rds <- rds[!finished]
        return(rds)
    }
    return(rds)
}


source_pkg_vers <- function(RGsetEpic){
    if (RGsetEpic@annotation[['array']] == "IlluminaHumanMethylation450k") {
        library("IlluminaHumanMethylation450kmanifest")
        library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    }
    if (RGsetEpic@annotation[['array']] == "IlluminaHumanMethylationEPIC") {
        library("IlluminaHumanMethylationEPICmanifest")
        library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
    }
    if (RGsetEpic@annotation[['array']] == "IlluminaHumanMethylationEPICv2") {
        epicV2script <- "/Volumes/CBioinformatics/Methylation/Rscripts/install_epic_v2_classifier.R"
        message("Installing package from source:\n", epicV2script)
        source(epicV2script)
        library("conumee2.0")
        library("IlluminaHumanMethylationEPICv2manifest")
        library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
        library("mnp.v12epicv2")
    }
}


get_new_pact_cnv <- function(samplename_data, sentrix.ids, i, idatPath = NULL, chrNum = NULL, doXY = F) {
    if (is.null(idatPath)) {
        idatPath <- getwd()
    }
    samName <- samplename_data[i]
    sampleEpic <- sentrix.ids[i]
    sam_out_png <- file.path(fs::path_home(), "Desktop", paste0(samName, "_cnv.png"))
    pathEpic <- file.path(idatPath, sampleEpic)
    RGsetEpic <- minfi::read.metharray(pathEpic, verbose = T, force = T)
    if (i == 1){
        source_pkg_vers(RGsetEpic)
    }
    MsetEpic <- minfi::preprocessIllumina(RGsetEpic, bg.correct = T, normalize = "controls")
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


loop_save_cnvs <- function(myDt) {
    samplename_data <- as.character(myDt[, 1])
    sentrix.ids <- as.character(myDt$SentrixID_Pos)
    for (i in 1:length(sentrix.ids)) {
        get_new_pact_cnv(samplename_data, sentrix.ids, i)
    }
}


copy_output_png <- function(outFolder = NULL) {
    if (is.null(outFolder)) {
        outFolder = "/Volumes/molecular/Molecular/MethylationClassifier/CNV_PNG"
    }
    desk <- file.path(fs::path_home(), "Desktop")
    the.cnvs <- dir(desk, "_cnv.png", full.names = T) %>% file.info() %>%
        rownames_to_column() %>% filter(as.Date(ctime) == Sys.Date()) %>% pull(rowname)

    if (length(the.cnvs) > 0) {
        savePath <- file.path(outFolder, basename(the.cnvs))
        message("\nCopying png files to Molecular folder:\n", outFolder, "\n")
        message(paste(capture.output(the.cnvs), collapse = '\n'))

        try(fs::file_copy(path = the.cnvs, new_path = savePath), T)

        if (any(!file.exists(savePath))) {
            message("The following failed to copy from the desktop:\n")
            print(basename(savePath[!file.exists(savePath)]))
        }

    } else{
        message("No CNV files found on Desktop to copy")
    }
}


generate_new_cnv <- function(myDt, asPNG = T) {
    mySentrix <- myDt[myDt[, "SentrixID_Pos"] %like% "_R0", ]

    if (nrow(mySentrix) > 0) {
        loop_save_cnvs(myDt)
    } else{
        message("The RD-number(s) do not have idat files in REDCap:/n")
        print(myDt)
    }
    msgCreated(mySentrix)
    while (!is.null(dev.list())) {dev.off()}
}


GetSampleList <- function(rds, sampleSheet="samplesheet.csv") {
    myDt <- read.csv(sampleSheet)
    toDrop <- myDt$Sample_Name %in% rds
    myDt <- myDt[toDrop,]
    rownames(myDt) <- 1:nrow(myDt)
    return(myDt)
}


try_cnv_make <- function(rds, token) {

    msg_rd_num(rds, token)
    sourceFuns2()
    get_pact_rds(rds, token)
    myDt <- GetSampleList(rds)

    tryCatch(
        expr = {
            generate_new_cnv(myDt)
        },
        error = function(e) {
            message("The following error occured:\n", e)
            message("\nTry checking the troubleshooting section on GitHub:\n")
            message("https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/PACT_scripts/README.md\n")
        }
    )
}


queue_cnv_maker <- function(output, token) {
    rds <- output$record_id[output$report_complete == "YES"]

    if (all(!is.null(rds)) == F | all(!is.na(rds)) == F | length(rds) == 0) {
        return(message(crayon::bgGreen("The PACT run has no cases with methylation.")))
    }

    rds <- CheckIfPngExists(rds)

    if (length(rds) > 0) {

        try_cnv_make(rds, token)
        try(copy_output_png(), silent = T)

    } else{
        message(crayon::bgGreen("No CNV png images to generate. Check the output directory."))
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
