#!/usr/bin/env Rscript
## Script name: pipelineHelper.R
## Purpose: source of global scripts to help execute Methylation Pipeline
## Date Created: August 5, 2021
## Version: 1.0.1
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025

gb <- globalenv(); assign("gb", gb)
dsh = "-----------"
bky <- function(...) {crayon::black$bgYellow$bold(paste(...))}
bkRed <- function(...) {crayon::bgRed$bold$white(paste(...))}
bkGrn <- function(...) {crayon::black$bgGreen$bold(paste(...))}
bkBlu <- function(...) {crayon::bgBlue$bold$white(paste(...))}

# Global Variables ----------------------------------
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
reportMd <- file.path(fs::path_home(),"report.Rmd")
QC_file <- file.path(fs::path_home(),"Methyl_QC.Rmd")
pipeLnk <-
    "https://github.com/NYU-Molecular-Pathology/Methylation/edit/main/pipelineHelper.R"
gb$UPLOAD_LOG_TSV <- "redcap_failed_file_uploads.tsv"

cbioLn <- switch(
    Sys.info()[['sysname']],
    "Darwin" = "/Volumes/CBioinformatics/Methylation/classifiers",
    "Linux" = "/gpfs/data/molecpathlab/production/Methylation/classifiers"
)

# List of three mount paths needed to run the pipleine
critialMnts <- c(
    "/Volumes/CBioinformatics/Methylation",
    "/Volumes/molecular/",
    "/Volumes/snudem01labspace/idats"
)
researchMount = "smb://research-cifs.nyumc.org/Research/"
molecularDrive = "smb://shares-cifs.nyumc.org/apps"

# Message Functions and Variables for Debugging -------------------------------
msgFunName <- function(pthLnk, funNam) {
    message("Executing function: ", bky(funNam), " from RScript in:", pthLnk)
}

msgParams <- function(...) {
    args <- list(...)
    names <- sapply(substitute(list(...))[-1], deparse)

    message(bkGrn("Params passed:"))
    for (i in seq_along(args)) {
        message(names[i], ": ", capture.output(print(args[[i]])), sep = "")
    }
}

# Helper function to return the index of priority selected samples first ------
reOrderRun <- function(selectRDs, sh="samplesheet.csv") {
    msgFunName(pipeLnk,"reOrderRun"); msgParams(selectRDs, sh)

    if (is.null(selectRDs)) {return(NULL)}
    allRd <- as.data.frame(read.csv(sh))
    runFirst <- which(allRd[,1] %in% selectRDs)
    runAfter <- which(!(allRd[,1] %in% selectRDs))
    return(c(runFirst,runAfter))
}


SetKnitProgress <- function() {
    library("knitr")
    library("progressr")
    library("cli")
    knitr::opts_knit$set(verbose = TRUE)
    options(knitr.package.verbose = FALSE)
    progressr::handlers(global = TRUE)
    progressr::handlers("cli")
    options(rmarkdown.render.message = FALSE)
    options(
        knitr.progress.fun = function(total, labels) {
            p <- progressr::progressor(
                steps = total,
                along = labels,
                trace = TRUE,
                auto_finish = FALSE,
                on_exit = FALSE
            )
            list(
                update = function(i) {
                    p(message = sprintf('Current Chunk: %s......', labels[i]),
                      class = "sticky")
                    newLab <- ifelse(labels[i] == "", "In-line Code", labels[i])
                    message("Chunk #", i, ": ", newLab, " - ",
                            round((i / total) * 100), "% complete")
                },
                done = function() p(type = 'finish')
            )
        }
    )
}


CopyRmdFile <- function(runID, rmdFile) {
    msgFunName(pipeLnk, "CopyRmdFile")
    message("runID: ", runID, " rmdFile: ",rmdFile)
    if (!file.exists(rmdFile)) {
        message(bkRed("rmdFile.rmd not found:"), "\n", rmdFile)
    }
    qcFileName = paste0(runID, "_QC.html") # output file name
    if (stringr::str_detect(rmdFile, pattern = "QC") == T) {
        if (file.exists(file.path(getwd(), qcFileName))) {
            warning(qcFileName, " Already Exists! Detete file to output new QC")
            return(NULL)
        }
    }else{
        if (file.exists(file.path(getwd(), rmdFile))) {
            warning(rmdFile, " Already Exists! Detete file to output new rmd")
            return(NULL)
        }
    }
    file.copy(
        from = rmdFile,
        to = file.path(getwd(), basename(rmdFile)),
        copy.mode = F,
        overwrite = F
    )
    return(qcFileName)
}


getRGset <- function(runPath, sentrix) {
    msgFunName(pipeLnk,"getRGset")
    barcode = stringr::str_split_fixed(sentrix, "_", 2)[1]
    RGsetEpic <- minfi::read.metharray(file.path(runPath, sentrix),
                                       verbose = T, force = T)
    aEpic = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19")
    a450k = c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19")
    arrayAnno <- RGsetEpic@annotation[['array']]
    if (arrayAnno == "IlluminaHumanMethylationEPICv2") {
        is_validation <- T
        is_validation <<- T
        return(RGsetEpic)
    }
    if (arrayAnno == "IlluminaHumanMethylationEPIC") {
        return(RGsetEpic)
    }
    if (arrayAnno == "IlluminaHumanMethylation450k") {
        RGsetEpic@annotation = a450k
        return(RGsetEpic)
    }
    if (as.numeric(barcode) <= as.numeric("204220033000")) {
        RGsetEpic@annotation = a450k
    } else{
        RGsetEpic@annotation = aEpic
    }
    return(RGsetEpic)
}


CheckSampleQCmetrics <- function(runID) {
    msgFunName(pipeLnk, "CheckSampleQCmetrics")
    qcValsFile <- file.path(getwd(), paste(runID, "qc_data.csv", sep = "_"))
    qc_cols <- c(
        "RD.number",
        "Log2sqrt.M.U.",
        "log2sqrt.R.G.",
        "BS_log2sqrt.R.G.",
        "log2sqrt.H.L.",
        "Pvalue"
    )
    qcVals <- as.data.frame(read.csv(qcValsFile)[qc_cols])
    qcVals <- gb$NameControl(qcVals, runID)
    qcVals[is.na(qcVals)] <- 0
    qcVals$Passed_SI <- as.integer(qcVals$Log2sqrt.M.U. <= 9.0)
    qcVals$Passed_BP <- as.integer(qcVals$log2sqrt.R.G. <= 11.0)
    qcVals$Passed_BS <- as.integer(qcVals$BS_log2sqrt.R.G. <= 10.0)
    qcVals$Passed_HC <- as.integer(qcVals$log2sqrt.H.L. <= 12.75)
    qcVals$Passed_NC <- as.integer(qcVals$Pvalue <= 0.95)
    qcVals$passed_qc <- with(qcVals, Passed_SI + Passed_BP + Passed_BS + Passed_HC + Passed_NC)
    qcVals$passed_qc <- ifelse(qcVals$passed_qc < 2, "yes", "no")
    final_qc <- data.frame(record_id = qcVals$RD.number,
                           qc_passed = qcVals$passed_qc)
    currVals <- final_qc$passed_qc
    if (any(is.na(currVals))) {
        final_qc$passed_qc[is.na(currVals)] <- "no"
    }
    return(final_qc)
}

# QC REPORT maker: knits the QC RMD file --------------------------------------
generateQCreport <- function(runID = NULL) {
    msgFunName(pipeLnk, "generateQCreport")
    if (is.null(runID)) {
        runID <- paste0(basename(getwd()))
    }
    qcFileName <- CopyRmdFile(runID, gb$QC_file)
    if (!is.null(qcFileName)) {
        outQCpath <- file.path(getwd(), qcFileName)
        rmdToKnit <- dir(getwd(), "*QC.Rmd", full.names = T)[1]
        deskRunDir <- file.path(fs::path_home(), "Desktop", runID)
        message("Knitting RMD file:\n", rmdToKnit, "\nWith params:")
        message(paste0(capture.output(
            data.frame(
                runID = runID,
                baseDir = deskRunDir,
                knitDir = getwd()
            )
        ), collapse = "\n"))
        rmarkdown::render(
            rmdToKnit,
            output_file = outQCpath,
            quiet = TRUE,
            params = list(runID = runID,
                          baseDir = deskRunDir,
                          knitDir = getwd())
        )
        qcCache <- stringr::str_replace_all(string = rmdToKnit, ".Rmd", "_cache")
        unlink(qcCache, recursive = T) #clear cache
        gb$uploadToRedcap(outQCpath, F)
    } else{
        message("Skipping QC render: QC file already exists in the directory!")
        message(bkRed("Delete QC rmd and QC.html to generate new output."))
    }
}


# Sends an email notification that the run is complete from redcap admin ------
launchEmailNotify <- function(runID) {
    msgFunName(pipeLnk, "launchEmailNotify")
    
    rcon <- redcapAPI::redcapConnection(gb$apiLink, gb$ApiToken)
    
    isMC = sjmisc::str_contains(runID, "MGDM") | sjmisc::str_contains(runID, "MC")
    is_validation <- sjmisc::str_contains(runID, "VAL")

    if (is_validation) {
        return(message("Run is Validation. No email will be generated"))
    }
    if (sjmisc::str_contains(runID, "EPICV1")) {
        return(message("Run is Validation. No email will be generated"))
    }

    com <- ifelse(isMC == T, "sample_qc", "sample_research")
    record_df = data.frame(record_id = paste0(runID, "_QC"), comments = com)
    record_data <- redcapAPI::castForImport(
        record_df, rcon, fields = c("record_id", "comments")
    )
    
    res <- redcapAPI::importRecords(rcon,
                             record_data,
                             returnContent = "ids",
                             overwriteBehavior = "overwrite")
    message(bkBlu("Check email to confirm notification created for", res))
}


#' Check for and Create a REDCap Record if It Does Not Exist
#'
#' This function checks whether a record with a given record_name exists in the REDCap database.
#' If the record is not found, it creates a new record with that record_name.
#'
#' @param record_name A character string representing the REDCap record identifier.
#' @param rcon A valid REDCap connection object created using redcapConnection().
#'
#' @return The result of the importRecords call if a new record is created,
#'         or a message indicating that the record already exists.
#'
#' @import redcapAPI
CreateRedcapRecord <- function(runID = NULL, recordWord = "QC") {
    rcon <- redcapAPI::redcapConnection(gb$apiLink, gb$ApiToken)
    record_name <-  paste0(runID, "_", recordWord)
    existing_record <- redcapAPI::exportRecordsTyped(rcon,
                                                     records = record_name)
    if (nrow(existing_record) == 0) {
        new_record_df <- data.frame(
            record_id = record_name,
            run_number = runID,
            stringsAsFactors = FALSE
        )
        new_record <- redcapAPI::castForImport(
            new_record_df, rcon, fields = c("record_id", "run_number")
            )
        import_result <- redcapAPI::importRecords(rcon, new_record)
        message(dsh, "Created REDCap Record: ", new_record_df[1], dsh, "\n")
    } else {
        message(paste(record_name, "already exists in the REDCap database."))
    }
}

# Checks if the QC File will be read ------------------------------------------
checkRunOutput <- function(runID) {
    msgFunName(pipeLnk,"checkRunOutput")
    redcapFi <- paste0(runID,"_Redcap.csv")
    csvLocation <- file.path(fs::path_home(),"Desktop", runID, redcapFi)

    if (!file.exists(csvLocation)) {
        message(bkRed("File not found:")," ", csvLocation)
        message(
            "QC Report may generate without the Summary Table which needs ",
            bkGrn(redcapFi))
    } else{
        outputPath <- getwd()
        if (file.exists(file.path(getwd(), redcapFi))) {
            outputPath <- file.path(getwd(), paste0(Sys.Date(), "_", redcapFi))
        }
        message("Copying file: ", outputPath)
        fs::file_copy(path = csvLocation,
                      new_path = outputPath,
                      overwrite = T)
    }
}

# Gets rid of Desktop files if run is successful ------------------------------
tidyUpFiles <- function(runID) {
    msgFunName(pipeLnk,"tidyUpFiles")
    deskDir <- file.path(fs::path_home(), "Desktop", runID)
    backupD <- file.path(gb$methDir, "csvRedcap")
    if (!dir.exists(backupD)) dir.create(backupD)
    file.copy(
        deskDir,
        backupD,
        overwrite = T,
        recursive = T,
        copy.mode = T
    )
    unlink(deskDir, T, T)
}


msgProgress <- function(msg, i, samList) {
    startEnd <- ifelse(msg == 1, "Now Running Index", "Completed Report")
    nOfTotal <- paste(i, "of", length(samList), dsh)
    cat("\n", bky(dsh, startEnd , nOfTotal),"\n")
}


msgRunUp <- function(sampleID, run_id, senLi) {
    message("\nsamplename_data = ", bkBlu(sampleID),
            "\nrun_id = ", bkBlu(run_id),
            "\npathEpic: ", bkBlu(file.path(getwd(), senLi)))
}


msgSamSheet <- function(samSh) {
    if (length(samSh) > 1) {
        message(bkRed("Multiple samplesheets found:"))
        cat(samSh, sep = "\n")
        samSh <- samSh[stringr::str_detect(samSh, pattern = "\\$", negate = T)]
        message(bkGrn("Reading the following .xlsm in current directory:"))
    } else{message(bkGrn("Reading the following .xlsm in current directory:"))}
    cat(samSh, sep = "\n")
    stopifnot(!length(samSh) > 1)
    return(samSh)
}


getRunData <- function(single_data, fi_end = ".html") {
    if (is.null(single_data)) {
        stop((bkRed("single_data is NULL: Check the SampleSheet.csv in the directory")))
    }

    tech1 <- nTech2 <- ""
    if (ncol(single_data) >= 10) {
        tech1 <- paste0(single_data[, 9])
        nTech2 <- paste0(single_data[, 10])
    }

    dat <- data.frame(
        sampleID  = paste0(single_data[, "Sample_Name"]),
        bnumber   = paste0(single_data[, "DNA_Number"]),
        senLi     = paste0(single_data[, "SentrixID_Pos"]),
        run_id    = paste0(single_data[, "RunID"]),
        mp_number = paste0(single_data[, "MP_num"]),
        tech      = tech1,
        tech2     = nTech2,
        outFi     = paste0(single_data[, 1], fi_end)
    )

    message(paste(capture.output(dat), collapse = "\n"))
    return(dat)
}


NameControl <- function(data, runId) {
    library("data.table")
    cntrl <- which(stringr::str_detect(
        data[, 1], pattern = stringr::regex('control', ignore_case = T)))
    if (length(cntrl) >= 1) {
        control_sams <- data[cntrl, 1]
        if (all(stringr::str_detect(control_sams, runId))) {
            return(data)
        }
        notNamed <- !stringr::str_detect(string = control_sams, pattern = runId)
        newCntrls <- paste(runId, control_sams[notNamed], sep = "_")
        control_sams[notNamed] <- newCntrls
        if (length(control_sams) > 1) {
            control_sams <- make.unique(control_sams, sep = "_")
        }
        data[cntrl, 1] <- control_sams
    } else{
        warning('No word "control" in RD-number found in samplesheet')
    }
    return(data)
}


ReadSamSheet <- function(samList) {
    msgFunName(pipeLnk, "ReadSamSheet")
    msgParams("samList")
    msgParams(samList)

    wb_path <- gb$GrabSampleSheet()
    if (is.null(wb_path)) {
        return(as.data.frame(read.csv("samplesheet.csv")))
    }

    xlSheets <- readxl::excel_sheets(wb_path)
    redSheet <- as.integer(which(grepl("REDCap", xlSheets) == T))

    message("Excel workbook sheet names:\n")
    message(paste(1:length(xlSheets), xlSheets, collapse = "\n"))

    if (length(redSheet) == 0) {
        warning('"REDCap" tab name is missing from workbook sheet names!')
        warning("Defaulting to reading sheet index #3 in:\n", wb_path)
        redSheet <- 3
    }

    message(bkGrn("Sheet index containing 'REDCap_Import':", redSheet), "\n")

    raw_data <- readxl::read_excel(wb_path, sheet = redSheet, range = "A1:M97",
                                   col_types = c("text"))

    samplesSheet <- as.data.frame(raw_data)[samList, 1:13]
    message(bkGrn("SampleSheet:"))
    print(samplesSheet)

    missing_rd <- samplesSheet$record_id == 0 | is.na(samplesSheet$record_id)

    if (any(missing_rd)) {
        warning(paste("Sample #", which(missing_rd), "is missing an RD-number!\n"))
        message(bkRed("Dropping sample #", which(missing_rd), "from SampleSheet!\n"))
        samplesSheet <- samplesSheet[!missing_rd, ]
    }

    return(samplesSheet)
}

# FUN: Parses the WetLab .xlsm sheet in the current directory -----------------
Check_sam_csv <- function(samList) {
    msgFunName(pipeLnk, "Check_sam_csv")
    msgParams("samList")
    msgParams(samList)

    require(rmarkdown)

    isMC <- sjmisc::str_contains(gb$runID, "MGDM") |
        sjmisc::str_contains(gb$runID, "MC")
    is_validation <- sjmisc::str_contains(gb$runID, "VAL")
    is_research <- grepl("MR", gb$runID)

    wksh <- ReadSamSheet(samList)
    if (is_validation) {
        reportMd <<-
            "/Volumes/CBioinformatics/Methylation/EPIC_V2_report_2.Rmd"
        CopyRmdFile(gb$runID, reportMd)
    }

    if (isMC == T & is_validation == F & is_research == F) {
        wksh <- NameControl(wksh, wksh$run_number[1])
    }

    stopifnot(!is.null(wksh))
    rownames(wksh) <- wksh[,1]
    return(wksh)
}


getRunList <- function(data, samList) {
    msgFunName(pipeLnk, "getRunList")
    all_html <- dir(pattern = ".html")
    if (length(all_html) == 0) {
        return(samList)
    }
    toRun <- unlist(lapply(samList, FUN = function(i) {
        currSam <- data[i, 1]
        if (any(stringr::str_detect(all_html, pattern = currSam))) {
            message(bky(basename(currSam), "already exists! Skipping..."),"\n")
        }else{return(i)}
    }))
    return(toRun)
}


make_knit_report <- function(dat, reportMd, params_init) {
    rmarkdown::render(
        input = reportMd, output_format = "html_document",
        output_file = file.path(getwd(), paste0(dat$outFi)),
        output_dir = getwd(), knit_root_dir = getwd(),
        clean = TRUE, quiet = FALSE,
        output_options = list(self_contained = TRUE, clean_supporting = TRUE),
        params = params_init
    )
}


handle_knit_error <- function(e, dat, params) {
    try(beepr::beep(1), silent = TRUE)
    message("Report generation failed for:", "\n", dat$outFi,
            "\nThe following error returned:\n", e)
    home_path <- Sys.getenv("HOME")
    tb <- capture.output(traceback())
    if (length(tb) > 0) {
        writeLines(tb, file.path(home_path, "error_log.txt"))
        message(tb)
    }
    saveRDS(params, file.path(home_path, "params.rds"))
    saveRDS(gb$chunk_env, file.path(home_path, "chunk_env.rds"))
    message("Saving logs to:", home_path)
    stop("Check error_log.txt, params.rds, and chunk_env.rds for details.")
}


# Knits the output HTML reports given single row of data from samplesheet -----
do_report <- function(single_data = NULL, genCn = FALSE) {
    msgFunName(pipeLnk, "do_report")
    msgParams("data")

    dat <- getRunData(single_data)

    RGsetEpic <- getRGset(runPath = getwd(), sentrix = dat$senLi)
    reportMd <- "/Volumes/CBioinformatics/Methylation/EPIC_V2_report_2.Rmd"

    if (genCn == T) {
        gb$generate_cnv_png(RGsetEpic, dat$sampleID)
    }

    msgRunUp(dat$sampleID, dat$run_id, dat$senLi)
    #message("Knitting report: ", reportMd)

    params_init <- list(
        token = gb$ApiToken,
        rundata  = dat,
        RGsetEpic = RGsetEpic,
        knitDir = getwd()
    )

    tryCatch(
        expr = make_knit_report(dat, reportMd, params_init),
        error = function(e) {handle_knit_error(e, dat, params_init)},
        finally = message("\nRunning next sample\n")
    )
}

# FUN: Knits teh V11 version html report --------------------------------------
do_report_v11 <- function(single_data = NULL, genCn = FALSE) {
    msgFunName(pipeLnk, "do_report_v11")
    msgParams("data")

    dat <- getRunData(single_data)
    new_filename <- sub("(\\.html)$", "_v11\\1", dat$outFi)
    dat$outFi <- new_filename
    RGsetEpic <- getRGset(runPath = getwd(), sentrix = dat$senLi)
    reportMd_v11 <- "/Volumes/CBioinformatics/Methylation/report.Rmd"

    msgRunUp(dat$sampleID, dat$run_id, dat$senLi)
    message("Knitting report: ", reportMd)

    params_init <- list(
        token = gb$ApiToken,
        rundata  = dat,
        RGsetEpic = RGsetEpic,
        knitDir = getwd()
    )

    tryCatch(
        expr = make_knit_report(dat, reportMd_v11, params_init),
        error = function(e) {handle_knit_error(e, dat, params_init)},
        finally = message("\nRunning next sample\n")
    )
}

# FUN: Generates the V11 version html report ----------------------------------
get_v11_reports <- function(your_csv){
    unloadNamespace("mnp.v12epicv2")
    library("mnp.v11b6")
    require("mnp.v11b6")
    library("utils")

    data <- utils::read.csv(your_csv, strip.white = T)
    samList <- 1:length(data$Sample_Name != 0)
    toRun <- gb$getRunList(data, samList)

    for (idx in toRun) {
        message("Running ", idx, " of ", length(toRun))
        sheetRunID <- data[idx, "RunID"]
        gb$runID <- sheetRunID
        assign("runID", sheetRunID)
        single_data = data[idx, ]
        do_report_v11(single_data)
    }
}

# FUN: Iterates over each sample in the csv file to generate a report ---------
# DEBUG: data <- read.csv("samplesheet.csv", strip.white=T)
loopRender <- function(samList = NULL, data, redcapUp = T) {
    msgFunName(pipeLnk, "loopRender")
    stopifnot(!is.null(data))
    if (is.null(samList)) {
        samList <- 1:length(data$Sample_Name != 0)
    }
    workbook_data <- Check_sam_csv(samList)
    toRun <- getRunList(data, samList)
    SetKnitProgress()
    requireNamespace("mnp.v12epicv2")
    currIdx = 1
    for (sam_idx in toRun) {
        totLeft <- length(toRun) - currIdx
        message(bkGrn(dsh, totLeft, "of", length(toRun), "samples remain to be run", dsh))
        remain_percent <- 100 - round((totLeft/length(toRun)) * 100 )
        perc_comp <- paste0(remain_percent, "%")
        message(bkGrn("Run is approximately", perc_comp, "complete"))

        msgProgress(1, sam_idx, samList)
        do_report(single_data = data[sam_idx, ], gb$genCn)
        msgProgress(2, sam_idx, samList)
        if (redcapUp == T) {
            sh_Dat <- workbook_data[sam_idx, ]
            gb$importSingle(sh_Dat)
        }
        currIdx = currIdx + 1
    }
    message(bkGrn(dsh, "RUN COMPLETE", dsh))
}


RenameFailed <- function(qcVals) {
    msgFunName(pipeLnk, "RenameFailed")

    if (!is.null(qcVals)) {
        message("qcVals")
        qcVals[is.na(qcVals)] <- "no"
        print(qcVals)

        if (any(qcVals$qc_passed == "no")) {
            file.list <- dir(getwd(), pattern = ".html", full.names = T)
            toRename <- qcVals$record_id[qcVals$qc_passed == "no"]
            for (rd_num in toRename) {
                findFile <- stringr::str_detect(file.list, pattern = rd_num)
                old_name <- file.list[findFile]
                is_renamed <- stringr::str_detect(old_name, pattern = "_QC_FAILED.html")
                if (!is_renamed) {
                    new_name <- stringr::str_replace_all(
                        old_name, ".html", "_QC_FAILED.html")
                    base::file.rename(old_name, new_name)
                }
            }
        }
    }
}


CreateControlRecords <- function(cntrl, runID, control_sams) {
    if (length(cntrl) == 1) {
        CreateRedcapRecord(runID, "control")
    } else {
        for (c_sam in 1:length(cntrl)) {
            curr_sam <- control_sams[c_sam]
            curr_splt <- stringr::str_split_fixed(curr_sam, "_", 2)
            CreateRedcapRecord(curr_splt[1,1], curr_splt[1,2])
        }
    }
}


check_control_sam1 <- function(data, cntrl, runID) {
    control_sams <- data[cntrl, 1]
    isNamed <- stringr::str_detect(string = control_sams, pattern = runID)
    if (!all(isNamed)) {
        notNamed <- !stringr::str_detect(string = control_sams, pattern = runID)
        newCntrls <- paste(runID, control_sams[notNamed], sep = "_")
        control_sams[notNamed] <- newCntrls
        if (length(control_sams) > 1) {
            control_sams <- make.unique(control_sams, sep = "_")
        }
        data[cntrl, 1] <- control_sams
        CreateControlRecords(cntrl, runID, control_sams)
    } else{
        CreateControlRecords(cntrl, runID, control_sams)
    }
    
    return(data)
}


check_csv_data <- function(sheetName = "samplesheet.csv") {
    data <- utils::read.csv(sheetName, strip.white = T)
    toKeep <- !(is.na(data[, 1]) | data[, 1] == 0 | data[, 1] == "")
    if (any(!toKeep)) {
        data <- data[toKeep, , drop = FALSE]
        rownames(data) <- NULL
    }

    sheetRunID <- paste0(data$RunID[1])
    if (sheetRunID != gb$runID) {
        message("Batch ID in sheet is: ", sheetRunID)
        message("Batch ID provided is: ", gb$runID)
        stopifnot(sheetRunID == gb$runID)
    }

    runID <- paste0(data$RunID[1])
    cntrl <- which(stringr::str_detect(
        data[, 1], pattern = stringr::regex('control', ignore_case = T)))

    if (length(cntrl) >= 1) {
        data <- check_control_sam1(data, cntrl, runID)
    } else {
        warning("No control found in the samplesheet csv file!")
    }

    return(data)
}


# MAIN: Generates Html reports with samplesheet.csv for V12_EPICV2 ----
makeHtmlReports <- function(runOrder = NULL,
                            skipQC = F,
                            email = T,
                            redcapUp = T) {
    msgFunName(pipeLnk, "makeHtmlReports")

    library("data.table")
    assign("genCn", FALSE, envir = gb)

    data <- check_csv_data()
    runID <- paste0(data$RunID[1])

    reportMd <- "/Volumes/CBioinformatics/Methylation/EPIC_V2_report_2.Rmd"
    CopyRmdFile(gb$runID, reportMd)

    library('mnp.v12epicv2')
    library("pander")
    library("htmltools")

    loopRender(runOrder, data, redcapUp)
    checkRunOutput(runID)

    qcVals <- NULL
    if (skipQC == F) {
        CreateRedcapRecord(runID)
        generateQCreport(runID)
        qcVals <- CheckSampleQCmetrics(runID)
        rcon <- redcapAPI::redcapConnection(gb$apiLink, gb$ApiToken)
        qcVals <- gb$NameControl(qcVals, runID)

        redcapAPI::importRecords(rcon, qcVals, "normal", "ids",
                                 logfile = "REDCapQCimports.txt")
    }

    if (grepl("TEST", runID)) {
        redcapUp <- email <- FALSE
    }

    # if (file.exists(gb$UPLOAD_LOG_TSV)) {
    #     file.list <- read.table(gb$UPLOAD_LOG_TSV)[,1]
    #     if (length(file.list) > 0) {
    #         gb$uploadToRedcap(file.list, T)
    #     }
    # }
    file.list <- dir(getwd(), pattern = ".html", full.names = TRUE)
    gb$uploadToRedcap(file.list, F)

    if (email == T) {
        RenameFailed(qcVals)
        gb$CombineClassAndQC(
            output_fi = paste0(runID, "_qc_data.csv"),
            gb$ApiToken, runDir = NULL, runID
        )
        launchEmailNotify(runID)
    }

    try(beepr::beep(2), T)
    tidyUpFiles(runID)
}


# FUN: Checks if all the paths are accessible to the Rscript location ---------
checkMounts <- function() {
    msgFunName(pipeLnk,"checkMounts")
    failMount <- lapply(critialMnts, function(driveMount) {
        ifelse(!dir.exists(driveMount), return(T), return(F))
    })
    if (any(failMount == T)) {
        toFix <- paste(critialMnts[which(failMount == T)])
        cat("PATH does not exist, ensure network drive is mounted:",
            bkRed(toFix),"\n")
        cat("You may need to re-mount one of the following paths:",
            bkBlu(researchMount,"CBioinformatics/"),
            bkGrn(researchMount,"snudem01lab/snudem01labspace"),
            bkBlu(molecularDrive,"/acc_pathology/molecular"), sep = "\n"
        )
        stopifnot(!any(failMount == T))
    }
}

# Checks and sets the base directory for executing the pipeline ---------------
CheckBaseDir <- function(baseFolder) {
    if (is.null(baseFolder)) {
        gb$baseDir <- gb$methDir <- gb$baseFolder <-
            "/Volumes/CBioinformatics/Methylation/Clinical_Runs"
    } else{
        gb$baseDir <- gb$methDir <- gb$baseFolder <- baseFolder
    }
    isDesktop <- stringr::str_detect(baseFolder, "Desktop")
    if (length(isDesktop) != 0)
        if (is.null(baseFolder) & isDesktop == T) {
            warning("The Methylation pipeline cannot be run from the Desktop")
            message("Try setting baseFolder to '~/Documents/' instead")
            stopifnot(isDesktop == F)
        }
    return(gb$baseFolder)
}

# Sets the working folder directory -------------------------------------------
SetBaseFolder <- function(token, baseFolder, runID) {
    baseFolder <- CheckBaseDir(baseFolder)
    methylPath <- gb$setRunDir(runID = gb$runID, workFolder = baseFolder)
    message("Working directory set to:\n", crayon::bgGreen(methylPath), "\n")
    gb$methDir <- gb$workFolder <- baseFolder
    gb$setVar("workFolder", baseFolder)
    gb$setVar("ApiToken", token)
    newPathRun <- file.path(baseFolder, runID)
    pathCdCmd <- paste("cd", newPathRun)
    message("Changing directory to:", pathCdCmd)
    system(pathCdCmd)
    setwd(file.path(baseFolder, runID))
    message(bkRed("Current Path"))
    message(getwd())
    return(baseFolder)
}


# Executes functions to check and create input and output files ---------------
PrepareRun <- function(token,
                       baseFolder = NULL,
                       runID,
                       runLocal = F,
                       rdInput = F) {
    mainHub <- "https://github.com/NYU-Molecular-Pathology/Methylation/main/R"
    gb$msgFunName(file.path(mainHub, "methylExpress.R"), "PrepareRun")
    if (runLocal == F) {
        gb$checkMounts()
        gb$checkValidRun(runID)
        gb$SetBaseFolder(token, baseFolder, runID)
        gb$copyWorksheetFile(runID = runID) # copies the xlsm file to cwd
        gb$readSheetWrite(runID = runID) # generates input .csv sheet from xlsm
        gb$get.idats() # Copy idat files to current folder from shared drives
        gb$moveSampleSheet(baseFolder, runID)  # copies outputs to Desktop QC
    } else{
        baseFolder <- SetBaseFolder(token, baseFolder, runID)
        assign("baseFolder", baseFolder)
        if (rdInput == F) {
            gb$RunLocalIdats(runID, token)
        }
    }
}


# Starts run knitting html reports loop ---------------------------------------
StartRun <- function(selectRDs = NULL, emailNotify = T, redcapUp = T) {
    mainHub <- "https://github.com/NYU-Molecular-Pathology/Methylation/main/R"
    gb$msgFunName(file.path(mainHub,"methylExpress.R"),"startRun")

    runOrder <- NULL
    if (!is.null(selectRDs)) {
        runOrder <- gb$reOrderRun(selectRDs)
    }
    gb$makeHtmlReports(
        runOrder = runOrder,   # Prioritize specific RD-numbers
        skipQC = F,            # Don't skip QC generation
        email = emailNotify,   # to email after Run complete
        redcapUp = redcapUp    # Flag to import files to REDCap
    )
}


GetPriorityCases <- function(selectRDs = NULL,
                             samSheet = "samplesheet.csv",
                             kwd = "BN0") {
    msgFunName(pipeLnk, "GetPriorityCases")
    csvFi <- read.csv(file.path(getwd(), samSheet))
    BN00 <- which(stringr::str_detect(csvFi$MP_num, kwd))
    if (length(BN00) > 0) {
        if (is.null(selectRDs)) {
            selectRDs <- c(csvFi$Sample_Name[BN00])
        }else{
            selectRDs <- c(selectRDs, csvFi$Sample_Name[BN00])
        }
        message("Prioritizing these cases first:\n",
                paste(capture.output(selectRDs), collapse = " "))
    }
    return(selectRDs)
}


AssignArgs <- function(runID,
                       baseFolder = NULL,
                       token = NULL,
                       selectRDs = NULL,
                       redcapUp,
                       gb) {
    msgFunName(pipeLnk,"AssignArgs")
    stopifnot(!is.null(token))
    stopifnot(!is.null(runID))
    assign("token", token, envir = gb)
    assign("ApiToken", token, envir = gb)
    assign("redcapUp", redcapUp, envir = gb)
    assign("redcapUpload", redcapUp, envir = gb)
    gb$defineParams(methDir = baseFolder, baseDir = baseFolder,
                    ApiToken = token)
    gb$setVar("runID", runID)
    if (!is.null(selectRDs)) {
        selectRDs <- stringr::str_split(selectRDs, ",")
    }
    return(selectRDs)
}


CheckIdatsCopied <- function() {
    msgFunName(pipeLnk,"CheckIdatsCopied")
    idatFiles <- dir(getwd(), pattern = ".idat", full.names = T)
    idatSize <- paste(round(file.info(idatFiles)$size/1024.0 ^ 2, 2),"MB")
    idatFiSi <- as.data.frame(table(idatSize))

    if (nrow(idatFiSi) > 1) {
        message(crayon::bgRed("Some idat files are not the same size!"))
        message("Make sure Red_ and Grn_.idat pairs are the same file size!")
        message(paste0(capture.output(idatFiSi), collapse = "\n"))
    }
}


StartCustomRun <- function(redcapUp = T) {
    msgFunName(pipeLnk,"StartCustomRun")
    gb$makeHtmlReports(
        skipQC = T,
        email = F,
        selectSams = NULL,
        redcapUp = redcapUp
    )
}
