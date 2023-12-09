#!/usr/bin/env Rscript
## ---------------------------
## Script name: pipelineHelper.R
## Purpose: source of global scripts to help execute Methylation Pipeline
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------
gb <- globalenv(); assign("gb", gb)
dsh = "-----------"
bky <- function(...) {crayon::black$bgYellow$bold(paste(...))}
bkRed <- function(...) {crayon::bgRed$bold$white(paste(...))}
bkGrn <- function(...) {crayon::black$bgGreen$bold(paste(...))}
bkBlu <- function(...) {crayon::bgBlue$bold$white(paste(...))}

# Global Variables ----------------------------------
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
reportMd <- file.path(fs::path_home(),"report.Rmd") # From curl github download
QC_file <- file.path(fs::path_home(),"Methyl_QC.Rmd") # From curl github download
pipeLnk <- "https://github.com/NYU-Molecular-Pathology/Methylation/edit/main/pipelineHelper.R"
cbioLn <- switch (Sys.info()[['sysname']],
                  "Darwin" = "/Volumes/CBioinformatics/Methylation/classifiers",
                  "Linux" = "~/molecpathlab/production/Methylation/classifiers"
)

predictionPath <- file.path(cbioLn, "mnp.v11b6/data/rfpred.v11b6.RData")

# List of three mount paths needed to run the pipleine
critialMnts <- c(
    "/Volumes/CBioinformatics/Methylation",
    "/Volumes/molecular/",
    "/Volumes/snudem01labspace/idats"
)
researchMount = "smb://research-cifs.nyumc.org/Research/"
molecularDrive = "smb://shares-cifs.nyumc.org/apps"

# Message Functions and Variables for Debugging ----------------------------------
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

# Helper function to return the index of priority selected samples first
reOrderRun <- function(selectRDs, sh=NULL) {
    msgFunName(pipeLnk,"reOrderRun"); msgParams(selectRDs, sh)

    if (is.null(selectRDs)) {return(NULL)}
    if (is.null(sh)) {sh<-"samplesheet.csv"}
    allRd <- as.data.frame(read.csv(sh))
    runFirst <- which(allRd[,1] %in% selectRDs)
    runAfter <- which(!(allRd[,1] %in% selectRDs))
    return(c(runFirst,runAfter))
}

# Saves the methyl CNV as a png file in the cwd
generateCNVpng <- function(RGsetEpic, sampleName) {
    msgFunName(pipeLnk,"generateCNVpng")

    imgName <- paste(sampleName, "cnv.png", sep="_")
    MsetEpic <- minfi::preprocessRaw(RGsetEpic)
    png(filename=imgName,width=1820, height=1040)
    mnp.v11b6::MNPcnvplot(MsetEpic)
    dev.off(); message("File saved:\n",imgName,"\n")
}


getRGset <- function(runPath, sentrix) {
    msgFunName(pipeLnk,"getRGset")
    barcode = stringr::str_split_fixed(sentrix, "_",2)[1]
    RGsetEpic <- minfi::read.metharray(file.path(runPath, sentrix), verbose = T, force = T)

    aEpic = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19")
    a450k = c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19")
    arrayAnno <- RGsetEpic@annotation[['array']]
    if (arrayAnno == "IlluminaHumanMethylationEPICv2") {
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


CheckSampleQCmetrics <- function(runID) {
    qcValsFile <- file.path(getwd(), paste(runID, "qc_data.csv", sep = "_"))
    qc_cols <- c("RD.number", "Log2sqrt.M.U.", "log2sqrt.R.G.", "BS_log2sqrt.R.G.", "log2sqrt.H.L.", "Pvalue")
    qcVals <- as.data.frame(read.csv(qcValsFile)[qc_cols])
    qcVals <- apply(qcVals, 2, function(x) ifelse(is.na(x), 0, x))
    qcVals$Passed_SI <- as.integer(qcVals$Log2sqrt.M.U. <= 9.0)
    qcVals$Passed_BP <- as.integer(qcVals$log2sqrt.R.G. <= 11.0)
    qcVals$Passed_BS <- as.integer(qcVals$BS_log2sqrt.R.G. <= 10.0)
    qcVals$Passed_HC <- as.integer(qcVals$log2sqrt.H.L. <= 12.75)
    qcVals$Passed_NC <- as.integer(qcVals$Pvalue <= 0.95)
    qcVals$passed_qc <- with(qcVals, Passed_SI + Passed_BP + Passed_BS + Passed_HC + Passed_NC)
    qcVals$passed_qc <- ifelse(qcVals$passed_qc < 2, "yes", "no")
    final_qc <- data.frame(record_id = qcVals$RD.number, qc_passed = qcVals$passed_qc)
    currVals <- final_qc$passed_qc
    if(any(is.na(currVals))){
        final_qc$passed_qc[is.na(currVals)] <- "no"
    }
    return(final_qc)
}

# QC REPORT maker: knits the QC RMD file
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
        rmarkdown::render(
            rmdToKnit,
            output_file = outQCpath,
            quiet = TRUE,
            params = list(runID = runID, baseDir = deskRunDir, knitDir = getwd())
        )
        qcCache <- stringr::str_replace_all(string = rmdToKnit, ".Rmd", "_cache")
        unlink(qcCache, recursive = T) #clear cache
        gb$uploadToRedcap(outQCpath, F)
    }else{
        message("Skipping QC render: QC file already exists in the directory.  Delete QC rmd and QC.html to generate new output.")
    }
}

# Sends an email notification that the run is complete from redcap admin
launchEmailNotify <- function(runID) {
    msgFunName(pipeLnk,"launchEmailNotify")
    isMC = sjmisc::str_contains(runID, "MGDM") | sjmisc::str_contains(runID, "MC")
    is_validation <- sjmisc::str_contains(runID, "VAL")
    if (is_validation) {
        return(message("Run is Validation. No email will be generated"))
    }
    com <- ifelse(isMC == T, "sample_qc", "sample_research") # research or clinical notification
    record = data.frame(record_id = paste0(runID, "_QC"), comments = com)
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox = T)
    RCurl::postForm(
        gb$apiLink,
        token = gb$ApiToken,
        content = 'record',
        format = 'json',
        type = 'flat',
        data = datarecord
    )
    cat(bkBlu("Check email to confirm Email Notification Created"), sep = "\n")
}

# FUN: Creates the QC record for the current run on redcap if it does not exist
CreateRedcapRecord <- function(runID = NULL, recordWord = "QC") {
    msgFunName(pipeLnk, "CreateRedcapRecord")
    if (is.null(runID)) {
        runID <- paste0(basename(gb$workDir))
    }
    record = c(record_id = paste0(runID, "_", recordWord), run_number = runID)
    cntrl <- jsonlite::toJSON(list(as.list(record)), auto_unbox = T)
    try(RCurl::postForm(
        gb$apiLink, token = gb$ApiToken, content = 'record', format = 'json',type = 'flat',data = cntrl),T)
    message(dsh, "Created REDCap Record: ", record[1], dsh, "\n")
}

# Check if the QC File will be read
checkRunOutput <- function(runID) {
    msgFunName(pipeLnk,"checkRunOutput")
    redcapFi <- paste0(runID,"_Redcap.csv")
    csvLocation <- file.path(fs::path_home(),"Desktop", runID, redcapFi)

    if (!file.exists(csvLocation)) {
        message(bkRed("File not found:")," ", csvLocation)
        message("QC Report may generate without the Summary Table which needs ", bkGrn(redcapFi))
    } else{
        outputPath <- getwd()
        if (file.exists(file.path(getwd(), redcapFi))) {
            outputPath <- file.path(getwd(), paste0(Sys.Date(), "_", redcapFi))
        }
        message("Copying file: ", outputPath)
        fs::file_copy(path=csvLocation, new_path=outputPath, overwrite = T)
    }
}

# gets rid of desktop files if run is successful
tidyUpFiles <- function(runID) {
    msgFunName(pipeLnk,"tidyUpFiles")
    deskDir <- file.path(fs::path_home(),"Desktop",runID)
    backupD <- file.path(gb$methDir,"csvRedcap")
    if (!dir.exists(backupD)) {dir.create(backupD)}
    file.copy(deskDir, backupD, overwrite=T, recursive = T, copy.mode = T)
    unlink(deskDir, T, T)
}


msgProgress <- function(msg, i, samList) {
    startEnd <- ifelse(msg == 1, "Now Running", "Completed Report")
    nOfTotal <- paste(i, "of", length(samList), dsh)
    cat("\n", bky(dsh, startEnd , nOfTotal),"\n")
}


msgRunUp <- function(sampleID, run_id, senLi) {
    message(
        "\nsamplename_data = ",
        bkBlu(sampleID),
        "\nrun_id = ",
        bkBlu(run_id),
        "\npathEpic: ",
        bkBlu(file.path(getwd(), senLi)),
        "\n"
    )
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


getRunData <- function(data) {
    if (ncol(data)<10) {
        tech1 <- nTech2 <- ""
    }else{
        tech1 <-paste0(data[, 9])
        nTech2 <- paste0(data[, 10])
    }
    runDt <- data.frame(
        sampleID = paste0(data[, 1]),
        bnumber = paste0(data[, 2]),
        senLi = paste0(data[, 5]),
        run_id = paste0(data[, 7]),
        mp_number = paste0(data[, 8]),
        tech = tech1,
        tech2 = nTech2,
        outFi = paste0(data[, 1], ".html")
    )
    return(runDt)
}


NameControl <- function(data, runId) {
    library("data.table")
    if (any(data[, 1] %like% 'control')) {
        cntrl <- which(data[, 1] %like% 'control') #DNA_Number
        data[cntrl, 1] <- paste0(runId, "_control")
    } else{
        if (any(data[, 2] %like% 'control')) {
            cntrl <- which(data[, 2] %like% 'control') #DNA_Number
            data[cntrl, 1] <- paste0(runId, "_control")
        } else{
            warning('No word "control" in RD-number found in samplesheet')
        }
        return(data)
    }
    return(data)
}


ReadSamSheet <- function(samList) {
    msgFunName(pipeLnk, "ReadSamSheet")
    msgParams("samList")
    msgParams(samList)

    samSh <- gb$GrabSampleSheet()
    if (is.null(samSh)) {
        return(as.data.frame(read.csv("samplesheet.csv")))
    }
    xlSheets <- readxl::excel_sheets(samSh)
    redSheet <- as.integer(which(grepl("REDCap",xlSheets)==T))
    message("Excel sheet names:")
    print(xlSheets)
    if (length(redSheet)==0) {
        warning("REDCap sheet not found in Sheet names, setting sheet #3")
        redSheet <- 3
    }
    message(bkGrn("Sheet Index containing 'REDCap_Import':"), " ", redSheet,"\n")

    samSh <- readxl::read_excel(samSh, sheet=redSheet, range = "A1:M97", col_types = c("text"))
    message(bkGrn("SampleSheet:"))

    samplesSheet <- as.data.frame(samSh)[samList, 1:13]
    print(samplesSheet)
    return(samplesSheet)
}

# FUN: Parses the WetLab .xlsm sheet in the current directory
checkSamSh <- function(samList) {
    msgFunName(pipeLnk, "checkSamSh")
    msgParams("samList")
    msgParams(samList)

    require(rmarkdown)
    wksh <- ReadSamSheet(samList)
    isMC = sjmisc::str_contains(gb$runID, "MGDM")|sjmisc::str_contains(gb$runID, "MC")
    is_validation <- sjmisc::str_contains(gb$runID, "VAL")
    if (isMC == T & is_validation == F) {
        wksh <- NameControl(wksh, wksh$run_number[1])
    }
    stopifnot(!is.null(wksh))
    rownames(wksh)<- wksh[,1]
    return(wksh)
}


getRunList <- function(data, samList) {
    msgFunName(pipeLnk, "getRunList")
    toRun <- unlist(lapply(samList, FUN=function(i) {
        nOutDir = file.path(getwd(), paste0(data[i, 1], ".html"))
        if (file.exists(nOutDir)) {
            message(bky(basename(nOutDir), "already exists! Skipping..."),"\n")
        }else{return(i)}
    }))
    return(toRun)
}


make_knit_report <- function(dat, reportMd, params_init) {
    rmarkdown::render(
        input = reportMd, output_format = "html_document",
        output_file = file.path(getwd(), paste0(dat$outFi)),
        output_dir = getwd(), knit_root_dir = getwd(),
        clean = TRUE, quiet = TRUE,
        output_options = list(self_contained = TRUE, clean_supporting = TRUE),
        params = params_init
    )
}


handle_knit_error <- function(e, dat, params) {
    try(beepr::beep(1), silent = TRUE)
    message("Report generation failed for:", "\n", dat$outFi, "\nThe following error returned:\n", e)
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


do_report <- function(single_data = NULL, genCn = FALSE) {
    msgFunName(pipeLnk, "do_report")
    msgParams("data")

    if (is.null(single_data)) {
        message(bkRed("single_data is NULL, check your SampleSheet.csv"))
        stopifnot(!is.null(single_data))
    }

    dat <- getRunData(single_data)
    RGsetEpic <- getRGset(runPath = getwd(), sentrix = dat$senLi)

    if (genCn == T) {
        generate_cnv_png(RGsetEpic, dat$sampleID)
    }

    msgRunUp(dat$sampleID, dat$run_id, dat$senLi)
    message("Knitting report: ", reportMd)

    params_init <- list(
        token = gb$ApiToken,
        rundata  = dat,
        RGsetEpic = RGsetEpic,
        knitDir = getwd()
    )

    message(paste0(capture.output(dat), collapse = "\n"))
    tryCatch(
        expr = make_knit_report(dat, reportMd, params_init),
        error = function(e) {handle_knit_error(e, dat, params_init)},
        finally = message("\nRunning next sample\n")
    )
}

# FUN: Iterates over each sample in the csv file to generate a report
loopRender <- function(samList = NULL, data, redcapUp = T) {
    msgFunName(pipeLnk, "loopRender")
    # Debug: data <- read.csv("samplesheet.csv", strip.white=T)
    stopifnot(!is.null(data))

    isMC = sjmisc::str_contains(gb$runID, "MGDM")|sjmisc::str_contains(gb$runID, "MC")
    is_validation <- sjmisc::str_contains(gb$runID, "VAL")
    if (isMC == T & is_validation == F) {
        data <- NameControl(data, data$RunID[1])
    }
    if (is.null(samList)) {
        samList <- 1:length(data$Sample_Name!=0)
    }
    wksh <- checkSamSh(samList)
    toRun <- getRunList(data, samList)

    is_validation <-  sjmisc::str_contains(gb$runID, "VAL")
    if (is_validation) {
        reportMd <<- "/Volumes/CBioinformatics/Methylation/EPIC_V2_report.Rmd"
        CopyRmdFile(gb$runID, reportMd)
    }

    currIdx = 1
    for (i in toRun) {
        message(bkGrn(dsh, currIdx, "of", length(toRun), "samples are complete", dsh))
        msgProgress(1, i, samList)
        single_data = data[i, ]
        do_report(single_data, gb$genCn)
        msgProgress(2, i, samList)
        if (redcapUp == T) {
            sh_Dat <- wksh[i, ]
            gb$importSingle(sh_Dat)
        }
        currIdx = currIdx + 1
    }
    message(bkGrn(dsh, "RUN COMPLETE", dsh))
}


RenameFailed <- function(qcVals) {
    if (!is.null(qcVals)) {
      message("qcVals")
      print(qcVals)
        if (any(qcVals$qc_passed == "no")) {
            file.list <- dir(getwd(), pattern = ".html", full.names = T)
            toRename <- qcVals$record_id[qcVals$qc_passed == "no"]
            for (rd_num in toRename) {
                findFile <- stringr::str_detect(file.list, pattern = rd_num)
                old_name <- file.list[findFile]
                new_name <- stringr::str_replace_all(old_name, ".html", "_QC_FAILED.html")
                base::file.rename(old_name, new_name)
            }
        }
    }
}


#' REPORT: Generates Html reports to cwd with samplesheet.csv
#' @param runPath The location of samplesheet.csv and idats
#' @param sheetName name of samplesheet if it is not "samplesheet.csv"
#' @param selectSams vector of integer index of specific samples to run from samplesheet
#' @param genCn also generate the CNV plot as a PNG file
#' @param skipQC default is false set as true to skip QC generation
#' @param email default is true, set to false to avoid email notification
#' @param cpReport default false, set TRUE to copy the reports to the Zdrive/research output Dir
#' @param redcapUp default is true, flag will upload output html files and dataframe to redcap
makeReports.v11b6 <- function(runPath = NULL,
                              sheetName = "samplesheet.csv",
                              selectSams = NULL,
                              genCn = F,
                              skipQC = F,
                              email = T,
                              cpReport = T,
                              redcapUp = T) {
    msgFunName(pipeLnk, "makeReports.v11b6")

    assign("genCn", genCn, envir = gb)
    data <- read.csv(sheetName, strip.white = T)
    runID <- paste0(data$RunID[1])

    isMC = sjmisc::str_contains(runID, "MGDM") | sjmisc::str_contains(runID, "MC")
    is_validation <- sjmisc::str_contains(runID, "VAL")
    if (isMC == T & is_validation == F) {
        CreateRedcapRecord(runID, "control")
    }

    loopRender(selectSams, data, redcapUp)

    checkRunOutput(runID)

    qcVals <- NULL
    if (skipQC == F) {
        CreateRedcapRecord(runID)
        generateQCreport(runID)
        qcVals <- CheckSampleQCmetrics(runID)
        rcon <- redcapAPI::redcapConnection(gb$apiLink, gb$ApiToken)
        redcapAPI::importRecords(rcon, qcVals, "normal", "ids", logfile = "REDCapQCimports.txt")
    }
    if (grepl("TEST", runID)) {
        cpReport = F; redcapUp = F; email = F
    }
    if (redcapUp == T) {
        file.list <- dir(pattern = ".html", full.names = T)
        gb$uploadToRedcap(file.list, T)
    }
    if (email == T) {
        RenameFailed(qcVals)
        gb$CombineClassAndQC(output_fi = paste0(runID, "_qc_data.csv"), gb$ApiToken, runDir = runPath, runID)
        launchEmailNotify(runID)
    }

    try(beepr::beep(2), T)
    tidyUpFiles(runID)
}

# FUN: Checks if all the paths are accessible to the Rscript location
checkMounts <- function() {
    msgFunName(pipeLnk,"checkMounts")
    failMount <- lapply(critialMnts, function(driveMount) {
        ifelse(!dir.exists(driveMount), return(T), return(F))
        })
    if (any(failMount == T)) {
        toFix <- paste(critialMnts[which(failMount == T)])
        cat("PATH does not exist, ensure network drive is mounted:", bkRed(toFix),"\n")
        cat("You may need to re-mount one of the following paths:",
            bkBlu(researchMount,"CBioinformatics/"),
            bkGrn(researchMount,"snudem01lab/snudem01labspace"),
            bkBlu(molecularDrive,"/acc_pathology/molecular"), sep = "\n"
        )
        stopifnot(!any(failMount == T))
    }
}

# Work Directory Functions ----------------------------------------------------
CheckBaseDir <- function(baseFolder) {
    if (is.null(baseFolder)) {
        gb$baseDir <- gb$methDir <- gb$baseFolder <- "/Volumes/CBioinformatics/Methylation/Clinical_Runs"
    }else{
        gb$baseDir <- gb$methDir <- gb$baseFolder <- baseFolder
        }
    isDesktop <- stringr::str_detect(baseFolder, "Desktop")
    if (length(isDesktop) != 0)
        if (is.null(baseFolder) & isDesktop == T) {
            warning("Trying to run methylation from Desktop working directory is not allowed")
            message("Try setting baseFolder to '~/Documents/' instead")
            stopifnot(isDesktop == F)
        }
    return(gb$baseFolder)
}

# Sets the working folder directory
SetBaseFolder <- function(token, baseFolder, runID) {
    baseFolder <- CheckBaseDir(baseFolder)
    methylPath <- gb$setRunDir(runID=gb$runID, workFolder = baseFolder)
    message("Working directory set to:\n", crayon::bgGreen(methylPath), "\n")
    gb$methDir <- gb$workFolder <- baseFolder
    gb$setVar("workFolder", baseFolder)
    gb$setVar("ApiToken", token) # assign the ApiToken & print params
    newPathRun <- file.path(baseFolder, runID)
    pathCdCmd <- paste("cd", newPathRun)
    message("Changing directory to:", pathCdCmd)
    system(pathCdCmd)
    setwd(file.path(baseFolder, runID))
    message(bkRed("Current Path"))
    message(getwd())
    return(baseFolder)
}

PrepareRun <- function(token, baseFolder=NULL, runID, runLocal=F, rdInput=F) {
    mainHub <- "https://github.com/NYU-Molecular-Pathology/Methylation/main/R"
    gb$msgFunName(file.path(mainHub,"methylExpress.R"),"PrepareRun")
    if (runLocal == F) {
        gb$checkMounts()
        gb$checkValidRun(runID)
        gb$SetBaseFolder(token, baseFolder, runID)
        gb$copyWorksheetFile(runID = runID) # copies the xlsm file
        gb$readSheetWrite(runID = runID) # reads xlsm and generates input .csv samplesheet
        gb$get.idats() # Copy idat files to current folder from molecular and snuderlabspace to cwd
        gb$moveSampleSheet(baseFolder, runID) #copies outputs temp to desktop for QC.Rmd
    } else{
        baseFolder <- SetBaseFolder(token, baseFolder, runID)
        assign("baseFolder", baseFolder)
        if (rdInput==F) {
            gb$RunLocalIdats(runID, token)
        }
    }
}

# Starts run knitting html reports loop
StartRun <- function(selectRDs = NULL, emailNotify = T, redcapUp = T) {
    mainHub <- "https://github.com/NYU-Molecular-Pathology/Methylation/main/R"
    gb$msgFunName(file.path(mainHub,"methylExpress.R"),"startRun")
    # Re-order sample report generation for priority
    if (!is.null(selectRDs)) {runOrder <- gb$reOrderRun(selectRDs) }else{runOrder <- NULL}
    gb$makeReports.v11b6(
        skipQC = F,            # Don't skip QC generation
        email = emailNotify,   # to email after Run complete
        cpReport = T,          # Flag to copy files to network drive
        selectSams = runOrder, # Prioritize specific RD-numbers
        redcapUp = redcapUp    # Flag to import files to REDCap
    )
}


GetPriorityCases <- function(selectRDs=NULL, samSheet = "samplesheet.csv", kwd="BN0") {
    msgFunName(pipeLnk,"GetPriorityCases")
    csvFi <- read.csv(file.path(getwd(), samSheet))
    BN00 <- which(stringr::str_detect(csvFi$MP_num, kwd))
    if (length(BN00) > 0) {
        if (is.null(selectRDs)) {
            selectRDs <- c(csvFi$Sample_Name[BN00])
        }else{
            selectRDs <- c(selectRDs, csvFi$Sample_Name[BN00])
        }
        message("Prioritizing these cases first:\n", paste(capture.output(selectRDs), collapse=" "))
    }
    return(selectRDs)
}


AssignArgs <- function(runID, baseFolder = NULL, token = NULL, selectRDs = NULL, redcapUp, gb) {
    msgFunName(pipeLnk,"AssignArgs")
    stopifnot(!is.null(token)); stopifnot(!is.null(runID))
    assign("token", token, envir = gb); assign("ApiToken", token, envir = gb)
    assign("redcapUp", redcapUp, envir = gb); assign("redcapUpload", redcapUp,  envir = gb)
    gb$defineParams(methDir = baseFolder, baseDir = baseFolder, ApiToken = token)
    gb$setVar("runID", runID)
    if (!is.null(selectRDs)) {
        selectRDs <- stringr::str_split(selectRDs, ",")
    }
    return(selectRDs)
}


CheckIdatsCopied <- function() {
    msgFunName(pipeLnk,"CheckIdatsCopied")
    idatFiles <- dir(getwd(), pattern=".idat", full.names = T)
    idatSize <- paste(round(file.info(idatFiles)$size/1024.0 ^ 2, 2),"MB")
    idatFiSi <- as.data.frame(table(idatSize))

    if (nrow(idatFiSi)>1) {
        message(crayon::bgRed("Some idat files are not the same size!"))
        message("Make sure Red and Grn .idat pairs are the same file size by re-copying any files again.")
        message(paste0(capture.output(idatFiSi), collapse="\n"))
    }
}


StartCustomRun <- function(redcapUp = T) {
    msgFunName(pipeLnk,"StartCustomRun")
    gb$makeReports.v11b6(skipQC = T, email = F, cpReport = F, selectSams = NULL, redcapUp = redcapUp)
}
