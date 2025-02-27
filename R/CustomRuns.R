#!/usr/bin/env Rscript
## Script name: CustomRuns.R
## Purpose: Source of functions to execute a customized methylation run
## Date Created:  August 29, 2022
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024

gb <- globalenv(); assign("gb", gb)
reportMd <- file.path(fs::path_home(),"report.Rmd") # From curl github download
sarcRmdFile = "/Volumes/CBioinformatics/Methylation/SarcReport13.Rmd"
cpInLnk4 = "https://github.com/NYU-Molecular-Pathology/Methylation/main/R/CustomRuns.R"

msgFunName <- function(pthLnk, funNam) {
    message("\nExecuting function: ", crayon::black$bgYellow(funNam),
            " from RScript in:\n", pthLnk, "\n")
}

GetLocalData <- function(rg){
    msgFunName(cpInLnk4,"GetLocalData")
    dat <- data.frame(
        sampleID = paste0(rg),
        bnumber = "NONE",
        senLi = paste0(rg),
        run_id = paste0(gb$runID),
        mp_number = "NONE",
        tech = "NONE",
        tech2 = "NONE",
        outFi = paste0(rg, ".html")
    )
    return(dat)
}

GetTargetData <- function(data) {
        msgFunName(cpInLnk4,"GetTargetData")
    runDt <- data.frame(
        sampleID = paste0(data[, "Sample_Name"]),
        bnumber = paste0(data[, "DNA_Number"]),
        senLi = paste0(data[, "SentrixID_Pos"]),
        run_id = paste0(data[, "RunID"]),
        mp_number = paste0(data[, "MP_num"]),
        tech = paste0("NONE"),
        tech2 = paste0("NONE"),
        outFi = paste0(data[, 1], ".html")
    )
    return(runDt)
}

KnitReportRmd <- function(dat, token, reportMd){
            msgFunName(cpInLnk4,"KnitReportRmd")
    message(paste0(capture.output(dat), collapse = "\n"))
    RGsetEpic <- suppressWarnings(gb$getRGset(getwd(), dat$senLi))
    RGset <- RGsetEpic[,1]
    assign("RGset", RGset)
    sample = 1
    rmarkdown::render(
        input = reportMd,
        output_format = "html_document",
        output_file = dat$outFi,
        output_dir = getwd(),
        quiet = F,
        output_options = c("self_contained = TRUE"),
        params = list(
            token = token,
            rundata = dat,
            RGsetEpic = RGsetEpic,
            knitDir = getwd()
            )
    )
}

loop_targets <- function(targets, reportMd){
    msgFunName(cpInLnk4,"loop_targets")
    reportMd <- file.path(fs::path_home(),"report.Rmd")
    require(compiler)
    enableJIT(3)
    mainPage = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/Report-Scripts/"
    scripts <- c("ClassTables.R", "MLH1_Functions.R", "PipeLineU.R",
                 "RedcapOutput.R", "TsneFunctions.R", "cnvggplotly.R")
    scripts <- paste0(mainPage, scripts)
    lapply(scripts, function(i) {message("Sourcing: ", i);devtools::source_url(i)})
    gb$moveSampleSheet(gb$baseFolder, gb$runID)
    gb$TryLoadUniD()
    try(silent=T, unloadNamespace("sarc.v12b6"))
    #gb$LoadMnpData(F)
    targets$RunID <- gb$runID
    for(i in 1:nrow(targets)){
        message("Sample ", i, " of ", nrow(targets))
        dat <- GetTargetData(data = targets[i,])
        KnitReportRmd(dat, token = gb$ApiToken, reportMd)
    }
}

loop_local <- function(RGSet){
                    msgFunName(cpInLnk4,"loop_local")
    reportMd <- file.path(fs::path_home(),"report.Rmd")
    for (rg in colnames(RGSet)) {
        thisSam <- RGSet[, rg]
        dat <- GetLocalData(rg)
        sentrix=dat$senLi
        RGsetEpic<-RGset<-thisSam
        rmarkdown::render(
            reportMd, "html_document", dat$outFi, getwd(), quiet = FALSE, output_options = c("self_contained = TRUE"),
            params = list(token = gb$ApiToken, rundata = dat, RGsetEpic = RGsetEpic, knitDir=getwd())
        )
    }
}

PromptInputCsv <- function(runID) {
    msgFunName(cpInLnk4,"PromptInputCsv")
    message('No idat files in current directory and no .xlsm with the runID"',
            runID,'"')
    message("Enter full path to csv file listing RD-numbers and press return")
    message("The CSV should have no header and list RD-numbers in Column A")
    message("Example: /Users/myName/Desktop/myFileList_rd_numbers.csv")
    csvFilePath <- 
        readline("Paste the full path to your csv file and hit return/Enter: ")
    csvFilePath <- as.character(csvFilePath)
    if (!file.exists(csvFilePath)) {
        warning("The file ",csvFilePath," does not exist!")
        cat("Try pasting the full path again and hit return/Enter: ")
        csvFilePath <- readLines(file("stdin"),1)
        csvFilePath <- as.character(csvFilePath)
        message(csvFilePath)
        stopifnot(file.exists(csvFilePath))
    }
    rd_numbers <- suppressWarnings(
        read.csv(file = csvFilePath, header = F, colClasses = "character"))[,1]
    rd_numbers <- as.data.frame(rd_numbers)
    message("The following RD-numbers were entered:")
    print(rd_numbers$rd_numbers)
    return(rd_numbers$rd_numbers)
}

ParseInputCsvPath <- function(samSheetIn) {
    csvFilePath <- as.character(samSheetIn)
    if (!file.exists(csvFilePath)) {
        warning("The file ", csvFilePath, " does not exist!")
        message(csvFilePath)
        stopifnot(file.exists(csvFilePath))
    }
    rd_numbers <- suppressWarnings(
        read.csv(file = csvFilePath, header = F, colClasses = "character"))[,1]
    rd_numbers <- as.data.frame(rd_numbers)
    message("The following RD-numbers were entered:")
    message(paste0(capture.output(rd_numbers), collapse = "\n"))
    return(rd_numbers$rd_numbers)
}


# Overwrite any existing samplesheet in current directory with new RD-number or file input and copy idats
MakeLocalSampleSheet <- function(runID,
                                 token,
                                 samSheetIn = NULL,
                                 rd_numbers = NULL,
                                 outputFi = NULL) {
    if (is.null(outputFi)) {
        outputFi <- "samplesheet_sarc.csv"
    }
    msgFunName(cpInLnk4, "MakeLocalSampleSheet")
    stopifnot(!is.null(token))
    token2 <- token
    idatScript <- "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/Research/pullRedcap_manual.R"
    if (is.null(rd_numbers)) {
        file_data <- NULL
        if (is.null(samSheetIn)) {
            file_data <- PromptInputCsv(runID)
        } else{
            if (grepl("\\.csv$", samSheetIn)) {
                file_data <- ParseInputCsvPath(samSheetIn)
            }
            if (grepl("\\.(xlsm|xlsx)$", samSheetIn)) {
                library("readxl")
                file_data <- read_excel(samSheetIn, sheet = 1, col_names = FALSE)[[1]]
            }
        }
        rd_numbers <- grep("^RD-", file_data, value = TRUE)
    }

    stopifnot(length(rd_numbers) >= 1)
    if (length(rd_numbers) == 1) {
        rd_numbers <- c(rd_numbers, rd_numbers)
    }
    stopifnot(all(grepl("^RD-", rd_numbers)))
    
    message("Sourcing: ", idatScript)
    devtools::source_url(idatScript)
    
    gb$token <- gb$ApiToken <- token <- token2
    
    gb$grabRDCopyIdat(rd_numbers,
                      token2,
                      copyIdats = T,
                      outputFi = outputFi)
}


RunFromSamplesheet <- function(samSh="samplesheet.csv"){
    msgFunName(cpInLnk4, "RunFromSamplesheet")
    targets <- read.csv(file=samSh, strip.white = T)
    targets[is.na(targets)] <- "NONE"
    write.csv(targets, samSh, quote = F,row.names = F)
    targets <- minfi::read.metharray.sheet(base = getwd(), pattern = basename(samSh))
    dupedSamples <- anyDuplicated(targets$Basename)
    if(dupedSamples!=0){
        warning("Removing Samples with identical (duplicated) Sentrix IDs")
        message("Some basenames are duplicated and will be dropped:")
        tgBn <- targets$Basename
        print(targets[tgBn==tgBn[anyDuplicated(tgBn)],])
        targets <- targets[-dupedSamples,]
    }
    rownames(targets) <- targets$Sample_Name
    return(targets)
}


RunLocalIdats <- function(runID, token, samSheet = "samplesheet.csv"){
  msgFunName(cpInLnk4,"RunLocalIdats"); stopifnot(!is.null(token))
  samSheet = file.path(getwd(), "samplesheet.csv")
  labWorkbook <- file.path(getwd(), paste0(runID,".xlsm"))
  if (file.exists(labWorkbook)) {
    gb$readSheetWrite(runID = runID)
    gb$moveSampleSheet(baseFolder, runID)
    return(message("Samplesheet Copied"))
  }
  idatFiles <- dir(path = getwd(), pattern = ".idat", full.names = T)
  if(!file.exists(samSheet) & length(idatFiles)>1) {
    idatBase <- unique(substring(idatFiles, 1, nchar(idatFiles) - 9))
    RGSet <- minfi::read.metharray(basenames = idatBase, force = T, verbose = T)
    loop_local(RGSet)
    return(message("Done!"))
  }
  if (!file.exists(samSheet)) {
    MakeLocalSampleSheet(runID, token)
    }
    targets <- RunFromSamplesheet(samSheet)
    targets$Basename <- paste0(file.path(getwd(), targets$SentrixID_Pos))
  gb$get.idats()
  loop_targets(targets)
  return(message("Done!"))
}


# Work Directory Functions ----------------------------------------------------
CheckBaseDir <- function(baseFolder){
    msgFunName(cpInLnk4,"CheckBaseDir")
    if(is.null(baseFolder)){
        gb$baseDir <- gb$methDir <- gb$baseFolder <- "/Volumes/CBioinformatics/Methylation/Clinical_Runs"
    }else{gb$baseDir <- gb$methDir <- gb$baseFolder <- baseFolder}
    isDesktop <- stringr::str_detect(baseFolder, "Desktop")
    if(length(isDesktop)!=0)
    if(is.null(baseFolder) & isDesktop==T) {
        warning("Trying to run methylation from Desktop working directory is not allowed")
        message("Try setting baseFolder to '~/Documents/' instead")
        stopifnot(isDesktop == F)
    }
    return(gb$baseFolder)
}

# Sets the working folder directory
SetBaseFolder <- function(token, baseFolder, runID){
        msgFunName(cpInLnk4,"SetBaseFolder")
    baseFolder <- CheckBaseDir(baseFolder)
    methylPath <- gb$setRunDir(runID=gb$runID, workFolder = baseFolder)
    message("Working directory set to:\n", crayon::bgGreen(methylPath), "\n")
    gb$methDir <- gb$workFolder <- baseFolder
    gb$setVar("workFolder", baseFolder)
    gb$setVar("ApiToken", token) # assign the ApiToken & print params
    setwd(file.path(baseFolder, runID))
    return(baseFolder)
}

PrepareRun <- function(token, baseFolder=NULL, runID, runLocal=F, rdInput=F){
            msgFunName(cpInLnk4,"PrepareRun")
    gb$checkMounts()
    if(runLocal==F){
        gb$checkValidRun(runID)
        gb$SetBaseFolder(token, baseFolder, runID)
        gb$copyWorksheetFile(runID = runID) # copies the xlsm file
        gb$readSheetWrite(runID = runID) # reads xlsm and generates input .csv samplesheet
        gb$get.idats() # Copy idat files to current folder from molecular and snuderlabspace to cwd
        gb$moveSampleSheet(baseFolder, runID) #copies outputs temp to desktop for QC.Rmd
    } else{
        gb$SetBaseFolder(token, baseFolder, runID)
        if(rdInput==F){
            gb$RunLocalIdats(runID, token)
        }
    }
}


SarcomaReport <- function(RGset, sampleID, output_dir = getwd()) {
                msgFunName(cpInLnk4,"SarcomaReport")
    output_file <- paste0(sampleID, "_sarc.html")
    RGset <- RGset[,1]
    rmarkdown::render(input= gb$sarcRmdFile, output_dir = output_dir, output_file = output_file)
}

PromptRDnumbers <- function(){
                    msgFunName(cpInLnk4,"PromptRDnumbers")
    message("Type your RD-numbers without quotes below")
    message("Separate with spaces & press Enter twice to exit")
    rd_numbers <- scan(what = " ")
    message("Values entered:\n", rd_numbers)
    rd_numbers <- rd_numbers[which(grepl(rd_numbers, pattern="^RD-"))]
    message("Valid RD-numbers entered:\n", paste0(capture.output(rd_numbers), collapse = '\n' ))
    if(length(rd_numbers)==0) {
        message("No valid RD-numbers entered! Try again.")
        PromptRDnumbers()
    }else(return(rd_numbers))
}

MakeBlankRun <- function(rd_numbers, token, outputFi="samplesheet.csv"){
                        msgFunName(cpInLnk4,"MakeBlankRun")
    stopifnot(length(rd_numbers)>=1 & length(rd_numbers)!=0 & stringr::str_detect(rd_numbers[1],"RD-"))
    if(!exists("grabRDCopyIdat")){
        idatScript <- "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/Research/pullRedcap_manual.R"
        message("Sourcing: ", idatScript)
        devtools::source_url(idatScript)
    }
    gb$grabRDCopyIdat(rd_numbers, token, outputFi=outputFi)
    return(as.data.frame(read.csv(outputFi)))
}

MakeSarcomaReport <- function(worksheet = "samplesheet.csv", targets = NULL) {
    msgFunName(cpInLnk4,"MakeSarcomaReport")
    if (is.null(targets)) {
        if (file.exists(worksheet)) {
            targets <- as.data.frame(read.csv(worksheet))
        } else{
            rd_numbers <- PromptRDnumbers()
            targets <- MakeBlankRun(rd_numbers, gb$token, worksheet)
            }
    }

    for (samIdx in 1:nrow(targets)) {
        message("Running ", samIdx, " of ", nrow(targets))
        RGset <- suppressWarnings(minfi::read.metharray(targets$Basename[samIdx]))
        if (RGset@annotation[["array"]] == "IlluminaHumanMethylationEPICv2") {
            sampleID <- targets[samIdx, 1]
            message("Sample ", sampleID, " is EPICV2!  Skipping...")
        } else{
            SarcomaReport(RGset, sampleID = targets[samIdx, 1])
        }
    }
}

StartRun <- function(selectRDs = NULL, emailNotify = T, redcapUp = T) {
    msgFunName(cpInLnk4,"StartRun")
    gb$msgFunName(paste0(mainHub,"methylExpress.R"),"startRun")
    # Re-order sample report generation for priority
    if (!is.null(selectRDs)) {runOrder <- gb$reOrderRun(selectRDs) }else{runOrder <- NULL}
    gb$makeHtmlReports(
        skipQC = F,            # Don't skip QC generation
        email = emailNotify,   # to email after Run complete
        selectSams = runOrder, # Prioritize specific RD-numbers
        redcapUp = redcapUp    # Flag to import files to REDCap
    )
}


LoadGitHubScripts <- function(ghRepo, scriptList){scripts = file.path(ghRepo, scriptList)
return(lapply(scripts, function(i){message("Sourcing: ", i); devtools::source_url(i)}))}

mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R"
rmdScripts <- c("ClassTables.R", "MLH1_Functions.R", "RedcapOutput.R")

LoadGitHubScripts(file.path(mainHub,"Report-Scripts"), rmdScripts)

