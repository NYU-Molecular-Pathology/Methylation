#!/usr/bin/env Rscript
## ---------------------------
## Script name: SaveCNVfromRD.R
## Purpose: Generate methylation CNV png plots from RD-numbers input
## Date Created: March 8, 2024
## Date Last Modified: March 8, 2024
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

args <- commandArgs(TRUE); library("base"); gb <- globalenv(); assign("gb", gb)
if(!requireNamespace("devtools")){install.packages("devtools", dependencies = T, ask = F)}
dsh = "\n================"
dsh2 = "================\n"

# Main arguments input in comandline
token <- args[1]
rd_numbers <- args[2]

# Displays the Input args -----------------------------------------------------------------------------
message(dsh,"Parameters input", dsh2)
message("token: ", token)
message("rd_numbers: ", rd_numbers,"\n")

stopifnot(!is.na(token))
stopifnot(!is.na(rd_numbers))

readFlag <- endsWith(rd_numbers, ".csv") == T
stopifnot(rlang::is_bool(readFlag))

# REDcap Heading Fields -----
flds = c("record_id", "b_number", "tm_number", "accession_number", "block", "diagnosis",
         "organ", "tissue_comments", "run_number", "nyu_mrn", "qc_passed")

# Check redcapAPI Package Version Compatibility -------------------------------------------------------
CheckREDCapVersion <- function(){
    if(!requireNamespace("redcapAPI")){
        params = list('nutterb/redcapAPI', dependencies = T, upgrade = "always", type = "source")
        do.call(devtools::install_github, c(params))
    }
    current_vers <- utils::packageVersion("redcapAPI")
    minimum_vers <- "2.7.4"
    if (is.na(current_vers) || utils::compareVersion(as.character(current_vers), minimum_vers) < 0) {
        install.packages("redcapAPI", ask = FALSE, update = TRUE, dependencies = TRUE)
    }
}


supM <- function(sobj){return(suppressMessages(suppressWarnings(sobj)))}

# FUN: Checks if z-drive is accessible to the Rscript
checkMounts <- function(){
    molecDrive = "/Volumes/molecular/MOLECULAR LAB ONLY"
    zDrive = "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular"
    failMount <- ifelse(dir.exists(molecDrive), T, F)
    if(failMount!=T){
        cat("PATH does not exist, ensure path is mounted:")
        cat(crayon::white$bgRed$bold(molecDrive))
        cat("You must mount the network Z-drive path:\n")
        cat(crayon::white$bgRed$bold(zDrive),"\n")
        stopifnot(!any(failMount==T))
    } else {message("\n", crayon::bgGreen("Z-drive path is accessible"), "\n")}
}


# Functions to load packages and get redcap info -----
loadPacks <- function(){
    pkgs = c("data.table", "foreach", "openxlsx", "jsonlite", "RCurl", "readxl",
             "stringr", "webshot", "tidyverse", "crayon", "plotly", "htmlwidgets", "tinytex")
    rlis = getOption("repos")
    rlis["CRAN"] = "http://cran.us.r-project.org"
    options(repos = rlis)
    invisible(lapply(pkgs, function(pk){
        if(supM(!requireNamespace(pk))){
            install.packages(pk, dependencies = T, verbose = T,
                             repos = "http://cran.us.r-project.org")
        }}))
    if(!requireNamespace("systemfonts")){devtools::install_github('r-lib/systemfonts')}
    if(!requireNamespace("redcapAPI")){install.packages("redcapAPI", dependencies = T, ask=F)}
    if(!requireNamespace("remotes")){install.packages("remotes", dependencies = T)}
    supM(library("redcapAPI"))
}


# FUN: Sets your directory and sources the helper functions
sourceFuns2 <- function(workingPath = NULL) {
    mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/"
    script.list <- c("R/SetRunParams.R", "R/CopyInputs.R",
                     "PACT_scripts/generateCNV.R", "Research/cnvFunctions.R")
    if (is.null(workingPath)) {workingPath = getwd()}
    scripts <- paste0(mainHub, script.list)
    invisible(lapply(scripts, function(i){
        suppressPackageStartupMessages(devtools::source_url(i))}))
    #supM(library("sest"))
    #supM(library("mnp.v11b6"))
    gb$setDirectory(workingPath)
    return(gb$defineParams())
}


msgRDs <- function(rds, token){
    message("\nRD-numbers with idats:\n", paste(rds, collapse="\n"))
    assign("rds", rds)
    message(dsh, crayon::bgMagenta("Starting CNV PNG Creation"), dsh2)
    ApiToken <- token
    assign("ApiToken", ApiToken)
}


grabRDs1 <- function(rd_numbers, token){
    result <- gb$search.redcap(rd_numbers, token)
    result <- result[!is.na(result$barcode_and_row_column),]
    samplesheet_ID = as.data.frame(stringr::str_split_fixed(result[,"barcode_and_row_column"],"_",2))
    gb$writeFromRedcap(result, samplesheet_ID) # writes API export as minfi dataframe sheet
    gb$get.idats()  # copies idat files from return to current directory
}

msgCreated <- function(mySentrix){
    pngFiles <- paste0(file.path(fs::path_home(),"Desktop",mySentrix[, 1]),"_cnv.png")
    cnvMade <- file.exists(pngFiles)
    if(any(cnvMade==F)){
        message("The following failed to be created:")
        print(pngFiles[!cnvMade])
        message("Try running again or check GitHub troubleshooting")
    }
    if(any(cnvMade==T)){
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


loopCNV <- function(mySentrix, asPNG){
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


CheckIfPngExists <- function(rds, outFolder = "/Volumes/molecular/Molecular/MethylationClassifier/CNV_PNG") {
    outpng <- paste0(rds, "_cnv.png")
    outFiles <- file.path(outFolder, outpng)
    finished <- file.exists(outFiles)
    if (any(finished)) {
        message(crayon::bgGreen("The following samples are completed and will be skipped:"))
        message(paste(capture.output(outFiles[finished]), collapse = '\n'))
        rds <- rds[!finished]
    }
    return(rds)
}


gb$SaveConumeePACT <-  function(x, sampleImg, doXY=F){
    chrAll <- paste0("chr", 1:22)
    if(doXY==T){
        chrAll <- "all"
    }
    message("Saving file to:\n", sampleImg)
    png(filename = sampleImg, width = 1820, height = 1040, res=150)
    conumee::CNV.genomeplot(x, chr = chrAll)
    invisible(dev.off())
}


gb$SaveCNVplotsPACT <-
    function(samplename_data,
             sentrix.ids,
             i,
             idatPath = NULL,
             chrNum = NULL,
             doXY = F) {

    if (is.null(idatPath)) {
        idatPath <- getwd()
    }
    samName <- samplename_data[i]
    sampleEpic <- sentrix.ids[i]
    sampleImg <- file.path(fs::path_home(), "Desktop", paste0(samName, "_cnv.png"))
    pathEpic <- file.path(idatPath, sampleEpic)
    RGsetEpic <- minfi::read.metharray(pathEpic, verbose = T, force = T)
    if(RGsetEpic@annotation[['array']] == "IlluminaHumanMethylation450k"){
        library("mnp.v11b4")
        MsetEpic <- mnp.v11b4::MNPpreprocessIllumina(
            RGsetEpic, bg.correct = T, normalize = "controls")
    }else{
        #library("mnp.v11b6")
        MsetEpic <- mnp.v11b6::MNPpreprocessIllumina(
            RGsetEpic, bg.correct = T, normalize = "controls")
    }
    x <- gb$customCNV(MsetEpic, samName, NULL)
    slot(x, 'detail', check = FALSE) <- NULL
    invisible(format(object.size(x), units = 'auto'))
    gb$SaveConumeePACT(x, sampleImg, F)
}


LoopSavePlainCNV3 <- function(targets) {
    samplename_data <- as.character(targets[,1])
    sentrix.ids <- as.character(targets$SentrixID_Pos)
    for (i in 1:length(sentrix.ids)) {
        gb$SaveCNVplotsPACT(samplename_data, sentrix.ids, i)
    }
}


TryCnvMaker <- function(myDt) {
    tryCatch(
        expr = {
            gb$makeCNV(myDt)
        },
        error = function(e) {
            message("The following error occured:\n", e)
            message("\n\nTry checking the troubleshooting section on GitHub:\n")
            message("https://github.com/NYU-Molecular-Pathology/Methylation/PACT_scripts/README.md\n"
            )
        },
        finally = {
            gb$copyOutputPng()
        }
    )
}


makeCNV <- function(myDt, asPNG = T) {
    mySentrix <- myDt[myDt[, "SentrixID_Pos"] %like% "_R0", ]
    if (nrow(mySentrix) > 0) {
        gb$LoopSavePlainCNV3(myDt)
    } else{
        message("The RD-number(s) do not have idat files in REDCap:/n")
        print(myDt)
    }
    msgCreated(mySentrix)
    while (!is.null(dev.list())) {dev.off()}
}


CheckIfPngExists <- function(rds, outFolder = "/Volumes/molecular/Molecular/MethylationClassifier/CNV_PNG") {
    outpng <- paste0(rds, "_cnv.png")
    outFiles <- file.path(outFolder, outpng)
    finished <- file.exists(outFiles)
    if (any(finished)) {
        message(
            crayon::bgGreen("The following samples are completed and will be skipped:"),
            "\n",
            paste(capture.output(outFiles[finished]), collapse = '\n')
        )
        rds <- rds[!finished]
    }
    return(rds)
}


TryCnvMaker <- function(myDt) {
    tryCatch(
        expr = {
            gb$makeCNV(myDt)
        },
        error = function(e) {
            message("The following error occured:\n", e)
            message("\n\nTry checking the troubleshooting section on GitHub:\n")
            message(
                "https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/PACT_scripts/README.md\n"
            )
        },
        finally = {
            gb$copyOutputPng()
        }
    )
}


GetSampleList <- function(rds, sampleSheet="samplesheet.csv"){
    myDt <- read.csv(sampleSheet)
    toDrop <- myDt$Sample_Name %in% rds
    myDt <- myDt[toDrop,]
    rownames(myDt) <- 1:length(rds)
    return(myDt)
}


Grn <- function(strMsg){
    message(crayon::bgGreen(strMsg))
}


QueCnvMaker <- function(output, token) {
    rds <- output$record_id[output$report_complete == "YES"]
    if (all(!is.null(rds)) == T & all(!is.na(rds)) == T & length(rds) > 0) {
        msgRDs(rds, token)
        rds <- CheckIfPngExists(rds)
        if (length(rds) > 0) {
            sourceFuns2()
            grabRDs1(rds, token)
            myDt <- GetSampleList(rds)
            TryCnvMaker(myDt)
            try(gb$copyOutputPng(), silent = T)
        } else{
            Grn("No CNV png images to generate.  Check the output directory.")
        }
    } else{
        Grn("The PACT run has no cases with methylation completed.")
        Grn("No CNV png images will generate.")
    }
}

# Load Required Packages ----------------------------------------------------------------
loadPacks()
CheckREDCapVersion()
checkMounts()

# CNV PNG Creation ----------------------------------------------------------------------
if(readFlag){
    rds <- as.data.frame(read.csv(rd_numbers, header = FALSE))[, 1]
    filtered_rds <- rds[grep("^RD-", rds)]
    output <- data.frame(record_id = filtered_rds)
}else{
    output <- data.frame(record_id = rd_numbers)
}

output$report_complete <- "YES"

if(nrow(output) > 0) {
    QueCnvMaker(output, token)
}
