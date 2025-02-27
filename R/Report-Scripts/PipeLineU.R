#!/usr/bin/env Rscript
## ---------------------------
## Script name: PipeLineU.R
## Purpose: source of global functions used for methylation UniD report.Rmd knitting
## Date Last Modified: January 19, 2024
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

TryLoadUniD <- function(rdsPath="/Volumes/CBioinformatics/Methylation/UniD/R/sysdata.rda"){
    library(verbose=F, warn.conflicts = F, quietly = T, package= "dplyr")
    library(verbose=F, warn.conflicts = F, quietly = T, package= "UniD")
    library(verbose=F, warn.conflicts = F, quietly = T, package= "minfi")
    library(verbose=F, warn.conflicts = F, quietly = T, package= "wateRmelon")
    try(library(verbose=F, warn.conflicts = F, quietly = T, package= "UniD"), silent = T)
    try(require("UniD"), silent = T)
    dataPath <- file.path(system.file(package = "UniD"), "R", "sysdata.rda")
    if (!file.exists(dataPath)) {fs::file_copy(path = rdsPath, new_path = dataPath)}
    load(dataPath)
}

UniD_dataqc <- function (loading, detP.cut = 0.05, bc.cut = 3, arrayType) {
    detP <- loading$detP
    sampleSet <- UniD:::GetsampleSet(loading$rgSet)
    Beta.raw <- minfi::getBeta(loading$Mset)
    Beta.raw <- as.data.frame(Beta.raw)
    bc <- wateRmelon::beadcount(loading$rgSet)
    bc[bc < bc.cut] <- NA
    Beta.raw[detP > detP.cut] <- NA
    Beta.raw[is.na(bc)] <- NA
    theFrac <- apply(detP, 2, function(x) length(which(x > detP.cut)))
    p.fail <- data.frame(
        Fail.Frac.detP = theFrac / nrow(detP),
        Fail.Frac.beadcount = colSums(is.na(bc)) / nrow(bc),
        Fail.Frac.NA = colSums(is.na(Beta.raw)) / nrow(Beta.raw)
        )
    cat(paste0(utils::capture.output(p.fail),collapse = "\n"))
    return(Beta.raw)
}

UniD_loadData <- function (sampleID, run_id=NULL) {
    if(grepl("control", sampleID)==T){sampleID <- "control"}
    if(is.null(run_id)|is.na(run_id)){run_id <- "NONE"}
    if(run_id == "NONE"){
      samSh <- "samplesheet.csv"
    }else{
      samSh <- paste0(run_id, "_samplesheet.csv")
    }
    
    inFile <- file.path(fs::path_home(),"Desktop", run_id, samSh)
    if(file.exists(inFile)){
        targets <- read.csv(inFile, strip.white = T)
    }else{
        message("Target samplesheet not found: ", inFile)
        samSh <- paste0("samplesheet.csv")
        inFile2 <- file.path(getwd(), samSh)
        message("Trying to find: ", inFile2)
        targets <- read.csv(inFile2, strip.white = T)
    }
    targRow <- targets[,1]==sampleID
    rgSet <- minfi::read.metharray.exp(targets = targets[targRow,], extended = T, force=T)
    detP <- minfi::detectionP(rgSet)
    Mset <- minfi::preprocessRaw(rgSet)
    loading <- list()
    loading[["rgSet"]] <- rgSet
    loading[["Mset"]] <- Mset
    loading[["detP"]] <- detP
    return(loading)
}
                     
PipelineU <- function(sampleID = "NONE", RGset, run_id = NULL) {
    is450k<-RGset@annotation[[1]] == "IlluminaHumanMethylation450k"
    if(is.null(run_id) | is.na(run_id)){run_id <- "NONE"}
    loading <- suppressWarnings(UniD_loadData(sampleID, run_id))
    outDir <- "."
    arrayType <- ifelse(is450k == F, "EPIC", "450k")
    Beta.raw <- UniD_dataqc(loading, arrayType = arrayType)
    Beta.BMIQ <- UniD::UniD.BMIQ(Beta.raw, arrayType, ".", write = F)
    Beta.clean <- UniD::UniD.probefilter(
        Beta.raw,
        outDir,
        filterXY = F,
        filterSNPHit = F,
        filterMultiHit = F,
        filterNonCG = F,
        filterNonEpic = F,
        arrayType = arrayType,
        filterSample = F ,
        0.1,
        F,
        write = F
    )
    
    Pred <- UniD::UniD.pred(
        inputdata = Beta.clean,
        inputdata.BMIQ = Beta.BMIQ,
        Pred.IDH = T,
        Pred.1p19q = T,
        Pred.ATRX = T,
        Pred.TERTp = T,
        Pred.ExpressSubtype = T,
        outDir = outDir,
        write = F
    )
    predU <- as.data.frame(Pred, row.names = NULL)
    predU$sample <- sampleID
    rownames(predU) <- NULL
    return(predU)
}


