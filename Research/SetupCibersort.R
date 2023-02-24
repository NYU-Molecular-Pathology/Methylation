#!/usr/bin/env R
gb <- globalenv(); assign("gb", gb)
options(install.packages.compile.from.source = "never")
options(install.packages.check.source = "no")

InstPkg <- function(pkg, character.only=T){install.packages(pkg, dependencies=T, verbose=T, ask=F, quiet = F)}
CheckReq <- function(pkgName){if(suppressWarnings(!require(pkgName, character.only = T, warn.conflicts = F))){InstPkg(pkgName)}}

isMac <- ifelse(Sys.info()[['sysname']] == "Darwin", T, F)
if(isMac==T){pkgPath <- "/Volumes/CBioinformatics/Methylation/Other_Packages"}

gitLink <- "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/Research/pullRedcap_manual.R"
path_to_file <- file.path(pkgPath, "MethylCIBERSORT_0.2.0.tar.gz")
path_to_ciber <- file.path(pkgPath, "CIBERSORT.R")

# Installation ----------------------
CheckReq("devtools")
CheckReq("BiocManager")
CheckReq("librarian")

# Install all Required Packages  ------------------------------------------------------------------------------------------------
package_list = c(
    "BiocGenerics",
    "IRanges",
    "S4Vectors",
    "ggpubr",
    "rstatix",
    "cellranger",
    "styler",
    "minfi",
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "IlluminaHumanMethylationEPICmanifest",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    "limma",
    "FlowSorted.Blood.450k",
    "caret",
    "glmnet",
    "NMF",
    "Rserve",
    "e1071",
    "parallel",
    "preprocessCore",
    "colorRamps",
    "ggplot2",
    "reshape2",
    "ggpubr",
    "readxl",
    "stringr",
    "tibble",
    "dplyr",
    "doParallel",
    "matrixStats",
    "sva"
)

librarian::shelf("tidyr", ask=F, update_all = F, quiet = F)
librarian::shelf(package_list, ask=F, update_all = F, quiet = F)
devtools::source_url(gitLink)

MergeSamSheets <- function(inputSheet, sampleSheet, colVariables){
    ogExcel <- as.data.frame(readxl::read_xlsx(inputSheet))
    if(any(duplicated(ogExcel))){
        toDrop <- which(duplicated(ogExcel))
        ogExcel <- ogExcel[-(toDrop),]
        rownames(ogExcel) <- 1:nrow(ogExcel)
    }
    redCsv <- as.data.frame(read.csv(sampleSheet))
    for (variableCol in colVariables) {
        redCsv[, variableCol] <- ""
    }
    for (varCol in colVariables) {
        for (sam in redCsv[,1]) {
            currSamIdx <- which(ogExcel[,1]==sam)
            data2add <- ogExcel[currSamIdx, varCol]
            newSamIdx <- which(redCsv[,1]==sam)
            redCsv[newSamIdx, varCol] <- data2add
        }}
    write.csv(redCsv, file = sampleSheet, quote = F, row.names = F)
    return(redCsv)
}


# FUN: Save object as Rdata file
SaveObj <- function(object, file.name){
    outfile <- file(file.name, "wb")
    serialize(object, outfile)
    close(outfile)
}

# FUN: Load object from Rdata file with progressbar
LoadRdatObj <- function(file.name, msgProg=T){
    library("foreach")
    library("utils")
    filesize <- file.info(file.name)$size
    chunksize <- ceiling(filesize / 100)
    pb <- txtProgressBar(min = 0, max = 100, style=3)
    infile <- file(file.name, "rb")
    if(msgProg==T){
        data <- foreach(it = icount(100), .combine = c) %do% {
            setTxtProgressBar(pb, it)
            readBin(infile, "raw", chunksize)
        }
    }else{
        data <- foreach(it = icount(100), .combine = c) %do% {readBin(infile, "raw", chunksize)}
    }
    close(infile)
    close(pb)
    return(unserialize(data))
}

if(!require("MethylCIBERSORT")){install.packages(path_to_file, repos = NULL, type="source")}
library("MethylCIBERSORT")
require("stats")
source(path_to_ciber)

# COMABAT batch correction between groups
batchCorrectBs <- function(betas, targets, batch_col) {
    targets$facilityBatch <- targets[, batch_col]
    batch <- targets$facilityBatch
    modcombat = model.matrix(~ 1, data = targets)
    gc(verbose = F)
    combat_beta <- sva::ComBat(
        dat = betas, batch = batch, mod = modcombat, par.prior = T, prior.plots = F
    )
    # fix if any beta >1 or < 0, so impute >1 is 1 and <0 is 0
    combat_beta <- ifelse(combat_beta < 0, 0, combat_beta)
    combat_beta <- ifelse(combat_beta > 1, 1, combat_beta)
    return(combat_beta)
}


GetCiberNoob <- function(targets, rgSetDataFi="rgsetData.Rdata", noobValData="noob_vals.Rdata"){
    if(!file.exists(rgSetDataFi)) {
        RGSet <- minfi::read.metharray.exp(targets = targets, force = T, verbose = T)
        gb$SaveObj(RGSet, file.name = rgSetDataFi)
    } else{
        RGSet <- gb$LoadRdatObj(rgSetDataFi)
    }
    # Apply noob normalization based on methods section of methylcibersort
    if(!file.exists(noobValData)) {
        MSet.noob <- minfi::preprocessNoob(RGSet, verbose = T)
        gb$SaveObj(MSet.noob, file.name = noobValData)
    } else{
        MSet.noob <- gb$LoadRdatObj(noobValData)
    }
    RSet <- minfi::ratioConvert(MSet.noob, what = "both", keepCN = T)
    GRset <- minfi::mapToGenome(RSet)
    return(GRset)
}


GetCiberBetas <- function(GRset, targets, batch_col=NULL, betaValData="ciberBeta.Rdata"){
    if(!file.exists(betaValData)) {
        beta <- minfi::getBeta(GRset)
        colnames(beta) <- targets[, 1]
        if (!is.null(batch_col)) {
            beta <- gb$batchCorrectBs(beta, targets, batch_col)
        }
        gb$SaveObj(beta, file.name = betaValData)
    } else{
        beta <- gb$LoadRdatObj(betaValData)
    }
    return(beta)
}


