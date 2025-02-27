#!/usr/bin/env Rscript
## ---------------------------
## Script name: EPICV2QC.R
## Purpose: source of global scripts imported for methylation QC for EPICV2
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------

if(!require("IlluminaHumanMethylationEPICv2manifest", mask.ok=T, warn.conflicts=F)){
    BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2manifest", update=F, ask=F)
}

if(!require("IlluminaHumanMethylationEPICv2anno.20a1.hg38", mask.ok=T, warn.conflicts=F)){
    BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38", update=F, ask=F)
}


summarizeWholeBunchV2 <- function (targets, rp.zero, verbose, ...) {
    RGset <- minfi::read.metharray.exp(targets = targets, ...)
    RGset@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
    if (rp.zero) {
        RGset <- MethylAid:::replaceZero(RGset)
    }
    RGset@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
    DP <- minfi::detectionP(RGset)
    DPfreq <- colSums(DP < 0.01, na.rm = TRUE)/nrow(DP)
    RG <-  MethylAid:::summarizeControls(RGset, "IlluminaHumanMethylationEPICv2")
    MU <-  MethylAid:::summarizeMUvalues(RGset)
    if (nrow(targets) > 1) {
        targets <- data.frame(apply(targets, 2, function(x) factor(as.character(x))), row.names = row.names(targets))
    }
    rownames(targets) <- colnames(MU)
    sData <- new("summarizedData", targets = targets, controls = RG$TypeControl, 
        Rcontrols = as.matrix(RG$R), Gcontrols = as.matrix(RG$G), DPfreq = DPfreq, MU = MU)
    return(sData)
}


summarizePerBatchV2 <- function (targets, batchSize, rp.zero, verbose, ...) {
    if (verbose) {message("Summarize data in batches...")}
    R <- G <- MU <- DPfreq <- c()
    tg <- targets
    while (nrow(tg) > 0) {
        ss <- 1:min(batchSize, nrow(tg))
        if (verbose) {message(paste("Summarizing", length(ss), "samples..."))}
        RGset <- minfi::read.metharray.exp(targets = tg[ss, ], ...)
        RGset@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
        if (rp.zero) {RGset <-  MethylAid:::replaceZero(RGset)}
        DP <- minfi::detectionP(RGset)
        DPfreq <- c(DPfreq, colSums(DP < 0.01, na.rm = TRUE)/nrow(DP))
        RG <-  MethylAid:::summarizeControls(RGset, "IlluminaHumanMethylationEPICv2")
        R <- cbind(R, RG$R)
        G <- cbind(G, RG$G)
        MU <- cbind(MU,  MethylAid:::summarizeMUvalues(RGset))
        tg <- tg[-ss, ]
    }
    if (nrow(targets) > 1) {
        targets <- data.frame(apply(targets, 2, function(x) factor(as.character(x))), row.names = row.names(targets))
    }
    rownames(targets) <- colnames(R) <- colnames(G) <- colnames(MU)
    sData <- new("summarizedData", targets = targets, controls = RG$TypeControl,
                 Rcontrols = R, Gcontrols = G, DPfreq = DPfreq, MU = MU)
    return(sData)
}

summarizeParallelV2 <- function (targets, batchSize, BPPARAM, rp.zero, verbose, ...) {
    nworkers <- BiocParallel::bpworkers(BPPARAM)
    y <- rep(1:nworkers, nrow(targets))
    y <- y[1:nrow(targets)]
    jobs <- split(targets, y)
    if (batchSize == -1) {
        sumParallel <- function(x) {
            summarizeWholeBunchV2(x, rp.zero = rp.zero, verbose = verbose, ...)
            }
    } else{
        sumParallel <- function(x) {
            summarizePerBatchV2(x, batchSize = batchSize, rp.zero = rp.zero, verbose = verbose, ...)
        }
    }
    res <- BiocParallel::bplapply(jobs, FUN = sumParallel, BPPARAM = BPPARAM)
    sData <- MethylAid:::reduce(res)
    return(sData)
}


summarize2 <- function (targets, batchSize = -1, BPPARAM = NULL, rp.zero = TRUE, verbose = TRUE, file = NULL, ...) {
    targets$None <- 1
    if (is.null(BPPARAM)) {
        if (batchSize == -1) {
            sData <- summarizeWholeBunchV2(targets, rp.zero, verbose, ...)
        } else {
            sData <- summarizePerBatchV2(targets, batchSize, rp.zero, verbose, ...)
        }
    } else {
        sData <- summarizeParallelV2(targets, batchSize, BPPARAM, rp.zero, verbose, ...)
    }
    sData@plotdata <- MethylAid:::prepareData(sData)
    if (!is.null(file)) {
        assign(basename(file), sData)
        save(list = basename(file), file = paste0(file, ".RData"))
        message(paste("Summarized data stored: ", paste0(file, ".RData")))
    }
    return(invisible(sData))
}


dropDataNa <- function(data_df){
    oldDataX <- data_df$final_data$x
    oldDataY <- data_df$final_data$y
    if(any(is.na(oldDataX))){
        oldDataX[is.na(oldDataX)] <- 0
        data_df$final_data$x <- oldDataX
    }
    if(any(is.na(oldDataY))){
        oldDataY[is.na(oldDataY)] <- 0
        data_df$final_data$y <- oldDataY
    }
    return(data_df)
}

