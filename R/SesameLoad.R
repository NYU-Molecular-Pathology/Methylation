## title: "sesameEthnicity"
## author: "Jonathan Serrano"
## date: "9/15/2021"
# https://bioconductor.org/packages/release/bioc/vignettes/sesame/inst/doc/sesame.html
# https://doi.org/10.1093/nar/gky691
gb <- globalenv(); assign("gb", gb)
sesaLnk = "https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/R/SesameLoad.R"

msgFunName <- function(pthLnk, funNam) {
    message("\nExecuting function: ", funNam, " from RScript in:\n", pthLnk, "\n")
}

rqPk <- function(pk) {if (!require(pk, character.only = T)) {
    install.packages(pk, warn.conflicts = F, ask = F, dependencies = T, Ncpus = 5)
    } else{library(pk, warn.conflicts = F, mask.ok = T, character.only = T, quietly = T)}
}

DirMake <- function(pathName){if(!dir.exists(pathName)){dir.create(pathName)}}

CheckCachePath <- function(){
    msgFunName(sesaLnk, "CheckCachePath")
    userName <- Sys.info()["user"][[1]]
    cachePath <- file.path("/Users", userName, "Library/Caches/org.R-project.R/R")
    DirMake(file.path(cachePath,"ExperimentHub"))
    DirMake(file.path(cachePath,"BiocFileCache"))
}

# Install Packages ------------------------------------------------------------------------
LoadPacman <- function(){
    msgFunName(sesaLnk, "LoadPacman")
    rqPk("devtools"); rqPk("pacman")
    if(!require("sesameData")){devtools::install_github("zwdzwd/sesameData")}
    if(!require("sesame")){devtools::install_github("zwdzwd/sesame")}
    pacman::p_load("parallel", "BiocParallel", "sesame", "sesameData", "fs", "SummarizedExperiment")
    CheckCachePath()
    tryCatch(sesameData::sesameDataCache(),
             error=function(e){message("Caught error:\n");message(e);sesameData::sesameDataCache()})
}


SwapSesNames <- function(sesameSet, targets) {
    msgFunName(sesaLnk, "SwapSesNames")
    idatName <- if(is.null(names(sesameSet))){colnames(sesameSet)} else{names(sesameSet)}
    matchIndex <- match(idatName, targets$SentrixID_Pos)
    newNamesId <- targets[matchIndex, "Sample_Name"]
    if (is.null(names(sesameSet))) {
        colnames(sesameSet) <- newNamesId
    } else{
        names(sesameSet) <- newNamesId
    }
    return(sesameSet)
}


SaveSesameTable <- function(sesameOutput, targets, outFiName){
    msgFunName(sesaLnk, "SaveSesameTable")
    sesameOutput$sesameID = rownames(sesameOutput)
    rownames(targets) <- targets$Sample_Name
    reOrdered <- targets[rownames(sesameOutput),]
    sesameOut <- cbind(sesameOutput, reOrdered)
    sesameOut <- as.data.frame(sesameOut,fix.empty.names=T)
    sesameOut <- apply(sesameOut,2,as.character)
    outputPath <- file.path(fs::path_home(), outFiName)
    write.table(sesameOut, file = outputPath, quote=F, sep=",", row.names=F)
    return(sesameOut)
}

GenerateSesameTab <- function(outFiName="samples", inputPath=NULL, sampleSheet="samplesheet.csv"){
    msgFunName(sesaLnk, "GenerateSesameTab")
    sesameData::sesameDataCacheAll()
    if(is.null(inputPath)){inputPath <- getwd()}
    outFiName <- paste(outFiName, "sesameOutput.csv", sep="_")

    targets <- read.csv(sampleSheet, row.names = NULL)
    allIdats <- file.path(inputPath, targets$SentrixID_Pos) #sesame::searchIDATprefixes(getwd())
    names(allIdats) <- targets$SentrixID_Pos
    mcc = parallel::detectCores() / 2

    betas <- sesame::openSesame(getwd(), BPPARAM = BiocParallel::MulticoreParam(mcc), func = getBetas)
    betas <- SwapSesNames(betas, targets)

    model <- sesame::sesameAnno_get("Anno/HM450/Clock_Horvath353.rds")
    tissueRef <- sesameData::sesameDataGet("EPIC.tissueSignature")

    samNames <- if(is.null(colnames(betas))){names(betas)} else{colnames(betas)}

    ages <- unlist(parallel::mclapply(samNames, FUN = function(n, betas, model) {
        nSam <- if(is.null(colnames(betas))){betas[[n]]}else{betas[, n]}
        theAge <- sesame::predictAge(nSam, model)
        names(theAge) <- n
        return(theAge)},
        mc.cores = mcc, betas, model))

    tissues <- unlist(parallel::mclapply(samNames, FUN = function(n, betas, tissueRef) {
        nSam <- if(is.null(colnames(betas))){betas[[n]]}else{betas[, n]}
        tissue <- if(length(nSam)==486427){"NOT EPIC"}else{sesame::inferTissue(nSam, reference = tissueRef)}
        names(tissue) <- n
        return(tissue)},
        mc.cores = mcc, betas, tissueRef))

    ssetsRaw = parallel::mclapply(allIdats, sesame::readIDATpair, mc.cores = mcc)
    ssets <- SwapSesNames(sesameSet = ssetsRaw, targets)

    sesame_sex <- parallel::mclapply(ssets, function(x){
        sesame::inferSex(x, platform = ifelse(nrow(x)==486427,"HM450", "EPIC"))}, mc.cores = mcc)

    karyotype_sex <- parallel::mclapply(ssets, sesame::inferSexKaryotypes, mc.cores = mcc)
    predicted_ethnicity <- parallel::mclapply(ssets, sesame::inferEthnicity, mc.cores = mcc)

    sesameOutDf <- data.frame(
        karyotype_sex = cbind(karyotype_sex),
        predicted_ethnicity = cbind(predicted_ethnicity),
        sesame_sex = cbind(sesame_sex),
        sesame_ages = cbind(ages),
        sesame_tissue = cbind(tissues)
    )

    sesameOut <- SaveSesameTable(sesameOutDf, targets, outFiName)
    return(sesameOut)
}

LoadPacman()

