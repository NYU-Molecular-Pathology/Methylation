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


rqPk <- function(pk) {
    if (!require(pk, character.only = T)) {
        install.packages(
            pk,
            warn.conflicts = F,
            ask = F,
            dependencies = T,
            Ncpus = 5
        )
    } else{
        library(
            pk,
            warn.conflicts = F,
            mask.ok = T,
            character.only = T,
            quietly = T
        )
    }
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
    sesameOutput$sesameID = rownames(sesameOutput)
    rownames(targets) <- targets$Sample_Name
    reOrdered <- targets[rownames(sesameOutput),]
    sesameOut <- cbind(sesameOutput, reOrdered)
    sesameOut <- as.data.frame(sesameOut,fix.empty.names=T)
    sesameOut <- apply(sesameOut,2,as.character)
    outputPath <- file.path(fs::path_home(),"Desktop", outFiName)
    write.table(sesameOut, file = outputPath, quote=F, sep=",", row.names=F)
    return(sesameOut)
}


GenerateSesameTab <- function(sampleSheet, outFiName="sesameOutput.csv"){
    msgFunName(sesaLnk, "GenerateSesameTab")
    sesameData::sesameDataCacheAll()
    mcc = parallel::detectCores()/2
    targets <- read.csv(sampleSheet, row.names=NULL)
    ssets = mclapply(sesame::searchIDATprefixes(getwd()), sesame::readIDATpair, mc.cores = mcc)
    ssets <- SwapSesNames(sesameSet = ssets, targets)
    betas = sesame::openSesame(getwd(), BPPARAM=BiocParallel::MulticoreParam(mcc), func = getBetas)
    betas <- SwapSesNames(sesameSet = betas, targets)

    sdfs = sesame::openSesame(getwd(), BPPARAM=BiocParallel::MulticoreParam(mcc), func = NULL) # return SigDF list
    allele_freqs = sesame::openSesame(getwd(), func = getAFs, BPPARAM=BiocParallel::MulticoreParam(mcc))
    sdfs = sesame::openSesame(sdfs, prep = "Q", func = NULL, BPPARAM=BiocParallel::MulticoreParam(mcc))
    model <- sesame::sesameAnno_get("Anno/HM450/Clock_Horvath353.rds")

    ages <- unlist(mclapply(colnames(betas), FUN = function(n){
        theAge <- sesame::predictAge(betas[,n], model)
        names(theAge) <- n
        return(theAge)
        }))
    sesame_sex <- mclapply(ssets, sesame::inferSex, platform="EPIC", mc.cores = mcc)
    karyotype_sex <- mclapply(ssets, sesame::inferSexKaryotypes, mc.cores = mcc)
    predicted_ethnicity <- mclapply(ssets, sesame::inferEthnicity, mc.cores = mcc)
    tissueRef <- sesameDataGet("EPIC.tissueSignature")

    tissues <- unlist(mclapply(
        colnames(betas), FUN = function(n){
            tissue <- sesame::inferTissue(betas[, n], reference = tissueRef)
            names(tissue) <- n
            return(tissue)},
        mc.cores = mcc))

    sesameOutput <- data.frame(
        karyotype_sex = cbind(karyotype_sex),
        predicted_ethnicity = cbind(predicted_ethnicity),
        sesame_sex = cbind(sesame_sex),
        sesame_ages = cbind(ages),
        sesame_tissue = cbind(tissues)
    )

    sesameOut <- SaveSesameTable(sesameOutput, targets, outFiName)
    return(sesameOut)
    #knitr::kable(sesameOut)
}
LoadPacman()
