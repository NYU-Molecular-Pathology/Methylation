library(verbose=F, warn.conflicts = F, quietly = T, package= "UniD")
library(verbose=F, warn.conflicts = F, quietly = T, package= "minfi")
library(verbose=F, warn.conflicts = F, quietly = T, package= "IlluminaHumanMethylationEPICmanifest")
data("IlluminaHumanMethylationEPICmanifest")

UniD_load <- function (sampleID, run_id) {
    samSh <- paste0(run_id,"_samplesheet.csv")
    inFile <- file.path("~","Desktop",run_id,samSh)
    if(file.exists(inFile)){
    targets <- read.csv(inFile, strip.white = T)
    }else{
        inFile <- "samplesheet.csv"
        targets <- read.csv(inFile, strip.white = T)
    }
    targRow <- targets[,1]==sampleID
    rgSet <- read.metharray.exp(targets = targets[targRow,], extended = T, force=T)
    detP <- minfi::detectionP(rgSet)
    Mset <- minfi::preprocessRaw(rgSet)
    loading <- list()
    loading[["rgSet"]] <- rgSet
    loading[["Mset"]] <- Mset
    loading[["detP"]] <- detP
    return(loading)
}

pipelineU <- function(sampleID, is450k=F, run_id) {
    require("UniD")
    dataPath <- file.path(system.file(package = "UniD"), "R", "sysdata.rda")
    if(!file.exists(dataPath)) {
        fs::file_copy(path = "/Volumes/CBioinformatics/Methylation/UniD/R/sysdata.rda", new_path = dataPath)
    }
    load(dataPath)
    loading <- suppressWarnings(UniD_load(sampleID, run_id))
    outDir <- file.path(getwd(),"UniD")
    if(!dir.exists(outDir)){dir.create(outDir)}
    arrayType <- ifelse(is450k==F, "EPIC", "450k")
    Beta.raw <- UniD::UniD.dataqc(loading, outDir, arrayType = arrayType, write = F)
    Beta.BMIQ <- UniD::UniD.BMIQ(Beta.raw, arrayType, outDir, write = F)
    Beta.clean <- UniD::UniD.probefilter(
        Beta.raw, outDir, is450k, is450k, is450k, is450k, is450k,
        arrayType, filterSample = T, 0.1, F, write = F
    )
    Pred <- UniD::UniD.pred(
        inputdata = Beta.clean, inputdata.BMIQ = Beta.BMIQ,
        Pred.IDH = T, Pred.1p19q = T, Pred.ATRX = T, Pred.TERTp = T,
        Pred.ExpressSubtype = T, outDir = outDir, write = F
    )
    pred <- as.data.frame(Pred,row.names = NULL)
    pred$sample <- sampleID
    rownames(pred) <- NULL
    outfi <- paste0(sampleID,"_prediction_uniD.csv")
    write.csv(pred, file.path(outDir, outfi))
    return(pred)
}
