gb <- globalenv(); assign("gb", gb)
# FUN: Obtain RGSet with Probes common to 450K Array ---------------------------
combine.EPIC.450K <- function(targets, batchFilter=NULL){
    if(length(targets$Batch)>0){
        message("Here are your targets$Batch options: ");print(unique(targets$Batch))
        targets_450k <- subset(targets, targets$Batch != batchFilter)
        targets_850k <- subset(targets, targets$Batch == batchFilter)
        if (!anyDuplicated(targets_450k$Basename) == F) {
            return(print(paste("450K duplicated basenames:", anyDuplicated(targets_450k$Basename))))}
        if (!anyDuplicated(targets_450k$Basename) == F) {
            return(print(paste("850K duplicated basenames:", anyDuplicated(targets_850k$Basename))))}
        else{cat("reading 450K arrays...")
            RGSet_450k <- minfi::read.metharray.exp(targets = targets_450k, force = T, verbose = T)
            cat("reading EPIC arrays...")
            RGSet_850k <- minfi::read.metharray.exp(targets = targets_850k, force = T, verbose = T)
            cat("Combining common probes...")
            RGSet <- minfi::combineArrays(RGSet_450k, RGSet_850k, outType = "IlluminaHumanMethylation450k")
        }
    }else{
        RGSet <- minfi::read.metharray.exp(targets = targets, force = T, verbose = T)
    }
    return(RGSet)
}

# Remove Low Quality and select probes using annotations
cleanUpProbes <- function(RGSet, saNames, targets, getfunorm=F){
    detP <- detectionP(RGSet)
    colnames(detP) <- saNames
    keep <- colMeans(detP) < 0.05
    RGSet <- RGSet[, keep]
    targets <- targets[keep,]
    detP <- detP[, keep]
    mSetSq <- preprocessQuantile(RGSet)
    detP <- detP[match(featureNames(mSetSq), rownames(detP)), ]
    keep <- rowSums(detP < 0.01) == ncol(mSetSq)
    gset.funnorm <- addSnpInfo(mSetSq[keep, ])
    gset.funnorm <- dropLociWithSnps(gset.funnorm,snps=c("SBE","CpG"),maf = 0)#drop loci with snps
    annot = getAnnotation(gset.funnorm) #getting annotation file
    sex_probes = annot$Name[annot$chr %in% c("chrX", "chrY")] # dropping sex probes
    gset.funnorm = gset.funnorm[!(rownames(gset.funnorm) %in% sex_probes),]
    if (getfunorm){return(gset.funnorm)} else {
        betas <- getBeta(gset.funnorm)
        colnames(betas) <- saNames
        return(betas)
    }
}

# Returns topVar beta probes instead of all probes
takeTopVariance <- function(betas, topVar){
    var_probes <- apply(betas, 1.0, var)
    select_var <- names(sort(var_probes[topVar], decreasing = T))
    top_var_beta <- betas[select_var, ]
    return(top_var_beta)
}

getSupervise <- function(the_beta, RGSet, topVar=1:10000){
    condition <- pData(RGSet)$Type
    dmp <- dmpFinder(the_beta, pheno = condition, type = "categorical")
    dmp <- cbind(dmp, ID = rownames(dmp))
    betas_df <- as.data.frame(the_beta)
    final_sam <- row.names((t(betas_df[row.names(dmp[topVar, ]), ]))) #topVar=1:10000
    betas_df <- betas_df[, colnames(betas_df) %in% final_sam]
    betas <- as.matrix(betas_df)
    return(betas)
}

# COMABAT batch correction between 2 groups
batchCorrectBs <- function(betas,RGSet,topVar=NULL, supervise = F) {
    batch = targets$Batch
    modcombat = model.matrix(~ 1, data = targets)
    combat_beta <- sva::ComBat(
        dat = betas,
        batch = batch,
        mod = modcombat,
        par.prior = T,
        prior.plots = F
    )
    if (supervise == T) {combat_beta <-getSupervise(combat_beta,RGSet)}
    if (!is.null(topVar)) {
        var_probes <- apply(combat_beta, 1.0, var)
        select_var <- names(sort(var_probes[topVar], decreasing = T))
        top_variable_beta <- combat_beta[select_var,]
    } else {top_variable_beta <- combat_beta}
    # fix if any beta >1 or < 0, so impute >1 is 1 and <0 is 0
    topVarBetas <- ifelse(top_variable_beta < 0, 0, top_variable_beta)
    topVarBetas <- ifelse(topVarBetas > 1, 1, topVarBetas)
    return(topVarBetas)
}

getRgset <- function(rgOut, targets){
if (!file.exists(rgOut)) {
  if (batchCorrect==T & !is.null(targets$Batch)) {
    RGSet <- gb$combine.EPIC.450K(targets = targets)
  } else{
    RGSet <- read.metharray.exp(base = getwd(), targets = targets, verbose = T, force = T)
  }
  saveRDS(RGSet, file = rgOut)
} else{
  RGSet <- readRDS(rgOut)
}
  return(RGSet)
}

getMdsPlot <-function(RGSet, samNames,samTypes, topN=1000) {
    mSetSq <- preprocessQuantile(RGSet)
    mSetSq <- addSnpInfo(mSetSq)
    mSetSq <- dropLociWithSnps(mSetSq, snps = c("SBE", "CpG"), maf = 0) # drop the loci which has snps
    annot = getAnnotation(mSetSq)
    sex_probes = annot$Name[annot$chr %in% c("chrX", "chrY")]
    mSetSq = mSetSq[!(rownames(mSetSq) %in% sex_probes), ]
    mSetSq.beta <- minfi::getBeta(mSetSq)
    colnames(mSetSq.beta) <- mSetSq$Sample_Name
    assign("mSetSq.beta",mSetSq.beta, envir=gb)
    plotNam <- paste0("top_",topN,"_msetBeta",".png")
    myColors <- RColorBrewer::brewer.pal(length(unique(samTypes)),"Set1")
    names(myColors) <- unique(samTypes)
    color_easy = c(myColors)[samTypes]
    png(filename=plotNam, width = 1400, height = 800, res=400) 
    plotMDS(
        mSetSq.beta,
        top = topN,
        gene.selection = "common",
        plot = T,
        #pch = c(1:length(unique(samTypes))),
        main = paste("Top", topN, "Common", "mSet Sq.beta", "MDS plot"),
        col =color_easy,
        labels=samNames
    )
    legend(
        "topright",
        legend = c(names(myColors)),
        col = paste(as.list(myColors)),
        pch = 15,
        cex = 1.5
    )
    dev.off()
    return(mSetSq.beta)
}

cleanRawProbes <- function(rawBetaDat, RGSet, samNames, targets) {
  if (!file.exists(rawBetaDat)) {
    betas <- gb$cleanUpProbes(RGSet, samNames, targets)
    saveRDS(betas, file = rawBetaDat)
  } else{
    betas <- readRDS(rawBetaDat)
  }
  return(betas)
}
