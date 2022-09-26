gb <- globalenv(); assign("gb", gb)
SilentLoadLib <- function(pkg){
    stopifnot(suppressPackageStartupMessages(library(
        verbose = F,
        warn.conflicts = F,
        quietly = T,
        character.only = T,
        package = pkg,
        logical.return = T
    )))
}
mealPkgs <- c("MEAL", "MultiDataSet",  "minfi", "ggplot2")
invisible(lapply(mealPkgs, SilentLoadLib))

GetMethSet <- function(RGSet){
    Mset <- minfi::preprocessIllumina(RGSet)
    meth <- mapToGenome(ratioConvert(Mset))
    rowData(meth) <- getAnnotation(meth)[, -c(1:3)]
    meth <- dropMethylationLoci(meth) ## Remove probes measuring SNPs
    meth <- dropLociWithSnps(meth) ## Remove probes with SNPs
    meth <- meth[!apply(getBeta(meth), 1, function(x) any(is.na(x))), ] ## Remove probes with NAs
    set.seed(0)  ## Select a subset of samples
    meth <- meth[sample(nrow(meth), 100000), ]
}

GetMealGene <- function(X, resAdj, listOfGenes){
    geneName <- paste0(listOfGenes[X])
    return(MEAL::getGeneVals(resAdj, geneName, genecol = "UCSC_RefGene_Name", fNames = c("chromosome", "start")))
}

GetNewGeneVals <- function(geneDataLi, resAdj) {
    listOfGenes <- geneDataLi[,"Gene"]
    X = 1:length(listOfGenes)
    
    geneValDf <- lapply(
        X,
        FUN = function(X, resAdj, listOfGenes) {
            geneName <- paste0(listOfGenes[X])
            theRow <-
                MEAL::getGeneVals(
                    resAdj,
                    geneName,
                    genecol = "UCSC_RefGene_Name",
                    fNames = c("chromosome", "start")
                )
            return(as.data.frame(theRow))
        },
        resAdj,
        listOfGenes
    )
    
    newValues <- data.frame()
    for (x in geneValDf) {
        newValues <- rbind(newValues, x)
    }
    
    for (gn in 1:length(newValues$UCSC_RefGene_Name)) {
        newName <-
            unlist(stringr::str_split(newValues$UCSC_RefGene_Name[gn], pattern = ";"))
        newName <- paste(unique(newName), collapse = ";")
        newValues[gn, "UCSC_RefGene_Name"] <- newName
    }
    newValues$Region <- ""
    for (rgn in 1:length(geneDataLi$Gene)) {
        theName <- geneDataLi$Gene[rgn]
        for (xName in 1:length(newValues$UCSC_RefGene_Name)) {
            theGene <-
                unlist(stringr::str_split(newValues$UCSC_RefGene_Name[xName], ";"))
            if (any(grepl(pattern = theName, theGene)) == T) {
                newValues[xName, "Region"] <- geneDataLi$Site[rgn]
            }
        }
    }
    return(newValues)
}
                        
# DMRcate -------------------
                        
PlotQCPdf <- function(targPairs, detP){
    pal <- brewer.pal(8,"Dark2")
    par(mfrow=c(1,2))
    barplot(colMeans(detP), col=pal[factor(targPairs$Sample_Group)], las=2,
            cex.names=0.8, ylab="Mean detection p-values")
    abline(h=0.05,col="red")
    legend("topleft", legend=levels(factor(targPairs$Sample_Group)), fill=pal,
           bg="white")
    
    barplot(colMeans(detP), col=pal[factor(targPairs$Sample_Group)], las=2,
            cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
    abline(h=0.05,col="red")
    legend("topleft", legend=levels(factor(targPairs$Sample_Group)), fill=pal,
           bg="white")
    qcReport(rgSet, sampNames=targPairs$ID, sampGroups=targPairs$Sample_Group,
             pdf="qcReport.pdf")
}
                        
PlotDensityMds <- 
function(targPairs,mSetSq){
    cat("\n\n")
    cat("### MDS")
    plotMDS(getM(mSetSq), top=1000, gene.selection="common",
            col=pal[factor(targPairs$Sample_Group)])
    legend("top", legend=levels(factor(targPairs$Sample_Group)), text.col=pal,
           bg="white", cex=0.7)
    cat("\n\n")
    cat("### MDS PCA 1 & 3")
    plotMDS(getM(mSetSq), top=1000, gene.selection="common",
            col=pal[factor(targPairs$Sample_Group)], dim=c(1,3))
    legend("top", legend=levels(factor(targPairs$Sample_Group)), text.col=pal,
           cex=0.7, bg="white")
    cat("\n\n")
    cat("### MDS PCA 2 & 3")
    plotMDS(getM(mSetSq), top=1000, gene.selection="common",
            col=pal[factor(targPairs$Sample_Group)], dim=c(2,3))
    legend("topleft", legend=levels(factor(targPairs$Sample_Group)), text.col=pal,
           cex=0.7, bg="white")
    cat("\n\n")
}
                        
PlotDimensions<-
function(mSetSqFlt,targPairs){
    par(mfrow=c(1,1))
      cat("### Common Genes Group")
    plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
            col=pal[factor(targPairs$Sample_Group)], cex=0.8)
    legend("right", legend=levels(factor(targPairs$Sample_Group)), text.col=pal,
           cex=0.65, bg="white")
    cat("\n\n")
    
     cat("### Common Genes Source")
    plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
            col=pal[factor(targPairs$Sample_Source)])
    legend("right", legend=levels(factor(targPairs$Sample_Source)), text.col=pal,
           cex=0.7, bg="white")
    cat("\n\n")
    par(mfrow=c(1,1))
    cat("### Dimensions PCA 1 & 3")
    # Examine higher dimensions to look at other sources of variation
    plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
            col=pal[factor(targPairs$Sample_Source)], dim=c(1,3))
    legend("right", legend=levels(factor(targPairs$Sample_Source)), text.col=pal,
           cex=0.7, bg="white")
    cat("\n\n")
    cat("### Dimensions PCA 2 & 3")
    plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
            col=pal[factor(targPairs$Sample_Source)], dim=c(2,3))
    legend("topright", legend=levels(factor(targPairs$Sample_Source)), text.col=pal,
           cex=0.7, bg="white")
    cat("\n\n")
    cat("### Dimensions PCA 3 & 4")
    plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
            col=pal[factor(targPairs$Sample_Source)], dim=c(3,4))
    legend("right", legend=levels(factor(targPairs$Sample_Source)), text.col=pal,
           cex=0.7, bg="white")
    cat("\n\n")
}
                        
GetSamPairs <- function(targets, RGSet, samplePairNum = NULL) {
    pal <- brewer.pal(8, "Dark2")
    targPairs <- targets
    targPairs$Sample_Source <- targPairs[,samplePairNum]
    toKeep <- !is.na(targPairs$Sample_Source)
    targPairs <- targPairs[toKeep, ]
    targPairs$ID <-
        paste(targPairs$Sample_Group, targPairs$Sample_Name, sep = ".")
    toDrop <- RGSet@colData@listData[["Sample_Name"]] %in% targPairs$Sample_Name
    rgSet <- RGSet[,toDrop]
    sampleNames(rgSet) <- targPairs$ID
    detP <- detectionP(rgSet)
    
    #PlotQCPdf(targPairs, detP)
    
    keep <- colMeans(detP) < 0.05
    rgSet <- rgSet[, keep]
    targPairs <- targPairs[keep, ]
    detP <- detP[, keep]
    mSetSq <- preprocessQuantile(rgSet)
    mSetRaw <- preprocessRaw(rgSet)
    cat("## Pairwise Analysis Pheochromocytoma {.tabset} ")
    gb$PlotDensityMds(targPairs, mSetSq)
    cat("\n\n")
    
    detP <- detP[match(featureNames(mSetSq), rownames(detP)), ]
    keep <- rowSums(detP < 0.01) == ncol(mSetSq)
    mSetSqFlt <- mSetSq[keep, ]
    annEPIC <- getAnnotation(rgSet)
    keep <-
        !(featureNames(mSetSqFlt) %in% annEPIC$Name[annEPIC$chr %in% c("chrX", "chrY")])
    mSetSqFlt <- mSetSqFlt[keep, ]
    mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
    
    return(list("targPairs" = targPairs,
                "mSetSqFlt" = mSetSqFlt))
}
                        
SubSetGroup <- function(strPatt, samPairs){
    idx <- which(grepl(samPairs, pattern = strPatt, ignore.case=T))
    return(samPairs[idx])
}
                        
