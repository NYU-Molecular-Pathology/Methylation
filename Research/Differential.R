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
    
    gb$PlotQCPdf(targPairs, detP)
    
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
                        
GetTracksPlot <- function(annEPICSub, bVals){
    gen <- "hg19"
    pal <- brewer.pal(8,"Dark2")
    # the index of the DMR that we will plot
    dmrIndex <- 1
    # extract chromosome number and location from DMR results
    chrom <- as.character(seqnames(results.ranges[dmrIndex]))
    start <- as.numeric(start(results.ranges[dmrIndex]))
    end <- as.numeric(end(results.ranges[dmrIndex]))
    # add 25% extra space to plot
    minbase <- start - (0.25*(end-start))
    maxbase <- end + (0.25*(end-start))
    dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
    islandHMM <- as.data.frame(read.csv(paste0(
        dataDirectory, "/model-based-cpg-islands-hg19-chr17.txt"),
        sep="\t", stringsAsFactors=FALSE, header=FALSE))
    colnames(islandHMM) <- c("Chr", "start", "end", "v1", "v2", "v3", "v4", "v5")
    islandHMM$start <- as.numeric(islandHMM$start)
    islandHMM$end <- as.numeric(islandHMM$end)
    #islandHMM <- makeGRangesFromDataFrame(islandHMM, keep.extra.columns = TRUE)
    islandData <- GRanges(seqnames=Rle(islandHMM[,1]),
                          ranges=IRanges(start=islandHMM[,2], end=islandHMM[,3]),
                          strand=Rle(strand("*" )))
    dnase <- read.csv(paste0(dataDirectory,"/wgEncodeRegDnaseClusteredV3chr17.bed"),
                      sep="\t",stringsAsFactors=FALSE,header=FALSE)

    dnaseData <- GRanges(seqnames=dnase[,1],
                         ranges=IRanges(start=dnase[,2], end=dnase[,3]),
                         strand=Rle("*"),
                         data=dnase[,5])

    iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name="")
    gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
    httr::set_config(httr::config(ssl_verifypeer = 0L))
    rTrack <- try(UcscTrack(genome=gen, chromosome=chrom, track="NCBI RefSeq",
                            from=minbase, to=maxbase, trackType="GeneRegionTrack",
                            rstarts="exonStarts", rends="exonEnds", gene="name",
                            symbol="name2", transcript="name", strand="strand",
                            fill="darkblue",stacking="squish", name="RefSeq",
                            showId=TRUE, geneSymbol=TRUE),silent = T)

    annEPICOrd <- annEPICSub[order(annEPICSub$chr,annEPICSub$pos),]
    bValsOrd <- bVals[match(annEPICOrd$Name,rownames(bVals)),]
    cpgData <- GRanges(seqnames=Rle(annEPICOrd$chr),
                       ranges=IRanges(start=annEPICOrd$pos, end=annEPICOrd$pos),
                       strand=Rle("*"),
                       betas=bValsOrd)
    # extract data on CpGs in DMR
    cpgData <- subsetByOverlaps(cpgData, results.ranges[dmrIndex])

    # methylation data track
    methTrack <- DataTrack(range=cpgData, groups=targPairs$Sample_Group,genome = gen,
                           chromosome=chrom, ylim=c(-0.05,1.05), col=pal,
                           type=c("a","p"), name="DNA Meth.\n(beta value)",
                           background.panel="white", legend=TRUE, cex.title=0.8,
                           cex.axis=0.8, cex.legend=0.8)
    # CpG island track
    islandTrack <- AnnotationTrack(range=islandData, genome=gen, name="CpG Is.",
                                   chromosome=chrom,fill="darkgreen")

    # DNaseI hypersensitive site data track
    dnaseTrack <- DataTrack(range=dnaseData, genome=gen, name="DNAseI",
                            type="gradient", chromosome=chrom)

    # DMR position data track
    dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR",
                                chromosome=chrom,fill="darkred")

    tracks <- list(iTrack, gTrack, methTrack, dmrTrack, islandTrack, dnaseTrack,
                   rTrack)
    return(list(
        "tracks"=tracks,
        "minbase"=minbase,
        "maxbase"=maxbase

    ))
}
                        
