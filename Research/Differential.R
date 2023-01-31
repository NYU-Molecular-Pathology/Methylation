gb <- globalenv(); assign("gb", gb)


SilentLoadLib <- function(pkg){
    stopifnot(suppressPackageStartupMessages(library(
        verbose = F,
        warn.conflicts = F,
        quietly = T,
        character.only = T,
        package = pkg,
        logical.return = T)))
}


LoadDifferentPacks <- function(outFigPath = NULL){
    if(is.null(outFigPath)){
        outFigPath <- file.path(getwd(),"figures","diffmean")
    }
    if(!dir.exists(outFigPath)){dir.create(outFigPath, recursive = T)}
    mealPkgs <- c("MEAL", "MultiDataSet",  "minfi", "ggplot2")
    otherPkg <- c(
        "RnBeads",
        "knitr",
        "limma",
        "minfi",
        "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        "IlluminaHumanMethylation450kmanifest",
        "RColorBrewer",
        "missMethyl",
        "minfiData",
        "Gviz",
        "DMRcate",
        "stringr",
        "methylationArrayAnalysis",
        "GenomicRanges"
    )
    librarian::shelf(mealPkgs, ask=F)
    #invisible(lapply(mealPkgs, SilentLoadLib))
    librarian::shelf(otherPkg, ask=F)
}


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
    listOfGenes <- geneDataLi[, "Gene"]
    X = 1:length(listOfGenes)
    geneValDf <- lapply(X, FUN = function(X, resAdj, listOfGenes) {
        geneName <- paste0(listOfGenes[X])
            theRow <- MEAL::getGeneVals(
                resAdj, geneName, genecol = "UCSC_RefGene_Name", fNames = c("chromosome", "start")
                )
            return(as.data.frame(theRow))
            }, resAdj, listOfGenes
        )
    
    newVals <- data.frame()
    for (x in geneValDf) {newVals <- rbind(newVals, x)}
    
    for (gn in 1:length(newVals$UCSC_RefGene_Name)) {
        newName <- unlist(stringr::str_split(newVals$UCSC_RefGene_Name[gn], pattern = ";"))
        newName <- paste(unique(newName), collapse = ";")
        newVals[gn, "UCSC_RefGene_Name"] <- newName
    }
    
    newVals$Region <- ""
    for (rgn in 1:length(geneDataLi$Gene)) {
        theName <- geneDataLi$Gene[rgn]
        for (xName in 1:length(newVals$UCSC_RefGene_Name)) {
            theGene <- unlist(stringr::str_split(newVals$UCSC_RefGene_Name[xName], ";"))
            if (any(grepl(pattern = theName, theGene)) == T) {
                newVals[xName, "Region"] <- geneDataLi$Site[rgn]
            }
        }
    }
    return(newVals)
}
                        
# DMRcate -------------------
                        
PlotQCPdf <- function(targPairs, detP,rgSet){
    pal <- brewer.pal(12,"Set3")
    legGroup <- levels(factor(targPairs$Sample_Group))
    par(mfrow=c(1,2))
    barplot(colMeans(detP), col=pal[factor(targPairs$Sample_Group)], las=2,
            cex.names=0.8, ylab="Mean detection p-values")
    abline(h=0.05,col="red")
    legend("topleft", legend=legGroup, fill=pal, bg="white")
    barplot(colMeans(detP), col=pal[factor(targPairs$Sample_Group)], las=2,
            cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
    abline(h=0.05,col="red")
    legend("topleft", legend=legGroup, fill=pal, bg="white")
    qcReport(rgSet, sampNames=targPairs$ID, sampGroups=targPairs$Sample_Group, pdf="qcReport.pdf")
}


GetGgMds <- function(mds,targPairs){
  scaleUp <- element_text(size = rel(1.2))
  toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Group = factor(targPairs$Sample_Group))
  theplot <- ggplot(toplot, size = 5, aes(Dim1, Dim2, colour = Group)) + geom_point(size = 3) + theme_bw() +
    theme( axis.title = scaleUp, legend.text = scaleUp, legend.title = scaleUp)
  return(theplot)
}
                        
                        
PlotDensityMds <- function(targPairs, mSetSq) {
  par(mfrow = c(1, 1))
  pal <- RColorBrewer::brewer.pal(8, "Dark2")
  legFact <- levels(factor(targPairs$Sample_Group))
  facPalCol <- pal[factor(targPairs$Sample_Group)]
  methylSet <- getM(mSetSq)
  cat("\n\n")
  cat("## MDS 1")
  cat("\n\n")
  mds1 <- plotMDS(methylSet, top = 1000, gene.selection = "common", col = facPalCol,
    dim = c(1, 2)
  )
  
  GetGgMds(mds1,targPairs)
  cat("\n\n")
  cat("## MDS PCA 1 & 3")
  cat("\n\n")
  mds2 <- plotMDS(methylSet, top = 1000, gene.selection = "common", col = facPalCol,
    dim = c(1, 3)
  )
  
  GetGgMds(mds2,targPairs)
  cat("\n\n")
  cat("## MDS PCA 2 & 3")
  cat("\n\n")
  mds3 <- plotMDS(methylSet, top = 1000, gene.selection = "common", col = facPalCol,
    dim = c(2, 3)
  )
  
  
  GetGgMds(mds3,targPairs)
  cat("\n\n")
}

                        
PlotDimensions <- function(mSetSqFlt,targPairs){
    par(mfrow=c(1,1))
    pal <- RColorBrewer::brewer.pal(8,"Dark2")
    legFact <- levels(factor(targPairs$Sample_Group))
    
    cat("\n\n")
    cat("## MDS Mset 1")
    cat("\n\n")
    
    plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
            col=pal[factor(targPairs$Sample_Group)], cex=0.8)
    legend("right", legend=legFact, text.col=pal, cex=0.65, bg="white")
    
    cat("\n\n")
    cat("## MDS Mset Sample Type")
    cat("\n\n")
    
    plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
            col=pal[factor(targPairs$Sample_Source)])
    legend("right", legend=levels(factor(targPairs$Sample_Source)), text.col=pal,
           cex=0.7, bg="white")
    
    cat("\n\n")
    cat("## mSetSqFlt Dim 1 & 3")
    cat("\n\n")
    
    # Examine higher dimensions to look at other sources of variation
    plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
            col=pal[factor(targPairs$Sample_Source)], dim=c(1,3))
    legend("right", legend=levels(factor(targPairs$Sample_Source)), text.col=pal,
           cex=0.7, bg="white")
    
    cat("\n\n")
    cat("## mSetSqFlt Dim 2 & 3")
    cat("\n\n")
    
    plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
            col=pal[factor(targPairs$Sample_Source)], dim=c(2,3))
    
    legend("topright", legend=levels(factor(targPairs$Sample_Source)), text.col=pal,
           cex=0.7, bg="white")
    
    cat("\n\n")
    cat("## mSetSqFlt Dim 3 & 4")
    cat("\n\n")
    
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
    targPairs$ID <- paste(targPairs$Sample_Group, targPairs$Sample_Name, sep = ".")
    toDrop <- RGSet@colData@listData[["Sample_Name"]] %in% targPairs$Sample_Name	
    rgSet <- RGSet[,toDrop]	
    sampleNames(rgSet) <- targPairs$ID
    detP <- detectionP(rgSet)
    PlotQCPdf(targPairs, detP, rgSet)
    keep <- colMeans(detP) < 0.05
    rgSet <- rgSet[, keep]
    targPairs <- targPairs[keep, ]
    detP <- detP[, keep]
    mSetSq <- preprocessQuantile(rgSet)
    mSetRaw <- preprocessRaw(rgSet)
    cat("\n\n")
    cat("# Pairwise Analysis MDS {.tabset}")
    cat("\n\n")
    gb$PlotDensityMds(targPairs, mSetSq)
    cat("\n\n")
    detP <- detP[match(featureNames(mSetSq), rownames(detP)), ]
    keep <- rowSums(detP < 0.01) == ncol(mSetSq)
    mSetSqFlt <- mSetSq[keep, ]
    annEPIC <- getAnnotation(rgSet)
    keep <- !(featureNames(mSetSqFlt) %in% annEPIC$Name[annEPIC$chr %in% c("chrX", "chrY")])
    mSetSqFlt <- mSetSqFlt[keep, ]
    mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
    return(list("targPairs" = targPairs, "mSetSqFlt" = mSetSqFlt))
}
                        
SubSetGroup <- function(strPatt, samPairs){
    idx <- which(grepl(samPairs, pattern = strPatt, ignore.case=T))
    return(samPairs[idx])
}
                        
GetTracksPlot <- function(annEPICSub, bVals, results.ranges, dmrIndex=1){
    library("GenomicRanges")
    library("knitr")
    library("limma")
    library("minfi")
    library("RColorBrewer")
    library("missMethyl")
    library("minfiData")
    library("Gviz")
    library("DMRcate")
    library("stringr")
    library("httr")
    library("methylationArrayAnalysis")
    data(cpgIslands)
    httr::set_config(httr::config(ssl_verifypeer = 0))
    httr::set_config(httr::config(ssl_verifyhost=0))
    #chr <- as.character(unique(seqnames(cpgIslands)))
    gen <- genome(cpgIslands)
    atrack <- AnnotationTrack(cpgIslands, name = "CpG")
    gen <- "hg19"
    pal <- brewer.pal(8,"Dark2")
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
     colnames(islandHMM) <- c("Chr", "start", "end", "length",	"CpGcount",	"GCcontent"	,"pctGC", "obsExp")
    islandHMM$start <- as.numeric(islandHMM$start)
    islandHMM$end <- as.numeric(islandHMM$end)
    #islandData<- makeGRangesFromDataFrame(islandHMM, keep.extra.columns = TRUE)
    
    islandData <- GRanges(seqnames=Rle(islandHMM[,1]),
                          ranges=IRanges(start=islandHMM[,2], end=islandHMM[,3]), strand=Rle(strand("*" )))
    
    dnase <- read.csv(paste0(dataDirectory,"/wgEncodeRegDnaseClusteredV3chr17.bed"),
                      sep="\t",stringsAsFactors=FALSE,header=FALSE)

    dnaseData <- GRanges(seqnames=dnase[,1],
                         ranges=IRanges(start=dnase[,2], end=dnase[,3]),
                         strand=Rle("*"),
                         data=dnase[,5])
    
    iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name="")
    gTrack <- GenomeAxisTrack(col="black", cex=1, name="GenomeAxis", fontcolor="black")
    rTrack <- gb$supM(UcscTrack(
        genome=gen, chromosome=chrom, track="UCSC Genes", from=minbase, to=maxbase, trackType="GeneRegionTrack", 
        start=minbase, end=maxbase, gene="name", symbol="name2", transcript="name", strand="strand", fill="darkblue",
        stacking="full", name="UCSC Genes", showId=TRUE, geneSymbol=TRUE)
                     )
    
    annEPICOrd <- annEPICSub[order(annEPICSub$chr,annEPICSub$pos),]
    bValsOrd <- bVals[match(annEPICOrd$Name,rownames(bVals)),]
    cpgData <- GRanges(seqnames=Rle(annEPICOrd$chr),
                       ranges=IRanges(start=annEPICOrd$pos, end=annEPICOrd$pos),
                       strand=Rle("*"),
                       betas=bValsOrd)
    # extract data on CpGs in DMR
    cpgData <- subsetByOverlaps(cpgData, results.ranges[dmrIndex])
    # methylation data track
    methTrack <- DataTrack(
        range=cpgData, groups=targPairs$Sample_Group,genome = gen, chromosome=chrom, 
        ylim=c(-0.05,1.05), col=pal, type=c("a","p"), name="DNA Meth.\n(beta value)",
        background.panel="white", legend=TRUE, cex.title=0.8, cex.axis=0.8, cex.legend=0.8
    )
    # CpG island track
    islandTrack <- AnnotationTrack(
        range=islandData, genome=gen, name="CpG Is.", chromosome=chrom, fill="darkgreen"
    )
    # DNaseI hypersensitive site data track
    dnaseTrack <- DataTrack(
        range=dnaseData, genome=gen, name="DNAseI", type="gradient", chromosome=chrom
    )
    # DMR position data track
    dmrTrack <- AnnotationTrack(
        start=start, end=end, genome=gen, name="DMR", chromosome=chrom, fill="darkred"
    )
    tracks <- list(iTrack, gTrack, methTrack, dmrTrack, islandTrack, dnaseTrack, rTrack)
    return(list("tracks"=tracks, "minbase"=minbase, "maxbase"=maxbase))
}


GetDesign <- function(targPairs) {
  cellType <- factor(targPairs$Sample_Group)
  individual <- factor(targPairs$Sample_Source)
  
  design <- model.matrix(~ 0 + cellType + individual, data = targPairs)
  colnames(design) <- c(levels(cellType), levels(individual)[-1])
  return(design)
}


ReadLoadDmps <- function(fit, contMat, annEPICSub, cateFile="./figures/diffmean/pairwise_DMPs.csv") {
  fit2 <- contrasts.fit(fit, contMat)
  fit2 <- eBayes(fit2)
  summary(decideTests(fit2))
  DMPs <- topTable(fit2, num = Inf, coef = 1, genelist = annEPICSub)
  
  if (!file.exists(cateFile)) {
      write.table(DMPs, file = cateFile, sep = ",", row.names = FALSE)
  }
  return(DMPs)
}

                                       
DrawFourCpg <- function(DMPs, bVals, targPhenotype){
    cat("## Top 4 most significantly differentially methylated CpG")
    cat("\n\n")
    par(mfrow = c(2, 2))
    par(mar=c(1,1,1,1))
    invisible(sapply(rownames(DMPs)[1:4], function(cpg) {
        minfi::plotCpg(bVals, cpg = cpg, pheno = targPhenotype, ylab = "Beta values")
    }))
    cat("\n\n")
}

GetDMRcolors <- function(targPairs){
    pal <- brewer.pal(8,"Dark2")
    groups <- pal[1:length(unique(targPairs$Sample_Group))]
    names(groups) <- levels(factor(targPairs$Sample_Group))
    cols <- groups[as.character(factor(targPairs$Sample_Group))]
    return(cols)
}

MakeDMRAnno <- function(targFit, contMat){
    targPairs <- targFit$targPairs
    mSetSqFlt <- targFit$mSetSqFlt
    design <- gb$GetDesign(targPairs)
    mVals <- minfi::getM(mSetSqFlt)
    myAnnotation <- DMRcate::cpg.annotate(
        object = mVals,
        datatype = "array",
        what = "M",
        analysis.type = "variability",
        design = design,
        contrasts = TRUE,
        cont.matrix = contMat,
        arraytype = "EPIC"
    )
    DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
    results.ranges <- DMRcate::extractRanges(DMRs)
    return(results.ranges)
}

DrawDMRAnno <- function(results.ranges, targPairs, bVals) {
    cols <- GetDMRcolors(targPairs)
    par(mfrow = c(1, 1))
    cat("## DMR.plot")
    cat("\n\n")
    DMRcate::DMR.plot(
        ranges = results.ranges,
        dmr = 2,
        CpGs = bVals,
        phen.col = cols,
        what = "Beta",
        arraytype = "EPIC",
        genome = "hg19"
    )
    cat("\n\n")
}
                        

GetEpicAnno <- function(RGSet, mVals){
  annEPIC <- getAnnotation(RGSet)
  annEPIC[match(rownames(mVals), annEPIC$Name), c(1:4, 12:19, 24:ncol(annEPIC))]
}
                        

DropDashes <- function(targCol){
   newValues <- stringr::str_replace_all(targCol, pattern="-", replacement = "_") 
   return(newValues)
}
                        
