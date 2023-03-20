gb <- globalenv(); assign("gb", gb)

generateSeg <- function(a, b, c) {
    x <- conumee::CNV.fit(a, b, c)
    x <- conumee::CNV.bin(x)
    x <- conumee::CNV.detail(x)
    x <- conumee::CNV.segment(x)
    return(x)
}

customCNV <- function (Mset, samName= NULL, sex = NULL) {
  if(is.null(samName)){
    samName <- colnames(Mset)[1]
  }
    path <- file.path(path.package("mnp.v11b6"), "ext")
    chiptype <- minfi::annotation(Mset)[[1]]
    if (is.null(sex)) {
        Rset <- minfi::ratioConvert(Mset, what = "both", keepCN = TRUE)
        sex <- ifelse(MNPgetSex(Rset)$predictedSex == "M", "Male", "Female")
    }
    if (chiptype == "IlluminaHumanMethylationEPIC") {
        load(file.path(path,"conumee_annotation_EPIC_B6.2019-11-29.RData"))
        cndata <- conumee::CNV.load(Mset, samName)
        if (sex == "Male") {
            load(file.path(path,"CNanalysis6_conumee_REF_M.2018-09-19.RData"))
            return(generateSeg(cndata, refM_epic, annoEPICxy))
        } else {
            load(file.path(path,"CNanalysis6_conumee_REF_F.2018-09-19.RData"))
            return(generateSeg(cndata, refF_epic, annoEPICxy))
            }
    } else {
        load(file.path(path,"CNanalysis4_conumee_ANNO.vh20150715.RData"))
        cndata <- conumee::CNV.load(Mset, samName)
        if (sex == "Male") {
            load(file.path(path,"CNanalysis4_conumee_REF-M.vh20150715.RData"))
             return(generateSeg(cndata, refM.data, annoXY))
        } else {
            load(file.path(path,"CNanalysis4_conumee_REF-F.vh20150715.RData"))
             return(generateSeg(cndata, refF.data, annoXY))
        }
    }
}


SaveCNVplots <- function(targets, idatPath = NULL) {
  if (is.null(idatPath)) {
    idatPath <- getwd()
  }
  samplename_data <- as.character(targets$Sample_ID)
  sentrix.ids <- as.character(targets$SentrixID_Pos)
  samGroup <- as.character(targets$Type)
  addCols = NULL
  for (i in 1:length(sentrix.ids)) {
    samName <- samplename_data[i]
    sampleEpic <- sentrix.ids[i]
    imgName <- paste(samName, sampleEpic, "CNV.png", sep="_")
    imgPath <- file.path(".", "figures", "cnv")
    imgFi <- file.path(imgPath, imgName)
    if (!file.exists(imgFi)) {
      pathEpic <- file.path(idatPath, sampleEpic)
      RGsetEpic <- read.metharray(pathEpic, verbose = T, force = T)
      MsetEpic <- mnp.v11b6::MNPpreprocessIllumina(RGsetEpic, bg.correct = TRUE, normalize = "controls")
      x <- gb$customCNV(MsetEpic, samName, NULL)
      slot(x, 'detail', check = FALSE) <- NULL
      invisible(format(object.size(x), units = 'auto'))
      png(filename = imgFi, width = 1820, height = 1040)
      cnvColors <- c("lightblue", "royalblue", "darkblue", "lightgrey", "#ffa19c", "red", "darkred")
      CNV.genomeplot(x, cols = cnvColors)
      invisible(dev.off())
    }
    knitr::include_graphics(imgName)
  }
}


# The GDC transforms copy number values into segment mean--equal to log2(copy-number/ 2).
# https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/
writeSegTab <- function(segFile, targets, idatPath= NULL) {
    if(is.null(idatPath)){idatPath <- getwd()}
  if (!file.exists(segFile)) {
      samplename_data <- as.character(targets$Sample_ID)
      sentrix.ids <- as.character(targets$SentrixID_Pos)
      samGroup <-as.character(targets$Type)
      addCols = NULL
      for (i in 1:length(sentrix.ids)) {
        samName <- samplename_data[i]
        sampleEpic <- sentrix.ids[i]
        pathEpic <- file.path(idatPath, sampleEpic)
        RGsetEpic <- read.metharray(pathEpic, verbose = T, force=T)
        MsetEpic <- mnp.v11b6::MNPpreprocessIllumina(RGsetEpic, bg.correct = TRUE, normalize = "controls")
        if(i==1){addCols <-T}else(addCols <- F)
        x <- gb$customCNV(MsetEpic, samName)
        q <- CNV.write(x)
        yy <- data.frame(x@detail$ratio)
        colnames(yy) <- paste(samplename_data[i])
        q$ID <- paste(samplename_data[i])
        q$group <- paste(samGroup[i])
        write.table(q, file = segFile, append = T, quote=F, sep=",", col.names=addCols, row.names=F)
    }
  }
}

writeDetailTab<-function(segFile, targets) {
  samNam <- as.character(targets$Sample_ID)
  sentrix.ids <- as.character(targets$SentrixID_Pos)
  samGroup <-as.character(targets$Type)
  if (!file.exists(segFile)) {
    addCols = NULL
    for (i in 1:length(sentrix.ids)) {
      sampleEpic <- sentrix.ids[i]
        pathEpic <- file.path(getwd(), sampleEpic)
        RGsetEpic <- read.metharray(pathEpic, verbose = T, force=T)
        MsetEpic <- mnp.v11b6::MNPpreprocessIllumina(RGsetEpic, bg.correct = TRUE, normalize = "controls")
        if(i==1){addCols <-T}else(addCols <- F)
        
        x <- gb$customCNV(MsetEpic)
        q <- CNV.write(x, what='detail')
        q$sample <- paste(samNam[i])
        zz <- CNV.write(x, what='probes')
        #Chromosome  Start    End    Feature 204776850101_R02C01
        colnames(zz) <- c("Chromosome", "Start", "End", "Probe", "Value")
        q$ID <- paste(samNam[i])
        q$group <- paste(samGroup[i])
        zz$ID <- paste(samNam[i])
        zz$group <- paste(samGroup[i])
        probeFi <- paste(samNam[i],"probeFile.csv",sep = "_")
        write.table(q, file = segFile, append = T, quote=F, sep=",", col.names=addCols, row.names=F)
        write.table(zz, file = probeFi, append = F, quote=F, sep=",", col.names=T, row.names=F)
    }
  }
}

savePlotPdf <- function(cnData, plotName, plotTitle) {
    freqPlot <- suppressMessages(
        GenVisR::cnFreq(
            cnData,
            genome = "hg19",
            plotType = "frequency",
            plot_title = plotTitle,
            CN_Loss_colour = "#00B0EB",
            CN_Gain_colour = "#00A36D"
        )
    )
    labelTitle <- paste("###", plotName, '\n\n')
    cat(labelTitle)
    print(freqPlot)
    cat('\n\n')
    pdf(
        file = paste0(plotName, "_cnv.pdf"),
        height = 5,
        width = 25
    )
    freqPlot
    invisible(dev.off())
}

filterGrp <- function(cnData, typeGroup){
    subCn <- cnData[cnData$group==typeGroup,]
    return(subCn[,1:(ncol(subCn)-1)])
}

SilentLoadLib <- function(pkg){
    suppressPackageStartupMessages(library(
        verbose = F,
        warn.conflicts = F,
        quietly = T,
        character.only = T,
        package = pkg
    ))
}

GetFreqPlot <- function(cnData, plotChr, plotTitle){
    plotChr <- unique(cnData$chromosome)
    freqPlot <- GenVisR::cnFreq(
        cnData,
        plotType = "frequency",
        plotChr = plotChr,
        plot_title = paste(plotTitle, "Samples Copy Number Frequency"),
        CN_Loss_colour = "#FF0000",
        CN_Gain_colour = "#00A36D",
        facet_lab_size = 9
    )
    freqPlot[["theme"]][["strip.text"]][["angle"]] <- -90
    freqPlot[["theme"]][["strip.text"]]$margin <- ggplot2::margin(t=10,r=10,b=10,l=10,unit="pt")
    freqPlot[["theme"]][["strip.placement"]] <- 'outside'
    freqPlot[["theme"]][["panel.spacing"]] <- unit(0, "points")
    freqPlot[["theme"]][["strip.placement"]] = "outside"
    freqPlot[["facet"]][["shrink"]] <- FALSE
    freqPlot[["coordinates"]][["clip"]] <- "off"
    return(freqPlot)
}

GetFreqData <- function(cnData,plotChr,plotTitle){
    freqDat <-
        suppressMessages(
            GenVisR::cnFreq(
                cnData,
                plotType = "frequency",
                plotChr = plotChr,
                out = "data",
                plot_title = paste(plotTitle, "Samples Copy Number Frequency"),
                CN_Loss_colour = "#FF0000",
                CN_Gain_colour = "#00A36D",
                facet_lab_size = 5
            )
        )
    return(freqDat)
}

SaveCnvData <- function(freqDat, plotName){
    cnvDir <- file.path(getwd(), "cnv")
    if (!dir.exists(cnvDir)) {dir.create(cnvDir, recursive=T)}
    suppressWarnings(write.csv(
        freqDat,
        file = file.path(cnvDir, paste(plotName, "cnv.csv", sep = "_")),
        row.names = F,
        col.names = T
    ))
    return(cnvDir)
}


SavePlotPng <- function(cnvDir, plotName, freqPlot, fileEnd="cnv.png"){
    freqFi <- file.path(cnvDir, paste(plotName, fileEnd, sep = "_"))
    png(
        file = freqFi,
        height = 8,
        width = 15,
        units = "in",
        res = 200
    )
    freqPlot
    invisible(dev.off())
}


SaveLoadCnvs <- function(cnData,
                           plotName,
                           plotTitle,
                           plotChr = c(paste0("chr", 1:22)),
                           saveImg = T) {
    SilentLoadLib("ggplot2")
    SilentLoadLib("GenVisR")
    SilentLoadLib("grDevices")

    freqPlot <- gb$GetFreqPlot(cnData, plotChr, plotTitle)
    freqDat <- gb$GetFreqData(cnData, plotChr, plotTitle)
    cnvDir <- gb$SaveCnvData(freqDat, plotName)
    labelTitle <- paste("###", plotName, "Samples" , '\n\n')
    cat(labelTitle)
    print(freqPlot)
    cat('\n\n')
    if (saveImg == T) {
        gb$SavePlotPng(cnvDir, plotName, freqPlot)
    }    
}


SaveClusters <- function(seg_clust_file, segFile){
  set1Nam<- c("chrom","loc.start","loc.end", "seg.mean","ID")
  set2Nam<- c("chr","start","end", "value","ID"   )
    detailVals <- as.data.frame(read.csv(segFile, row.names=NULL))
    detailVals <- detailVals[,set1Nam]
    #detailVals <- detailVals[,set2Nam]
    colnames(detailVals) <- c("chromosome", "start", "end", "segmean", "sample")
    detailVals$segmean <- (2**(detailVals$segmean))*2
    write.table(detailVals,file=seg_clust_file,sep = "\t",row.names = F)
    cnData <- read.delim(seg_clust_file,header = T,sep = "\t",row.names=NULL)
    return(cnData)
}


grabClusterDat <- function(seg_clust_file,segFile){
     if(!file.exists(seg_clust_file)){
        cnData <- SaveClusters(seg_clust_file,segFile)} else{
        cnData <- read.delim(seg_clust_file,header = T,sep = "\t",row.names=NULL)
        }
    cnData$start <- as.numeric(cnData$start)
    cnData$end <- as.numeric(cnData$end)
    return(cnData)
}

grabGsetFun <- function(gsetFile, RGSet, targets) {
  if (!file.exists(file.path(getwd(), gsetFile))) {
    gset.funnorm <- gb$cleanUpProbes(RGSet, targets, getfunorm = T)
    saveRDS(gset.funnorm, file = gsetFile)
  } else{
    gset.funnorm <- readRDS(gsetFile)
  }
  return(gset.funnorm)
}

grabGsetBeta <- function(gsetbeta, gset.funnorm) {
  if (!file.exists(file.path(getwd(), gsetbeta))) {
    gset.funnorm.beta <- supM(minfi::getBeta(gset.funnorm))
    saveRDS(gset.funnorm.beta, file = gsetbeta)
  } else{
    gset.funnorm.beta <- readRDS(gsetbeta)
  }
  return(gset.funnorm.beta)
}


#colnames(gset.funnorm.beta) <- gset.funnorm$Sample_Group
getDmpData <- function(ClusfiNam, gset.funnorm, condition){
if(!file.exists(file.path(getwd(),ClusfiNam))){
  annot = minfi::getAnnotation(gset.funnorm)
  gset.funnorm.beta <- gb$grabGsetBeta(gBetaFile , gset.funnorm)
  dmp <- minfi::dmpFinder(gset.funnorm.beta, pheno = condition, type = "categorical")
  dmp <- cbind(dmp, ID = rownames(dmp))
  dmp_annot_combined <- cbind(annot[row.names(dmp), ], dmp)
  dmp_annot_combined_df <- as.data.frame(dmp_annot_combined)
  dmp_annot_significant <- subset(dmp_annot_combined_df, dmp_annot_combined_df$qval < 0.05)
  write.csv(dmp_annot_significant,file=ClusfiNam)
  theDmpData <- read.csv(ClusfiNam,header = T,sep = ",")
  }else{
    theDmpData <- read.csv(ClusfiNam,header = T,sep = ",")
    }
  gene_char_split <- strsplit(as.character(theDmpData$UCSC_RefGene_Name),';')
  gene_char_unlist <- unlist(gene_char_split)
  gene_char_unique <- as.data.frame(unique(gene_char_unlist))
  colnames(gene_char_unique) <- 'Genes_By_Sample'
  return(gene_char_unique)
}


##Convert to entrz ids for kegg enrichment analysis##
entrz2kegg <- function(gene_char_unique){
  ids <- supM(clusterProfiler::bitr(gene_char_unique$Genes_By_Sample, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db"))
  kk <- supM(clusterProfiler::enrichKEGG(gene = ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05, universe=NULL))
  enrichplot::dotplot(kk)
  kk_final <- DOSE::setReadable(kk, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID")
  return(kk_final)
}


# Removes repeated gene names in Manifest annotations
removeDupeGene <- function(newOvGene) {
  if(!file.exists("newOvGene.rds")){
  newOvGene2 <- unlist(mclapply(
    X = 1:nrow(newOvGene),
    FUN = function(rw) {
      rowGene <- stringr::str_split(newOvGene[rw, 1], ";")
      if (rowGene == "") {rowGene = "NO.GENE.NAME"}
      newRow <- paste(unique(rowGene[[1]]), collapse = "/")
      return(newRow)
    }
  ))
  saveRDS(newOvGene2, "newOvGene.rds")
  }else{ 
    newOvGene2 <- readRDS("newOvGene.rds")
    }
  return(newOvGene2)
}


writeDetailTab<-function(segFile, targets) {
  samNam <- as.character(targets$Sample_ID)
  sentrix.ids <- as.character(targets$SentrixID_Pos)
  samGroup <-as.character(targets$Type)
  if (!file.exists(segFile)) {
    addCols = NULL
    for (i in 1:length(sentrix.ids)) {
      sampleEpic <- sentrix.ids[i]
        pathEpic <- file.path(getwd(), sampleEpic)
        RGsetEpic <- read.metharray(pathEpic, verbose = T, force=T)
        MsetEpic <- mnp.v11b6::MNPpreprocessIllumina(
          RGsetEpic, bg.correct = TRUE, normalize = "controls")
        if(i==1){addCols <-T}else(addCols <- F)
        x <- gb$customCNV(MsetEpic)
        #cnv_df<-getCnvPlot(x)
        q <- CNV.write(x, what='detail')
        q$sample <- paste(samNam[i])
        zz <- CNV.write(x, what='probes')
        #Chromosome  Start    End    Feature 204776850101_R02C01
        #c("chromosome", "start", "end",  "sample","segmean")
        colnames(zz) <- c("Chromosome", "Start", "End", "Probe", "Value")
        q$ID <- paste(samNam[i])
        q$group <- paste(samGroup[i])
        zz$ID <- paste(samNam[i])
        zz$group <- paste(samGroup[i])
        probeFi <- paste(samNam[i],"probeFile.csv",sep = "_")
        write.table(q, file = segFile, append = T, quote=F, sep=",",
                    col.names=addCols, row.names=F)
        write.table(zz, file = probeFi, append = F, quote=F, sep=",", col.names=T, row.names=F)
    }
  }
}


rmDupeAnnotation <- function(dupeRows) {
    dupeRows2 <- unlist(mclapply(
    X = 1:length(dupeRows),
    FUN = function(rw) {
      rowGene <- stringr::str_split(dupeRows[rw], ";")[[1]]
      if (rowGene == "") {rowGene = "NO.NAME"}
      newRow <- paste(unique(rowGene), collapse = "/")
      return(newRow)}))
  return(dupeRows2)
}


GenCNVdataGroupSave <- function(cnData, targets, col_samGrp, plotChr = NULL) {
  cnData$group <- NULL
  targets$Type <- targets[col_samGrp, ]
  for (sn in unique(cnData$sample)) {
    selec <- cnData$sample == sn
    cnData$group[selec] <- targets$Type[targets$Sample_Name == sn]
  }
  
  for (samgrp in unique(cnData$group)) {
    plotTitle <- paste0(samgrp)
    cnDat2 <- gb$filterGrp(cnData, samgrp)
    gb$SaveLoadCnvs(cnDat2, plotName = plotTitle, plotTitle, plotChr=plotChr)
  }
}

