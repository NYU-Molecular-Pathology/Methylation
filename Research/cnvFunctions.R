#!/usr/bin/env Rscript
## ---------------------------
## Script name: Params_output.R
## Purpose: Source global functions used for research report output copy number variation
## Date Created: May 17, 2022
## Date Last Modified: February 29, 2024
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

gb <- globalenv(); assign("gb", gb)

generateSeg <- function(a, b, c) {
    x <- conumee::CNV.fit(a, b, c)
    x <- conumee::CNV.bin(x)
    x <- conumee::CNV.detail(x)
    x <- conumee::CNV.segment(x)
    return(x)
}


customCNV <- function (Mset, samName = NULL, sex = NULL, customAnno = NULL) {
    gb <- globalenv()
    if(is.null(samName)){samName <- colnames(Mset)[1]}

    chiptype <- minfi::annotation(Mset)[[1]]
    Rset <- minfi::ratioConvert(Mset, what = "both", keepCN = TRUE)
    cndata <- conumee::CNV.load(Mset, samName)

    if (chiptype == "IlluminaHumanMethylationEPIC") {
        #require("mnp.v11b6")

        if (is.null(sex)) {
            sex <- ifelse(mnp.v11b6::MNPgetSex(Rset)$predictedSex == "M", "Male", "Female")
        }

        path <- file.path(path.package("mnp.v11b6"), "ext")
        load(file.path(path,"conumee_annotation_EPIC_B6.2019-11-29.RData"), envir = gb)
        load(file.path(path,"CNanalysis6_conumee_REF_M.2018-09-19.RData"), envir = gb)
        load(file.path(path,"CNanalysis6_conumee_REF_F.2018-09-19.RData"), envir = gb)
        sexRefData <- if(sex == "Male"){gb$refM_epic} else{gb$refF_epic}
        mainAnno <- gb$annoEPICxy
    } else {
        require("mnp.v11b4")

        if (is.null(sex)) {
            sex <- ifelse(mnp.v11b4::MNPgetSex(Rset)$predictedSex == "M", "Male", "Female")
        }

        path <- file.path(path.package("mnp.v11b4"), "ext")
        load(file.path(path,"CNanalysis4_conumee_ANNO.vh20150715.RData"), envir = gb)
        load(file.path(path,"CNanalysis4_conumee_REF-M.vh20150715.RData"), envir = gb)
        load(file.path(path,"CNanalysis4_conumee_REF-F.vh20150715.RData"), envir = gb)
        sexRefData <- if(sex == "Male"){gb$refM.data}else{gb$refF.data}
        mainAnno <- gb$annoXY
    }

    if(!is.null(customAnno)){
        return(generateSeg(cndata, sexRefData, customAnno))
    } else{
        return(generateSeg(cndata, sexRefData, mainAnno))
    }

}


SaveCNVplots <- function(samplename_data, sentrix.ids, i, idatPath = NULL) {
  if (is.null(idatPath)) {idatPath <- getwd()}
  sample_id <- samplename_data[i]
  idatBase <- sentrix.ids[i]
  imgName <- paste(sample_id, idatBase, "CNV.png", sep = "_")
  imgFi <- file.path(file.path(".", "figures", "cnv"), imgName)
  if (!file.exists(imgFi)) {
    samRgset <- minfi::read.metharray(basenames = file.path(idatPath, idatBase), verbose = F, force = T)
    MsetEpic <- mnp.v11b6::MNPpreprocessIllumina(samRgset, bg.correct = TRUE, normalize = "controls")
    x <- gb$customCNV(MsetEpic, sample_id, NULL)
    slot(x, 'detail', check = FALSE) <- NULL
    invisible(format(object.size(x), units = 'auto'))
    png(filename = imgFi, width = 1820, height = 1040)
    cnvColors <- c("lightblue", "royalblue", "darkblue", "lightgrey", "#ffa19c", "red", "darkred")
    conumee::CNV.genomeplot(x, cols = cnvColors)
    invisible(dev.off())
  }
  cat(paste("###", sample_id, "\n\n"))
  postImg <- paste0("![", sample_id, "](", imgFi, ")")
  cat("\n\n")
  cat(postImg)
  cat("\n\n")
}


SetPrintGrid <- function(sentrix.ids) {
  if (length(sentrix.ids) %% 2 == 0) {
    xRow <- length(sentrix.ids)
  } else{
    xRow <- length(sentrix.ids) + 1
  }
  return(par(mfrow = c(xRow, 2)))
}

# The GDC transforms copy number values into segment mean--equal to log2(copy-number/ 2).
# Any copy number gain from 10 and above is an amplification
# https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/
GetGeneCopyNumb <- function(cnv_data, samName){
    cnvGeneRatio <- data.frame(cnv_data@detail$ratio)
    colnames(cnvGeneRatio) <- "log2_ratio"
    cnvGeneRatio$Gene <- rownames(cnvGeneRatio)
    cnvGeneRatio$CopyNumber <- round((2**(cnvGeneRatio[,1]))*2, 2)
    cnvGeneRatio$GainLoss <- "NONE"
    cnvGeneRatio$GainLoss[cnvGeneRatio$CopyNumber > 2 ] <- "Gain"
    cnvGeneRatio$GainLoss[cnvGeneRatio$CopyNumber <= 1] <- "Loss"
    cnvGeneRatio$GainLoss[cnvGeneRatio$CopyNumber >= 10 ] <- "Amplification"
    cnvGeneRatio$ID <- samName
    cnvGeneRatio <- cnvGeneRatio[, c("ID", "Gene", "log2_ratio", "CopyNumber", "GainLoss")]
    return(cnvGeneRatio)
}


SaveHtmlCnv <- function(samName, cnvPath, cnv_data) {
    htmlFi <- paste(samName, "cnv.html", sep = "_")
    htmlOut <- file.path(cnvPath, htmlFi)
    plotObj <- suppressMessages(mnp.v11b6::MNPcnvggplotly(cnv_data, getTables = F))
    suppressWarnings(htmlwidgets::saveWidget(plotObj, htmlFi, selfcontained = T))
    file.rename(htmlFi, htmlOut)
}

GrabMsetData <- function(idatBasename){
    rgSet <- minfi::read.metharray(idatBasename, verbose = T, force = T)
    if(rgSet@annotation[["array"]] == "IlluminaHumanMethylationEPIC"){
        #require("mnp.v11b6")
        Mset <- minfi::preprocessIllumina(rgSet, normalize = "controls")
    }else{
        require("mnp.v11b4")
        Mset <- minfi::preprocessIllumina(rgSet, normalize = "controls")
    }
    return(Mset)
}


GetCNVObj <- function(segFile, targets, i, idatPath, custom_anno) {
    segFile <- basename(segFile)
    samplename_data <- as.character(targets[,1])
    sentrix.ids <- as.character(targets$SentrixID_Pos)
    
    samGroup <- as.character(targets$Type)
    
    segPath <- file.path(gb$runDir, "data", "CNV_segments")
    cnvPath <- file.path(gb$runDir, "figures", "cnv", "interactive_cnv")
    
    #if(!dir.exists(segPath)){dir.create(segPath)}
    #if(!dir.exists(cnvPath)){dir.create(cnvPath)}
    
    samName <- samplename_data[i]
    sample_id <- sentrix.ids[i]
    idatBasename <- file.path(idatPath, sample_id)
    currFi <- file.path(segPath, paste0(samName, "_", segFile))
    geneFi <- file.path(segPath, paste("genes", samName, segFile, sep = "_"))
    
    if (file.exists(currFi) & file.exists(geneFi)) {
        message("Files already exist! Skipping:\n", currFi, "\n", geneFi)
    } else{
        Mset <- GrabMsetData(idatBasename)
        cnv_data <- gb$customCNV(Mset, samName, customAnno = custom_anno)
        cnvObjOut <- conumee::CNV.write(cnv_data)
        cnvObjOut$ID <- paste(samName)
        cnvObjOut$group <- paste(samGroup[i])
        cnvGeneRatio <- GetGeneCopyNumb(cnv_data, samName)
        write.table(cnvObjOut, file = currFi, quote=F, sep=",", row.names=F)
        write.table(cnvGeneRatio, file = geneFi, quote=F, sep=",", row.names=F)
        SaveHtmlCnv(samName, cnvPath, cnv_data)
    }
}


writeSegTab <- function(segFile = NULL, targets = NULL, idatPath = NULL, custom_anno = NULL) {
    if(is.null(segFile)){segFile <- paste0(format(Sys.Date(),"%b%d"), "_segfile.csv")} #_segmentVals.csv
    if(is.null(idatPath)){idatPath <- file.path(gb$runDir, "idats")}
    if(is.null(targets)){targets <- as.data.frame(read.csv("samplesheet.csv"))}
    stopifnot(any(stringr::str_detect(colnames(targets),"SentrixID_Pos")))
    
    if(!any(stringr::str_detect(colnames(targets),"Type"))){
        targets$Type <- "Sample_Group_1"
    }
    
    sentrix.ids <- as.character(targets$SentrixID_Pos)
    
    for (i in 1:length(sentrix.ids)) {
        GetCNVObj(segFile, targets, i, idatPath, custom_anno)
    }
    
    cnvSegDir <- file.path(gb$runDir, "data", "CNV_segments")
    segPath <- cnvSegDir
    segFileAll <- file.path(segPath, basename(segFile))
    
    if(!file.exists(segFileAll)){
        allSegFiles <- dir(path = segPath, pattern = "_segfile.csv", full.names = TRUE)
        cnvSegFiles <- allSegFiles[which(stringr::str_detect(basename(allSegFiles), "genes_", TRUE))]
        allSegs <- vroom::vroom(cnvSegFiles, col_types = paste(rep("c", 10), collapse = ""))
        write.csv(allSegs, file = segFileAll, quote = F, row.names = F)
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
    cnvDir <- file.path(gb$runDir, "figures", "cnv", "data")
    if (!dir.exists(cnvDir)) {dir.create(cnvDir, recursive=T)}
    suppressWarnings(write.csv(
        freqDat,
        file = file.path(cnvDir, paste(plotName, "cnv.csv", sep = "_")),
        row.names = F, col.names = T))
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


SaveLoadCnvs <- function(cnData, plotName = "All", plotTitle = "All", plotChr = c(paste0("chr", 1:22)), saveImg = T, isGrp=F) {
    SilentLoadLib("ggplot2")
    SilentLoadLib("GenVisR")
    SilentLoadLib("grDevices")
    if(isGrp==FALSE){
      cat("\n\n")
      cat(paste("## Sample", as.character(gb$col_samGrp), "Groups CNV Frequency {.tabset}"))
      cat("\n\n")
    }
    freqPlot <- suppressMessages(gb$GetFreqPlot(cnData, plotChr, plotTitle))
    freqDat <- suppressMessages(gb$GetFreqData(cnData, plotChr, plotTitle))
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
    set1Nam <- c("chrom", "loc.start", "loc.end", "seg.mean", "ID")
    set2Nam <- c("chr", "start", "end", "value", "ID")
    detailVals <- as.data.frame(read.csv(segFile, row.names=NULL))
    detailVals <- detailVals[, set1Nam]
    #detailVals <- detailVals[,set2Nam]
    colnames(detailVals) <- c("chromosome", "start", "end", "segmean", "sample")
    detailVals$CopyNumber <- round((2**(detailVals$segmean))*2, 1)
    write.table(detailVals, seg_clust_file, sep = "\t", row.names = F)
    return(cnData)
}


grabClusterDat <- function(seg_clust_file, segFile, targets, idatPath){
    writeSegTab(segFile, targets, idatPath)
    if(!file.exists(seg_clust_file)) {
        cnData <- SaveClusters(seg_clust_file, segFile)
    } else{
        cnData <- read.delim(seg_clust_file, header = T, sep = "\t", row.names=NULL)
    }
    cnData$start <- as.numeric(cnData$start)
    cnData$end <- as.numeric(cnData$end)
    cnData <- cnData[,-6]
    return(cnData)
}


grabGsetFun <- function(RGSet, targets, gb) {
  gSetFunData <- gb$gsetFile
  if (!file.exists(gSetFunData)) {
    gset.funnorm <- gb$cleanUpProbes(RGSet, targets, gb, getfunorm = T)
    gb$SaveObj(gset.funnorm, file.name = gSetFunData)
    } else{
      gset.funnorm <- gb$LoadRdatObj(gSetFunData, F)
      }
  return(gset.funnorm)
}


grabGsetBeta <- function(gsetbeta, gset.funnorm) {
    gsetbetaFi <- gb$gBetaFile  #file.path(getwd(), gsetbeta)
  if (!file.exists(gsetbetaFi)) {
    gset.funnorm.beta <- supM(minfi::getBeta(gset.funnorm))
    gb$SaveObj(gset.funnorm.beta, file.name = gsetbetaFi)
  } else{
    gset.funnorm.beta <- gb$LoadRdatObj(gsetbetaFi, F)
  }
  return(gset.funnorm.beta)
}


DedupeUniq <- function(StrObj){
  return(as.data.frame(unique(unlist(strsplit(as.character(StrObj),';')))))
}

getDmpData <- function(ClusfiNam, gset.funnorm, condition, gb) {
    annot = minfi::getAnnotation(gset.funnorm)
    gset.funnorm.beta <- gb$grabGsetBeta(gb$gBetaFile , gset.funnorm)
    dmp <- minfi::dmpFinder(gset.funnorm.beta, pheno = condition, type = "categorical")
    dmp <- cbind(dmp, ID = rownames(dmp))
    dmp_annot_combined <- cbind(annot[row.names(dmp),], dmp)
    dmp_annot_combined_df <- as.data.frame(dmp_annot_combined)
    dmpAnnoSign <- subset(dmp_annot_combined_df, dmp_annot_combined_df$qval < 0.2)
    dmpAnnoSign <- subset(dmpAnnoSign, dmpAnnoSign$pval < 0.05)
    write.csv(dmpAnnoSign, file = ClusfiNam, quote = F)
    theDmpData <- read.csv(ClusfiNam, header = T, sep = ",")
    gene_char_unique1 <- DedupeUniq(theDmpData$GencodeCompV12_Accession)
    gene_char_unique2 <- DedupeUniq(theDmpData$UCSC_RefGene_Accession)
    gene_char_unique3 <- DedupeUniq(theDmpData$GencodeCompV12_NAME)
    gene_char_unique4 <- DedupeUniq(theDmpData$UCSC_RefGene_Name)
    gene_char_unique3 <- c(gene_char_unique3[, 1], gene_char_unique4[, 1])
    gene_char_unique <- list(
        GencodeAccession = gene_char_unique1[, 1],
        UCSCrefseq = gene_char_unique2[, 1],
        GeneNames = gene_char_unique3
        )
  return(gene_char_unique)
}



##Convert to entrz ids for kegg enrichment analysis##
# entrz2kegg <- function(gene_char_unique){
#   ids <- supM(clusterProfiler::bitr(gene_char_unique$Genes_By_Sample, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db"))
#   kk <- supM(clusterProfiler::enrichKEGG(gene = ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05, universe=NULL))
#   enrichplot::dotplot(kk)
#   kk_final <- DOSE::setReadable(kk, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID")
#   return(kk_final)
# }


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


writeDetailTab <- function(segFile, targets) {
  samNam <- as.character(targets$Sample_ID)
  sentrix.ids <- as.character(targets$SentrixID_Pos)
  samGroup <-as.character(targets$Type)
  if (!file.exists(segFile)) {
    addCols = NULL
    for (i in 1:length(sentrix.ids)) {
      sampleEpic <- sentrix.ids[i]
        pathEpic <- file.path(gb$runDir, sampleEpic)
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
  targets$Type <- targets[, col_samGrp]
  cat('\n\n')
  cat(paste("## Sample", col_samGrp, "Groups CNV Frequency {.tabset}"))
  cat('\n\n')
  for (sn in unique(cnData$sample)) {
    selec <- cnData$sample == sn
    cnData$group[selec] <- targets$Type[targets$Sample_Name == sn]
  }
  for (samgrp in unique(cnData$group)) {
    plotTitle <- paste0(samgrp)
    cnDat2 <- gb$filterGrp(cnData, samgrp)
    gb$SaveLoadCnvs(cnDat2, plotName = plotTitle, plotTitle, plotChr = plotChr, isGrp = T)
  }
    cat('\n\n')
}


LoopSavePlainCNV <- function(targets) {
    cat("\n\n")
    cat('## Sample CNV Plots {.tabset}')
    samplename_data <- as.character(targets$Sample_ID)
    sentrix.ids <- as.character(targets$SentrixID_Pos)
    cat("\n\n")
    for (i in 1:length(sentrix.ids)) {
        gb$SaveCNVplots(samplename_data, sentrix.ids, i, gb$idatPath)
        cat("\n\n")
    }
}


SaveChromoCNV <- function(sampleImg, x, chrNum) {
    imgName <- paste(sampleImg, "CNV",chrNum, ".png", sep = "_")
    imgFi <- file.path(".", "figures", "cnv", imgName)
    png(filename = imgFi, width = 1820, height = 1040, res=150)
    #cnvColors <- c("lightblue", "royalblue", "darkblue", "lightgrey", "#ffa19c", "red", "darkred")
    conumee::CNV.genomeplot(x, chr = chrNum)
    invisible(dev.off())
}

SaveAllChromoCNV <- function(x, sampleImg, doXY=F){
    chrAll <- paste0("chr", 1:22)
    if(doXY==T){
        chrAll <- "all"
    }
    imgName <- paste(sampleImg, "CNV_AllChr.png", sep = "_")
    imgFi <- file.path(".", "figures", "cnv", imgName)
    png(filename = imgFi, width = 1820, height = 1040, res=150)
    #cnvColors <- c("lightblue", "royalblue", "darkblue", "lightgrey", "#ffa19c", "red", "darkred")
    conumee::CNV.genomeplot(x, chr = chrAll)
    invisible(dev.off()) 
}

# https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/Research/cnvFunctions.R
SaveCNVplotsCustom <- function(samplename_data, sentrix.ids, i, idatPath = NULL, chrNum=NULL, doXY=F) {
    if (is.null(idatPath)) {
        idatPath <- getwd()
    }
    samName <- samplename_data[i]
    sampleEpic <- sentrix.ids[i]
    sampleImg <- paste(samName, sampleEpic, sep = "_")
    if(!file.exists(sampleImg)){
        pathEpic <- file.path(idatPath, sampleEpic)
        RGsetEpic <- read.metharray(pathEpic, verbose = T, force = T)
        MsetEpic <- mnp.v11b6::MNPpreprocessIllumina(RGsetEpic, bg.correct = TRUE, normalize = "controls")
        x <- gb$customCNV(MsetEpic, samName, NULL)
        slot(x, 'detail', check = FALSE) <- NULL
        invisible(format(object.size(x), units = 'auto'))
        if(!is.null(chrNum)){
            for(nChromo in chrNum){
                SaveChromoCNV(sampleImg, x, chrNum= paste("chr", nChromo)) 
            }
        }
        SaveAllChromoCNV(x, sampleImg, doXY)
    }
}


LoopSavePlainCNV2 <- function(targets) {
    cat("\n\n")
    cat('## Sample CNV Plots {.tabset}')
    samplename_data <- as.character(targets[,1])
    sentrix.ids <- as.character(targets$SentrixID_Pos)
    cat("\n\n")
    for (i in 1:length(sentrix.ids)) {
        SaveCNVplotsCustom(samplename_data, sentrix.ids, i, gb$idatPath)
        cat("\n\n")
    }
}

