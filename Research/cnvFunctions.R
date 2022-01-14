gb <- globalenv(); assign("gb", gb)

generateSeg <- function(a, b, c) {
    x <- conumee::CNV.fit(a, b, c)
    x <- conumee::CNV.bin(x)
    x <- conumee::CNV.detail(x)
    x <- conumee::CNV.segment(x)
    return(x)
}

customCNV <- function (Mset, sex = NULL) {
    path <- file.path(path.package("mnp.v11b6"), "ext")
    chiptype <- minfi::annotation(Mset)[[1]]
    if (is.null(sex)) {
        Rset <- minfi::ratioConvert(Mset, what = "both", keepCN = TRUE)
        sex <- ifelse(MNPgetSex(Rset)$predictedSex == "M", "Male", "Female")
    }
    if (chiptype == "IlluminaHumanMethylationEPIC") {
        load(file.path(path,"conumee_annotation_EPIC_B6.2019-11-29.RData"))
        cndata <- conumee::CNV.load(Mset)
        if (sex == "Male") {
            load(file.path(path,"CNanalysis6_conumee_REF_M.2018-09-19.RData"))
            return(generateSeg(cndata, refM_epic, annoEPICxy))
        } else {
            load(file.path(path,"CNanalysis6_conumee_REF_F.2018-09-19.RData"))
            return(generateSeg(cndata, refF_epic, annoEPICxy))
            }
    } else {
        load(file.path(path,"CNanalysis4_conumee_ANNO.vh20150715.RData"))
        cndata <- conumee::CNV.load(Mset)
        if (sex == "Male") {
            load(file.path(path,"CNanalysis4_conumee_REF-M.vh20150715.RData"))
             return(generateSeg(cndata, refM.data, annoXY))
        } else {
            load(file.path(path,"CNanalysis4_conumee_REF-F.vh20150715.RData"))
             return(generateSeg(cndata, refF.data, annoXY))
        }
    }
}

writeSegTab <- function(segFile, sentrix.ids) {
  if (!file.exists(segFile)) {
    addCols = NULL
    for (i in 1:length(sentrix.ids)) {
      sampleEpic <- sentrix.ids[i]
        pathEpic <- file.path(getwd(), sampleEpic)
        RGsetEpic <- read.metharray(pathEpic, verbose = T, force=T)
        MsetEpic <- mnp.v11b6::MNPpreprocessIllumina(RGsetEpic, bg.correct = TRUE, normalize = "controls")
        if(i==1){addCols <-T}else(addCols <- F)
        x <- gb$customCNV(MsetEpic)
        q <- CNV.write(x)
        yy <- data.frame(x@detail$ratio)
        colnames(yy) <- paste(samplename_data[i])
        q$ID <- paste(samplename_data[i])
        q$group <- paste(samGroup[i])
        write.table(q, file = segFile, append = T, quote=F, sep=",", col.names=addCols, row.names=F)
    }
  }
}

saveClusters <- function(seg_clust_file,segFile){
    detailVals <- as.data.frame(read.csv(segFile, row.names=NULL))[,c("chrom","loc.start","loc.end", "seg.mean","ID")]
    colnames(detailVals) <- c("chromosome", "start", "end", "segmean", "sample")
    detailVals$segmean <- (2**(detailVals$segmean))*2
    write.table(detailVals,file=seg_clust_file,sep = "\t",row.names = F)
    cnData <- read.delim(seg_clust_file,header = T,sep = "\t",row.names=NULL)
    return(cnData)
}


savePlotPdf <- function(cnData, plotName, plotTitle){
    freqPlot <- suppressMessages(
      GenVisR::cnFreq(cnData,genome = "hg19",plotType = "frequency", 
                      plot_title=plotTitle, CN_Loss_colour = "#00B0EB", 
                      CN_Gain_colour = "#00A36D"))
    labelTitle <- paste("###", plotName,'\n\n')
    cat(labelTitle)
    print(freqPlot)
    cat('\n\n')
    pdf(file=paste0(plotName,"_cnv.pdf"),height = 5,width = 25)
    freqPlot
    dev.off()
}


filterGrp <- function(cnData, typeGroup){
    subCn <- cnData[cnData$group==typeGroup,]
    return(subCn[,1:(ncol(subCn)-1)])
}

savePlotPng <- function(cnData, plotName, plotTitle, plotChr = c(paste0("chr",1:22)), saveImg=F){
    freqPlot <- suppressMessages(GenVisR::cnFreq(cnData, plotType = "frequency", plotChr = plotChr,
                                                 plot_title= paste(plotTitle,"Samples Copy Number Frequency"),
                                                 CN_Loss_colour = "#00B0EB", CN_Gain_colour = "#00A36D"))
    freqDat <- suppressMessages(GenVisR::cnFreq(cnData, plotType = "frequency", plotChr = plotChr, out = "data",
                                                plot_title= paste(plotTitle,"Samples Copy Number Frequency"),
                                                CN_Loss_colour = "#00B0EB", CN_Gain_colour = "#00A36D"))
    write.csv(freqDat, file = file.path(getwd(), paste(plotName, "cnv.csv", sep = "_")), row.names=F, col.names=T)
    if(saveImg==F){
        labelTitle <- paste("###", plotName, "Samples" ,'\n\n')
        cat(labelTitle)
        print(freqPlot)
        cat('\n\n')
    }else{
        png(file=file.path(getwd(), paste(plotName,"cnv.png", sep="_")), height = 8, width = 15, units = "in", res = 200)
        freqPlot
        dev.off()
    }
}

grabClusterDat <- function(seg_clust_file,segFile){
cnData <- if(!file.exists(seg_clust_file)){
  gb$saveClusters(seg_clust_file,segFile)} else{
    read.delim(seg_clust_file,header = T,sep = "\t",row.names=NULL)
    }
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


