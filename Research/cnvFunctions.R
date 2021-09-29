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

writeSegTab <- function(segFile) {
  if (!file.exists(segFile)) {
    addCols = NULL
    for (i in 1:length(sentrix.ids)) {
      sampleEpic <- sentrix.ids[i]
      pathEpic <- file.path(getwd(), sampleEpic)
      RGsetEpic <- minfi::read.metharray(pathEpic, verbose = T, force = T)
      MsetEpic <- mnp.v11b6::MNPpreprocessIllumina(RGsetEpic, bg.correct = TRUE, normalize = "controls")
      if (i == 1) {
        addCols <- T
      } else
        (addCols <- F)
      x <- customCNV(MsetEpic)
      q <- conumee::CNV.write(x)
      yy <- data.frame(x@detail$ratio)
      colnames(yy) <- paste(samplename_data[i])
      q$ID <- paste(samplename_data[i])
      q$group <- paste(samGroup[i])
      write.table(q, file = segFile, append = T, quote=F, sep=",", col.names=addCols, row.names=F)
    }
  }
}

saveClusters <- function(seg_clust_file,segFile){
    detailVals <- as.data.frame(read.csv(segFile, row.names=NULL))[,c("chrom","loc.start","loc.end", "seg.mean","ID","group")]
    colnames(detailVals) <- c("chromosome", "start", "end", "segmean", "sample","group")
    
    detailVals$segmean <- (2**(detailVals$segmean))*2
    sexChr <- detailVals$chromosome=="chrX" | detailVals$chromosome== "chrY"
    detailVals <- detailVals[!sexChr,]
    write.table(detailVals,file=seg_clust_file,sep = "\t",row.names = F)
    cnData <- read.delim(seg_clust_file,header = T,sep = "\t",row.names=NULL)
    return(cnData)
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
if(!file.exists(seg_clust_file)){
    cnData <- gb$saveClusters(seg_clust_file,segFile)
}else{
    cnData <- read.delim(seg_clust_file,header = T,sep = "\t",row.names=NULL)
}
return(cnData)
}
