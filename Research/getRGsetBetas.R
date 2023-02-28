gb <- globalenv(); assign("gb", gb)
library("utils")


# FUN: Save object as Rdata file ---------------------------------------------------------------
SaveObj <- function(object, file.name){
    outfile <- file(file.name, "wb")
    serialize(object, outfile)
    close(outfile)
}


# FUN: Load object from Rdata file with progressbar ---------------------------------------------
LoadRdatObj <- function(file.name, msgProg=T){
    library("foreach")
    library("utils")
    filesize <- file.info(file.name)$size
    chunksize <- ceiling(filesize / 100)
    pb <- txtProgressBar(min = 0, max = 100, style=3)
    infile <- file(file.name, "rb")
    if(msgProg==T){
    data <- foreach(it = icount(100), .combine = c) %do% {
        setTxtProgressBar(pb, it)
        readBin(infile, "raw", chunksize)
    }
    }else{
        data <- foreach(it = icount(100), .combine = c) %do% {readBin(infile, "raw", chunksize)}
    }
    close(infile)
    close(pb)
    return(unserialize(data))
}

# FUN: Checks Minfi Targets for Duplicate Arrays ------------------------------------------------
CheckDupeArrays <- function(arrayBase, arrayType){
    if (anyDuplicated(arrayBase) != 0) {
        duped <- paste(anyDuplicated(arrayBase), collapse=" ")
        message(paste(arrayType, "duplicated basenames:", duped))
        stopifnot(anyDuplicated(arrayBase) != 0)
    }
}


# FUN: Obtain RGSet with Probes common to 450K Array --------------------------------------------
combine.EPIC.450K <- function(targets, gb, batchFilter="EPIC") {
    stopifnot(any(gb$col_arrayType %in% colnames(targets)))
    targets$Batch <- targets[, gb$col_arrayType]
    message("Here are your targets$Batch options:\n", paste(unique(targets$Batch), collapse = " & "))
    stopifnot(any(batchFilter %in% targets$Batch) & length(targets$Batch) > 0)
    targets_450k <- subset(targets, targets$Batch != batchFilter)
    targets_EPIC <- subset(targets, targets$Batch == batchFilter)
    CheckDupeArrays(targets_450k, "450K")
    CheckDupeArrays(targets_EPIC, "EPIC")
    message("reading 450K arrays...")
    RGSet_450k <- minfi::read.metharray.exp(targets = targets_450k, force = T, verbose = T)
    message("reading EPIC arrays...")
    RGSet_EPIC <- minfi::read.metharray.exp(targets = targets_EPIC, force = T, verbose = T)
    message("Combining common probes...")
    RGSet <- minfi::combineArrays(RGSet_450k, RGSet_EPIC, outType = "IlluminaHumanMethylation450k")
    return(RGSet)
}

# Remove Low Quality and select probes using annotations
cleanUpProbes <- function(RGSet, targets, getfunorm=F){
    library("minfi")  
    td <- 
    dateFile <- paste(Sys.Date(),"DetPvals.Rdata", sep = "_")
    detPfile <- paste0("data", .Platform$file.sep, dateFile)
    if(!file.exists(detPfile)){
      detP <- minfi::detectionP(RGSet)
      gb$SaveObj(detP, file.name = detPfile)
    }else{
      detP <- gb$LoadRdatObj(detPfile)
    }
    colnames(detP) <- RGSet@colData@listData[["Sample_Name"]]
    keep <- colMeans(detP) < 0.05
    RGSet <- RGSet[, keep]
    targets <- targets[keep,]
    detP <- detP[, keep]
    dropping <- table(keep)["FALSE"]>0
    if(!is.na(dropping)){if(dropping==TRUE){message(
        "Dropping probes: ", paste0(capture.output(as.data.frame(table(keep))), collapse = "\n"))}}
    mSetSq <- suppressWarnings(preprocessQuantile(RGSet))
    detP <- detP[match(featureNames(mSetSq), rownames(detP)), ]
    keep <- rowSums(detP < 0.01) == ncol(mSetSq)
    gset.funnorm <- addSnpInfo(mSetSq[keep, ])
    gset.funnorm <- dropLociWithSnps(gset.funnorm,snps=c("SBE","CpG"),maf = 0)#drop loci with snps
    annot = getAnnotation(gset.funnorm) #getting annotation file
    sex_probes = annot$Name[annot$chr %in% c("chrX", "chrY")] # dropping sex probes
    gset.funnorm = gset.funnorm[!(rownames(gset.funnorm) %in% sex_probes),]
    if (getfunorm) {
        return(gset.funnorm)
    } else {
        betas <- minfi::getBeta(gset.funnorm)
        colnames(betas) <- RGSet@colData@listData[["Sample_Name"]]
        return(betas)
    }
}


# Returns topVar beta probes instead of all probes
takeTopVariance <- function(betas, topVar){
    var_probes <- apply(betas, 1.0, var)
    select_var <- names(sort(var_probes, decreasing = T))[topVar]
    top_var_beta <- betas[select_var, ]
    return(top_var_beta)
}


getSupervise <- function(the_beta, RGSet, topVar=1:10000, cutOff=0.05, dmpTyp = "categorical" ){
    condition <- pData(RGSet)$Type
    dmp <- dmpFinder(the_beta, pheno = condition, type = dmpTyp)
    dmp <- cbind(dmp, ID = rownames(dmp))
    betas_df <- as.data.frame(the_beta)
    qVals <- dmp$qval < cutOff
    if(table(qVals)[["TRUE"]]>=max(topVar)) {dmp <- dmp[qVals, ]}
    dmp <- dmp[dmp$pval < cutOff, ]
    final_sam <- row.names((t(betas_df[row.names(dmp[topVar, ]), ]))) #topVar=1:10000
    betas_df <- betas_df[, colnames(betas_df) %in% final_sam]
    write.csv(betas_df, "dmp_methylation.csv", quote = F)
    betas <- as.matrix(betas_df)
    return(betas)
}


# Checks if supervised rds data exists to load, else calculates it 
loadSupervise <- function(RGSet, betas, supbetaOut, varProbes, col_sentrix="SentrixID_Pos", dmpTyp = "categorical") {
  if (!file.exists(supbetaOut)) {
    rgRows <- RGSet@colData@rownames # ensure poor samples are dropped
    rgLiDat <- RGSet@colData@listData
    dropFilter <- rgLiDat[["Sample_ID"]] %in% colnames(betas)
    samFltr <- rgLiDat[[col_sentrix]][dropFilter]
    newRg <- RGSet[, rgRows %in% samFltr]
    topVar = 1:max(varProbes)
    # supervised dmp top Variance with getSupervise
    superbetas <- gb$getSupervise(betas, newRg, topVar, dmpTyp = dmpTyp) 
    SaveObj(superbetas, file.name = supbetaOut)
  } else{
    superbetas <- LoadRdatObj(supbetaOut)
  }
  return(superbetas)
}


# COMABAT batch correction between groups
batchCorrectBs <- function(betas, targets, batch_col) {
    targets$facilityBatch <- targets[, batch_col]
    batch <- targets$facilityBatch
    modcombat = model.matrix(~ 1, data = targets)
    gc(verbose = F)
    combat_beta <- sva::ComBat(
        dat = betas, batch = batch, mod = modcombat, par.prior = T, prior.plots = F
    )
    # fix if any beta >1 or < 0, so impute >1 is 1 and <0 is 0
    combat_beta <- ifelse(combat_beta < 0, 0, combat_beta)
    combat_beta <- ifelse(combat_beta > 1, 1, combat_beta)
    return(combat_beta)
}


RemoveBatchEffect <- function(betas, targets, gb) {
    if (file.exists(gb$combatOut)) {
      betas <- gb$LoadRdatObj(gb$combatOut)
      return(betas)
    }
    betas <- gb$batchCorrectBs(betas, targets, gb$batch_col)
    gb$SaveObj(betas, file.name = gb$combatOut)
    return(betas)
}


CleanUpSheetRows <- function(sheet, idatPath, targets){
    sheet <- as.data.frame(sheet)
    sheet$Barcode <-  targets[,gb$col_sentrix]
    sheet[,gb$col_sentrix] <- NA
    sheet[,gb$col_sentrix] <- sheet$Barcode
    splitSentrix <- stringr::str_split_fixed(sheet$Barcode, "_", 2)
    sheet$Slide <- splitSentrix[,1]
    sheet$Array <- splitSentrix[,2]
    sheet$Basename <- file.path(idatPath, paste0(sheet$Slide, "_" , sheet$Array))
    sheet$Sample_Name <- targets[, gb$col_samNames]
    shFiles <- c(paste0(sheet$Basename, "_Grn.idat"), paste0(sheet$Basename, "_Red.idat"))
    if(any(!file.exists(shFiles))){
      warning("Some idat files are missing and samples will be dropped!")
      message("The following cases are missing:")
      message(paste0(capture.output(sheet[!file.exists(shFiles),]), collapse="\n"))
      }
    sheet <- sheet[file.exists(shFiles),]
    toDrop <- anyDuplicated(sheet)
    sheet <- sheet[-toDrop,]
    sheet <- sheet[!is.na(sheet$Sample_Name), ]
    stopifnot(nrow(sheet) > 0)
    return(sheet)
}


getRgset <- function(rgOut, targets, mergeProbes = F, csvPath = "samplesheet.csv", idatPath = NULL){
  require("minfi")
  if(is.null(idatPath)){idatPath <- getwd()}
    gc(verbose = F)
    if (file.exists(rgOut)) {
      RGSet <- gb$LoadRdatObj(rgOut)
      }else{
        if (mergeProbes == T & !is.null(targets$Batch)) {
            RGSet <- gb$combine.EPIC.450K(targets = targets)
        } else{
          sheet <- NULL
          isPathway <- stringr::str_detect(csvPath, .Platform$file.sep)
          if (isPathway == T) {
              sheet <- minfi::read.metharray.sheet(dirname(csvPath), pattern = basename(csvPath))
          }else{
              sheet <- minfi::read.metharray.sheet(getwd(), pattern = csvPath)
          }
          sheet <- as.data.frame(sheet)
          sheet <- CleanUpSheetRows(sheet, idatPath, targets)
          sheet <- sheet[!is.na(sheet$Sample_Name),]
          RGSet <- minfi::read.metharray.exp(base = idatPath, targets = sheet, verbose = T, force = T)
          }
        gb$SaveObj(RGSet, file.name = rgOut)
      }
    return(RGSet)
}


GrabMinfiSheet <- function(idatPath, csvPath){
    isPathway <- stringr::str_detect(csvPath, .Platform$file.sep)
    if (isPathway == T) {
        sheet <- minfi::read.metharray.sheet(idatPath, pattern =  basename(csvPath))
    }else{
      if(!file.exists(file.path(idatPath, csvPath))){file.copy(csvPath, idatPath)}  
      sheet <- minfi::read.metharray.sheet(idatPath, pattern = csvPath)
    }
    sheet <- gb$CleanUpSheetRows(sheet, gb$idatPath, targets)
    return(sheet)
}



GetRgsetDat <- function(csvPath = "samplesheet.csv", gb){
    require("minfi")
    gc(verbose = F)
     if (file.exists(gb$rgOut)) {
        RGSet <- gb$LoadRdatObj(gb$rgOut)
        return(RGSet)
     }
    if(is.null(gb$idatPath)){gb$idatPath <- getwd()}
    sheet <- GrabMinfiSheet(gb$idatPath, csvPath)
    if (gb$mergeProbes == T & !is.null(gb$col_arrayType)) {
        RGSet <- gb$combine.EPIC.450K(targets = sheet, gb)
        gb$SaveObj(RGSet, file.name = gb$rgOut)
        return(RGSet)
    }
    RGSet <- minfi::read.metharray.exp(base = gb$idatPath, targets = sheet, verbose = T, force = T)
    gb$SaveObj(RGSet, file.name = gb$rgOut)
    return(RGSet)
}


cleanRawProbes <- function(RGSet, targets, gb) {
  gc(verbose = F)
  if (!file.exists(gb$rawBetaFi)) {
    betas <- gb$cleanUpProbes(RGSet = RGSet, targets = targets)
    SaveObj(betas, file.name = gb$rawBetaFi)
  } else{
    betas <- LoadRdatObj(gb$rawBetaFi)
  }
  if (!is.null(gb$batch_col)) {
    betas <- gb$RemoveBatchEffect(betas, targets, gb)
  }
  return(betas)
}


# Check if MDS data already loaded
getMdsPlot <- function(RGSet, samNames, topN=1000) {
    mSetSq <- preprocessQuantile(RGSet)
    mSetSq <- addSnpInfo(mSetSq)
    mSetSq <- dropLociWithSnps(mSetSq, snps = c("SBE", "CpG"), maf = 0) # drop the loci which has snps
    annot = getAnnotation(mSetSq)
    sex_probes = annot$Name[annot$chr %in% c("chrX", "chrY")]
    mSetSq = mSetSq[!(rownames(mSetSq) %in% sex_probes), ]
    mSetSq.beta <- minfi::getBeta(mSetSq)
    colnames(mSetSq.beta) <- mSetSq$Sample_Name
    assign("mSetSq.beta",mSetSq.beta, envir=gb)
    return(mSetSq.beta)
}


plot.mds <- function(mSetSq.beta, targets, varName, topN) {
    varColorCol <- paste0(varName, "_color")
    myColors <- c(targets[,varColorCol])
    names(myColors) <- targets$Sample_Name
    plotNam <- paste0(varName, "_top_", topN, "_msetBeta", ".png")
    outPlotNam <- file.path(getwd(),"figures","mds")
    if(!dir.exists(outPlotNam)){dir.create(outPlotNam, recursive = T)}
    outPlotFi <- file.path(outPlotNam, plotNam)
    png(filename = outPlotFi, width = 12, height = 8, res = 200, units = "in")
    limma::plotMDS(
      mSetSq.beta, top = topN, gene.selection = "common",
      plot = T, col = myColors, pch=19, labels=colnames(mSetSq.beta),
      main = paste("Top", topN, "Common", "mSet Sq.beta", "MDS plot")
    )
    names(myColors) <- targets[,varName]
    legend("topright", legend = c(unique(names(myColors))), 
           col = paste(as.list(unique(myColors))), pch = 15, cex = 0.8)
    invisible(dev.off())
    thepng <- outPlotFi #paste0('./', plotNam)
    return(knitr::include_graphics(thepng))
}


dropGroup <- function(targets, filterCol=NULL, group2rm=NULL) {
    remov <- ifelse(!is.null(group2rm), c(targets[,filterCol] == group2rm), targets[,filterCol] != group2rm)
    return(targets[!remov, ])
}


dropBadQc <- function(targets, betas, csvFi="samplesheet.csv") {
    dropping <- targets$Sample_Name %in% colnames(betas)
    if(table(dropping)[1]>=1){
        theMissing <- targets[!dropping, ]
        targets <- targets[dropping, ]
    }else{
        theMissing <- NULL
    }
    rownames(targets) <- 1:nrow(targets)
    write.csv(targets, file = csvFi, quote = F, row.names = F)
    return(theMissing)
}


# Matches RGset to any dropped samples of cleaned Beta Values
cleanRgset <- function(RGSet, col_sentrix, betas) {
  rgDat <- RGSet@colData@listData
  rgDat[[gb$col_sentrix]] <- rgDat[["Barcode"]]
  senDat <- rgDat[[gb$col_sentrix]]
  dropSam <- rgDat[["Sample_Name"]] %in% colnames(betas)
  rgRows <- RGSet@colData@rownames
  return(RGSet[, rgRows %in% senDat[dropSam]])
}


checkMdsRds <- function(mbfile, runDir, RGSet, targets) {
  if (file.exists(file.path(runDir, mbfile))) {
    mSetSq.beta <- LoadRdatObj(mbfile)
  } else{
    mSetSq.beta <- gb$supM(gb$getMdsPlot(RGSet, targets$Sample_ID, targets$Type))
    SaveObj(mSetSq.beta, file.name = file.path(gb$runDir, gb$mbfile))
  }
  return(mSetSq.beta)
}


GetMsetSq <- function(RGSet, targets, mbfile) {
  if (file.exists(mbfile)) {
    mSetSq.beta <- LoadRdatObj(mbfile)
  } else{
    mSetSq.beta <- gb$supM(gb$getMdsPlot(RGSet, targets$Sample_ID))
    SaveObj(mSetSq.beta, file = mbfile)
  }
  return(mSetSq.beta)
}


getTargCsv <- function(csvFi = "samplesheet.csv", betas=NULL, samName_col="Sample_Name") {
  gc(verbose = F)
  targets <- read.csv(csvFi, strip.white = T, row.names = NULL)
  targets <- dfTargets(targets)
  if(!is.null(betas)){
      targets <- targets[targets[,samName_col] %in% colnames(betas),]
  }
  stopifnot(nrow(targets) > 0)
  return(targets)
}


MatchRGtargets <- function(RGSet, targets, sampleSheet=NULL){
  if(is.null(sampleSheet)){sampleSheet <- "samplesheet.csv"}
  #file.path(getwd(),"csv","samplesheet.csv")}
  toDrop <- targets[,gb$col_sentrix] %in% RGSet@colData@rownames
  targets <- targets[toDrop,]
  write.csv(targets, sampleSheet, quote=F, row.names=F)
  targets <- gb$getTargCsv(sampleSheet)
  return(targets)
}


