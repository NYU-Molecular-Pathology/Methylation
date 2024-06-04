gb <- globalenv(); assign("gb", gb)
library("utils")


PasteCaptureData <- function(objDat) {
    return(paste0(capture.output(as.data.frame(objDat)), collapse = "\n"))
}


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
    
    infile <- file(file.name, "rb")
    if(msgProg==T){
        pb <- txtProgressBar(min = 0, max = 100, style=3)
        data <- foreach(it = icount(100), .combine = c) %do% {
            setTxtProgressBar(pb, it)
            readBin(infile, "raw", chunksize)
        }
        close(infile)
        close(pb)
    } else{
        data <- foreach(it = icount(100), .combine = c) %do% {readBin(infile, "raw", chunksize)}
        close(infile)
    }
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
    row.names(targets_450k) <- 1:nrow(targets_450k)
    row.names(targets_EPIC) <- 1:nrow(targets_EPIC)
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


MsgDropping <- function(keep){
    dropping <- table(keep)["FALSE"] > 0
    if (!is.na(dropping)) {
        if (dropping == TRUE) {
            drpMsg <- PasteCaptureData(table(keep))
            message("Dropping probes:\n", drpMsg)
        }
    }
}


DropSexSnpProbes <- function(detP, mSetSq, sexSnps=T){
    detP <- detP[match(featureNames(mSetSq), rownames(detP)), ]
    keep <- rowSums(detP < 0.01) == ncol(mSetSq)
    gset.funnorm <- addSnpInfo(mSetSq[keep, ])
    if(sexSnps==T){
      gset.funnorm <- minfi::dropLociWithSnps(gset.funnorm, snps = c("SBE", "CpG"), maf = 0)
      annot = minfi::getAnnotation(gset.funnorm)
      sex_probes = annot$Name[annot$chr %in% c("chrX", "chrY")]
      gset.funnorm = gset.funnorm[!(rownames(gset.funnorm) %in% sex_probes),]
      }
    return(gset.funnorm)
}


cleanUpProbes <- function(RGSet, targets, gb, getfunorm = F, getNoob=F, sexSnps = T){
    library("minfi")  
    detP <- minfi::detectionP(RGSet)
    samNameCol <- colnames(targets)[1]
    colnames(detP) <- RGSet@colData@listData[[samNameCol]]
    keep <- colMeans(detP) < 0.05
    RGSet <- RGSet[, keep]
    targets <- targets[keep,]
    detP <- detP[, keep]
    MsgDropping(keep)
    if(getNoob==T){
      mSetSq <- suppressWarnings(minfi::preprocessFunnorm(RGSet))
    }else{
      mSetSq <- suppressWarnings(minfi::preprocessQuantile(RGSet))
    }
    gset.funnorm <- DropSexSnpProbes(detP, mSetSq, sexSnps)
    if (getfunorm) {
        return(gset.funnorm)
    } 
    betas <- minfi::getBeta(gset.funnorm)
    colnames(betas) <- RGSet@colData@listData[[samNameCol]]
    return(betas)
}


# Returns topVar beta probes instead of all probes
takeTopVariance <- function(betas, topVar){
    var_probes <- apply(betas, 1.0, var)
    select_var <- names(sort(var_probes, decreasing = T))[topVar]
    top_var_beta <- betas[select_var, ]
    return(top_var_beta)
}


getSupervise <- function(the_beta, RGSet, topVar = 1:10000, cutOff = 0.05, dmpTyp = "categorical", superVar = NULL){
    condition <- pData(RGSet)[, superVar]
    stopifnot(length(condition) == ncol(the_beta))
    dmp <- minfi::dmpFinder(the_beta, pheno = condition, type = dmpTyp)
    dmp <- cbind(dmp, ID = rownames(dmp))
    betas_df <- as.data.frame(the_beta)
    qVals <- dmp$qval < cutOff
    if (any(qVals)) {
      if (table(qVals)[["TRUE"]] >= max(topVar)) {
        dmp <- dmp[qVals, ]
      }
    }
    dmp <- dmp[order(dmp$pval),]
    topDmp <- dmp[topVar, ]
    topProbes <- rownames(topDmp)
    final_sam <- row.names(t(betas_df[row.names(topDmp), ])) #topVar=1:10000
    betasDmp <- betas_df[topProbes, colnames(betas_df) %in% final_sam]
    dmpOutFi <- paste("top", max(topVar), superVar, "dmp_values.csv", sep = "_")
    dmpOutDir <- file.path(".", "figures", "csv")
    if(!dir.exists(dmpOutDir)){dir.create(dmpOutDir)}
    write.csv(topDmp, file.path(dmpOutDir, dmpOutFi) , quote = F)
    write.csv(betasDmp, file.path(dmpOutDir,paste0(superVar,"_dmp_betas.csv")), quote = F)
    betasDmp <- as.matrix(betasDmp)
    return(betasDmp)
}


# Checks if supervised rds data exists to load, else calculates it 
loadSupervise <- function(RGSet, betas, gb, superVar = NULL, dmpTyp = "categorical") {
    if (gb$supervisedRun == F) {
      return(NULL)
    }
    supbetaOutFi <- paste0(gb$supbetaOut, "_", superVar, ".Rdata")
    if (!file.exists(supbetaOutFi)) {
      rgRows <- RGSet@colData@rownames # ensure poor samples are dropped
      rgLiDat <- RGSet@colData@listData
      if(any(is.na(rgLiDat[["Sample_ID"]]))){
          rgLiDat[["Sample_ID"]] <- rgLiDat[["Sample_Name"]]
          dropFilter <- rgLiDat[["Sample_Name"]] %in% colnames(betas)
      }else{
          dropFilter <- rgLiDat[["Sample_ID"]] %in% colnames(betas)
      }
      samFltr <- rgLiDat[["SentrixID_Pos"]][dropFilter]
      newRg <- RGSet[, rgRows %in% samFltr]
      topVar = 1:max(gb$varProbes)
      # supervised dmp top Variance with getSupervise
      superbetas <- gb$getSupervise(betas, newRg, topVar, dmpTyp = dmpTyp, superVar = superVar)
      SaveObj(superbetas, file.name = supbetaOutFi)
    } else{
      superbetas <- LoadRdatObj(supbetaOutFi)
    }
    return(superbetas)
}


LoopSuperviseCategory <- function(gb, RGSet, betas) {
  if(gb$supervisedRun){
      for(i in 1:length(gb$selectedVars)){
          superbetas <- gb$loadSupervise(RGSet, betas, gb, gb$selectedVars[i], dmpTyp = "categorical")
          varName <- paste0("superbetas", i)
          assign(x = varName, value = superbetas, envir = gb)
      }
  }
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
    rownames(targets) <- targets[,1]
    if(is.null(gb$col_sentrix)){
        gb$col_sentrix <- "SentrixID_Pos"
    }
    if(is.null(gb$col_samNames)){
        if(!any(stringr::str_detect(colnames(targets),"Sample_Name"))){
            targets$Sample_Name <- targets[,1]
        }
        gb$col_samNames <- "Sample_Name"
    }
    sheetNames <- sheet[,1]
    targets <- targets[sheetNames,]
    sheet$Barcode <- targets[, gb$col_sentrix]
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


getRgset <- function(rgOut = NULL, targets=NULL, mergeProbes = F, csvPath = "samplesheet.csv", idatPath = NULL){
    require("minfi")
    if(is.null(targets)){targets <- as.data.frame(read.csv("samplesheet.csv"))}
    if(is.null(rgOut)){rgOut <- paste0(format(Sys.Date(),"%b%d"), "_RGset.Rdata")}
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
            sheet <- CleanUpSheetRows(as.data.frame(sheet), idatPath, targets)
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
        file.copy(csvPath, idatPath, overwrite = T)
        sheet <- minfi::read.metharray.sheet(idatPath, pattern =  basename(csvPath))
    }else{
        file.copy(csvPath, idatPath, overwrite = T)
        sheet <- minfi::read.metharray.sheet(idatPath, pattern = csvPath)
    }
    if(any(duplicated(sheet$Sample_Name))){
        message("DUPLICATE SAMPLES DROPPED")
        message(paste0(capture.output(sheet[duplicated(sheet$Sample_Name),]), collapse="\n"))
        toKeep <- !duplicated(sheet$Sample_Name)
        sheet <- sheet[toKeep,]
        rownames(sheet) <- 1:nrow(sheet)
    }
    sheet <- gb$CleanUpSheetRows(sheet, gb$idatPath, targets)
    return(sheet)
}


GetRgsetDat <- function(csvPath = "samplesheet.csv", gb) {
    require("minfi")
    gc(verbose = F)
    if (file.exists(gb$rgOut)) {
        RGSet <- gb$LoadRdatObj(gb$rgOut)
        return(RGSet)
    }
    if (is.null(gb$idatPath)) {
        gb$idatPath <- getwd()
    }
    sheet <- GrabMinfiSheet(gb$idatPath, csvPath)
    
    if (gb$mergeProbes == T & !is.null(gb$col_arrayType)) {
        message("Merging 450K and EPIC probes...")
        RGSet <- gb$combine.EPIC.450K(targets = sheet, gb)
        gb$SaveObj(RGSet, file.name = gb$rgOut)
        return(RGSet)
    } else{
        RGSet <- minfi::read.metharray.exp(base = gb$idatPath, targets = sheet, verbose = T, force = T)
        gb$SaveObj(RGSet, file.name = gb$rgOut)
        return(RGSet)
    }
}


MatchBetas2Targets <- function(targets, betas){
  targets <- targets[targets[,gb$col_samNames] %in% colnames(betas),]
  row.names(targets) <- 1:nrow(targets)
  return(targets)
}


cleanRawProbes <- function(RGSet, targets, gb) {
    gc(verbose = F)
    if (!file.exists(gb$rawBetaFi)) {
        betas <- gb$cleanUpProbes(RGSet, targets, gb, getNoob = gb$getNoob)
        SaveObj(betas, file.name = gb$rawBetaFi)
    } else{
        betas <- LoadRdatObj(gb$rawBetaFi)
    }
    if (!is.null(gb$batch_col)) {
        targets <<- MatchBetas2Targets(targets, betas)
        betas <- gb$RemoveBatchEffect(betas, targets, gb)
    }
    targets <<- MatchBetas2Targets(targets, betas)
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
    myColors <- c(targets[, varColorCol])
    names(myColors) <- targets[,varName]
    mds <- limma::plotMDS(mSetSq.beta, top = topN, labels = colnames(mSetSq.beta), 
                          gene.selection = "common", plot = F)
    toplot <- data.frame(PCA1 = mds$x, PCA2 = mds$y, Group = factor(names(myColors)))
    outMDSplot <- 
      ggplot2::ggplot(toplot, aes(PCA1, PCA2, colour = Group)) + 
      ggplot2::geom_point(aes(color = Group), size = 5, alpha = 0.65) +
      ggplot2::scale_color_manual(values = myColors) + ggplot2::theme_bw() + 
      labs(main = paste("Top", topN, "Common", "mSet Sq.beta", "MDS plot")) +
      theme(plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
            legend.position = "right", legend.box = "vertical")
    return(outMDSplot)
}


LoopPlotMds <- function(mSetSq.beta, targets, selectedVars, topN) {
  for (varNum in 1:length(unique(selectedVars))) {
    cat("\n\n")
    cat(paste("## MDS plot of Variable:", selectedVars[varNum]))
    cat("\n\n")
    mdsFi <- plot.mds(mSetSq.beta, targets, selectedVars[varNum], topN)
    mdsFi[["theme"]][["legend.position"]] <- "bottom"
    mdsFi[["theme"]][["legend.text"]][["size"]] <- 12
    mdsFi[["theme"]][["axis.text"]][["size"]] <- 12
    mdsFi[["theme"]][["axis.title.x"]][["size"]] <- 14
    mdsFi[["theme"]][["axis.title.y"]][["size"]] <- 14
    print(mdsFi)
    cat("\n\n")
  }
    cat("\n\n")
    cat(":::")
    cat("\n\n")
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
    targets <<- gb$getTargCsv(csvFi)
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


GetArrayTypes <- function(targets,
                          arrayColumn = NULL,
                          outputFi = "samplesheet_annotated_epic.csv",
                          sentrixID_col = "SentrixID_Pos",
                          idatPath = NULL) {
    if (is.null(arrayColumn)) {
        targets$ArrayType <- NULL
        arrayColumn <- "ArrayType"
    }
    if (is.null(idatPath)) {
        idatPath <- getwd()
    }
    targets$Basename <- file.path(idatPath, targets[, sentrixID_col])
    sentrix.ids <- as.character(targets[, sentrixID_col])
    outCsvDir <- file.path(getwd(), "csv", outputFi)
    if (!file.exists(outCsvDir)) {
        message("Adding ArrayTypes (EPIC vs 450k) to Samplesheet saving in:\n", outCsvDir)
        for (i in 1:length(sentrix.ids)) {
            sample_idat <- sentrix.ids[i]
            idat_path <- file.path(idatPath, sample_idat)
            rgset <- minfi::read.metharray(idat_path, verbose = T, force = T)
            is450K <- rgset@annotation[["array"]] == "IlluminaHumanMethylation450k"
            if (is450K) {
                targets[i, arrayColumn] <- "450k"
            } else{
                targets[i, arrayColumn] <- "EPIC"
            }
        }
        write.csv(targets, file = outCsvDir, row.names = F, quote = F)
    }
    targets <- as.data.frame(read.csv(outCsvDir))
    return(targets)
}


Get_Unfiltered_Betas <- function(idatPath, samsheet, dataFi="data/unfilteredBetas.Rdata"){
    if (file.exists(dataFi)) {
        rawBetas <- gb$LoadRdatObj(dataFi)
    }else{
         sheet <- GrabMinfiSheet(gb$idatPath, gb$samsheet)
    rawRGSet <- minfi::read.metharray.exp(base = idatPath, targets = sheet, verbose = T, force = T)
    rawBetas <- minfi::getBeta(rawRGSet)
    colnames(rawBetas) <- sheet[,1]
    gb$SaveObj(rawBetas, file.name = dataFi)
    }
    return(rawBetas)
}

