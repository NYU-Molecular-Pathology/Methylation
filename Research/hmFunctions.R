#!/usr/bin/env Rscript
## ---------------------------
## Script name: hmFunctions.R
## Purpose: source of global scripts for methylation heatmap analysis
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

gb <- globalenv(); assign("gb", gb)

# FUN: Loads required libraries into the environment
LoadHeatMapLibs <- function(){
  suppressPackageStartupMessages(library("ComplexHeatmap"))
  suppressMessages(require("ComplexHeatmap"))
  ComplexHeatmap::ht_opt(message = FALSE)
  #suppressPackageStartupMessages(library('magick'))
  suppressMessages(require('Cairo'))
  library("grid")
  suppressMessages(require("grid"))
  suppressMessages(options(bitmapType='cairo'))
  #ht_opt$message = FALSE
}

create_hm_colors <- function(targets1, varColumns) {
    if (length(varColumns) <= 1) {
      if ("Type" %in% varColumns == F) {
        targets$Type <- targets[, varColumns[1]]
        varColumns <- c(varColumns, "Type")
      }
      
      if (length(unique(varColumns)) == 1) {
        targets$NewCol <- targets[, varColumns[1]]
        varColumns <- c(varColumns, "NewCol")
      }
    }
    dat <- targets1[, varColumns]
    anno_df <- data.frame(dat)
    colorValues <- list(varColumns)
    names(colorValues[[1]]) <- varColumns
    
    for (curr_col in colnames(dat)) {
        color_col <- paste0(curr_col, "_color")
        curr_vals <- dat[, curr_col]
        curr_colors <- targets1[, color_col]
        names(curr_colors) <- curr_vals
        colorValues[[curr_col]] <- curr_colors
    }
  
    ha <- suppressWarnings(gb$getHeatAnno(colorValues, anno_df))
    return(ha)
}


# Returns heatmap annotations to colors on variables
getHeatAnno <- function(colorValues, anno_df) {
  ha <- ComplexHeatmap::HeatmapAnnotation(
    df = anno_df,
    col = colorValues,
    na_col = "white",
    annotation_legend_param = list(
      title_gp = gpar(
        fontsize = 12,
        lineheight = 2,
        fontface = "bold"
      ),
      labels_gp = gpar(fontsize = 12, lineheight = 2),
      grid_height = unit(0.8, "cm"),
      grid_width = unit(0.8, "cm"),
      gt_render = list(padding = unit(c(2, 2, 2, 2), "cm"))
    ),
    annotation_name_gp = gpar(
      fontsize = 12,
      fontface = "bold",
      lineheight = 2,
      gt_render = list(padding = unit(c(2, 2, 2, 2), "cm"))
    ),
    annotation_name_side = "left"
  )
  par(c(1.5, 1.5, 1.5, 2.5) + 0.5)
  return(ha)
}

#' assignColors returns a ComplexHeatmap annotation object with a list of variable columns you have
#'
#' @param targets your dataframe of sample IDs, and variables i.e. Age, Sex, Tumor, etc.
#' @param varColumns string name of the columns you want to create colors for in the complexheatmap within targets
#' @param col_vect NULL is default, these are all the colors you want to use in you heatmap annotations
assignColors <- function(targets, varColumns = c("Type", "Grade"), col_vect = NULL) {
    if (is.null(col_vect)) {col_vect <- pals::glasbey()}
    col_vect <- col_vect[sample(1:length(col_vect))]
    #dimnames(targets)[[2]]
    dat <- targets[,varColumns]
    anno_df <- data.frame(dat)
    vars2Color <- as.list(lapply(dat, unique))
    colorValues <-lapply(vars2Color, function(x) {x = (col_vect)[1:(length(x))]})
    for (x in 1:length(vars2Color)) {
        for (varNum in 1:length(vars2Color[x])) {
          names(colorValues[x][[1]]) = c(vars2Color[x][[1]])
        }
    }
    if(all(names(colorValues)==names(anno_df))){
        return(getHeatAnno(colorValues,anno_df))
    }
}

## Define the Heatmap ----------------------
drawHeatMap <- function(yourHeatMap) {
    return(
        ComplexHeatmap::draw(
            yourHeatMap,
            merge_legend = F,
            ht_gap = unit(3, "cm"),
            heatmap_legend_side = "left",
            annotation_legend_side = "right",
            column_title_gp = grid::gpar(fontsize = 10),
            padding = unit(c(1, 1, 1, 1), "cm")
        )
    )
}


GetHmDimensions <- function(hmTopNumbers){
    hm_width <- hmTopNumbers@matrix_param[["width"]]
    hm_ht <- hmTopNumbers@matrix_param[["height"]]
    hm_width <- round(grid::convertUnit(hm_width, "inches"), 1)
    hm_ht <- round(grid::convertUnit(hm_ht, "inches"), 1)
    hm_width <- as.numeric(hm_width) + 5
    hm_ht <- as.numeric(hm_ht) + 3
    hm_asp <- hm_ht/hm_width
    
    if (hm_ht > 24 | hm_width > 24) {
      hm_ht <- hm_ht / 10
      hm_width <- hm_width / 10
    }
    message(paste("Height is:", hm_ht, "Width is:", hm_width))
    
    knitr::opts_chunk$set(
      fig.dim = NULL,
      fig.width = hm_width,
      fig.height = hm_ht,
      dpi = 350,
      out.width = '100%',
      fig.asp = hm_asp
    )
    
    return(
      knitr::opts_current$set(
        fig.dim = NULL,
        fig.width = hm_width,
        fig.height = hm_ht,
        dpi = 350,
        out.width = '100%',
        fig.asp = hm_asp
      )
    )
}


calc_ht_size = function(ht, unit = "inch") {
    pdf(NULL)
    ht = draw(ht)
    w = ComplexHeatmap:::width(ht)
    w = convertX(w, unit, valueOnly = TRUE)
    h = ComplexHeatmap:::height(ht)
    h = convertY(h, unit, valueOnly = TRUE)
    dev.off()
    c(w, h)
}

getHeatMap <- function(betaRanges,
                       titleValue,
                       ha,
                       geneNams = F,
                       colSplt = NULL,
                       rwsplt = NULL,
                       hideTopAnno = F,
                       hideSamName = F,
                       rowClus = T,
                       colClus = T) {
    if (is.null(colSplt)) {
        hmLevels <- ha@anno_list[[1]]@color_mapping@levels
        colSplt <- length(unique(hmLevels))
    }
    col_fun2 <- circlize::colorRamp2(
        c(0, 0.25, 0.5, 0.75, 1),
        c("darkblue", "deepskyblue", "white", "tomato", "red")
    )
    titleOfPlot <- paste("Heatmap of", titleValue, sep = " ")
    
    # Determine whether to label rows based on the number of rows
    toLabRows <- nrow(betaRanges) <= 60
    geneNams <- geneNams && toLabRows
    
    # Adjust row height based on the number of rows
    rowTall <- if (toLabRows) 5 else 2
    shrinkRows <- nrow(betaRanges) > 800
    if (shrinkRows == T) {
        rowTall <- ifelse(nrow(betaRanges) > 2000, 0.06, 0.6) 
    }
    
    try(knitr::opts_chunk$set(out.width = '100%'), silent = TRUE)
    
    # Create the heat map
    hmTopNumbers <- ComplexHeatmap::Heatmap(
        betaRanges,
        col = col_fun2,
        cluster_columns = colClus,
        cluster_rows = rowClus,
        show_column_names = hideSamName,
        column_names_gp = gpar(fontsize = 12),
        show_row_names = geneNams,
        row_names_side = "left",
        row_title_side = "left",
        row_names_gp = gpar(fontsize = 10),
        row_title_gp = gpar(fontsize = 12, fontface = "bold"),
        show_row_dend = FALSE,
        show_column_dend = TRUE,
        width = ncol(betaRanges) * unit(5, "mm"),
        height = nrow(betaRanges) * unit(rowTall, "mm"),
        show_heatmap_legend = TRUE,
        top_annotation = ha,
        column_title = titleOfPlot,
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        heatmap_legend_param = list(
            title = "Beta Value",
            labels_gp = gpar(fontsize = 14),
            title_gp = gpar(fontsize = 14, fontface = "bold"),
            legend_direction = "vertical",
            legend_padding = unit(5, "mm"),
            annotation_legend_side = "right",
            legend_height = unit(2.5, "in")
        ),
        column_split = colSplt,
        row_split = rwsplt,
        use_raster = shrinkRows,
        raster_device = "png",
        raster_quality = 1,
        raster_resize_mat = T
    )
    
    if (hideTopAnno) {
        hmTopNumbers@top_annotation <- NULL
    }
    if (hideSamName) {
        hmTopNumbers@column_names_param[["labels"]] <- NULL
        hmTopNumbers@column_names_param[["show"]] <- FALSE
    }
    
    GetHmDimensions(hmTopNumbers)
    return(drawHeatMap(hmTopNumbers))
}

assignColors2 <- function(targets, varColumns = c("Type","Origin"), col_vect = NULL) {
    if (is.null(col_vect)) {col_vect <- pals::glasbey()}
    #col_vect <- col_vect[sample(1:length(col_vect))]
    col_vect[6] = "#eb7d34" #changing dark forest to orange color
    col_vect[4] = "#ADD8E6"
    dimnames(targets)[[2]]
    dat <- targets[,varColumns] # varColumns
    #colours <- col_vect[1:(length(unique(dat)))]
    anno_df <- data.frame(dat)
    vars2Color <- as.list(lapply(dat, unique))
    colorValues <-lapply(vars2Color, function(x) {x = (col_vect)[1:(length(x))]})
    for (x in 1:length(vars2Color)) {
        for (varNum in 1:length(vars2Color[x])) {names(colorValues[x][[1]]) = c(vars2Color[x][[1]])}
    }
    if(all(names(colorValues)==names(anno_df))){return(gb$getHeatAnno(colorValues,anno_df))}
}

modifyHaColors <- function(varColumns, targets, ha) {
  targets <- gb$colorTargets(targets, varColumns)
  for (vCol in varColumns) {
    colorColName <- paste0(vCol, "_color")
    colours <- unique(targets[, colorColName])
    varNames <- sort(unique(targets[, vCol]))
    colours <- colours[1:length(varNames)]
    names(colours) <- varNames
    colMap <- ha@anno_list[[vCol]]@color_mapping
    colMap@colors <- colMap@full_col <- colours
    colMap@levels <- names(colours)
    ha@anno_list[[vCol]]@color_mapping <- colMap
  }
  return(ha)
}

assignColors3 <- function(targets, varColumns = c("Type", "Origin"), col_vect = NULL, manualNames = NULL) {
    message("\n~~~~Targets Dimnames:\n", paste(dimnames(targets)[[2]], collapse = "\n"))
    
    col_vect <- gb$GetDefaultColors(col_vect)
    targets <- gb$colorTargets(targets, varColumns)
    targets <- gb$FixNullNaVars(targets, varColumns)
    
    if (length(varColumns) <= 1) {
      if ("Type" %in% varColumns == F) {
        targets$Type <- targets[, varColumns[1]]
        varColumns <- c(varColumns, "Type")
      }
      
      if (length(unique(varColumns)) == 1) {
        targets$NewCol <- targets[, varColumns[1]]
        varColumns <- c(varColumns, "NewCol")
      }
    }
    
    dat <- targets[, varColumns] # varColumns
    anno_df <- data.frame(dat)
    vars2Color <- as.list(lapply(dat, function(x) {return(sort(unique(x)))}))
    
    colorValues <- lapply(vars2Color, function(x) {
      x = (col_vect)[1:(length(x))]
    })
    
    colorValues <- lapply(X = names(colorValues), function(X) {
      currCol <- colorValues[X][[1]]
      curNam <- vars2Color[[X]]
      names(currCol) <- curNam
      colorValues[X][[1]] <- currCol
      names(colorValues[X]) <- X
      return(colorValues[X][[1]])
    })
    
    names(colorValues) <- varColumns
    stopifnot(all(names(colorValues) == names(anno_df)))
  
    for(x in names(colorValues)){
      targCol <- paste0(x, "_color")
      currColumn <- colorValues[[x]]
      for(y in names(currColumn)){
        matchCol <- targets[,targCol]
        matchSam <- which(targets[,x] == y)[1]
        newColSub <- matchCol[matchSam]
        colorValues[[x]][[y]] <- newColSub
      }
    }
  
    ha <- suppressWarnings(gb$getHeatAnno(colorValues, anno_df))
    return(ha)
  }


addGeneName <- function(RGSet, oldBeta, geneNameHm=F) {
    if(geneNameHm==F){
        return(oldBeta)
    }
    annot <- minfi::getAnnotation(RGSet)
    geneNameLi <- annot[rownames(oldBeta), "UCSC_RefGene_Name"]
    for (idx in 1:length(geneNameLi)) {
        gn <- geneNameLi[idx]
        if(gn == ""){
            gn <- annot@listData[["GencodeCompV12_NAME"]][idx]
        }
        if(!is.na(gn)){
            if (gn != "") {
                new_gn <- paste(stringr::str_split(gn, ";", simplify = T)[1, ])
                new_gn <-
                    paste(paste(unique(new_gn)), sep = "_", collapse = "_")
                new_gn <- paste(new_gn, rownames(oldBeta)[idx], sep = "_")
                rownames(oldBeta)[idx] <- new_gn
            }
        }
    }
    return(oldBeta)
}


# FUNCTION: Removes duplicated gene annotations from the heatmap rownames if desired after addGeneName is called
removeDupeAnnot <- function(geneBetas,out.fi="top_variance_genes_list.csv") {
  geneRows <- rownames(geneBetas)
  for (rw in 1:length(geneRows)) {
    spltRw <- stringr::str_split_fixed(geneRows[rw], pattern = "_", 2)
    rowGene <- spltRw[1, 1]
    if (rowGene == "") {rowGene = "NO.GENE.NAME"}
    rowProb <- spltRw[1, 2]
    rowGene <- stringr::str_split(rowGene, ";")
    rowGene <- paste(unique(rowGene[[1]]), collapse = ".")
    newRow <- paste(rowGene, rowProb, sep = "_")
    geneRows[rw] <- newRow
  }
  toWrite <- as.data.frame(stringr::str_split_fixed(geneRows,"_",2))
  colnames(toWrite)<- c("Gene_Names", "Probe_Name")
  write.csv(toWrite, out.fi, quote=F,row.names=F)
  geneRows <- lapply(geneRows, function(x){stringr::str_replace(x,"_cg",".")})
  rownames(geneBetas)<- geneRows
  return(geneBetas)
}


#' grabProbes function that lists probes annotated by matching gene name and gene region
#' @param your_genes  a character list of genes i.e. from a file c("ALK","ETMR","HDAC")
#' @param RGSet  the name or path to your minfi RGset data
#' @param region string vector of gene selections like c("TSS200", "Body") can one gene region "TSS200"
# FUNCTION: lists probes annotated by matching gene name and gene region
grabProbes <- function(your_genes, RGSet, region){
    stopifnot(is.character(your_genes))
    listRows=10000
    bp <- ifelse((length(region)>1), paste(region, collapse="|"), region)
    annot <- minfi::getAnnotation(RGSet)
    g=NULL
    geneAnno <- paste(annot$UCSC_RefGene_Name, annot$UCSC_RefGene_Group, sep = "_")
    pick <- function(x, fltr){return(do.call(grepl,c(list(pattern=x), fltr)))}
    z <- foreach::foreach(g = your_genes, .combine = "cbind") %dopar% {
        fltr <- list(x = geneAnno, ignore.case = T)
        y <- cbind(rownames(annot)[pick(g, fltr) & pick(bp, fltr)])
        naFill <- rbind(y, cbind(rep(NA, (listRows - nrow(y)))))
        colnames(naFill) <- g
        return(naFill)
    }
    naDrops <- lapply(1:ncol(z), function(i){which(!is.na(z[,i]))})
    cutt <- length(naDrops[[which.max(lapply(naDrops, function(x) sum(lengths(x))))]])
    return(z[1:cutt,])
}
                                           

SaveHmPng <- function(imgFile, hm) {
    wd_original <- as.numeric(hm@ht_list_param[["width"]])
    ht_original <- as.numeric(hm@ht_list_param[["height"]])
    resolution <- 150
    wd <- round(wd_original * 8)
    ht <- round(ht_original * 8)
    if (ht > 10000 | wd > 10000) {
      wd <- round(wd / 2)
      ht <- round(ht / 2)
      resolution <- 100
    }
    grDevices::png(
      file = imgFile,
      width = wd,
      height = ht,
      units = "px",
      res = resolution
    )
    gb$drawHeatMap(hm)
    invisible(dev.off())
}

                                           
# FUN: Alternate heatmap function                                           
getHeatMap2 <- function(betaRanges, titleValue, ha, geneNamesHeatMap = F, colSplt = NULL, rwsplt = NULL, rwOrder = NULL, clusRows = F, clusCols = F){
  col_fun3 <- circlize::colorRamp2(c(0, 0.20, 0.25, 0.5, 0.75, 1), c("black","darkblue","deepskyblue", "white", "tomato","red"))
  titleOfPlot <- paste("Heatmap of", titleValue, sep = " ")
   hmTopNumbers <- ComplexHeatmap::Heatmap(betaRanges, 
        col = gb$col_fun3,  ## Define the color scale
        cluster_rows = clusRows,
        cluster_row_slices = F,
        cluster_column_slices = F,
        cluster_columns = clusCols,  ## Cluster the columns
        show_column_names = F,  ## Show the Column Names (which is sample #)
        column_names_gp = gpar(fontsize = 12),  ## Column Name Size
        show_row_names = geneNamesHeatMap,  ## Show Row names (which is probes)
        row_order=rwOrder,    
        row_names_side = "left",
        row_title_side = "left",
        row_names_gp = gpar(fontsize = 10),
        row_title_gp = gpar(fontsize = 12, fontface = "bold"),
        show_row_dend = F,
        show_column_dend = F,
        use_raster=T,
        show_heatmap_legend = T,
        top_annotation = ha,
        column_title = titleOfPlot,
        column_title_gp = gpar(fontsize = 14,fontface = "bold"),
        raster_device = "CairoPNG",
        raster_quality = 3,
        heatmap_legend_param = list(
            title = "Beta Value",
            labels_gp = gpar( fontsize = 12),
            title_gp = gpar(fontsize = 12, fontface = "bold"),
            legend_direction = "vertical",
            heatmap_legend_side = "right", 
            annotation_legend_side = "right",
            legend_height =  unit(2.5, "in")
        ),
        column_split = colSplt,
        row_split= rwsplt,
        heatmap_width = unit(10, "in"),
        heatmap_height = unit(20, "in")
    )
    return(gb$drawHeatMap(hmTopNumbers))
}


# FUN: Returns blank data frame to fill                                            
GetEmptyDf <- function(csvColumns,betas){
    temp <- matrix(nrow= ncol(csvColumns), ncol=dim(betas)[2])
    avgBetas <- data.frame(temp)
    colnames(avgBetas) <- colnames(betas)
    rownames(avgBetas) <- colnames(csvColumns)
    return(avgBetas)
}


GetProbeAverage <- function(csvColumns, betas, pathwayName){
    avgBetas <- GetEmptyDf(csvColumns,betas)
    for (geneCol in 1:ncol(csvColumns)) {
        probeList <- csvColumns[,geneCol]
        for (sam in 1:dim(betas)[2]) {
            probeBetas <- betas[rownames(betas) %in% probeList, sam]
            sampleNam <- colnames(betas)[sam]
            geneNam <- colnames(csvColumns)[geneCol]
            probeAvg <-mean(probeBetas)
            avgBetas[geneNam,sampleNam] <- probeAvg
        }
    }
    outDir <- file.path(getwd(), "csv")
    pathwayName <- paste(stringr::str_split(pathwayName, " ", simplify = T), collapse="_")
    outFile <- paste(pathwayName, "avgBetas_per_gene.csv" , sep = "_")
    outPath <- file.path(outDir, outFile)
    avgBetas$XXX <- rownames(avgBetas)
    write.csv(avgBetas, file = outPath, row.names=F, na="", quote=F)
    return(avgBetas)
}


#' grabProbes function that lists probes annotated by matching gene name and gene region
#' @param your_genes  a character list of genes i.e. from a file c("ALK","ETMR","HDAC")
#' @param RGSet  the name or path to your minfi RGset data
#' @param region string vector of gene selections like c("TSS200", "Body") can one gene region "TSS200"
GetProbesGenes <- function(your_genes, RGSet, region="Island"){
  stopifnot(is.character(your_genes))
  doParallel::registerDoParallel(cores=6)
  listRows=15000
  bp <- ifelse((length(region)>1), paste(region, collapse="|"), region)
  annot <- minfi::getAnnotation(RGSet)
  g=NULL
  geneAnno <- paste(annot$UCSC_RefGene_Name, annot@listData[["Relation_to_Island"]], sep = "_") #annot$UCSC_RefGene_Group
  pick <- function(x, fltr){return(do.call(grepl,c(list(pattern=x), fltr)))}
  z <- foreach::foreach(g = your_genes, .combine = "cbind") %dopar% {
    fltr <- list(x = geneAnno, ignore.case = T)
    y <- cbind(rownames(annot)[pick(g, fltr) & pick(bp, fltr)])
    naFill <- rbind(y, cbind(rep(NA, (listRows - nrow(y)))))
    colnames(naFill) <- g
    return(naFill)
  }
  naDrops <- lapply(1:ncol(z), function(i){which(!is.na(z[,i]))})
  cutt <- length(naDrops[[which.max(lapply(naDrops, function(x) sum(lengths(x))))]])
  return(z[1:cutt,])
}


AddGeneProbeChrName <- function(RGSet, oldBeta) {
  annot <- minfi::getAnnotation(RGSet)
  geneNameLi <- annot[rownames(oldBeta), "UCSC_RefGene_Name"]
  promotorLi <- annot[rownames(oldBeta), "GencodeBasicV12_Group"]
  chromoLi <- annot[rownames(oldBeta), "chr"]
  
  for(x in 1:length(promotorLi)){
    theReg <- promotorLi[x]
    if(theReg!=""){
      regSplit <- paste(stringr::str_split(theReg, pattern = ";", simplify = T))
      newReg <- paste(unique(regSplit), sep = ";", collapse = ";")
      promotorLi[x] <- newReg
    }
  }
  for (idx in 1:length(geneNameLi)) {
    gn <- geneNameLi[idx]
    
    if(!is.na(gn)){
      if (gn != "") {
        new_gn <- paste(stringr::str_split(gn, ";", simplify = T)[1, ])
        new_gn <-
          paste(paste(unique(new_gn)), sep = ";", collapse = ";")
        new_gn <- paste(new_gn, rownames(oldBeta)[idx], promotorLi[idx], chromoLi[idx], sep = "_")
        rownames(oldBeta)[idx] <- new_gn
      }
    }
  }
  return(oldBeta)
}
                                           
GetCsvGeneColumns <- function(pathwayName, z){
    outFile1 <- paste(pathwayName,"genes_columns_withMissing.csv", sep = "_")
    outDir <- file.path(getwd(),"csv")
    if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}
    write.csv(z, file = file.path(outDir, outFile1), row.names=F, na="")
    z <- as.data.frame(z)
    toDrop <- c(is.na(z[1,]))
    z <- z[,!toDrop]
    outFile2 <- paste(pathwayName,"genes_columns.csv", sep = "_")
    write.csv(z, file = file.path(outDir, outFile2), row.names=F, na="")
    csvColumns <- as.data.frame(read.csv(file.path(outDir, outFile2)))
    return(csvColumns)
}

GetHeatMapGenes <-  function(betaRanges, titleValue, ha, geneNamesHeatMap = F, colSplt = NULL, rwsplt = NULL){
  col_fun2 <- circlize::colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("darkblue","deepskyblue", "white", "tomato","red"))  
  titleOfPlot <- paste("Heatmap of", titleValue, sep = " ")
    hmTopNumbers <- ComplexHeatmap::Heatmap(
        betaRanges,
        col = gb$col_fun2,  ## Define the color scale
        cluster_columns = T,  ## Cluster the columns
        #cluster_rows = rowcluster,
        #raster_resize_mat = TRUE,
        show_column_names = T,  ## Show the Column Names (which is sample #)
        column_names_gp = gpar(fontsize = 12),  ## Column Name Size
        show_row_names = geneNamesHeatMap,  ## Show Row names (which is probes)
        row_names_side = "left",
        row_title_side = "left",
        row_names_gp = gpar(fontsize = 10),
        row_title_gp = gpar(fontsize = 10, fontface = "bold"),
        show_row_dend = F,
        show_column_dend = T,
        use_raster=T,
        show_heatmap_legend = T,
        top_annotation = ha,
        column_title = titleOfPlot,
        column_title_gp = gpar(fontsize = 12,fontface = "bold"),
        raster_device = "CairoPNG",
        raster_quality = 3,
        heatmap_legend_param = list(
            title = "Beta Value",
            #legend_height = unit(5, "in"),
            labels_gp = gpar( fontsize = 12),
            title_gp = gpar(fontsize = 12, fontface = "bold"),
            legend_direction = "vertical",
            heatmap_legend_side = "right", annotation_legend_side = "right",
            legend_height =  unit(2.5, "in")
        ),
        column_split = colSplt,
        row_split= rwsplt,
        heatmap_width = unit(10, "in"),
        heatmap_height = unit(30, "in")
    )
    #size = gb$calc_ht_size(hmTopNumbers)
    #size
    return(gb$drawHeatMap(hmTopNumbers))
}


CheckGeneOutput <- function(pathwayName) {
    outDir <- file.path(getwd(), "csv")
    if(!dir.exists(outDir)){dir.create(outDir)}
    pathwayName <- paste(stringr::str_split(pathwayName, " ", simplify = T), collapse="_")
    outFile <- paste(pathwayName, "avgBetas_per_gene.csv" , sep = "_")
    outPath <- file.path(outDir, outFile)
    if(file.exists(outPath)){
        return(outPath)
    }else{
        return(FALSE)    
    }
}

CatImg <- function(subTitle, imgFile){
    cat("\n\n")
    imgCat <- paste0("![", subTitle, "](", imgFile, ")")
    cat(imgCat)
    cat("\n\n")
    cat("\n\n")
}

GetAvgBetaAnno <- function(targets, avgBetas, selectedVars){
    toKeep <- which(targets[, 1] %in% colnames(avgBetas))
    targets1 <- targets[toKeep, ]
    rownames(targets1) <- 1:nrow(targets1)
    rownames(targets1) <- targets1[, 1]
    targets1 <- targets1[colnames(avgBetas), ]
    rownames(targets1) <- 1:nrow(targets1)
    targets1 <- gb$colorTargets(targets1, varColumns = selectedVars)
    ha <- gb$AnnotateHmVars(targets1, varColumns = selectedVars)
    ha <- gb$FilterHmAnno(ha, selectedVars) # drop any unwanted columns
    ha <- gb$MatchHaLegend(ha, selectedVars, targets1)
    return(ha)
}

GetAvgGeneHeatMap <- function(avgBetas, titleValue, ha, geneNamesHeatMap=T, colSplt = NULL, rwsplt=NULL){
    col_fun2 <- circlize::colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("darkblue","deepskyblue", "white", "tomato","red"))
    titleOfPlot <- paste("Heatmap of", titleValue)
    toLabRows <- ifelse(nrow(avgBetas)<=50, T, F)
    if(geneNamesHeatMap == T){
        geneNamesHeatMap <- toLabRows
    }
    colTall <- ncol(avgBetas)*unit(5, "mm")
    if(ncol(avgBetas)*5 <= 20){
        colTall <- colTall + unit(45, "mm")
    }
    rowTall <- ifelse(toLabRows==T, 5, 2)
    hmTopNumbers <- ComplexHeatmap::Heatmap(
        avgBetas,
        col = gb$col_fun2,  ## Define the color scale
        cluster_columns = T,  ## Cluster the columns
        cluster_rows = cluster::diana,
        show_column_names = T,  ## Show the Column Names (which is sample #)
        column_names_gp = gpar(fontsize = 12),  ## Column Name Size
        show_row_names = geneNamesHeatMap,  ## Show Row names (which is probes)
        row_names_side = "left",
        row_title_side = "left",
        row_names_gp = gpar(fontsize = 10),
        row_title_gp = gpar(fontsize = 12, fontface = "bold"),
        show_row_dend = F,
        show_column_dend = T,
        width = colTall,
        height = nrow(avgBetas)*unit(rowTall, "mm"),
        use_raster=T,
        show_heatmap_legend = T,
        top_annotation = ha,
        column_title = titleOfPlot,
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        raster_device = "CairoPNG",
        raster_quality = 1,
        raster_resize_mat = TRUE,
        heatmap_legend_param = list(
            title = "Beta Value",
            labels_gp = gpar( fontsize = 14),
            title_gp = gpar(fontsize = 14, fontface = "bold"),
            legend_direction = "vertical",
            heatmap_legend_side = "right", annotation_legend_side = "right",
            legend_height =  unit(1, "in")
        ),
        column_split = colSplt,
        row_split = rwsplt
    )
    GetHmDimensions(hmTopNumbers)
    return(gb$drawHeatMap(hmTopNumbers))
}


LoopPathwayHeatMap <- function(pathWayGenes, RGSet, betas, targets) {
  cat("\n\n")
  cat("\n\n")
  cat("## HeatMaps of Genes In Pathways {.tabset}")
  cat("\n\n")
  hmOutPath <- getwd()
  for (pathRow in 1:nrow(pathWayGenes)) {
    currPathway <- pathWayGenes[pathRow, ]
    pathwayName <- paste0(gsub(" ", "_", currPathway$Description))
    message("Looping Pathway Creation for: ", pathwayName)
    avgExist <- gb$CheckGeneOutput(pathwayName)
    avgBetas <- NULL
    if (avgExist == F) {
      your_genes <- stringr::str_split(currPathway$geneID, pattern = "/")[[1]]
      z <- gb$GetProbesGenes(your_genes, RGSet, region = c("Island"))
      csvColumns <- gb$GetCsvGeneColumns(pathwayName, z)
      avgBetas <- gb$GetProbeAverage(csvColumns, betas, pathwayName)
    } else{
      avgBetas <- read.csv(avgExist)
    }
    rwNm = NULL
    rwNm <- which(colnames(avgBetas) == "XXX")
    row.names(avgBetas) <- avgBetas[, rwNm]
    avgBetas <- avgBetas[, -rwNm]
    titleValue <- paste("Average Probe Beta Values for", currPathway$Description, "Genes")
    avgBetas <- as.matrix(avgBetas)
    avgBetas <- na.omit(avgBetas)
    msgTitle <- paste("###", paste0(currPathway$Description), "\n\n")
    cat("\n\n")
    cat(msgTitle)
    cat("\n\n")
    outDir <- file.path(".", "figures", "heatmaps")
    imgFile <- file.path(outDir, paste0("hm_genes_", pathwayName, ".png"))
    if (!file.exists(imgFile)) {
      colnames(avgBetas) <- gsub(".", "-",  colnames(avgBetas), fixed = TRUE)
      ha <- GetAvgBetaAnno(targets, avgBetas, gb$selectedVars)
      hm <- gb$GetAvgGeneHeatMap(avgBetas, titleValue, ha, geneNamesHeatMap = T)
      hm
      cat("\n\n")
      gb$SaveHmPng(fi_prefix = "hm_genes_", fi_suffix = ".png", hm, topvar = pathwayName, outDir = NULL)
    }
    CatImg(paste0(currPathway$Description), imgFile)
    if (nrow(avgBetas) >= 50) {
      imgFile2 <- file.path(outDir, paste0("hm_genes_top_50_", pathwayName, ".png"))
      msgTitle2 <- paste("###", paste0(currPathway$Description), "Top 50", "\n\n")
      cat("\n\n")
      cat(msgTitle2)
      cat("\n\n")
      if (!file.exists(imgFile2)) {
        avgBetas2 <- gb$takeTopVariance(avgBetas, 1:50)
        colnames(avgBetas) <- gsub(".", "-",  colnames(avgBetas), fixed = TRUE)
        ha <- GetAvgBetaAnno(targets, avgBetas, gb$selectedVars)
        hm <- gb$GetAvgGeneHeatMap(avgBetas2, titleValue, ha, geneNamesHeatMap = T)
        hm
        cat("\n\n")
        gb$SaveHmPng(fi_prefix = "hm_genes_top_50_", fi_suffix = ".png", hm, topvar = pathwayName, outDir = NULL)
      }
      CatImg(paste0(currPathway$Description), imgFile2)
    }
  }
  cat("\n\n")
}

                                           
  LoopSaveHm <- function(hm, topNum, fi_prefix = "unsuper_hm_top", fi_suffix = "notAnnot"){
    supervise <- ifelse(stringr::str_detect(fi_prefix, pattern = "unsuper"), "unsupervised", "supervised")
    hmOutPath <- file.path(".", "figures", "heatmaps", supervise)
    if (!dir.exists(hmOutPath)) {dir.create(hmOutPath, recursive = T)}
    imgFile <- file.path(hmOutPath, paste0(fi_prefix, "_", topNum, "_", fi_suffix, ".png"))
    if (!file.exists(imgFile)) {SaveHmPng(imgFile, hm)}
    cat("\n\n")
    txtLink <- paste0("[", "Click_to_View_", topNum,"_", supervise, "](", imgFile, ")")
    cat(txtLink)
    cat("\n\n")
}
                                           

LoopPrintHeatMap <- function(gb, unBetas, ha, fi_prefix = "unsuper_hm_top") {
  geneNams <- gb$addGenesHm
  colSplt <- gb$colSplitHm
  hideTopAnno <- gb$hideTopAnno
  hideSamName <- gb$hideSamName
  fi_suffix <- "Annotated"
  if (hideTopAnno == T) {
    fi_suffix <- "notAnnot"
  } 
  
  ComplexHeatmap::ht_opt("message" = FALSE)
  for (topNum in gb$varProbes) {
    invisible(gc(verbose = F))
    cat(paste("###", "Top", topNum, "\n\n"))
    hmTitle <- paste("Top", topNum, "Variance Probes Beta Values")
    bv <- unBetas[1:topNum,]
    currHm <- invisible(
        gb$getHeatMap(
            bv,
            hmTitle,
            ha,
            geneNams,
            colSplt = colSplt,
            hideTopAnno = hideTopAnno,
            hideSamName = hideSamName
        )
    )
    LoopSaveHm(hm = currHm, topNum, fi_prefix, fi_suffix)
    cat('\n\n')
  }
}


AnnotateHmVars <- function(targets1, varColumns){
  ha <- gb$assignColors3(targets1, varColumns)
  ha <- gb$modifyHaColors(varColumns, targets1, ha)
  return(ha)
}


# Dropping any unwanted column annotation
FilterHmAnno <- function(ha, varToPlot) {
    toKeep <- names(ha@anno_list) %in% varToPlot
    others <- unlist(lapply(names(ha@anno_list), function(anno_name) {
        length(ha@anno_list[[anno_name]]@color_mapping@levels) > 1
    }))
    toKeep <- toKeep & others
    ha@anno_list <- ha@anno_list[toKeep]
    ha@height <- grid::unit(1, "mm")
    return(ha)
}


PrintSingleHeatMap <- function(bv, ha, geneNams, colSplt = 3, topvar = 1, hmTitle = "Probe-Specific Beta Values") {
    cat(paste("###", "Probe-Specific Beta Value HeatMap", topvar, "\n\n"))
    hm <- getHeatMap(bv, hmTitle, ha, geneNams, colSplt = colSplt)
    hm
    cat('\n\n')
    hmOutPath <- file.path(".", "figures", "heatmaps")
    if (!dir.exists(hmOutPath)) {dir.create(hmOutPath, recursive = T)}
    SaveHmPng("hm_custom", "_probes.png", hm, topvar, outDir = hmOutPath)
}
                                           
                                           
MatchHaLegend <- function(ha, selectedVars, targets1){
    for (x in gb$selectedVars) {
        colorCol <- paste0(x, "_color")
        varColor <- unique(targets1[, colorCol])
        varTypes <- unique(targets1[, x])
        varList <- NULL
        for (y in 1:length(varColor)){
            colorIdx <- which(targets1[,colorCol] == varColor[y])[1]
            newColor <- targets1[colorIdx,x]
            if(!is.null(varList)){
                varList <- c(varList, newColor)
            }else{
                varList <- newColor
            }
        }
        names(varColor) <- varList
        ha@anno_list[[x]]@color_mapping@colors <- varColor
        ha@anno_list[[x]]@color_mapping@levels <- names(varColor)
        ha@anno_list[[x]]@color_mapping@full_col <- varColor
    }
    return(ha)
}
                                           
GetHeatMapData <-  function(targets, betas, RGSet, gb, getAll = F, varToFilter = NULL){
    varProbes <- gb$hmVarProbes
    varProbes <<- gb$hmVarProbes  # Change top variance for heatmaps
    gb$MsgNullVar(gb, "col_sentrix", "The Sentrix IDs for the samples", '"Sentrix_Pos"')
    gb$MsgNullVar(gb ,"col_samNames", "The unique sample names or IDs", '"RD_number"')
    gb$MsgNullVar(gb, "selectedVars", "The Phenotype or Metadata for the samples", 'c("Cell_type", "Tissue_Type")')

    RGSet@colData@listData$Sample_ID <- RGSet@colData@listData[[gb$col_samNames]]
    targets1 <- gb$SubsetTargets(targets, varToFilter)
    betas1 <- betas[, targets1[, gb$col_samNames]]
    unBetas <- gb$tierBetas(betas1, gb$col_sentrix, RGSet, getAll = T) # unsupervised betas
    #unBetas <- gb$addGeneName(RGSet, unBetas, gb$addGenesHm)
    ha <- create_hm_colors(targets1, varColumns = gb$selectedVars)
    ha <- gb$FilterHmAnno(ha, varToPlot = gb$selectedVars) # drop any unwanted columns
    ha <- gb$MatchHaLegend(ha, gb$selectedVars, targets1)
    hmPlotData <- list("unBetas" = unBetas, "ha" = ha)
    return(hmPlotData)
}
                                           
SubsetTargets <- function(targets, varToFilter = NULL){
  if(is.null(varToFilter)){return(targets)}
  targets1 <- targets[varToFilter,] # exclude of targets are subsetted
  rownames(targets1) <- 1:nrow(targets1) # ensures rows are chronological
  return(targets1)
}


GetSuperHmData <- function(gb, targets, RGSet) {
  library("grDevices")
  if (gb$supervisedRun) {
    for (i in 1:length(gb$selectedVars)) {
          superVar <- gb$selectedVars[i]
          varName <- paste0("hmPlotData", i)
          superFile <- paste0(gb$supbetaOut, "_", superVar, ".Rdata")
          superbetas1 <- gb$LoadRdatObj(superFile, msgProg = F)
          superHmFi <- paste0(gb$supbetaOut, "_hm_", superVar, ".Rdata")
          if (file.exists(superHmFi)) {
              hmPlotData <- gb$LoadRdatObj(superHmFi, msgProg = F)
          } else {
              hmPlotData <- gb$GetHeatMapData(targets, superbetas1, RGSet, gb)
              gb$SaveObj(hmPlotData, file.name = superHmFi)
          }
          assign(x = varName, value = hmPlotData, envir = gb)
      }
  }
}

                                         
LoopSupervisedHm <- function(gb) {
    if(gb$supervisedRun){
        for(i in 1:length(gb$selectedVars)){
          invisible(gc(verbose = F))
          cat("\n\n")
          hmPlotData <- eval(parse(text = paste0("gb$hmPlotData", i)))
          gb$LoopPrintHeatMap(gb, hmPlotData$unBetas, hmPlotData$ha,
                              fi_prefix = paste0("super_", gb$selectedVars[i], "_hm_top"))
          cat("\n\n")
          invisible(gc(verbose = F))
        }
    } else {
        cat("\n\nNo supervised HeatMap analysis output\n\n")
    }
}

