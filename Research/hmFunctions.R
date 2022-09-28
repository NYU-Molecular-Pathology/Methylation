#' getHeatAnno returns a heatmap annotation object where you input the list of variables you would like
#'
#' @param colorValues the variable name paired with color for each sample by variable
#' @param anno_df the dataframe that ou are annotating, these are the filtered columns of the varColumns
gb <- globalenv(); assign("gb", gb)

suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressMessages(require("ComplexHeatmap"))
suppressPackageStartupMessages(library('magick'))
suppressMessages(require('Cairo'))
library("grid")
suppressMessages(require("grid"))
ht_opt$message = FALSE
suppressMessages(options(bitmapType='cairo'))

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
col_fun2 <- circlize::colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("darkblue","deepskyblue", "white", "tomato","red"))

drawHeatMap <- function(yourHeatMap) {
    return(
        ComplexHeatmap::draw(
            yourHeatMap,
            merge_legend = T,
            ht_gap = unit(3, "cm"),
            heatmap_legend_side = "left",
            annotation_legend_side = "left",
            padding = unit(c(20, 20, 20, 20), "mm")
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

getHeatMap <- function(betaRanges, titleValue, ha, geneNamesHeatMap=F, colSplt = NULL, rwsplt=NULL){
    titleOfPlot <- paste("Heatmap of",titleValue,sep = " ")
    hmTopNumbers <- ComplexHeatmap::Heatmap(
        betaRanges,
        col = gb$col_fun2,  ## Define the color scale
        cluster_columns = T,  ## Cluster the columns
        #cluster_rows = rowcluster,
        #raster_resize_mat = TRUE,
        show_column_names = F,  ## Show the Column Names (which is sample #)
        column_names_gp = gpar(fontsize = 12),  ## Column Name Size
        show_row_names = geneNamesHeatMap,  ## Show Row names (which is probes)
        row_names_side = "left",
        row_title_side = "left",
        row_names_gp = gpar(fontsize = 10),
        row_title_gp = gpar(fontsize = 12, fontface = "bold"),
        show_row_dend = F,
        show_column_dend = T,
        use_raster=T,
        show_heatmap_legend = T,
        top_annotation = ha,
        column_title = titleOfPlot,
        column_title_gp = gpar(fontsize = 14,fontface = "bold"),
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
        heatmap_height = unit(20, "in")
    )
    #size = gb$calc_ht_size(hmTopNumbers)
    #size
    return(gb$drawHeatMap(hmTopNumbers))
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

modifyHaColors <- function(varColumns,targets, ha){
  for (vCol in varColumns) {
      colours <- unique(targets$color)
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

assignColors3 <- function(targets, varColumns = c("Type","Origin"), col_vect = NULL, manualNames=NULL) {
    dimnames(targets)[[2]]
    if ("Type" %in% varColumns == F) {
        targets$Type <- targets[, varColumns[1]]
        varColumns <- c(varColumns, "Type")
    }
    if (length(unique(varColumns)) == 1) {
        targets$NewCol <- targets[, varColumns[1]]
        varColumns <- c(varColumns, "NewCol")
    }  
    dat <- targets[,varColumns] # varColumns
    anno_df <- data.frame(dat)
    col_vect <-unique(targets$color)
    vars2Color <- as.list(lapply(dat, function(x){return(sort(unique(x)))}))
    colorValues <-lapply(vars2Color, function(x) {x = (col_vect)[1:(length(x))]})
    colorValues <-lapply(X=names(colorValues), function(X) {
        currCol <- colorValues[X][[1]]
        curNam <- vars2Color[[X]]
        names(currCol)<-curNam
        colorValues[X][[1]] <- currCol
        names(colorValues[X])<-X
        return(colorValues[X][[1]])}
    )
    names(colorValues)<- varColumns
    stopifnot(all(names(colorValues)==names(anno_df)))
    ha <- suppressWarnings(gb$getHeatAnno(colorValues,anno_df))
    return(ha)
}

addGeneName <- function(annot, oldBeta) {
  geneName <- annot[rownames(oldBeta), "UCSC_RefGene_Name"]
  rownames(oldBeta) <- paste(geneName, rownames(oldBeta), sep = "_")
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
grabProbes <- function(your_genes, RGSet, region){
  stopifnot(is.character(your_genes)); listRows=10000
  bp <- ifelse((length(region)>1), paste(region, collapse="|"), region); annot <- minfi::getAnnotation(RGSet); g=NULL
  geneAnno <- paste(annot$UCSC_RefGene_Name, annot$UCSC_RefGene_Group, sep = "_")
  pick <- function(x, fltr){return(do.call(grepl,c(list(pattern=x), fltr)))}
  z <- foreach::foreach(g = your_genes, .combine = "cbind") %dopar% {
    fltr <- list(x = geneAnno, ignore.case = T); y <- cbind(rownames(annot)[pick(g, fltr) & pick(bp, fltr)])
    naFill <- rbind(y, cbind(rep(NA, (listRows - nrow(y)))));colnames(naFill) <- g
    return(naFill)
  }
  naDrops <- lapply(1:ncol(z), function(i){which(!is.na(z[,i]))})
  cutt <- length(naDrops[[which.max(lapply(naDrops, function(x) sum(lengths(x))))]])
  return(z[1:cutt,])
}
                                           
col_fun3 <- circlize::colorRamp2(c(0, 0.20, 0.25, 0.5, 0.75, 1), c("black","darkblue","deepskyblue", "white", "tomato","red"))
                                           
saveHmPng <- function(fi_prefix, fi_suffix, hm, topvar = "", outDir=NULL) {
  if(is.null(outDir)){outDir<-gb$runDir}
  hmOutDir <- file.path(outDir, "HeatMaps")
  if(!dir.exists(hmOutDir)){dir.create(hmOutDir)}
  imgFile <- file.path(hmOutDir, paste0(fi_prefix, topvar, fi_suffix))
  wd <- as.numeric(hm@ht_list_param[["width"]]) + 5
  ht <- as.numeric(hm@ht_list_param[["height"]]) + 5
  png(
    file = imgFile,
    width = wd,
    height = ht,
    units = "mm",
    res = 200
  )
  ComplexHeatmap::draw(hm)
  invisible(dev.off())
}  
                                           
                                           
gb$getHeatMap2 <-
  function(betaRanges,
           titleValue,
           ha,
           geneNamesHeatMap = F,
           colSplt = NULL,
           rwsplt = NULL,
           rwOrder = NULL, 
           clusRows = F,
           clusCols = F) {
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
    
    outDir <- file.path(getwd(), "data")
    outFile <- paste(pathwayName, "avgBetas_per_gene.csv" , sep = "_")
    outPath <- file.path(outDir, outFile)
    write.csv(avgBetas, file = outPath, row.names=F, na="")
    return(avgBetas)
}

#' grabProbes function that lists probes annotated by matching gene name and gene region
#' @param your_genes  a character list of genes i.e. from a file c("ALK","ETMR","HDAC")
#' @param RGSet  the name or path to your minfi RGset data
#' @param region string vector of gene selections like c("TSS200", "Body") can one gene region "TSS200"
GetProbesGenes <- function(your_genes, RGSet, region){
  stopifnot(is.character(your_genes)); listRows=10000
  bp <- ifelse((length(region)>1), paste(region, collapse="|"), region); annot <- minfi::getAnnotation(RGSet); g=NULL
  geneAnno <- paste(annot$UCSC_RefGene_Name, annot$UCSC_RefGene_Group, sep = "_")
  pick <- function(x, fltr){return(do.call(grepl,c(list(pattern=x), fltr)))}
  z <- foreach::foreach(g = your_genes, .combine = "cbind") %dopar% {
    fltr <- list(x = geneAnno, ignore.case = T); y <- cbind(rownames(annot)[pick(g, fltr) & pick(bp, fltr)])
    naFill <- rbind(y, cbind(rep(NA, (listRows - nrow(y)))));colnames(naFill) <- g
    return(naFill)
  }
  naDrops <- lapply(1:ncol(z), function(i){which(!is.na(z[,i]))})
  cutt <- length(naDrops[[which.max(lapply(naDrops, function(x) sum(lengths(x))))]])
  return(z[1:cutt,])
}
                                           
GetCsvGeneColumns <- function(pathwayName, z){
    outFile1 <- paste(pathwayName,"genes_columns_withMissing.csv", sep = "_")
    outDir <- file.path(getwd(),"data")
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

GetHeatMapGenes <-function(betaRanges, titleValue, ha, geneNamesHeatMap=F, colSplt = NULL, rwsplt=NULL){
    titleOfPlot <- paste("Heatmap of",titleValue,sep = " ")
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
        row_title_gp = gpar(fontsize = 12, fontface = "bold"),
        show_row_dend = F,
        show_column_dend = T,
        use_raster=T,
        show_heatmap_legend = T,
        top_annotation = ha,
        column_title = titleOfPlot,
        column_title_gp = gpar(fontsize = 14,fontface = "bold"),
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
        heatmap_height = unit(25, "in")
    )
    #size = gb$calc_ht_size(hmTopNumbers)
    #size
    return(gb$drawHeatMap(hmTopNumbers))
}
