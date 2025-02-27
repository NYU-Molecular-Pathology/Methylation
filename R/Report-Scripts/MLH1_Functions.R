#!/usr/bin/env Rscript
## ---------------------------
## Script name: MLH1_Functions.R
## Purpose: source of global functions used for calculation of methylation MLH1 values in report.Rmd
## Date Last Modified: January 19, 2024
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

GetMLH1Data <- function(ratioSet, RGset){
    #Output is DF with sample information, betas for relevant MLH1 loci
    MLH1_CGIDS_og = c("cg23658326","cg11600697","cg21490561","cg00893636")
    foundProbes <- MLH1_CGIDS_og %in% rownames(minfi::getBeta(ratioSet))
    MLH1_CGIDS <- MLH1_CGIDS_og[foundProbes]
    mlh.b = minfi::getBeta(ratioSet)[MLH1_CGIDS,]
    mlh.b <- as.data.frame(mlh.b)
    for (bal in 1:length(MLH1_CGIDS)) {
      if (is.nan(mlh.b[bal, 1])) {
        mlh.b[bal, 1] = 0.00
      }
      mlh.b[bal, 1] <- round(mlh.b[bal, 1], digits = 3)
      if (mlh.b[bal, 1] < 0.09) {
        mlh.b[bal, 1] = 0.00
      }
      if (mlh.b[bal, 1] > 1.00) {
        mlh.b[bal, 1] = 1.00
      } else{
        (mlh.b[bal, 1] = mlh.b[bal, 1])
      }
    }
    #Calculate detection of Pvalues using minfi method #c(MGMT_CGIDS,MLH1_CGIDS),]
    detection.pvals = minfi::detectionP(RGset)
    l.pvals = detection.pvals[c(MLH1_CGIDS),]
    l.pvals <- as.data.frame(l.pvals)
    for (lpv in 1:length(MLH1_CGIDS)) {
      if (is.nan(l.pvals[lpv, 1])) {
        l.pvals[lpv, 1] = 0.00
      }
      l.pvals[lpv, 1] <- round(l.pvals[lpv, 1], digits = 3)
      if (l.pvals[lpv, 1] < 0.09) {
        l.pvals[lpv, 1] = 0.00
      }
      if (l.pvals[lpv, 1] > 1.00) {
        l.pvals[lpv, 1] = 1.00
      } else {
        l.pvals[lpv, 1] = l.pvals[lpv, 1]
      }
    }
    #calculate number of methylated loci with beta values above cutoff for mlh1
    pos.loci = 0
    for (mlhProbe in rownames(mlh.b)) {
      if("cg23658326" == mlhProbe){
        pos.loci <- pos.loci + (mlh.b["cg23658326",] >= 0.18)
      }
      if("cg11600697" == mlhProbe){
        pos.loci <- pos.loci + (mlh.b["cg11600697",] >= 0.27)
      }
      if("cg21490561" == mlhProbe){
        pos.loci <- pos.loci + (mlh.b["cg21490561",] >= 0.13)
      }
      if("cg00893636" == mlhProbe){
        pos.loci <- pos.loci + (mlh.b["cg00893636",] >= 0.09)
      }
    }
    #pos.loci = (mlh.b[1,] >= 0.18) + (mlh.b[2,] >= 0.27) + (mlh.b[3,] >= 0.13) + (mlh.b[4,] >= 0.09)
    m.reslt = rep("", length(pos.loci))
    m.reslt[pos.loci == 4] = "POSITIVE"
    m.reslt[pos.loci >= 0 & pos.loci < 3] = "NEGATIVE"
    m.reslt[pos.loci == 3] = "INDETERMINATE"
    
    all.data = data.frame(
        Sample_Name = 1,
        cg23658326.MLH1 = mlh.b["cg23658326", ], cg11600697.MLH1 = mlh.b["cg11600697", ],
        cg21490561.MLH1 = mlh.b["cg21490561", ], cg00893636.MLH1 = mlh.b["cg00893636", ],
        MLH1.pos.loci = pos.loci, m.reslt = m.reslt,
        Pvalue.cg23658326 = l.pvals["cg23658326", ], Pvalue.cg11600697 = l.pvals["cg11600697", ],
        Pvalue.cg21490561 = l.pvals["cg21490561", ], Pvalue.cg00893636 = l.pvals["cg00893636", ], stringsAsFactors = F
    )
    return(all.data)
}

## Plot Beta densities for samples and plot ----------------------------------------------------------
getPlotTable = function(clin.data) {
    title.string = c(
        "MLH1 Result" = clin.data$m.reslt[1],
        "Total Positive Loci" = clin.data$MLH1.pos.loci[1]
    )
    title.string <- as.data.frame(title.string)
    colnames(title.string) = NULL
    return(title.string)
}

## clean up beta.matrix values ------------
sanitizeDense <- function(beta.matrix){
    beta.matrix <- as.data.frame(beta.matrix)
    beta.matrix[is.na(beta.matrix)] <- 0
    beta.matrix[beta.matrix < 0] <- 0
    beta.matrix[beta.matrix > 1] <- 1
    beta.matrix <- round(beta.matrix, digits = 3)
    density <- density(beta.matrix[, 1])
    density_df <- data.frame(x = density$x, y = density$y)
    return(density_df)
}

## Print out the Plotly for Mlh1 ----------------------------------------------------------
renderPlot <- function(beta.matrix, clin.res) {
    na_vals <- which(is.na(clin.res))
    if(length(na_vals) > 0){
      new_names <- colnames(clin.res)
      new_names[na_vals] <- paste(new_names[na_vals], "PROBE_MISSING", sep="_")
      colnames(clin.res) <- new_names
    }
    clin.res[is.na(clin.res)] <- 0
    clinVals <- as.data.frame(t(clin.res[2:5]))
    vals_df <- data.frame(
      x = clinVals[, 1],
      y = c(1, 0.9, 0.8, 0.7),
      Probe.Names = rownames(clinVals)
    )
    density_df <- sanitizeDense(beta.matrix)
    fig <- ggplot2::ggplot(vals_df) + geom_bar(aes(x, y, fill = Probe.Names), width = 0.05, stat = 'identity') +
      geom_line(density_df, mapping = aes(x = x, y = y)) +
      xlab('Beta Values') + ylab('Density') + theme(panel.background = element_blank())
    return(suppressWarnings(fig))
}

drawPlotTab <- function(mlhP){
    btso = c("bordered")
    kgb <- c("striped",font_size = 12, bootstrap_options = btso, position = "left")
    txtc = "text-align:center;"
    be = c(booktabs = T, escape = F, linesep = "")
    theMlhTab <- 
        mlhP %>% knitr::kable("html", be, align = 'clc') %>%
        kableExtra::kable_styling(kgb, full_width = F, position = "left") %>%
        kableExtra::column_spec(column = 2, background = "rgb(204, 255, 204)", extra_css = txtc) %>%
        kableExtra::column_spec(column = 2, bold = T, border_right = T)
    return(theMlhTab)
}

Mlh1Pipeline <- function(RGset) {
    Mset = minfi::preprocessIllumina(RGset, bg.correct = TRUE, normalize = "controls")
    ratioSet = ratioConvert(Mset, what = "both", keepCN = FALSE)
    clin.res <- GetMLH1Data(ratioSet, RGset)
    par(mar = c(5, 6, 4, 1) + .1)
    beta.matrix = getBeta(ratioSet)
    theMlhplot <- renderPlot(beta.matrix, clin.res)
    mlhP <- getPlotTable(clin.data = clin.res)
    return(list(
        "theMlhplot" = theMlhplot,
        "mlhP" = drawPlotTab(mlhP),
        "theValue" = clin.res
    ))
}

