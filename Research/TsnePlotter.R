gb <- globalenv(); assign("gb", gb)
# FUN: Returns optimal number of PCs given X ---------------------------
library(dplyr)
require('foreach')
library(ggplot2)
getOptPC <- function(expr, getPval = F, n= 60){
    PC <- prcomp(t(log10(expr + 1)), center = T, scale = F)
    expl_var <- PC$sdev ^ 2 / sum(PC$sdev ^ 2)
    expl_var_perm <- matrix(NA, ncol = length(PC$sdev), nrow = 10)
    expr_perm <- apply(expr, 2, sample)
    PC_perm <- prcomp(t(log10(expr_perm + 1)), center = T, scale = F)
    for (k in 1:10) {expl_var_perm[k,] <- (PC_perm$sdev ^ 2) / sum((PC_perm$sdev) ^ 2)}
    pval <- apply(t(expl_var_perm) >= expl_var, 1, sum) / 10
    if (getPval) {return(pval)} else
        optPC <- head(which(pval >= 0.05), 1) - 1
    return(optPC)
}

# FUN: Returns optimal perplexity given X
getPerplexity <- function(expr) {
    N_cells <- dim(expr)[2]
    return(round(sqrt(N_cells), 0))
}

generateTvals <- function(betas) {
    betas <- as.data.frame(betas)
    X <- t(betas)
    optPC <- getOptPC(X)
    opPer <- getPerplexity(betas)
    
    tryCatch(
        expr = {
            TSNE <- Rtsne::Rtsne(X, dims = 3, theta = 0, perplexity = opPer,
                    initial_dims = optPC, verbose = F, partial_pca=T, max_iter = 10000)
            return(TSNE)
        },
        error = function(e) {
            TSNE <- Rtsne::Rtsne(X, dims = 3, theta = 0, perplexity = opPer-1,
                    initial_dims = optPC, verbose = F, partial_pca=T, max_iter = 10000)
            return(TSNE)
        }
    )
}

# Returns TSNE x,y,z coordinates and groups labels colors
getTsneVal <- function(TSNE, saNames, samGrp, colorGrp, symGrp) {
        tsne_plot <- data.frame(
            x = TSNE$Y[, 1],
            y = TSNE$Y[, 2],
            z = TSNE$Y[, 3],
            GROUPS = samGrp, col = colorGrp,
            samples = saNames, symbol = symGrp)
        return(tsne_plot)
 }

# Final TSNE plot title will combine: clusterTitle + titleMain + varProbes
generateTitles <- function(clusterTitle, topTitle, titlemain) {
    library('foreach'); require('foreach')
    tsne_titles <- as.data.frame(foreach::foreach(
        ti = 1:length(clusterTitle), .combine = "rbind", .packages = c('foreach')
    ) %do% {
        foreach::foreach(
            topV = 1:length(topTitle), .combine = "rbind", .packages = c('foreach')
        ) %do% {
            return(paste0(clusterTitle[ti], " ", titlemain, " ", topTitle[topV], " Variance Probes"))
        }
    })[, 1]
    return(tsne_titles)
}

# FUNCTION: Appends ggrepel labels to T-sne plot
addPlotLabels <- function(groupTsne, tsneData){
    groupTsne <- groupTsne + ggrepel::geom_label_repel(
        data = tsneData,
        aes(x = x,  y = y, label = samples, size = 2),
        alpha = 0.85, segment.alpha = 0.70,
        nudge_x = -30, nudge_y = 4, direction = "both", fontface = "bold",
        box.padding = 0.5, fill = "pink", min.segment.length = 0.01,
        color = "black", label.size = 1.0, size = 3,
        label.padding = unit(0.5, "lines"), label.r = unit(0.5, "lines"),
        force = 8, show.legend = F)
    
    return(groupTsne)
    
}


AddPlotSymbols <- function(tsne_plot, groupTsne, symFlags){
  if(symFlags==F) {
    symShape <- shapeVals <- shapeLabels <- NULL
    groupTsne <- groupTsne +
      geom_point(aes(x, y, color = tsne_plot$GROUPS, shape = symShape), size = 5, alpha = 0.85)
    return(groupTsne)
  }else{
    shapeVals <- c(19, 17, 15, 7, 8, 9, 1, 3, 4, 5)
    plotSymLen <- 1:length(unique(tsne_plot$symbol))
    if (length(plotSymLen) > 10) {
      shapeVals <- c(1:25, 33:127)
    }
    sv <- shapeVals[plotSymLen]
    shapeLabels <- levels(as.factor(tsne_plot$symbol))
    if (any(is.na(shapeLabels))) {
      shapeLabels <- colorLabel
    }
    symShape <- as.factor(tsne_plot$symbol)
    groupTsne <- groupTsne +
      geom_point(aes(x, y, color = tsne_plot$GROUPS, shape = symShape), size = 5, alpha = 0.85) +
      scale_shape_manual(name = "Sample Type", values = sv, labels = shapeLabels)
    return(groupTsne)
    }
}


# FUNCTION: Generates a T-sne plot with or without additional labels
genTsnePlot <- function(tsne_plot, titleLabel, groupToLabel = NULL,
                        symbolsLabel = NULL, colorLabel = NULL, names2Label = NULL)
    {
    require('ggplot2')
    colours <- unique(tsne_plot$col)
    symFlags <- !is.null(symbolsLabel)
    devAskNewPage(ask=F) #options("device.ask.default"=F)
    # Parameters for text & geom_label_repel
    et <- ggplot2::element_text(size = 12)
    eBlank <- ggplot2::element_blank()
    # Creating Main ggplot Object
    groupTsne <- ggplot(tsne_plot, aes(x=tsne_plot$x,y=tsne_plot$y,group=tsne_plot$GROUPS)) 
    # Adding Symbols if provided
    groupTsne <- AddPlotSymbols(tsne_plot, groupTsne, symFlags)
    # Adding Plot Color, Theme, and Axis Labels
    groupTsne <- groupTsne +
      scale_color_manual(values = colours, name = "Sample Label") +
      labs(color = colorLabel, size = 5, x = "TSNE 1", y = "TSNE 2") + theme_bw(base_size = 16) +
      theme(panel.grid = eBlank, panel.background = eBlank, text = et, plot.title = et, 
            legend.text = et, legend.position = "bottom", axis.text = et) +
      ggtitle(label = titleLabel) +
      guides(fill = guide_legend(title = "Sample Legend", ncol = 2, nrow = 30, byrow = TRUE))
    # Below only runs to label sample IDs if label group is provided
    if (!is.null(groupToLabel)) {
        tsneData = subset(tsne_plot, tsne_plot$samples == groupToLabel)
        groupTsne <- addPlotLabels(groupTsne, tsneData)
    }
    # Below runs to label specific sample IDs if provided
    if (!is.null(names2Label)) {
        groupTsne <- addPlotLabels(groupTsne, tsne_plot)
    }
    return(groupTsne)
}


gb$tierBetas <- function(betas, col_sentrix, RGSet,
                         batchCorrect = F, getSuper = F, topVar = 1:10000) {
    rgLiDat <- RGSet@colData@listData
    selectSams <- rgLiDat[[col_sentrix]][rgLiDat[["Sample_ID"]] %in% colnames(betas)]
    rgColRows <- RGSet@colData@rownames
    newRgset <- RGSet[, rgColRows %in% selectSams]
    if (batchCorrect == T) {
      if (getSuper == T) {
        superbetas <- gb$batchCorrectBs(betas, newRgset, topVar, T)
        return(superbetas)
      } else{
        unBetas <- gb$batchCorrectBs(betas, newRgset , topVar)
        return(unBetas)
      }
    } else{
      if (getSuper == T) {
        superbetas <- gb$getSupervise(betas, newRgset, topVar)
        return(superbetas)
      } else{
        unBetas <- gb$takeTopVariance(betas, topVar)
        return(unBetas)
      }
    }
}

gb$selectPlots <- function(doPlotly=F,tplots,ty,tps,outDirs){
  if (doPlotly == F) {
      for (zz in 1:nrow(outDirs)) {
        xx <- paste('###', tps[ty, 1], outDirs[zz, 3], '\n\n')
        cat(xx)
        fig <- tplots[[zz]]
        supM(print(fig))
        cat('\n\n')
      }
    } else{
      for(zz in 1:length(tplots))  {
        xx <- paste('###', tps[ty, 1], outDirs[zz, 3], '\n\n')
        cat(xx)
        fig <- tplots[[zz]]
        op <- gb$makePlotly(fig)
        supM(print(htmltools::tagList(ggplotly(op))))
        cat('\n\n')
      }
    }
  return(assign("diagPlot",tplots[[1]]))
}

doMultiple <- function(allBetas1,tsne_titles, outDirs, targets1, tps,ty,custom){
  plotN = NULL
  
  tsneList <- lapply(X = 1:length(allBetas1), FUN = function(X) {
      return(suppressMessages(gb$generateTvals(allBetas1[[X]])))
    })
  
  plotList <-
    list(foreach::foreach(plotN = 1:length(tsneList), .packages = "foreach") %do%
           {
             ###################### TO CHANGE ########################
             gc(verbose = F)
             tsne_plot <- gb$getTsneVal(
               TSNE = tsneList[[plotN]],
               saNames = targets1$SampleFilter,
               samGrp = targets1$PointColors,
               colorGrp = targets1$color, #targets1$Type,#
               symGrp = targets1$Sym_Shape#targets1[, gb$col_shapes]
             ) #targets1[,col_samGroup]
             ###################### TO CHANGE ########################
             return(tsne_plot)
           })[[1]]
  return(plotList)
}

getTopPlot <- function(samNames){
    mds <- limma::plotMDS(
        gb$mSetSq.beta,
        top = 1000,
        gene.selection = "common",
        plot = T,
        ndim = 3
    )
    toplot_Histo <- data.frame(
        Dim1 = mds$cmdscale.out[, 1],
        Dim2 = mds$cmdscale.out[, 2],
        Dim3 = mds$cmdscale.out[, 3],
        Sample = samNames
    )
    return(toplot_Histo)
}

GenerateUnsuperTsne <- function(targets1, betas, gb, colorVariable = NULL, shapeVariable = NULL, sampleGrouping = "Sample_Group", isSuper=NULL){
    if(is.null(colorVariable)){
        colorVariable <- gb$col_samTypes
    }
    if(is.null(shapeVariable)){
        shapeVariable <-  gb$col_shapes
    }
    if(is.null(isSuper)){
        isSuper <- F
    }
    if(isSuper == F){
        tsne_titles = gb$tsne_titles[1:3]
    }else{
        tsne_titles = gb$tsne_titles[4:6]
        betas <- gb$superbetas
    }
    gb$subsetBetas(
        targFilter = sampleGrouping,
        samGroup = colorVariable,
        betas = betas,
        targets = targets1,
        samShapes = shapeVariable,
        samNames = gb$col_samNames,
        tsne_titles = tsne_titles,
        names2Label = gb$names2Label,
        supervised = isSuper
    )
    invisible(gc(verbose = F))
}

