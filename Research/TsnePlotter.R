gb <- globalenv(); assign("gb", gb)
# FUN: Returns optimal number of PCs given X ---------------------------
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


# FUN: Generate T-sne values given betas --------
generateTvals <- function(betas) {
    betas <- as.data.frame(betas)
    X <- t(betas)
    optPC <- getOptPC(X)
    opPer <- getPerplexity(betas)
    
    tryCatch(
        expr = {
            TSNE <-
                Rtsne::Rtsne(
                    X,
                    dims = 3,
                    theta = 0,
                    perplexity = opPer,
                    initial_dims = optPC,
                    verbose = T,
                    max_iter = 10000
                )
            message(
                "\nOptimal Perplexity for T-sne is: ",
                opPer,
                "\n",
                "Optimal #PCs are: ",
                optPC,
                "\n"
            )
            return(TSNE)
        },
        error = function(e) {
            
            TSNE <-
                Rtsne::Rtsne(
                    X,
                    dims = 3,
                    theta = 0,
                    perplexity = opPer-1,
                    initial_dims = optPC,
                    verbose = T,
                    max_iter = 10000
                )
            message(
                "\nOptimal Perplexity for T-sne is: ",
                opPer,
                "\n",
                "Optimal #PCs are: ",
                optPC,
                "\n"
            )
            return(TSNE)
        }
    )
}

genTsnePlot <- function(tsne_plot, titleLabel, groupToLabel = NULL, symbolsLabel=NULL, colorLabel=NULL,names2Label=NULL){
    col_vect <- pals::glasbey()
    colours <- col_vect[1:(length(unique(tsne_plot$GROUPS)))]
    colours[6] = "#eb7d34" #changing dark forest to orange color
    symFlags <- !is.null(symbolsLabel)
    #options("device.ask.default"=F)
    devAskNewPage(ask=F)
    if(symFlags==T){
        shapeVals <- c(17, 19, 15, 7, 8, 9, 1, 3, 4, 5)
        sv <- shapeVals[1:length(unique(tsne_plot$symbol))]
        shapeLabels <- "Sample Type" #levels(as.factor(tsne_plot$symbol))
        symShape <- as.factor(tsne_plot$symbol)
    }else{symShape <- shapeVals <- shapeLabels <- NULL}
    et <- element_text(size = 24)
    groupTsne <- ggplot(tsne_plot,aes(x=tsne_plot$x,y=tsne_plot$y,group=tsne_plot$GROUPS)) +
        geom_point(aes(x,y,color=tsne_plot$GROUPS, shape= symShape), size=4,alpha=0.85)
    if(symFlags==T){
        groupTsne <- groupTsne + scale_shape_manual(name=symbolsLabel,values=sv,labels=shapeLabels)
    }
    groupTsne <- groupTsne +
        scale_color_manual(values=colours) +
        labs(color=colorLabel,size=4) + theme_bw(base_size=24) +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
              panel.background=element_blank(),text=et) +
        ggtitle(label=titleLabel) +
        theme(plot.title=et,legend.text=et,text=et,legend.position="right",axis.text.y=et,axis.text.x=et) +
        labs(x="TSNE 1",y="TSNE 2")
    # Below only runs to label sample IDs if label group is provided
    if (!is.null(groupToLabel)) {
        groupTsne <- groupTsne + ggrepel::geom_label_repel(
            data = subset(tsne_plot, tsne_plot$samples == groupToLabel),
            aes(x = x,  y = y, label = samples, size = 2),
            alpha = 0.85, segment.alpha = 0.70, nudge_x = -30, nudge_y = 4,
            direction = "both", fontface = "bold",
            box.padding = 0.5, fill = "pink", min.segment.length = 0.01, color = "black",
            label.size = 1.0, size = 4, label.padding = unit(0.5, "lines"),
            label.r = unit(0.5, "lines"), force = 8, show.legend = F)
    }
    
    if (!is.null(names2Label)) {
        groupTsne <- groupTsne + ggrepel::geom_label_repel(
            data = tsne_plot,
            aes(x = x,  y = y, label = samples, size = 2),
            alpha = 0.85, segment.alpha = 0.70, nudge_x = -30, nudge_y = 4,
            direction = "both", fontface = "bold",
            box.padding = 0.5, fill = "pink", min.segment.length = 0.01, color = "black",
            label.size = 1.0, size = 4, label.padding = unit(0.5, "lines"),
            label.r = unit(0.5, "lines"), force = 8, show.legend = F)
    }
    return(groupTsne)
}

 getTsneVal <- function(TSNE, saNames, samGrp, colorGrp, symGrp) {
        tsne_plot <- data.frame(
            x = TSNE$Y[, 1],
            y = TSNE$Y[, 2],
            z = TSNE$Y[, 3],
            GROUPS = samGrp, col = colorGrp,
            samples = saNames, symbol = symGrp)
        return(tsne_plot)
 }

generateTitles <- function(clusterTitle, topTitle, titlemain) {
  tsne_titles <- as.data.frame(foreach::foreach(
    ti = 1:length(clusterTitle),
    .combine = "rbind",
    .packages = 'foreach'
  ) %do% {
    foreach::foreach(
      topV = 1:length(topTitle),
      .combine = "rbind",
      .packages = 'foreach'
    ) %do% {
      return(paste0(
        clusterTitle[ti],
        " ",
        titlemain,
        " ",
        topTitle[topV],
        " Variance Probes"
      ))
    }
  })[, 1]
  print(tsne_titles)
  return(tsne_titles)
}

getTopPlot <- function(samNames){
          mds <-
            limma::plotMDS(
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
# FUN: Generate T-sne plot given TSNE values --------
plotSaver <- function(outDirs,tsne_titles,tps,ty,plotList,custom) {
  plotN=NULL
  options("device.ask.default"=F)
  invisible(foreach::foreach(plotN = 1:nrow(tps)) %do% {
    pL<- plotList[[plotN]]
    message("Plot to Render:")
    print(plotList[[plotN]])
    thePlot <- gb$genTsnePlot(
      tsne_plot=plotList[[plotN]],
      titleLabel=tsne_titles[plotN],
      symbolsLabel = pL$GROUPS,
      colorLabel = pL$symbol,
      names2Label = pL$samples
      )
    outDir <- paste(outDirs[plotN, ], collapse = "_")
    fn <- paste(outDir, colnames(tps), tps[ty, 1], "tnse.png", sep = "_")
    devAskNewPage(ask=F)
    ggplot2::ggsave(
      filename = fn,
      device = "png",
      path = getwd(),
      plot = thePlot,
      width = 15,
      height = 15,
      units = "in"
    )
    options("device.ask.default"=F)
    message("Plot saved:")
    cat(fn, sep = "\n")
  })
}
