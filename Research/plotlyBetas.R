gb <- globalenv(); assign("gb", gb)
require("ggplot2")
printPlotTab <- function(thePlot, theTabName){
  cat(theTabName)
  plot.new()
  return(supM(print(thePlot)))
}

grabPngNames <- function(tsne_titles=NULL, keywrd="Top"){
  outDirs <- stringr::str_split_fixed(tsne_titles, " ", 16)
  colNumb <- grep(keywrd, outDirs[1,], ignore.case = T)
  outDirs <- outDirs[, c(1, colNumb, colNumb + 1)]
  return(outDirs)
}


FormatLegendText <- function(fig){
  #gb$SaveObj(fig, "test_fig.Rdata")
  #fig <- gb$LoadRdatObj("test_fig.Rdata")
  #stopifnot(FALSE)
  fig2 <- fig + geom_point(aes(x, y, color = "Shape", shape = symbol), color = fig$data$col, stroke = 2, size = 10)
  otherPlot <- gb$supM(plotly::ggplotly(
      fig2, dynamicTicks = T, width = 1200, height = 800, source = "A", layerData = 1))
    figGrps <- paste(fig$data$GROUPS, fig$data$symbol, sep = ",") #fig$data$x, fig$data$y,
    otherPlot[["x"]][["layout"]][["shapes"]][[1]][["line"]][["width"]] <- 2
    for (nSam in 1:length(otherPlot[["x"]][["data"]])) {
        lgndGrp <- otherPlot[["x"]][["data"]][[nSam]][["legendgroup"]]
        nlgndGrp <- stringr::str_remove_all(lgndGrp, "[()]")
        lgndSplt <- stringr::str_split_fixed(nlgndGrp, ",", 2)
        isShape <- lgndSplt[, 2] == "NA"
        otherPlot[["x"]][["data"]][[nSam]][["marker"]][["opacity"]] <- 1
        otherPlot[["x"]][["data"]][[nSam]][["marker"]][["size"]] <- 18
        otherPlot[["x"]][["data"]][[nSam]][["marker"]][["line"]][["width"]] <- 2
        otherPlot[["x"]][["data"]][[nSam]][["marker"]][["line"]][["opacity"]] <- 1
        if (isShape) {
            nGroupSplit <- stringr::str_split_fixed(nlgndGrp, ",", 2)[1, 1]
            theSamLabs <- which(lgndSplt[1,1]==fig$data$symbol)
            otherPlot[["x"]][["data"]][[nSam]]$text <- paste("Sample:", fig$data$samples[theSamLabs])
            otherPlot[["x"]][["data"]][[nSam]]$visible <- 'legendonly'
            otherPlot[["x"]][["data"]][[nSam]][["marker"]][["line"]][["color"]] <- "grey"
            otherPlot[["x"]][["data"]][[nSam]][["marker"]][["color"]] <- "darkgrey"
        } else{
            nGroupSplit <- stringr::str_split_fixed(nlgndGrp, ",", 2)[1, 1]
            theSamLabs <- which(nlgndGrp==figGrps)
            otherPlot[["x"]][["data"]][[nSam]]$text <- paste("Sample:", fig$data$samples[theSamLabs])
            if (otherPlot[["x"]][["data"]][[nSam]][["marker"]][["symbol"]] != "circle") {
                otherPlot[["x"]][["data"]][[nSam]][["showlegend"]] <- F
                otherPlot[["x"]][["data"]][[nSam]][["marker"]][["line"]] <- list(width = 2) 
                fillColor <- otherPlot[["x"]][["data"]][[nSam]][["marker"]][["color"]]
                colSplit <- stringr::str_split_fixed(fillColor, ",", 4)
                colSplit[1,4] <- "1)"
                otherPlot[["x"]][["data"]][[nSam]][["marker"]][["line"]][["color"]] <- paste(colSplit, collapse = ",")
            }
        }
        otherPlot[["x"]][["data"]][[nSam]][["name"]] <-
            otherPlot[["x"]][["data"]][[nSam]][["legendgroup"]] <- nGroupSplit
    }
    return(otherPlot)
}


FormatPlotlyLayout <- function(otherPlot) {
    otherPlot$x[["layout"]][["annotations"]] <- NULL
    opLayout <- otherPlot[["x"]][["layout"]]
    opLayout[["font"]][["size"]] <- 12
    opLayout[["title"]][["font"]][["size"]] <- 20
    opLayout[["xaxis"]][["title"]][["font"]][["size"]] <- 14
    opLayout[["yaxis"]][["title"]][["font"]][["size"]] <- 14
    opLayout[["xaxis"]][["tickfont"]][["size"]] <- 12
    opLayout[["yaxis"]][["tickfont"]][["size"]] <- 12
    otherPlot[["x"]][["layout"]] <- opLayout
    otherPlot[["x"]][["layout"]][["legend"]][["title"]][["text"]] <- "Sample Color & Shape Legend"
    pltTtl <- otherPlot[["x"]][["layout"]][["title"]][["text"]]
    otherPlot[["x"]][["layout"]][["margin"]]$t <- 100
    return(otherPlot)
}


GetPlotlySymbols <- function(fig, uniGrp){
    markerSyms <- c(
        "circle", "square", "diamond", "triangle-up", "cross", "X", "triangle-down", "triangle-left", 
        "triangle-right", "triangle-ne", "triangle-se", "triangle-sw", "triangle-nw", "pentagon",
        "hexagon", "hexagon2", "octagon", "star", "hexagram", "star-triangle-up", "star-triangle-down",
        "star-square", "star-diamond", "diamond-tall", "diamond-wide", "hourglass", "bowtie", "circle-cross",
        "circle-x", "square-cross", "square-x", "diamond-cross"
    )
    markerSyms <- c(markerSyms, paste0(markerSyms,"-open"))
    if(length(unique(fig[["data"]]$symbol))>1){
        markerSyms <- markerSyms[1:length(uniGrp)]
    }else{
        markerSyms <- rep("circle", length(uniGrp))
    }
    return(markerSyms)
}


FormatHoverInfo <- function(otherPlot){
    opInfo <- otherPlot[["x"]][["data"]]
    for (sam in 1:length(opInfo)) {
        opInfo[[sam]][["hoverinfo"]] <- "none"
        opInfo[[sam]][["marker"]][["size"]] <- 18
        opInfo[[sam]][["hoverinfo"]] <- c("text")
    }
    otherPlot[["x"]][["data"]] <- opInfo
    return(otherPlot)
}


FormatPlotLabels <- function(fig, otherPlot, uniGrp, markerSyms){
  figDat <- fig[["data"]]
  for (grpIdx in 1:length(otherPlot$x$data)) {
    grpNam <- otherPlot$x$data[[grpIdx]]$legendgroup
    samNames <-
      stringr::str_split(otherPlot[["x"]][["data"]][[grpIdx]]$text,
                         "Sample: ", simplify = T)[, 2]
    theIdx <- match(samNames, fig$data$samples)
    currGrp <- figDat$symbol == grpNam
    newText <- paste(
      paste(otherPlot[["x"]][["data"]][[grpIdx]]$text),
      paste("x:", round(figDat$x[theIdx], 2)),
      paste("y:", round(figDat$y[theIdx], 2)),
      paste("Group:", figDat$GROUPS[theIdx]),
      paste("Symbol:", figDat$symbol[theIdx]),
      paste("Color:", figDat$col[theIdx]),
      #paste("NAME:", figDat$samples[theIdx]),
      sep = "<br />"
    )
    if (grpNam == otherPlot[["x"]][["data"]][[grpIdx]]$name) {
      otherPlot[["x"]][["data"]][[grpIdx]]$text <- newText
    }
  }
  return(otherPlot)
}


FormatPlotlyLegend <- function(otherPlot) {
  mrg <- list(
    l = 50,
    r = 50,
    b = 100,
    t = 100,
    pad = 4
  )
  ttl <-
    list(text = "<b>Shapes and Colors Legend</b><br>", font = list(size = 14))
  lgnd <- list(title = ttl, font = list(size = 14))
  otherPlot <- otherPlot %>% plotly::layout(margin = mrg, legend = lgnd)
  return(otherPlot)
}


makePlotly <- function(fig) {
    otherPlot <- NULL
    otherPlot <- FormatLegendText(fig)
    otherPlot <- FormatPlotlyLayout(otherPlot)
    otherPlot <- FormatHoverInfo(otherPlot)
    uniGrp <- unlist(lapply(
        X = 1:length(otherPlot[["x"]][["data"]]), FUN = function(X) {return(otherPlot[["x"]][["data"]][[X]]$name)}))
    markerSyms <- GetPlotlySymbols(fig, uniGrp)
    otherPlot <- FormatPlotLabels(fig, otherPlot, uniGrp, markerSyms)
    otherPlot <- FormatPlotlyLegend(otherPlot)
    return(otherPlot)
}


GetFlatPlots <- function(fig){
    fig <- fig +
        theme(
            legend.direction = "vertical",
            legend.margin = ggplot2::margin(t = -25),
            legend.box.margin = ggplot2::margin(0, 0, 0, 0),
            legend.justification = "right",
            legend.position = "right"
        )
    leg <- supM(cowplot::get_legend(fig))
    if(length(leg)>10){
        fig <- fig + theme(legend.position = "none")
        supM(print(fig))
        grid::grid.newpage()
        grid::grid.draw(leg)
    }else{
        fig <- fig +
            theme(legend.direction = "vertical", legend.margin = ggplot2::margin(t = 0))
        supM(print(fig))
    }
}


selectPlots <- function(doPlotly = F, tplots, ty, tps, outDirs) {
    for (zz in 1:nrow(outDirs)) {
        invisible(gc())
        fig <- tplots[[zz]]
        tabStart <- paste('### Top', tps[ty, 1], outDirs[zz, 3])
        
        cat(paste(tabStart, '(Interactive)','\n\n'))
        gc(verbose = F)
        op <- gb$makePlotly(fig)
        supM(print(htmltools::tagList(ggplotly(op))))
        cat('\n\n')
        
        cat(paste(tabStart, '\n\n'))
        gc(verbose = F)
        GetFlatPlots(fig)
        cat('\n\n')
    }
    #return(assign("diagPlot", tplots[[1]]))
}


gb$grabAllBeta <- grabAllBeta <- function(targets1, betas, supervised = F) {
    if (supervised == T) {
        supBets <- betas[, targets1$SampleFilter] # filtering betas
        allBetas1 <- list(supBets[1:100,], supBets[1:1000,], supBets)
    } else{
        betas1 <- betas[, targets1$SampleFilter] # filtering betas
        if (file.exists(file.path(gb$runDir, gb$unbetaVariance))) {
            unBets <- gb$LoadRdatObj(gb$unbetaVariance, msgProg=F)
        } else{
            unBets <- gb$takeTopVariance(betas1, topVar = 1:10000)
            gb$SaveObj(unBets, file = file.path(gb$runDir, gb$unbetaVariance))
        }
        unBets <- unBets[, targets1$SampleFilter]
        allBetas1 <- list(unBets[1:100, ], unBets[1:1000, ], unBets)
    }
    return(allBetas1)
}


doMultiple <- function(allBetas1, targets1){
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
               symGrp = targets1$Sym_Shape #targets1[, gb$col_shapes]
             ) #targets1[,col_samGroup]
             ###################### TO CHANGE ########################
             return(tsne_plot)
           })[[1]]
  return(plotList)
}


plotSaver <- function(tsne_titles, plotList, names2Label = NULL) {
        plotN = NULL
        options("device.ask.default" = F)
        pltList <- foreach::foreach(plotN = 1:length(plotList), .packages = "foreach") %do% {
          pL <- plotList[[plotN]]
                ###################### MODIFY TO CUSTOMIZE ########################
                gc(verbose = F)
                return(
                    gb$genTsnePlot(
                        tsne_plot = plotList[[plotN]],
                        titleLabel = tsne_titles[plotN],
                        symbolsLabel = pL$symbol,
                        colorLabel = pL$GROUPS,
                        names2Label = names2Label # any specific sample names to label on the plot
                    )
                )
                ###################### MODIFY TO CUSTOMIZE ########################
        }
        return(pltList)
}


MessageTsneLoop <- function(targets, custom) {
    message("Current Sample Group TSNE: ",
            custom,
            "\nAll Target Group Filters: ",
            paste(unique(targets$SamGroups)))
}


FormatPlotTargets <- function(targets, samGroup, samShapes, samNames, targFilter){
    targets <- gb$colorTargets(targets, samGroup)
    targets$SamGroups <- targets$SampleFilter <- NULL
    targets$SampleFilter <- targets[, samNames] # creating new column
    targets$SamGroups <- targets[, targFilter]
    targets$Sym_Shape <- targets[, samShapes]
    targets$PointColors <- targets[, samGroup]
    colorColName <- paste0(samGroup, "_color")
    targets$color <- targets[, colorColName]
    return(targets)
}


subsetBetas <- function(targFilter, samGroup, betas, targets, samShapes, samNames,
                        tsne_titles, doPlotly = F, supervised = F, names2Label = NULL)
{
    targets <- FormatPlotTargets(targets, samGroup, samShapes, samNames, targFilter)
    tps <- as.data.frame(unique(targets$SamGroups))
    ty = NULL
    for (ty in 1:nrow(tps)) {
        gc(verbose = F)
        custom = tps[ty, 1]
        MessageTsneLoop(targets, custom)
        # Filter The Beta Values ---------
        targets1 <- targets[targets$SamGroups == custom, ]
        allBetas1 <- grabAllBeta(targets1, betas, supervised)
        outDirs <- gb$grabPngNames(tsne_titles)
        tplots <- plotList <- NULL
        # Get T-sne Values ---------
        plotList <- gb$doMultiple(allBetas1, targets1)
        gc(verbose = F)
        tplots <- gb$plotSaver(tsne_titles, plotList, names2Label)
        gb$selectPlots(doPlotly, tplots, ty, tps, outDirs)
    }
}


takeTopVariance <- function(betas, topVar = 1:10000){
    var_probes <- apply(betas, 1, var) # vars <- apply(gset.funnorm.beta,1,var)
    select_var <- names(sort(var_probes, decreasing = T)) # select_var <- names(sort(vars,decreasing = TRUE))
    sorted_betas <- betas[select_var, ] # top_variable_beta <- gset.funnorm.beta[select_var,]
    top_var_beta <- sorted_betas[topVar, ] # top_variable_beta[1:10000,]
    return(top_var_beta)
}

gb$tierBetas <- function(betas, col_sentrix, RGSet, batchCorrect = F, getSuper = F, topVar = 1:10000) {
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

SubsetTargets <- function(targets, varToFilter = NULL){
  if(is.null(varToFilter)){return(targets)}
  targets1 <- targets[varToFilter,] # exclude of targets are subsetted
  rownames(targets1) <- 1:nrow(targets1) # ensures rows are chronological
  return(targets1)
}

assign("subsetBetas", subsetBetas)
