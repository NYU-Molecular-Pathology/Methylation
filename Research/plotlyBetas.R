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


FormatMarkerData <- function(otherPlot) {
    numList <- 1:length(otherPlot[["x"]][["data"]])
    for (nSam in numList) {
        fillColor <- otherPlot[["x"]][["data"]][[nSam]][["marker"]][["color"]]
        if (any(fillColor == 'darkgrey') & length(fillColor) == 1) {
            lgndGrp <- otherPlot[["x"]][["data"]][[nSam]][["legendgroup"]]
            otherPlot[["x"]][["data"]][[nSam]][["showlegend"]] <- F
            for (idx in numList[-nSam]) {
                lgndGrp2 <- otherPlot[["x"]][["data"]][[idx]][["legendgroup"]]
                if (lgndGrp == lgndGrp2) {
                    otherPlot[["x"]][["data"]][[nSam]][["marker"]][["color"]] <-
                        otherPlot[["x"]][["data"]][[idx]][["marker"]][["color"]]
                    otherPlot[["x"]][["data"]][[nSam]][["marker"]][["line"]][["color"]] <-
                        otherPlot[["x"]][["data"]][[idx]][["marker"]][["line"]][["color"]]
                }
            }
        }
    }
    return(otherPlot)
}


FormatShapeData <- function(isSingle, otherPlot, nSam) {
    fillColor <- otherPlot[["x"]][["data"]][[nSam]][["marker"]][["color"]]
    otherPlot[["x"]][["data"]][[nSam]][["marker"]][["line"]] <- list(width = 2)
    colSplit <- stringr::str_split_fixed(fillColor, ",", 4)
    colSplit[1, 4] <- "1)"
    otherPlot[["x"]][["data"]][[nSam]][["marker"]][["line"]][["color"]] <- paste(colSplit, collapse = ",")
    if (any(fillColor == 'darkgrey') & length(fillColor) == 1) {
        otherPlot[["x"]][["data"]][[nSam]][["showlegend"]] <- F
    }
    return(otherPlot)
}


FormatShapeColors <- function(isSingle, nlgndGrp, lgndSplt, fig, otherPlot, nSam) {
    if (isSingle == F) {
        nGroupSplit <- stringr::str_split_fixed(nlgndGrp, ",", 2)[1, 1]
        theSamLabs <- which(lgndSplt[1,1] == fig$data$symbol)
        otherPlot[["x"]][["data"]][[nSam]]$text <- paste("Sample:", fig$data$samples[theSamLabs])
        otherPlot[["x"]][["data"]][[nSam]]$visible <- 'legendonly'
        otherPlot[["x"]][["data"]][[nSam]][["marker"]][["line"]][["color"]] <- "black"
        otherPlot[["x"]][["data"]][[nSam]][["marker"]][["color"]] <- "darkgrey"
        # Hides the color legend
        for (i in seq_along(otherPlot[["x"]][["data"]])) {
            if (!is.null(otherPlot[["x"]][["data"]][[i]][["name"]])) {
                fillColor <- otherPlot[["x"]][["data"]][[i]][["marker"]][["color"]]
                if (fillColor != otherPlot[["x"]][["data"]][[nSam]][["marker"]][["color"]]) {
                    otherPlot[["x"]][["data"]][[i]][["showlegend"]] <- F
                }
            }
        }
        return(otherPlot)
    } else {
        otherPlot[["x"]][["data"]][[nSam]]$visible <- F
        return(otherPlot)
    }
}


#gb$SaveObj(fig, "test_fig.Rdata")
#fig <- gb$LoadRdatObj("test_fig.Rdata")
#stopifnot(FALSE)

FormatLegendText <- function(fig){

    fig2 <- fig + 
        geom_point(aes(
            fig$data$x, fig$data$y, color = "Shape", shape = fig$data$symbol),
            color = fig$data$col, stroke = 2, size = 10
            )
    fig2 <- fig2 + guides(shape = "none")
    
    if (length(fig2[["layers"]]) > 1) {
        fig2$layers <- fig2$layers[-2]
    }

    otherPlot <- gb$supM(plotly::ggplotly(fig2, dynamicTicks = T, width = 1200, height = 800,
                                          source = "A", layerData = 1))
    figGrps <- paste(fig$data$GROUPS, fig$data$symbol, sep = ",")

    otherPlot[["x"]][["layout"]][["shapes"]][[1]][["line"]][["width"]] <- 2

    for (nSam in 1:length(otherPlot[["x"]][["data"]])) {
        lgndGrp <- otherPlot[["x"]][["data"]][[nSam]][["legendgroup"]]
        if (!is.null(lgndGrp)) {
            nlgndGrp <- stringr::str_remove_all(lgndGrp, "[()]")
            lgndSplt <- stringr::str_split_fixed(nlgndGrp, ",", 2)
            isShape <- lgndSplt[, 2] == 1
            isSingle <- F
            if (lgndSplt[, 1] == "" | lgndSplt[, 2] == "") {
                nlgndGrp <- paste0(nlgndGrp, ",", nlgndGrp)
                isSingle <- T
            }
            otherPlot[["x"]][["data"]][[nSam]][["marker"]][["opacity"]] <- 1
            otherPlot[["x"]][["data"]][[nSam]][["marker"]][["size"]] <- 18
            otherPlot[["x"]][["data"]][[nSam]][["marker"]][["line"]][["width"]] <- 2
            otherPlot[["x"]][["data"]][[nSam]][["marker"]][["line"]][["opacity"]] <- 1
            nGroupSplit <- stringr::str_split_fixed(nlgndGrp, ",", 2)[1, 1]
            if (isShape) {
                otherPlot <- FormatShapeColors(isSingle, nlgndGrp, lgndSplt, fig, otherPlot, nSam)
            } else{
                theSamLabs <- which(nlgndGrp == figGrps)
                otherPlot[["x"]][["data"]][[nSam]]$text <- paste("Sample:", fig$data$samples[theSamLabs])
                curr_shape_type <- otherPlot[["x"]][["data"]][[nSam]][["marker"]][["symbol"]]
                otherPlot <- FormatShapeData(isSingle, otherPlot, nSam)
            }
            otherPlot[["x"]][["data"]][[nSam]][["name"]] <-
                otherPlot[["x"]][["data"]][[nSam]][["legendgroup"]] <- nGroupSplit
        }
    }

    otherPlot <- FormatMarkerData(otherPlot)
    return(otherPlot)
}


FormatPlotlyLayout <- function(otherPlot, gb) {
    otherPlot$x[["layout"]][["annotations"]] <- NULL
    opLayout <- otherPlot[["x"]][["layout"]]
    opLayout[["font"]][["size"]] <- 12
    opLayout[["title"]][["font"]][["size"]] <- 20
    opLayout[["xaxis"]][["title"]][["font"]][["size"]] <- 14
    opLayout[["yaxis"]][["title"]][["font"]][["size"]] <- 14
    opLayout[["xaxis"]][["tickfont"]][["size"]] <- 12
    opLayout[["yaxis"]][["tickfont"]][["size"]] <- 12
    otherPlot[["x"]][["layout"]] <- opLayout
    pltTtl <- otherPlot[["x"]][["layout"]][["title"]][["text"]]
    otherPlot[["x"]][["layout"]][["margin"]]$t <- 100
    numList <- 1:length(otherPlot[["x"]][["data"]])
    shapeText <- paste("Symbols:", R.utils::capitalize(gb$samShapes))
    #colorText <- paste("Colors:", R.utils::capitalize(gb$samGroup))
    for(nSam in numList) {
        if (!is.null(otherPlot[["x"]][["data"]][[nSam]][["visible"]])) {
            otherPlot[["x"]][["data"]][[nSam]]$visible <- 'legendonly'
            otherPlot[["x"]][["data"]][[nSam]]$showlegend <- TRUE
            otherPlot[["x"]][["data"]][[nSam]][["legendgroup"]] <- paste("Shape Legend", otherPlot[["x"]][["data"]][[nSam]][["legendgroup"]])
            otherPlot[["x"]][["data"]][[nSam]][["marker"]][["line"]][["color"]] <- "black"
            otherPlot[["x"]][["data"]][[nSam]][["marker"]][["color"]] <- "darkgrey"
            otherPlot[["x"]][["data"]][[nSam]][["legendgrouptitle"]] <- list(text = shapeText, font = list(color="black", size=14))
        }
        #else{
            #otherPlot[["x"]][["data"]][[nSam]][["legendgroup"]] <- "Color Legend"
            #otherPlot[["x"]][["data"]][[nSam]][["legendgrouptitle"]] <- list(text = colorText, font = list(color="black", size=14))
        #}
    }
    otherPlot[["x"]][["layout"]][["legend"]][["title"]][["text"]] <- "Shape & Color Legend"
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
  isNoShapes <- all(markerSyms == "circle")
  if (isNoShapes) {
      for (grpIdx in 1:length(otherPlot$x$data)) {
        currGrp <- otherPlot$x$data[[grpIdx]]$legendgroup
        matchedRows <- which(figDat$GROUPS == currGrp)
        currData <- figDat[matchedRows,]
        otherPlot$x$data[[grpIdx]]$x <- currData$x
        otherPlot$x$data[[grpIdx]]$y <- currData$y
        newText <- paste(
            paste("Sample:", currData$samples),
            paste("x:", round(currData$x, 2)),
            paste("y:", round(currData$y, 2)),
            paste("Group:", currGrp),
           paste("Color:", currData$col),
            sep = "<br />"
        )
        otherPlot$x$data[[grpIdx]]$text <- newText
      }
    return(otherPlot)
  }
  
  for (grpIdx in 1:length(otherPlot$x$data)) {
    grpNam <- otherPlot$x$data[[grpIdx]]$legendgroup
    if (!is.null(grpNam)) {
      samNames <- stringr::str_split(otherPlot[["x"]][["data"]][[grpIdx]]$text,
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
    list(text = "<b>Sample Groups Legend</b><br>", font = list(size = 14))
  lgnd <- list(title = ttl, font = list(size = 14))
  otherPlot <- otherPlot %>% plotly::layout(margin = mrg, legend = lgnd)
  return(otherPlot)
}


makePlotly <- function(fig, gb) {
  otherPlot <- NULL
  otherPlot <- FormatLegendText(fig)
  otherPlot <- FormatPlotlyLayout(otherPlot, gb)
  otherPlot <- FormatHoverInfo(otherPlot)
  uniGrp <- unlist(lapply(
    X = 1:length(otherPlot[["x"]][["data"]]),
    FUN = function(X) {return(otherPlot[["x"]][["data"]][[X]]$name)}
    ))
  markerSyms <- GetPlotlySymbols(fig, uniGrp)
  otherPlot <- FormatPlotLabels(fig, otherPlot, uniGrp, markerSyms)
  otherPlot <- FormatPlotlyLegend(otherPlot)
  return(otherPlot)
}


GetFlatPlots <- function(fig){
  options(
    repr.plot.width = 19,
    repr.plot.height = 12,
    repr.plot.res = 350
  )
    fig <- fig +
        ggplot2::theme(
            legend.direction = "vertical",
            legend.margin = ggplot2::margin(t = -25),
            legend.box.margin = ggplot2::margin(0, 0, 0, 0),
            legend.justification = "right",
            legend.position = "right"
        )
    leg <- supM(cowplot::get_legend(fig))
    
    if (length(leg) > 10) {
        fig <- fig + theme(legend.position = "none")
        supM(print(fig))
        grid::grid.newpage()
        grid::grid.draw(leg)
    } else{
        fig <- fig +
            ggplot2::theme(
                legend.direction = "vertical",
                legend.margin = ggplot2::margin(t = 0),
                legend.position = "right"
            )
            
        #supM(print(fig))
    }
}


GetFlatPlotOnly <- function(fig){
  options(repr.plot.width = 19, repr.plot.height = 12, repr.plot.res = 350)
    fig <- fig +
        theme(
            legend.direction = "vertical",
            legend.margin = ggplot2::margin(t = -25),
            legend.box.margin = ggplot2::margin(0, 0, 0, 0),
            legend.justification = "right",
            legend.position = "bottom"
        )
    leg <- supM(cowplot::get_legend(fig))
    if (length(leg) > 10) {
        fig <- fig + theme(legend.position = "none")
        supM(print(fig))
        grid::grid.newpage()
        grid::grid.draw(leg)
    }else{
        fig <- fig + theme(legend.direction = "vertical",
                           legend.margin = ggplot2::margin(t = 0))
        options(repr.plot.width = 19, repr.plot.height = 12, repr.plot.res = 350)
    }
    return(fig)
}


# Function to clean and format labels
cleanLabel <- function(label) {
  label <- stringr::str_replace_all(label, c("# " = "", "#" = "", "\n" = "", " " = "_"))
  label <- stringr::str_replace_all(label, c("\\(" = "", "\\)" = "", "/" = "_"))
  label <- stringr::str_replace_all(label, 'Interactive', "")
  gc(verbose = F)
  return(label)
}


renamePlotlyPng <- function(label, path) {
  plotlyFile <- file.path(path, paste0(label, ".png"))
  currPngs <- dir(path, pattern = ".png", full.names = T)
  info <- file.info(currPngs)
  newest_png <- which.max(info$mtime)
  currPng <- paste(currPngs[newest_png])
  cat(paste("\ncurrPng:", currPng))
  original_png <- file.path(path, currPng)
  cat(paste("from =",original_png,"\n", "to =", plotlyFile,"\n", "dirname", dirname(plotlyFile)))
  file.rename(from = original_png, to = plotlyFile)
}

# Function to save ggplot plot
saveGgplotImage <- function(plot, label, path) {
  ggplotFile <- file.path(path, paste0(label, ".png"))
  ggplot2::ggsave(filename = ggplotFile, plot = plot, width = 19, height = 12, dpi = 350)
}


GrabPngPath <- function(op) {
  currLab <- op[["x"]][["layout"]][["title"]][["text"]]
  supervise <- ifelse(stringr::str_detect(currLab, pattern = "Unsuper"), 
                      "unsupervised", "supervised")
  pngFiPath <- file.path(gb$runDir, "figures", "tsne", supervise)
  return(pngFiPath)
}


SaveGgplotPng <- function(newLab, pngFiPath, fig) {
    newLab2 <- stringr::str_replace_all(newLab, 'Interactive', "")
    pngFi2 <- file.path(pngFiPath, paste0(newLab2, ".png"))
    fig2 <- GetFlatPlotOnly(fig)
    fig2 <- fig2 + ggplot2::theme(legend.position = "right", legend.box = "vertical")
    ggplot2::ggsave(filename = pngFi2, plot = fig2, units = "in", width = 10, height = 6, dpi = 350)
    
    fi_path <- stringr::str_split_fixed(pngFi2, "/figures/", 2)[1,2]
    fi_path <- file.path(".", "figures", fi_path)
    altTxt <- paste("![File", newLab2,"Plot](")
    cat(paste0(altTxt, fi_path, "){height=600px}"))
    cat("\n\n")
}
        
        
selectPlots <- function(doPlotly = F, tplots, ty, tps, outDirs) {
  for (zz in 1:nrow(outDirs)) {
    invisible(gc())
    fig <- tplots[[zz]]
    tabStart <- paste('### Top', tps[ty, 1], outDirs[zz, 3])
    
    cat(paste(tabStart, '(Interactive)','\n\n'))
    newLab <- cleanLabel(tabStart)
    op <- gb$makePlotly(fig, gb)
    supM(print(htmltools::tagList(plotly::ggplotly(op))))
    cat('\n\n')
    
    pngFiPath <- GrabPngPath(op)
    #renamePlotlyPng(newLab, pngFiPath)
    cat(paste(tabStart, '\n\n'))
    gc(verbose = F)
    options(repr.plot.width = 19, repr.plot.height = 12, repr.plot.res = 350)
    GetFlatPlots(fig)
    cat('\n\n')
    SaveGgplotPng(newLab, pngFiPath, fig)
  }
}


gb$grabAllBeta <- grabAllBeta <- function(targets1, betas, supervised = F) {
    if (supervised == T) {
        supBets <- betas[, targets1$SampleFilter] # filtering betas
        allBetas1 <- list(supBets[1:100,], supBets[1:1000,], supBets)
    } else{
        betas1 <- betas[, targets1$SampleFilter] # filtering betas
        if (file.exists(gb$unbetaVariance)) {
            unBets <- gb$LoadRdatObj(gb$unbetaVariance, msgProg=F)
        } else{
            unBets <- gb$takeTopVariance(betas1, topVar = 1:10000)
            gb$SaveObj(unBets, gb$unbetaVariance)
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
                        tsne_titles, doPlotly = F, supervised = F, names2Label = NULL){
    targets <- FormatPlotTargets(targets, samGroup, samShapes, samNames, targFilter)
    tps <- as.data.frame(unique(targets$SamGroups))
    ty = NULL
    samShapes <<- samShapes
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


SubsetTargets <- function(targets, varToFilter = NULL){
  if(is.null(varToFilter)){return(targets)}
  targets1 <- targets[varToFilter,] # exclude of targets are subsetted
  rownames(targets1) <- 1:nrow(targets1) # ensures rows are chronological
  return(targets1)
}


assign("subsetBetas", subsetBetas)

GeneratePlotList <- function(targets1, allBetas1, gb) {
  lapply(seq_along(allBetas1), function(X) {
    tsne_plot <-
      gb$getTsneVal(
        suppressMessages(gb$generateTvals(allBetas1[[X]])),
        saNames = targets1$SampleFilter,
        samGrp = targets1$PointColors,
        colorGrp = targets1$color,
        symGrp = targets1$Sym_Shape
      )
    return(tsne_plot)
  })
}

ProcessAndSavePlots <- function(plotList, custom, outDirs, gb, tsne_titles) {
    for (plotIndex in seq_along(plotList)) {
        invisible(gc())
      fig <- gb$genTsnePlot(
        tsne_plot = plotList[[plotIndex]],
        titleLabel = tsne_titles[plotIndex],
        symbolsLabel = plotList[[plotIndex]]$symbol,
        colorLabel = plotList[[plotIndex]]$GROUPS,
        names2Label = gb$names2Label
      )
        tabStart <- paste('### Top', custom, outDirs[plotIndex, 3])
        cat(paste(tabStart, '(Interactive)', '\n\n'))
        newLab <- gb$cleanLabel(tabStart)
        op <- gb$makePlotly(fig, gb)
        print(htmltools::tagList(plotly::ggplotly(op)))
        cat('\n\n')
        pngFiPath <- gb$GrabPngPath(op)
        cat(paste(tabStart, '\n\n'))
        gc(verbose = F)
        options(repr.plot.width = 19, repr.plot.height = 12, repr.plot.res = 350)
        gb$GetFlatPlots(fig)
        cat('\n\n')
        gb$SaveGgplotPng(newLab, pngFiPath, fig)
    }
}

run_tsne_plots <- function(
    targets,
    betas,
    gb,
    colorVariable = "Type",
    shapeVariable = "Shape",
    sampleGrouping = "Sample_Group",
    isSuper = FALSE) {
  gb$CatShapeColor(colorVariable, shapeVariable, 
                   preFix = ifelse(isSuper, "Supervised", "Unsupervised"))
    tsne_titles <- if (isSuper) {gb$tsne_titles[4:6]} else {gb$tsne_titles[1:3]}
    for (custom in unique(targets[[sampleGrouping]])) {
        targets1 <- targets[targets[[sampleGrouping]] == custom, ]
        targets1 <- gb$FormatPlotTargets(targets1, colorVariable, shapeVariable,
                                         gb$col_samNames, sampleGrouping)
        allBetas1 <- gb$grabAllBeta(targets1, betas, isSuper)
        outDirs <- gb$grabPngNames(tsne_titles)
        plotList <- GeneratePlotList(targets1, allBetas1, gb)
        ProcessAndSavePlots(plotList, custom, outDirs, gb, tsne_titles)
    }
}

