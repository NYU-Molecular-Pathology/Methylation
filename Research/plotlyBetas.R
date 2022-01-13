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

makePlotly <- function(fig) {
    otherPlot <-
        supM(plotly::ggplotly(
            fig,
            dynamicTicks = T,
            height = 800,
            width = 1200
        ))
    otherPlot$x[["layout"]][["annotations"]] <- NULL
    opLayout <- otherPlot[["x"]][["layout"]]
    opLayout[["font"]][["size"]] <- 12
    opLayout[["title"]][["font"]][["size"]] <- 24
    opLayout[["xaxis"]][["title"]][["font"]][["size"]] <- 16
    opLayout[["yaxis"]][["title"]][["font"]][["size"]] <- 16
    opLayout[["xaxis"]][["tickfont"]][["size"]] <- 16
    opLayout[["yaxis"]][["tickfont"]][["size"]] <- 16
    otherPlot[["x"]][["layout"]] <- opLayout
    hoverinfo <- paste0("Sample: ", fig[["data"]][["samples"]], "</br></br>")
    opInfo <- otherPlot[["x"]][["data"]]
    
    for (sam in 1:length(opInfo)) {
        opInfo[[sam]][["hoverinfo"]] <- "none"
        opInfo[[sam]][["marker"]][["symbol"]] <- 'circle'
        opInfo[[sam]][["marker"]][["size"]] <- 10
        opInfo[[sam]][["hoverinfo"]] <- c("text") #"text",
    }
    
    otherPlot[["x"]][["data"]] <- opInfo
    uniGrp <- unique(fig[["data"]]$GROUPS)
    
    for (grpT in 1:length(uniGrp)) {
        grpNam <- uniGrp[grpT]
        currGrp <- fig[["data"]]$GROUPS == grpNam
        samLabs <- fig[["data"]]$samples
        
        for (grpLabNam in 1:length(otherPlot[["x"]][["data"]])) {
             currLabNam <- otherPlot[["x"]][["data"]][[grpLabNam]]$name
            currLabNam <- stringr::str_split_fixed(currLabNam,",",2)[,2]
            currLabNam <- stringr::str_remove_all(currLabNam, ")")
            message(currLabNam, " and ", grpNam)
            if (currLabNam == grpNam) {
              hoverinfo <- paste0("Sample: ", samLabs[currGrp], " (", grpNam, ")", "</br></br>")
              otherPlot[["x"]][["data"]][[grpLabNam]][["text"]] <- hoverinfo
            }
        }
    }
    otherPlot <-
        otherPlot %>% plotly::layout(legend = list(
            title = list(text = "<b>Legend</b><br>", font = list(size = 24)),
            font = list(size = 14)
        ))
    
    return(otherPlot)
}

selectPlots <- function(doPlotly=F,tplots,ty,tps,outDirs){
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
        op <- makePlotly(fig)
        
        supM(print(htmltools::tagList(ggplotly(op))))
        cat('\n\n')
      }
    }
}

grabAllBeta <- function(targets1, betas) {
  betas1 <- betas[ ,targets1$SampleFilter] # filtering betas
  unBets <- gb$takeTopVariance(betas1, topVar = 1:10000)
  allBetas1 <-list(unBets[1:100,], unBets[1:1000,], unBets)
  return(allBetas1)
}

subsetBetas <- function(targFilter,samGroup, betas, targets, samNames,tsne_titles, doPlotly=F) {
  tps <- unique(targFilter)
  targets$SamGroups <- targets$SampleFilter <- NULL
  targets$SampleFilter <- samNames # creating new column
  targets$SamGroups <- targets[, colnames(targFilter)]
  ty = NULL
  for (ty in 1:nrow(tps)) {
    custom = tps[ty, 1]
    message("Current: ", custom)
    targets1 <- targets[targFilter == custom,]
    allBetas1 <- grabAllBeta(targets1,betas)
    outDirs <- grabPngNames(tsne_titles)[1:3, ]
    plotList = NULL
    plotList <-  gb$doMultiple(allBetas1,tsne_titles, outDirs, targets1, tps,ty, custom)
    tplots <- NULL
    
    tplots <- plotSaver2(outDirs, tsne_titles, tps, ty, plotList, custom)
    selectPlots(doPlotly,tplots,ty,tps,outDirs)
    }
}

tierBetas <- function(betas, batchCorrect = F, getSuper = F) {
    selectSams <- RGSet@colData@listData[["SentrixID_Position"]][RGSet@colData@listData[["Sample_ID"]] %in%  colnames(betas)]
    if (batchCorrect == T) {
        unBetas <- gb$batchCorrectBs(betas, RGSet[,RGSet@colData@rownames==selectSams], topVar = 1:10000)
        superbetas <- gb$batchCorrectBs(betas, RGSet[,RGSet@colData@rownames==selectSams], topVar = 1:10000, T)
    } else{ unBetas <- gb$takeTopVariance(betas, topVar = 1:10000)}
    if (getSuper == T) {
        superbetas <- gb$getSupervise(betas, RGSet[,RGSet@colData@rownames==selectSams], topVar = 1:10000)
        return(superbetas)
    } else{return(unBetas)}
}

assign("subsetBetas", subsetBetas)
