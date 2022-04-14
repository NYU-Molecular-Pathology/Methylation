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

makePlotly<-function(fig) {
    otherPlot <-
        supM(plotly::ggplotly(
            fig,
            dynamicTicks = T,
            height = 800, #800
            width = 1200 #1200
        ))
    otherPlot$x[["layout"]][["annotations"]] <- NULL
    opLayout <- otherPlot[["x"]][["layout"]]
    opLayout[["font"]][["size"]] <- 12
    opLayout[["title"]][["font"]][["size"]] <- 20
    opLayout[["xaxis"]][["title"]][["font"]][["size"]] <- 14
    opLayout[["yaxis"]][["title"]][["font"]][["size"]] <- 14
    opLayout[["xaxis"]][["tickfont"]][["size"]] <- 12
    opLayout[["yaxis"]][["tickfont"]][["size"]] <- 12
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
            currLabNam <- stringr::str_remove_all(currLabNam, "[()]")
            #message(currLabNam, " and ", grpNam)
            if (currLabNam == grpNam) {
              hoverinfo <- paste0("Sample: ", samLabs[currGrp] , " (", grpNam, ")", "</br></br>")
              otherPlot[["x"]][["data"]][[grpLabNam]][["text"]] <- hoverinfo
            }
        }
    }
    otherPlot <-
        otherPlot %>% plotly::layout(legend = list(
            title = list(text = "<b>Legend</b><br>", font = list(size = 14)),
            font = list(size = 12)
        ))

    return(otherPlot)
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

#grabAllBeta <-function(targets1, betas) {
#  betas1 <- betas[ ,targets1$SampleFilter] # filtering betas
#  unBets <- gb$takeTopVariance(betas1, topVar = 1:10000)
#  allBetas1 <-list(unBets[1:100,], unBets[1:1000,], unBets)
#  return(allBetas1)
#}

grabAllBeta <- function(targets1, betas) {
    betas1 <- betas[ ,targets1$SampleFilter] # filtering betas
    if(file.exists(file.path(gb$runDir,gb$unbetaVariance))){
      unBets <- readRDS(unbetaVariance)
    }else{
    unBets <- gb$takeTopVariance(betas1, topVar = 1:10000)
    saveRDS(unBets, file=file.path(gb$runDir,gb$unbetaVariance))
    }
    allBetas1 <-list(unBets[1:100,], unBets[1:1000,], unBets)
    return(allBetas1)
}

doMultiple <- function(allBetas1,tsne_titles, outDirs, targets1, tps,ty,custom){
    plotN = NULL
    tsneList <-lapply(X = 1:length(allBetas1), FUN=function(X){
      return(suppressMessages(gb$generateTvals(allBetas1[[X]])))
      })
    plotList <- list(foreach::foreach(plotN = 1:length(tsneList),.packages="foreach")%do%{
###################### TO CHANGE ########################      
      tsne_plot <- gb$getTsneVal(
        TSNE = tsneList[[plotN]],
        saNames = targets1$SampleFilter,
        samGrp = targets1$Type,
        colorGrp = targets1$color,
        symGrp = targets1[,gb$col_samGrp]
        ) #targets1[,col_samGroup]
###################### TO CHANGE ########################
      return(tsne_plot)
      })[[1]]
    return(plotList)
}

plotSaver <- function(outDirs,tsne_titles,tps,ty,plotList,custom) {
  plotN=NULL
  options("device.ask.default"=F)
 pltList <- foreach::foreach(plotN = 1:length(plotList),.packages="foreach") %do% {
    pL<- plotList[[plotN]]
###################### TO CHANGE ########################
    return(gb$genTsnePlot(
      tsne_plot=plotList[[plotN]],
      titleLabel=tsne_titles[plotN],
      symbolsLabel = pL$symbol,
      colorLabel = pL$GROUPS,
      names2Label = NULL)) # any specific sample names to label on the plot
###################### TO CHANGE ########################
 }
 return(pltList)
}

gb$subsetBetas <- 
  function(targFilter,samGroup, betas, targets, samNames,tsne_titles, doPlotly=F, supervised=F) {
    tps <- unique(targFilter)
    targets$SamGroups <- targets$SampleFilter <- NULL
    targets$SampleFilter <- samNames # creating new column
    targets$SamGroups <- targets[, colnames(targFilter)]
    ty = NULL
    for (ty in 1:nrow(tps)) {
      custom = tps[ty, 1]
      message("Current Sample Group TSNE: ", custom)
      targets1 <- targets[targFilter == custom, ]
      if (supervised == T) {
        allBetas1 <- list(betas[1:100, ], betas[1:1000, ], betas)
      } else{
        allBetas1 <- gb$grabAllBeta(targets1, betas)
      }
      outDirs <- gb$grabPngNames(tsne_titles)[,]
      plotList = NULL
      plotList <- gb$doMultiple(allBetas1, tsne_titles, outDirs, targets1, tps, ty, custom)
      tplots <- NULL
      tplots <- gb$plotSaver(outDirs, tsne_titles, tps, ty, plotList, custom)
      gb$selectPlots(doPlotly, tplots, ty, tps, outDirs)
    }
}

takeTopVariance <- function(betas, topVar){
    var_probes <- apply(betas, 1.0, var)
    select_var <- names(sort(var_probes[topVar], decreasing = T))
    top_var_beta <- betas[select_var, ]
    return(top_var_beta)
}

gb$tierBetas <-
  function(betas,
           col_sentrix,
           RGSet,
           batchCorrect = F,
           getSuper = F,
           topVar = 1:10000) {
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

assign("subsetBetas", subsetBetas)
