gb <- globalenv(); assign("gb", gb)
setDirectory<-function(foldr) {
    bsDir = paste("cd", foldr); mm2 = crayon::white$bgRed("Location Not Found:", foldr)
    if (dir.exists(foldr)) {system(bsDir);setwd(foldr);assign("runDir", foldr)} else{warning(mm2)}
}

#knitr::opts_chunk$set(fig.width=15, fig.height=10)

if(!require("librarian")){install.packages("librarian", dependencies=T, verbose=T, Ncpus = 6, quiet=T)}

# Load/install Packages
pkgs <- c(
  'animation',
  'knitr',
  'readr',
  'ggplot2',
  'devtools',
  "BiocManager",
  'BiocVersion',
  'GenVisR',
  'parallel',
  'iterators',
  'doParallel',
  "foreach",
  "dplyr",
  "DT",
  "kableExtra",
  "foreach",
  "minfi",
  "magick",
  "magrittr",
  "conumee",
  "mnp.v11b6",
  "enrichplot",
  "limma",
  "stringr",
  "tidyverse",
  "tidyr",
  "DOSE",
  "ggnewscale",
  "enrichplot",
  "ComplexHeatmap",
  "assertr",
  "Biobase",
  "magick",
  "irlba",
  "magick",
  "grid",
  "RColorBrewer",
  "Cairo",
  "gridExtra"
)
librarian::shelf(pkgs, ask=F, warn.conflicts=F)

require("gridExtra")
require("Cairo")
require("minfi")
require("ComplexHeatmap")
if(!require("GenVisR")){BiocManager::install("GenVisR")}

makeDt <- 
  function(targets) {
    dtOpts <-
      list(
        columnDefs = list(list(className = 'dt-center', targets = "_all")),
        scrollX = T,
        scrollY = T,
        info = F,
        pageLength = 20,
        autoWidth = F,
        rownames = F,
        lengthChange = F,
        searchable = T
      )
    theDt <- DT::datatable(
      targets,
      selection = "single",
      autoHideNavigation = F,
      options = dtOpts
    )
    return(theDt)
}

smallTab <- function(dtObj) {
  dtTable <- knitr::kable(dtObj, row.names = F, "html")
  dtTable <- kableExtra::kable_styling(
      dtTable,
      bootstrap_options = c("striped", "condensed"),
      full_width = F,
      position = "left"
  )
    
 dtTable <- kableExtra::column_spec(dtTable, 1:ncol(dtObj), width = "4cm")
  return(dtTable)
}

# Adds colors to csv targets file to maintain same color scheme between plots
colorTargets <- function(targets, varColumns = c("Type","Origin"), col_vect = NULL) {
    if (is.null(col_vect)) {col_vect <- pals::glasbey()}
    col_vect[6] = "#eb7d34" #changing dark forest to orange color
    col_vect[4] = "#ADD8E6"
    col_vect <- unique(col_vect)
    message("Dimnames:\n",paste(dimnames(targets)[[2]], collapse = " | "))
    stopifnot(all(varColumns %in% dimnames(targets)[[2]]))
    dat <- targets[,varColumns] # varColumns
    anno_df <- data.frame(dat)
    vars2Color <- as.list(lapply(dat, unique))
    colorValues <-lapply(vars2Color, function(x) {x = (col_vect)[1:(length(x))]})
    for (x in 1:length(vars2Color)) {
        for (varNum in 1:length(vars2Color[x])) {
          names(colorValues[x][[1]]) = c(vars2Color[x][[1]])}
      }
      targets$color <- NULL
      for (colNam in varColumns) {
        for (samNam in names(colorValues[colNam][[1]])) {
          currColor <- targets$Type == samNam
          targets$color[currColor] <- paste0(colorValues$Type[samNam])
        }
      }
      return(targets)
}

getColors <- function(samTypes) {
  library('RColorBrewer')
  qualCols = brewer.pal.info[brewer.pal.info$category == 'qual', ]
  col_vector = unlist(mapply(brewer.pal, qualCols$maxcolors, rownames(qualCols)))
  col_vector <- unique(col_vector)
  myColors = col_vector[1:length(unique(samTypes))]
  names(myColors) <- unique(samTypes)
  return(myColors)
}

sanitizeSheet <- function(inputFi, samsheet) {
    library("magrittr")
    library("dplyr")
  samSh <- readxl::read_excel(inputFi)
  samSh <- samSh %>% dplyr::mutate_all(stringr::str_replace, ",", "")
  write.csv(samSh, samsheet, quote = F, row.names = F)
  targets <- read.csv(samsheet, strip.white = T)
  return(targets)
}

animation::ani.options(autobrowse = FALSE); options(width=1200)
library("mnp.v11b6")
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)  
require("tidyverse")
require("plotly")
require("ggplot2")
library("DT")

printMissing <- function(theMissing) {
  if (!is.null(theMissing) & nrow(theMissing) > 0) {
    selCols <-
      c(
        "Sample_Name",
        gb$col_sentrix,
        gb$col_samTypes,
        gb$col_samGrp,
        "Sample_ID",
        "Sample_Group"
      )
    return(gb$smallTab(theMissing[, selCols]))
  } else{
    cat("\n\nNone\n\n")
  }
}


loadHtmlTag <- function(){
    require("tidyverse")
    require("plotly")
    require("ggplot2")
    library("DT")
    htmltools::tagList(DT::datatable(cars, options=list(pageLength = 15, width="100%")))
    par(ask=F); devAskNewPage(ask = F)
    return(htmltools::tagList(plotly::ggplotly(ggplot2::ggplot())))
}
