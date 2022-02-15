gb <- globalenv(); assign("gb", gb)
setDirectory<-function(foldr) {
    bsDir = paste("cd", foldr); mm2 = crayon::white$bgRed("Location Not Found:", foldr)
    if (dir.exists(foldr)) {system(bsDir);setwd(foldr);assign("runDir", foldr)} else{warning(mm2)}
}

#knitr::opts_chunk$set(fig.width=15, fig.height=10)
library("foreach")
require("foreach")
library("dplyr")
library("kableExtra")
if(!require("GenVisR")){BiocManager::install("GenVisR")}
library("minfi")
require("minfi")
library("conumee")
library("mnp.v11b6")
library("enrichplot")
library("limma")
library("stringr")
library(tidyverse)
require(tidyr)
library("DOSE")
library(ggnewscale)
library("enrichplot")
require("enrichplot")
if(!require("assertr")){install.packages("assertr")}
suppressPackageStartupMessages(library(ComplexHeatmap))
require("ComplexHeatmap")
library("Biobase")
if(!require("magick")){install.packages("magick", quiet =T)}
if(!require("irlba")){install.packages("irlba", quiet =T)}
library(magick)
library(grid)
require(grid)
library(gridExtra)
require(gridExtra)
library(Cairo)
require(Cairo)

makeDt <- function(targets) {
  dtOpts <- list(scrollX = T, scrollY = T, info = F, pageLength=15,autoWidth = T, rownames = F,lengthChange = T,searchable = T)
  theDt <- DT::datatable(
    targets, selection = "single", autoHideNavigation = T, options = dtOpts
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
