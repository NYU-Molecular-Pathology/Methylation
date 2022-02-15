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
  "Cairo",
  "gridExtra"
)
librarian::shelf(pkgs, ask=F, warn.conflicts=F)

require("gridExtra")
require("Cairo")
require("minfi")
require("ComplexHeatmap")
if(!require("GenVisR")){BiocManager::install("GenVisR")}

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
