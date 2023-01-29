gb <- globalenv(); assign("gb", gb)
# TODO: add args to load Differential + parrallel envir

# https://bookdown.org/yihui/rmarkdown-cookbook/custom-knit.html
supM <- function(sobj){return(suppressMessages(suppressWarnings(sobj)))}
supPk <- function(sobj){return(suppressPackageStartupMessages(sobj))}
supSrt <- function(pk){return(suppressPackageStartupMessages(suppressWarnings(pk)))}

mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/Research/"
rFiles = c(
    #"all_installer.R",
    "getRGsetBetas.R",
    "TsnePlotter.R",
    "getIdatFiles.R",
    "hmFunctions.R",
    "plotlyBetas.R",
    "cnvFunctions.R",
    "tableMaker.R",
    "pathviews.R",
    "pullRedcap_manual.R",
    "cleanSamples.R"
#    "Differential.R"
)

pkgs <-
  c(
    "MEAL",
    "SummarizedExperiment",
    "S4Vectors",
    "minfi",
    "maxprobes",
    "limma",
    "IRanges",
    "impute",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "GenomicRanges",
    "conumee",
    "BiocGenerics",
    "Biobase",
    "cowplot"
  )

scripts = paste0(mainHub, rFiles)
suppressWarnings(lapply(scripts, function(i){message("Sourcing: ", i);devtools::source_url(i)}))
supSrt(librarian::shelf(pkgs, ask = F, update_all = F, quiet = FALSE))
extraDeps <- c("minfiData")
supSrt(librarian::shelf(extraDeps, ask = F, update_all = F, quiet = FALSE))
if(!require("maxprobes")){ devtools::install_github("markgene/maxprobes", dependencies = T)}
if(!require("cnv.methyl")){devtools::install_github("https://github.com/ijcBIT/cnv.methyl.git", dependencies = T)}

knitOpt <- list(
    echo = FALSE,
    message = FALSE,
    warning = FALSE,
    cache = FALSE,
    error = TRUE,
    rows.print = 15,
    fig.keep = "all",
    fig.show = "asis",
    class.source = "bg-success",
    results = 'asis',
    fig.path = "figures/"
)

animation::ani.options(autobrowse = FALSE)
options(width = 1500)
library("mnp.v11b6")
library("magrittr")
library("dplyr")

setKnitDir <- function(runDir) {
  syscmd <- paste("cd", runDir)
  system(syscmd); setwd(runDir)
  knitr::opts_knit$set(root.dir = runDir)
}

CheckDirCreate <- function(pathLocation){
    dataOutDir <- file.path(getwd(), pathLocation)
    if (!dir.exists(dataOutDir)) {dir.create(dataOutDir)}
}

set.seed(1234)
Sys.setenv('R_MAX_VSIZE'=32000000000)
