gb <- globalenv(); assign("gb", gb)
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
    "cleanSamples.R",
    "Differential.R"
)
scripts <- paste0(mainHub, rFiles)
pkgs <- c("MEAL","SummarizedExperiment", "S4Vectors", "minfi", "limma", "IRanges", "impute", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "GenomicRanges", "conumee", "BiocGenerics", "Biobase", "cowplot")

supSrt(lapply(scripts, function(i){message("Sourcing: ", i);devtools::source_url(i)}))
supSrt(librarian::shelf(pkgs, ask = F, update_all = F, quiet = FALSE))
#if(!require("maxprobes")){ devtools::install_github("markgene/maxprobes", dependencies = T)}
#if(!require("cnv.methyl")){devtools::install_github("https://github.com/ijcBIT/cnv.methyl.git", dependencies = T)}

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

#knitr::opts_knit$set(knitOpt); knitr::opts_chunk$set(knitOpt)
animation::ani.options(autobrowse = FALSE)
options(width=1800)

library("mnp.v11b6")
library("magrittr")
library("dplyr")

setKnitDir <- function(runDir) {
    syscmd <- paste("cd", runDir)
    system(syscmd); setwd(runDir)
    knitr::opts_knit$set(root.dir = runDir)
}

set.seed(1234)