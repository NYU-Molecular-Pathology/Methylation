gb <- globalenv(); assign("gb", gb)
# TODO: add args to load Differential + parrallel envir
if(!require("devtools", warn.conflicts = F)){install.packages("devtools", dependencies=T)}
# https://bookdown.org/yihui/rmarkdown-cookbook/custom-knit.html
supM <- function(sobj){return(suppressMessages(suppressWarnings(sobj)))}
supPk <- function(sobj){return(suppressPackageStartupMessages(sobj))}
supSrt <- function(pk){return(suppressPackageStartupMessages(suppressWarnings(pk)))}

if(Sys.info()[['sysname']]=="Linux") {
    if(!require("rprofile")){devtools::install_github("csgillespie/rprofile", dependencies = T)}
    rprofile::set_startup_options(show.signif.stars = FALSE, useFancyQuotes = FALSE, Ncpus = parallel::detectCores()-2)
    Sys.setenv("R_PROFILE"=file.path(Sys.getenv("HOME"), "Rprofile.site"))
    Sys.setenv(IMAGEMAGICK_V6_HOME="/gpfs/share/apps/imagemagick/6.9.10/bin/convert")
}

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
scripts = paste("https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/LoadInstallPackages.R", scripts)
suppressWarnings(lapply(scripts, function(i){message("Sourcing: ", i);devtools::source_url(i)}))
supSrt(librarian::shelf(pkgs, ask = F, update_all = F, quiet = FALSE))

#supSrt(librarian::shelf(ask = F, update_all = F, quiet = F, dependencies = T, libs=extraDeps))
#easypackages::packages(extraDeps, prompt = F, Ncpus = 6)

if(!require("minfiData")){BiocManager::install("minfiData", update=F, ask=F, dependencies=T)}
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

try(animation::ani.options(autobrowse = FALSE), silent=T)
options(width = 1500)
library("magrittr")
library("dplyr")
require("minfi")


setKnitDir <- function(runDir) {
    syscmd <- paste("cd", runDir)
    system(syscmd); setwd(runDir)
    knitr::opts_knit$set(root.dir = runDir)
}


CheckDirCreate <- function(pathLocation){
    dataOutDir <- file.path(getwd(), pathLocation)
    if (!dir.exists(dataOutDir)) {dir.create(dataOutDir)}
}


GetRunDirPath <- function(projectName){
    userNam <- paste0(Sys.info()[["user"]])
    if (Sys.info()[['sysname']] == "Darwin") {
        runDir <- file.path("","Users",userNam,"Documents","Rprojects_local", projectName)
    }else{
        runDir <- file.path("/gpfs/home", userNam,"methylation_projects", projectName)
    }
    return(runDir)
}


#Sys.setenv('R_MAX_VSIZE'=32000000000)
set.seed(1234)

