gb <- globalenv(); assign("gb", gb)
formals(library)$quietly <- T
formals(library)$warn.conflicts <- F
formals(require)$warn.conflicts <- F
formals(install.packages)$dependencies <- T
formals(install.packages)$verbose <- T
formals(install.packages)$ask <- F 

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
mainLnk = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/LoadInstallPackages.R"
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
scripts = c(mainLnk, scripts)
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
options(knitr.package.verbose = TRUE)

gb$tempOpts <- list(
heatmaps = list(
    echo = FALSE,
    fig.width = 8,
    fig.height = 11,
    dpi = 300,
    fig.width='80%',
    fig.keep = 'all',
    results = 'asis',
    fig.path = "figures/heatmaps/"
),

tnsePlot = list(
    echo = FALSE,
    error = FALSE,
    dpi = 300,
    out.height = "650px",
    out.width = "3850px",
    results = 'asis',
    fig.path = "figures/tsne/",
    fig.keep = 'all'
)
,
samcnv = list(
    echo = FALSE,
    fig.height = 8,
    fig.width = 15,
    fig.asp = 0.75,
    fig.keep = 'all',
    fig.path = "figures/cnv/",
    results = 'asis',
    include = TRUE
)
,
cnvplots = list(
    echo = FALSE,
    fig.keep = 'all',
    fig.path = 'figures/cnv/',
    fig.show = 'asis',
    results = 'asis',
    fig.height = 10,
    fig.width = 16,
    dpi = 300,
    include = TRUE,
    out.width = '100%'
)
,
mdsPlot = list(
    echo = FALSE,
    fig.height = 8,
    fig.width = 12,
    out.width = '100%',
    dpi = 300,
    fig.keep = 'all',
    include = TRUE
)
)
