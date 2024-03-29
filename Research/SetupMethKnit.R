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
if(!require("knitrProgressBar")){devtools::install_github("rmflight/knitrProgressBar")}
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
if(!require("maxprobes")){ devtools::install_github("markgene/maxprobes", dependencies = T, upgrade="never")}
if(!require("cnv.methyl")){devtools::install_github("https://github.com/ijcBIT/cnv.methyl.git", dependencies = T, upgrade="never")}

SetKnitOpts <- function(){
    knitOpt <- list(
    echo = FALSE,
    message = FALSE,
    warning = FALSE,
    cache = FALSE,
    error = TRUE,
    rows.print = 15,
    fig.keep = "all",
    fig.show = "asis",
    out.width = '100%',
    class.source = "bg-success",
    results = 'asis',
    fig.path = "figures/",
    opts.label = TRUE,
    progress = TRUE, verbose = TRUE
)
    return(knitr::opts_knit$set(knitOpt))
}

SetChunkOpts <- function(){
    knitOpt <- list(
    echo = FALSE,
    message = FALSE,
    warning = FALSE,
    cache = FALSE,
    error = TRUE,
    rows.print = 15,
    fig.keep = "all",
    fig.show = "asis",
    out.width = '100%',
    class.source = "bg-success",
    results = 'asis',
    fig.path = "figures/",
    opts.label = TRUE,
    progress = TRUE, verbose = TRUE
)
    return(knitr::opts_chunk$set(knitOpt))
}

try(animation::ani.options(autobrowse = FALSE), silent=T)
#options(width = 3500)
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

GrabKnitTemplates <- function(){
    tempOpts <- list(
mdsPlot = list(
            echo = FALSE,
            fig.height = 6,
            fig.width = 12,
            out.width = '100%',
            dpi = 350,
            fig.keep = 'all',
            include = TRUE
        )
        ,
heatmaps = list(
    echo = FALSE,
    fig.width = 18, fig.height = 22,
    dpi = 350,
    fig.asp = 0.8,
    out.width='100%',
    results = 'asis',
    fig.keep = 'all',
    fig.path = "figures/heatmaps/"
)
,
tnsePlot = list(
    echo = FALSE,
    error = FALSE,
    dpi = 350,
    results = 'asis',
    fig.height=7,
    fig.width=9,
    fig.keep = 'all',
    out.height="100%", out.width="100%",
    fig.path = "figures/tsne/"
)
,
samcnv = list(
    echo = FALSE,
    fig.keep = 'all',
    fig.height = 8,
    fig.width = 15,
    fig.asp = 0.75,
    results = 'asis',
    fig.path = "figures/cnv/",
    include = TRUE
)
,
cnvplots = list(
    echo = FALSE,
    fig.keep = 'all',
    fig.show = 'asis',
    results = 'asis',
    fig.height = 8, #10,
    fig.width = 15, #16,
    out.width = '100%',
    dpi = 350,
    fig.path = 'figures/cnv/',
    include = TRUE
)
,
clustProf = list(echo=FALSE, dpi=350, fig.height=7, fig.width=9, out.width='800px', fig.path = "figures/cluster/", include=TRUE)
)
#return(invisible(lapply(X = names(tempOpts),FUN=function(X){knitr::opts_template$set(tempOpts[[X]])})))
return(knitr::opts_template$set(
    cnvplots = tempOpts[["cnvplots"]],
    mdsPlot = tempOpts[["mdsPlot"]],
    heatmaps = tempOpts[["heatmaps"]],
    tnsePlot = tempOpts[["tnsePlot"]],
    samcnv = tempOpts[["samcnv"]],
    clustProf = tempOpts[["clustProf"]]
))
}


htmlClose <- function(){return(cat("\n\n:::\n\n"))}

# library(reticulate)

# # Attempt to import the kaleido module
# try_import <- function() {
#   tryCatch({
#     import("kaleido")
#     TRUE  # Return TRUE if import succeeded
#   }, error = function(e) {
#     FALSE  # Return FALSE if import failed
#   })
# }
# reticulate::use_python("/usr/local/bin/python", required = TRUE)
# # Check if kaleido is installed; if not, install it
# if (!try_import()) {
#   reticulate::py_install("kaleido", method = "virtualenv")
# }
