#!/usr/bin/env Rscript
## Script name: SetupMethKnit.R
## Purpose: source of global scripts and parameters to knit Methyl research projects
## Date Created: Jan 16, 2023
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024

gb <- globalenv()
assign("gb", gb)

formals(library)$quietly <- T
formals(library)$warn.conflicts <- F
formals(require)$warn.conflicts <- F
formals(requireNamespace)$quietly <- T
formals(install.packages)$dependencies <- T
formals(install.packages)$ask <- F

# TODO: add args to load Differential + parrallel envir

if(!require("devtools", warn.conflicts = F)) {
    install.packages("devtools", dependencies = T, ask = F)
}

if (!requireNamespace("knitr", quietly = T)) {
    install.packages("knitr", ask = F, dependencies = T, type = "binary")
}

library("knitr")

# https://bookdown.org/yihui/rmarkdown-cookbook/custom-knit.html
supM <- function(sobj){return(suppressMessages(suppressWarnings(sobj)))}
supPk <- function(sobj){return(suppressPackageStartupMessages(sobj))}
supSrt <- function(pk){return(suppressPackageStartupMessages(suppressWarnings(pk)))}

download_pkg_unzip <- function(git_repo, zip_name = "main.zip") {
    repo_url <- file.path("https://github.com", git_repo, "archive/refs/heads")
    url_path <- file.path(repo_url, zip_name)
    local_dir <- file.path(fs::path_home(), "github_pkgs")
    out_file <- file.path(local_dir, paste0(basename(git_repo), "-", zip_name))
    utils::download.file(url = url_path, method = "auto", destfile = out_file)
    utils::unzip(out_file, exdir = local_dir)
}


# FUN: Downloads Github repo locally then installs ----------------------------
local_github_pkg_install <- function(git_repo) {
    repo_url <- file.path("https://github.com", git_repo, "archive/refs/heads")
    local_dir <- file.path(fs::path_home(), "github_pkgs")
    
    if (!dir.exists(local_dir)) dir.create(local_dir)
    tryCatch(
        expr = {
            download_pkg_unzip(git_repo, zip_name = "main.zip")
        },
        warning = function(e){
            download_pkg_unzip(git_repo, zip_name = "master.zip")
        },
        error = function(e){
            download_pkg_unzip(git_repo, zip_name = "master.zip")
        }
    )
    
    unzipped_dir <- list.dirs(local_dir, full.names = TRUE, recursive = FALSE)
    pkg_dir <- unzipped_dir[grepl(basename(git_repo), unzipped_dir)]
    
    install.packages(pkg_dir, repo = NULL, type = "source", dependencies = T)
}


# FUN: Installs package from a Github repository ------------------------------
install_repo <- function(git_repo){
    tryCatch(
        expr = {
            devtools::install_github(
                git_repo, dependencies = T, upgrade = "always", type = "source",
                auth_token = NULL, subdir = basename(git_repo))
        },
        error = function(e) {
            devtools::install_github(
                git_repo, dependencies = T, upgrade = "always", type = "source",
                auth_token = NULL)
        }
    )
}


try_github_inst <- function(git_repo){
    message("Installing from repo: ", git_repo)
    tryCatch(
        expr = {
            install_repo(git_repo)
        },
        error = function(e) {
            local_github_pkg_install(git_repo)
        }
    )
}


check_git_install <- function(git_repo) {
    pkg <- basename(git_repo)
    if (!requireNamespace(pkg, quietly = T)) {
        try_github_inst(git_repo)
    }
    library(pkg, character.only = T, logical.return = T)
}

load_pkg <- function(pkg){
    if (!requireNamespace(pkg, quietly = T)) {
        tryCatch(
            install.packages(pkg, dependencies = T, ask = F),
            error = BiocManager::install(pkg, update = F, ask = F)
        )
    }
    pkg_success <- library(pkg, character.only = T, logical.return = T)
    return(paste(pkg, "loaded:", pkg_success))
}

cnv_deps <- c("textshaping", "fstcore", "fst", "fstcore")

lapply(cnv_deps, load_pkg)
#try(devtools::install_github("https://github.com/ijcBIT/cnv.methyl.git"), T)
#check_git_install("ijcBIT/cnv.methyl")

if (Sys.info()[['sysname']] == "Linux") {
    load_pkg("rprofile")
    rprofile::set_startup_options(
        show.signif.stars = FALSE,
        useFancyQuotes = FALSE,
        Ncpus = parallel::detectCores() - 2
    )
    Sys.setenv("R_PROFILE" = file.path(Sys.getenv("HOME"), "Rprofile.site"))
    Sys.setenv(IMAGEMAGICK_V6_HOME = "/gpfs/share/apps/imagemagick/6.9.10/bin/convert")
}

mainHub <-
    "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main"
research_url <- file.path(mainHub, "Research")
mainLnk <- file.path(mainHub, "R/LoadInstallPackages.R")

rFiles = c(
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
)

pkgs <- c(
    "MEAL",
    "SummarizedExperiment",
    "S4Vectors",
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
    "cowplot",
    "fst",
    "itertools",
    "Cairo",
    "tinytex"
)

scripts <- file.path(research_url, rFiles)
suppressWarnings(lapply(scripts, devtools::source_url))


if (!requireNamespace("minfi", T)) {
    devtools::install_github(
        "mwsill/minfi",
        upgrade = "never",
        force = T,
        dependencies = T
    )
    devtools::install_github(
        "mwsill/IlluminaHumanMethylationEPICv2manifest",
        upgrade = "always",
        force = T,
        dependencies = T
    )
}

minfiVers <- as.character(utils::packageVersion("minfi"))
if (minfiVers != "1.43.1") {
    devtools::install_github("mwsill/minfi", upgrade = "never", force = T, dependencies = T)
    devtools::install_github("mwsill/IlluminaHumanMethylationEPICv2manifest", upgrade = "always", force = T, dependencies = T)
}

supSrt(librarian::shelf(pkgs, ask = F, update_all = F, quiet = FALSE))

check_git_install("rmflight/knitrProgressBar")
check_git_install("Ryo-N7/tvthemes")
if (!requireNamespace("mdthemes", quietly = T)) {
    remotes::install_github("thomas-neitmann/mdthemes", upgrade = "never")
    #check_git_install("thomas-neitmann/mdthemes")
}

load_pkg("minfiData")
check_git_install("markgene/maxprobes")

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
library("minfi")

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

gb$knitOpt <- SetKnitOpts()

htmlClose <- function(){return(cat("\n\n:::\n\n"))}

if(Sys.info()[['sysname']] == "Linux") {
    library("reticulate")
    # Attempt to import the kaleido module
    try_import <- function() {
        tryCatch({
            import("kaleido")
            TRUE
        }, error = function(e) {
            FALSE
        })
    }
    reticulate::use_python("/usr/local/bin/python", required = TRUE)
    # Check if kaleido is installed; if not, install it
    if (!try_import()) {
        reticulate::py_install("kaleido", method = "virtualenv")
    }
}

#try(tinytex::tlmgr_update(self = TRUE, all = TRUE), T)
gb$SetKnitOpts()
gb$sourceParams()
gb$assignOpts()
