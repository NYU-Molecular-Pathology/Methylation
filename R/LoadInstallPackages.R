#!/usr/bin/env Rscript
## Script name: LoadInstallPackages.R
## Purpose: Check if the required dependencies are installed and loaded
## Date Created: August 5, 2022
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025


# Define the specific versions and mirrors
snapshot_date <- "2025-05-01"
bioc_version  <- "3.22"

options(
    download.file.method = "curl",
    repos = c(
        CRAN = sprintf("https://packagemanager.posit.co/cran/%s", snapshot_date),
        BiocManager::repositories(version = bioc_version)
    )
)

cbioLn <- "/Volumes/CBioinformatics/Methylation/classifiers"
epicV2script <-
    "/Volumes/CBioinformatics/Methylation/Rscripts/install_epic_v2_classifier.R"

# Function: Returns the directory where brew is installed if not then NULL ----
getBrewDir <- function(module_name) tryCatch(
    system(sprintf("brew --prefix %s", module_name), intern = TRUE, ignore.stderr = TRUE),
    error   = function(e) NULL,
    warning = function(e) NULL
)

proj_path <- getBrewDir("proj")
sqli_path <- getBrewDir("sqlite")
uniD_path <- file.path(cbioLn, "UniD")

Sys.setenv(PROJ_LIBS = file.path(proj_path, "lib"))
Sys.setenv(SQLITE3_LIBS = file.path(sqli_path, "lib"))
typeSrc <- "binary"
classPacks <- c(
    mgmtstp27 = "https://github.com/badozor/mgmtstp27/raw/master/archive/mgmtstp27_0.6-3.tar.gz",
    mnpqc = paste0(file.path(cbioLn, "mnpqc_0.1.0.tar.gz"))
)

if (Sys.info()[['sysname']] != "Darwin") {typeSrc <- "source"}

# Cran Packages ---------------------------------------------------------------
cranPkgs <- c(
    # Fundamental utilities and core dependencies
    "Rcpp",
    "rlang",
    "digest",
    "RCurl",
    "R.utils",
    "future.apply",
    "fastmatch",

    # Development tools
    "devtools",
    "usethis",

    # Data manipulation and transformation
    "data.table",
    "reshape2",
    "stringr",

    # Plotting and graphics
    "ggplot2",
    "scales",
    "RColorBrewer",
    "gridExtra",
    "ggfortify",
    "ggrepel",
    "gplots",
    "plotly",

    # Reporting and document generation
    "knitr",
    "kableExtra",
    "xtable",
    "pander",

    # Web and interactive tools
    "htmltools",
    "webshot",
    "DT",

    # Miscellaneous tools and specialized packages
    "sjmisc",
    "qdapTools",
    "beepr",
    "magick",
    "ade4",
    "targets",
    "readxl",
    "tinytex",
    "httpuv",
    "shiny",
    "pals",
    "Polychrome"
)


# GitHub Packages -------------------------------------------------------------
gHubPkgs <- c(
    'rstudio/rmarkdown',
    'dasonk/docstring',
    'nutterb/redcapAPI',
    "r-lib/crayon"
)

# BioConductor Packages  ------------------------------------------------------
biocPkgs <- c(
    # Core Bioconductor infrastructure and utilities
    "BiocParallel",
    "Biobase",
    "GenomeInfoDb",
    "MatrixGenerics",
    "rngtools",

    # Data acquisition and storage
    "GEOquery",
    "HDF5Array",

    # Illumina methylation annotations
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "IlluminaHumanMethylationEPICmanifest",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",

    # Methylation data processing
    "lumi",
    "methylumi",
    "MethylAid",
    "bumphunter",
    "limma",
    "conumee",
    "randomForest",
    "glmnet",
    "Rtsne"
)


# Extra Libraries -------------------------------------------------------------
extraPks <- c(
    "doParallel",
    "foreach",
    "tidyverse",
    "plyr",
    "dplyr",
    "sqldf",
    "mlr",
    "wateRmelon",
    "RPMM",
    "impute",
    "cli",
    "progressr"
)

# Helper Functions ------------------------------------------------------------

# Loads a library if not in current NameSpace
load_lib <- function(libName) {
    lib.opts <- list(package = libName, character.only = T, verbose = F)
    if (!libName %in% loadedNamespaces()) {
        suppressPackageStartupMessages(do.call(library, c(lib.opts)))
    }
    message(libName, " ...load successful")
}

# Checks if a package is installed
not_installed <- function(pkgName) return(!pkgName %in% rownames(installed.packages()))

# FUNC: Installs package from github link
gh.inst <- function(git_repo, ...) {
    pkg <- basename(git_repo)
    message("Checking ", pkg, "...")
    gh.opt <- list(pkg, dependencies = T, upgrade = "never", ...)
    if (not_installed(pkg)) {
        params <- c(gh.opt, list(type = "source"))
        tryCatch(
            do.call(devtools::install_github, params),
            error = function(cond) {
                params <- c(gh.opt, list(type = typeSrc))
                do.call(devtools::install_github, params)
            }
        )
    }
    load_lib(pkg)
}

# FUNC: Installs package from Source link
source_inst <- function(pkg) {
    message("Checking ", names(pkg), "...")
    if (not_installed(names(pkg))) {
        params <- list(pkgs = pkg[[1]], repos = NULL, type = "source", dependencies = T)
        tryCatch(
            do.call(install.packages, c(params, Ncpus = 6)),
            error = function(cond) {
                do.call(install.packages, c(params, method = "libcurl"))
            },
            warning = function(cond) {
                try(do.call(install.packages, c(params, method = "auto")), TRUE)
            }
        )
    }
    load_lib(names(pkg))
}

# FUNC: Installs package from Bioconductor
bio_install <- function(pknm) {
    message("Checking ", pknm, "...")
    bio.opt <- list(pkgs = pknm, update = F, ask = FALSE, dependencies = T)
    if (not_installed(pknm)) {
        tryCatch(
            do.call(BiocManager::install, c(bio.opt)),
            error = function(cond) {
                message(cond)
                pak::pkg_install(pknm, ask = FALSE)
            }
        )
    }
    load_lib(pknm)
}

# FUNC: Installs and updates autoloading of the needs package
manage_needs <- function() {
    # Install Needs
    tryCatch(
        install.packages("needs", dependencies = TRUE, verbose = TRUE,
                         Ncpus = 4, ask = FALSE, type = "binary"),
        error = function(cond) {
            devtools::install_github("joshkatz/needs", ref = "development",
                                     dependencies = TRUE, verbose = TRUE)
        }
    )
    # Configure Needs
    sysfile <- system.file("extdata", "promptUser", package = "needs")
    options(needs.promptUser = FALSE)
    invisible(needs:::autoload(TRUE))
    if (file.exists(sysfile) & length(readLines(sysfile, warn = FALSE)) == 0) {
        try(write(0, file = sysfile), silent = TRUE)
    }
    # Update .Rprofile
    txt1 <- "^[:blank:]*autoload\\(\"needs\", \"needs\"\\)"
    txt2 <- "\n\nautoload(\"needs\", \"needs\")\n\n"
    profPath <- Sys.getenv("R_PROFILE", unset = NA)
    siteProf <- ifelse(is.na(profPath) || profPath == "",
                       file.path(Sys.getenv("HOME"), "Rprofile.site"), profPath
                       )
    if (!file.exists(siteProf)) try(file.create(siteProf), silent = TRUE)
    if (file.exists(siteProf)) {
        lines <- readLines(siteProf, warn = FALSE)
        if (!any(grepl(txt1, lines))) {
            try(write(txt2, file = siteProf, append = TRUE), silent = TRUE)
        }
    }
}

# FUNC: ensures the main package installers are loaded
loadMainPkgs <- function() {
    cpus = as.integer(system("sysctl -n hw.logicalcpu", intern = TRUE))
    options(Ncpus = cpus - 2)

    if (not_installed("devtools")) install.packages("devtools", dependencies = TRUE, ask = FALSE)
    if (not_installed("BiocManager")) install.packages("BiocManager", dependencies = TRUE, ask = FALSE)
    load_lib("devtools")
    load_lib("BiocManager")

    if (not_installed("zip")) install.packages("zip", dependencies = TRUE, ask = FALSE)
    load_lib("zip")

    if (not_installed("pak")) install.packages("pak", type = "binary", ask = FALSE, dependencies = T)
    library("pak")
    # CRAN from PPM snapshot (compact spec), Bioconductor official URLs
    pak::repo_add(CRAN = sprintf("PPM@%s", snapshot_date))
    
    # Supply Bioconductor repos explicitly
    bioc_repos <- BiocManager::repositories(version = bioc_version)
    bioc_repos <- bioc_repos[!names(bioc_repos) %in% c("CRAN")]
    pak::repo_add(.list = as.list(bioc_repos))

    if (not_installed("needs")) manage_needs()
    options(needs.promptUser = FALSE)
    options(promptUser = FALSE)
    options(device.ask.default = FALSE)
    options(install.packages.check.source = "yes")
    gh.inst("jakesherman/easypackages")
    if (not_installed("tidyverse")) pak::pkg_install("tidyverse", ask = FALSE)
}

# FUNC: Applies each function install type (Github, CRAN, BioConductor)
load_install_pkgs <- function(pkgList, instFun) {
    invisible(lapply(X = 1:length(pkgList), function(X) instFun(pkgList[X])))
}

# FUNC: Ensures the correct minfi and manifest are installed for the classifier
install_minfi <- function() {
    if (isNamespaceLoaded("minfi")) try(unloadNamespace("minfi"), silent = TRUE)
    if (not_installed("IlluminaHumanMethylationEPICv2manifest"))
        devtools::install_github("mwsill/IlluminaHumanMethylationEPICv2manifest", upgrade = "always", force = T, dependencies = T)
    devtools::install_github("mwsill/minfi", upgrade = "never", force = T, dependencies = T)
}

# FUNC: Checks if the package in the input list is already installed
checkNeeded <- function(pkgList) {
    neededPkgs <- pkgList[!pkgList %in% rownames(installed.packages())]
    return(neededPkgs)
}

# FUNC: Ensures correct version of Plotly to avoid issues loading ggplot2
check_plotly <- function(){
    plotly_version <- package_version("4.10.4")

    if (packageVersion("plotly") < plotly_version) {
        if ("package:plotly" %in% search()) {
            try(detach("package:plotly", unload = TRUE, character.only = TRUE), TRUE)
        }
        if ("plotly" %in% loadedNamespaces()) {
            try(unloadNamespace("plotly"), TRUE)
        }
        library(pak)
        pak::pkg_install("plotly", upgrade = TRUE, ask = FALSE, dependencies = TRUE)
    }
}

# FUNC: Loads the list of packages from CRAN, GitHub, and BioConductor needed
loadPacks <- function() {
    loadMainPkgs()

    neededPkgs <- checkNeeded(cranPkgs)
    if (length(neededPkgs) > 0) pak::pkg_install(neededPkgs, ask = FALSE)

    if (not_installed("mdthemes")) {
        remotes::install_github("thomas-neitmann/mdthemes", upgrade = "never")
    }
  
    if (not_installed("minfi")) install_minfi()

    load_install_pkgs(gHubPkgs, gh.inst)

    if (not_installed("MethylAid")) {
        tryCatch(
            BiocManager::install("MethylAid", dependencies = T, ask = FALSE, update = F),
            error = function(e) pak::pkg_install("MethylAid", ask = FALSE)
        )
    }
    check_plotly()
    neededPkgs <- checkNeeded(biocPkgs)
    if (length(neededPkgs) > 0) pak::pkg_install(neededPkgs, ask = FALSE)

    suppressWarnings(load_install_pkgs(biocPkgs, bio_install))

    if (as.character(utils::packageVersion("minfi")) != "1.43.1") install_minfi()
    if (not_installed("librarian")) install.packages("librarian", dependencies = T, ask = FALSE)

    neededPkgs <- checkNeeded(extraPks)
    if (length(neededPkgs) > 0) pak::pkg_install(neededPkgs, ask = FALSE)

    load_lib('grid')
    load_lib("pals")
    load_lib("stringr")
}

# FUNC: Ensures the network drives are mounted
checkMounts <- function() {
    paths <- c("/Volumes/CBioinformatics/" = "CBioinformatics", "/Volumes/molecular/Molecular" = "Molecular")
    missingPaths <- names(paths)[!dir.exists(names(paths))]
    if (length(missingPaths) > 0)
        stop("You do not have the following drives mounted:\n",
                paste(missingPaths, collapse = "\n"))
}

# FUNC: Ensures the latest classifier version is installed
checkEpicV2 <- function(pkg, epicV2script) {
    if (pkg %in% rownames(installed.packages())) {
        message(paste("Package", pkg, "is installed with version:", utils::packageVersion(pkg)))
    } else{
        message("Installing package from source:\n", epicV2script)
        source(epicV2script)
    }
}

# FUNC: Checks if the package UniD and it's requirements are installed
check_uniD_pkg <- function() {
    uniDpkgs <- c("ade4", "methylumi", "mlr")
    librarian::shelf(uniDpkgs, ask = FALSE, update_all = F, quiet = F)
    try(bio_install("impute"), silent = T)
    try(bio_install("wateRmelon"), silent = T)
    if (not_installed("UniD")) {
        try(install.packages(uniD_path, type = "source", dependencies = T, repo = NULL), silent = T)
    }
}

# MAIN: Loads all packages and functions --------------------------------------
startLoadingAll <- function(typeSrc, epicV2script) {
    message("Your R library path(s):\n", paste(.libPaths(), collapse = "\n"))
    checkMounts()
    library("base")
    Sys.setenv(TORCH_INSTALL = "1")
    options(needs.auto = TRUE)
    loadPacks()
    suppressWarnings(load_install_pkgs(classPacks, source_inst))
    #check_uniD_pkg()
    checkEpicV2("mnp.v12epicv2", epicV2script)
}

startLoadingAll(typeSrc, epicV2script)
