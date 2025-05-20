#!/usr/bin/env Rscript
## Script name: LoadInstallPackages.R
## Purpose: Check if the required dependencies are installed and loaded
## Date Created: August 5, 2022
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025

isMacOS <- Sys.info()[['sysname']] == "Darwin"
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

# Additional Classifier package Dependencies:
if (isMacOS) {
    Sys.setenv(PROJ_LIBS = file.path(proj_path, "lib"))
    Sys.setenv(SQLITE3_LIBS = file.path(sqli_path, "lib"))
    typeSrc <- "binary"
    bothType <- "both"
    classPacks <- c(
        mgmtstp27 = "https://github.com/badozor/mgmtstp27/raw/master/archive/mgmtstp27_0.6-3.tar.gz",
        mnpqc = paste0(file.path(cbioLn, "mnpqc_0.1.0.tar.gz"))
    )
} else{
    options(BioC_mirror = "https://packagemanager.rstudio.com/bioconductor")
    options(repos = c(CRAN = "https://packagemanager.rstudio.com/cran/__linux__/jammy/latest"))
    bothType <- typeSrc <- "source"
    classPacks <- c(
        sest = "https://github.com/jungch/sest/raw/master/sest.tar",
        mgmtstp27 = "https://github.com/badozor/mgmtstp27/raw/master/archive/mgmtstp27_0.6-3.tar.gz"
    )
}

# Cran Packages ---------------------------------------------------------------
cranPkgs <- c(
    # Fundamental utilities and core dependencies
    "Rcpp",
    "rlang",
    "digest",
    "RCurl",
    "R.utils",
    "optparse",
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
    #'rstudio/rstudioapi',
    'nutterb/redcapAPI',
    "r-lib/crayon",
    "thomas-neitmann/mdthemes"
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

    # Illumina methylation annotation (450k)
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",

    # Illumina methylation annotation (EPIC)
    "IlluminaHumanMethylationEPICmanifest",
    #"IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",

    # Methylation data processing
    "lumi",
    "methylumi",
    "MethylAid",

    # Analysis and machine learning methods
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
    if (!libName %in% loadedNamespaces()) suppressPackageStartupMessages(do.call(library, c(lib.opts)))
    message(libName, " ...load successful")
}

# Checks if a package is installed
not_installed <- function(pkgName) return(!pkgName %in% rownames(installed.packages()))

# Messages the current package being checked for installation
msgCheck <- function(pkg, warn = F) {
    if (warn == F) message("Checking ", pkg, "...") else message("\nError caught for package:\n", pkg, "\n-----------")
}


# FUNC: Installs package from github link
gh.inst <- function(git_repo, ...) {
    pkg <- basename(git_repo)
    msgCheck(pkg)
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
source_inst <- function(fn) {
    msgCheck(names(fn))
    if (not_installed(names(fn))) {
        params <- list(pkgs = fn[[1]], repos = NULL, type = "source", dependencies = T)
        tryCatch(
            do.call(install.packages, c(params, Ncpus = 6)),
            error = function(cond) {
                msgCheck(cond, T)
                do.call(install.packages, c(params, method = "libcurl"))
            },
            warning = function(cond) {
                msgCheck(cond, T)
                try(do.call(install.packages, c(params, method = "auto")), TRUE)
            }
        )
    }
    load_lib(names(fn))
}

# FUNC: Installs package from Bioconductor
bc.inst <- function(pknm) {
    message("Checking ", pknm, "...")
    bio.opt <- list(pkgs = pknm, update = F, ask = F, dependencies = T)
    if (not_installed(pknm)) {
        tryCatch(
            do.call(BiocManager::install, c(bio.opt)),
            error = function(cond) {
                message("Error with Package install:\n", pknm)
                message(cond)
                pak::pkg_install(pknm, ask = F)
            }
        )
    }
    load_lib(pknm)
}

# FUNC: Installs and updates autoloading of the needs package
manageNeeds <- function() {

    updateRProfile <- function() {
        txt1 <- "^[:blank:]*autoload\\(\"needs\", \"needs\"\\)"
        txt2 <- "\n\nautoload(\"needs\", \"needs\")\n\n"
        profPath <- Sys.getenv("R_PROFILE", unset = NA)
        siteProf <- ifelse(
            is.na(profPath) || profPath == "",
            file.path(Sys.getenv("HOME"), "Rprofile.site"),
            profPath
        )

        if (!file.exists(siteProf)) try(file.create(siteProf), silent = TRUE)

        if (file.exists(siteProf)) {
            lines <- readLines(siteProf, warn = FALSE)
            if (!any(grepl(txt1, lines))) {
                try(write(txt2, file = siteProf, append = TRUE), silent = TRUE)
            }
        }
    }

    configureNeeds <- function() {
        sysfile <- system.file("extdata", "promptUser", package = "needs")
        options(needs.promptUser = FALSE)
        invisible(needs:::autoload(TRUE))
        if (file.exists(sysfile)) {
            if (length(readLines(sysfile, warn = FALSE)) == 0) {
                try(write(0, file = sysfile), silent = TRUE)
            }
        }
    }

    installNeeds <- function() {
        tryCatch(
            install.packages("needs", dependencies = TRUE, verbose = TRUE, Ncpus = 4, ask = FALSE, type = "binary"),
            error = function(cond) {
                devtools::install_github("joshkatz/needs", ref = "development", dependencies = TRUE, verbose = TRUE, upgrade = "always")
            }
        )
    }

    if (not_installed("needs")) {
        installNeeds()
        configureNeeds()
        updateRProfile()
    }

}

# FUNC: ensures the main package installers are loaded
loadMainPkgs <- function() {
    if (not_installed("devtools")) install.packages("devtools", dependencies = TRUE, ask = F)
    if (not_installed("BiocManager")) install.packages("BiocManager", dependencies = TRUE, ask = F)
    load_lib("devtools")
    load_lib("BiocManager")

    if (not_installed("zip")) install.packages("zip", dependencies = TRUE, ask = F)
    load_lib("zip")

    if (not_installed("pak")) install.packages("pak", type = "binary", ask = F, dependencies = T)
    library("pak")
    pak::repo_add(PPM = "PPM@2025-04-01")

    manageNeeds()
    options(needs.promptUser = FALSE)
    options(promptUser = FALSE)
    options(device.ask.default = FALSE)
    options(install.packages.check.source = "yes")

    gh.inst("jakesherman/easypackages")
    if (not_installed("tidyverse")) pak::pkg_install("tidyverse", ask = F)
}

# FUNC: Applies each function install type (Github, CRAN, BioConductor) to Pkgs
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
    if (length(neededPkgs) > 0) pak::pkg_install(neededPkgs, ask = F)

    if (not_installed("minfi")) install_minfi()

    load_install_pkgs(gHubPkgs, gh.inst)

    if (not_installed("MethylAid")) {
        tryCatch(
            BiocManager::install("MethylAid", dependencies = T, ask = F, update = F),
            error = function(e) pak::pkg_install("MethylAid", ask = F)
            )
    }
    check_plotly()
    neededPkgs <- checkNeeded(biocPkgs)
    if (length(neededPkgs) > 0) pak::pkg_install(neededPkgs, ask = F)

    suppressWarnings(load_install_pkgs(biocPkgs, bc.inst))

    if (as.character(utils::packageVersion("minfi")) != "1.43.1") install_minfi()
    if (not_installed("librarian")) install.packages("librarian", dependencies = T, ask = F)

    neededPkgs <- checkNeeded(extraPks)
    if (length(neededPkgs) > 0) pak::pkg_install(neededPkgs, ask = F)

    load_lib('grid')
    load_lib("pals")
    load_lib("stringr")
}

# FUNC: Ensures the network drives are mounted
checkMounts <- function() {
    paths <- c("/Volumes/CBioinformatics/" = "CBioinformatics", "/Volumes/molecular/Molecular" = "Molecular")
    message("You have the following drives mounted:")
    system("ls /Volumes")
    missingPaths <- names(paths)[!dir.exists(names(paths))]
    if (length(missingPaths) > 0)
        message("You do not have the following drives mounted:\n", paste(missingPaths, collapse = "\n"))
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
    librarian::shelf(uniDpkgs, ask = F, update_all = F, quiet = F)
    try(bc.inst("impute"), silent = T)
    try(bc.inst("wateRmelon"), silent = T)
    if (not_installed("UniD")) {
        try(install.packages(uniD_path, type = "source", dependencies = T, repo = NULL), silent = T)
    }
}

# MAIN: Loads all packages and functions --------------------------------------
startLoadingAll <- function(isMacOS, typeSrc, epicV2script) {
    message("Your R library path(s):\n", paste(.libPaths(), collapse = "\n"))
    options(Ncpus = 6)
    checkMounts()
    library("base")
    options(repos = c(CRAN = "http://cran.us.r-project.org"))
    loadPacks()
    suppressWarnings(load_install_pkgs(classPacks, source_inst))
    if (isMacOS) {
        #check_uniD_pkg()
        checkEpicV2("mnp.v12epicv2", epicV2script)
    }
}

startLoadingAll(isMacOS, typeSrc, epicV2script)
