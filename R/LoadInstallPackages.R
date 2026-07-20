#!/usr/bin/env Rscript
## Script name: LoadInstallPackages.R
## Purpose: Check if the required dependencies are installed and loaded
## Date Created: August 5, 2022
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2026

ncpus <- suppressWarnings(try(parallel::detectCores(logical = TRUE), silent = TRUE))
if (inherits(ncpus, "try-error") || length(ncpus) != 1L || is.na(ncpus) || ncpus < 1L) ncpus <- 1L

core_deps <- c("Depends", "Imports", "LinkingTo")
install_params <- list(dependencies = core_deps, clean = TRUE, ask = FALSE,
                       INSTALL_opts = "--no-multiarch")
interactive <- function() FALSE
http_user_agent <- sprintf("R/%s R (%s; %s; %s)", getRversion(), R.version$platform,
                           R.version$arch, R.version$os)

options(
    repos = c(CRAN = "https://posit.co"),
    BioC_mirror = "https://bioconductor.org",
    Ncpus = max(1L, ncpus - 1L),
    timeout = 1200L,
    askYesNo = function(...) { return(FALSE) },
    download.file.method = "libcurl",
    HTTPUserAgent = http_user_agent,
    device.ask.default = FALSE,
    install.packages.compile.from.source = "never",
    install.packages.check.source = "no",
    needs.promptUser = FALSE,
    promptUser = FALSE,
    needs.auto = TRUE,
    warn = 1L
)

Sys.setenv(
    R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true",
    BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS = "TRUE",
    PAK_NO_UPDATE_MESSAGE = "true",
    LIBARROW_BINARY = "FALSE",
    LIBARROW_BUILD = "FALSE",
    ARROW_USE_PKG_CONFIG = "TRUE",
    NOT_CRAN = "FALSE",
    ARROW_R_DEV = "TRUE"
)

if (!interactive()) Sys.setenv(RGL_USE_NULL = "true", TORCH_INSTALL = "1")
if (identical(R.version$os, "linux-gnu")) Sys.setenv(DEBIAN_FRONTEND = "noninteractive")

# Checks if a package is installed
not_installed <- function(pkgName) return(!pkgName %in% rownames(installed.packages()))

# Installs package if missing
check_install <- function(pkgs) {
    for (pkg in pkgs) if (not_installed(pkg)) do.call(install.packages, c(pkg, install_params))
}

# Loads a library if not in current NameSpace
load_lib <- function(pkg, loaded = FALSE) {
    lib_opts <- list(package = pkg, character.only = TRUE, verbose = FALSE, logical.return = TRUE)
    if (!pkg %in% loadedNamespaces()) {
        loaded <- suppressPackageStartupMessages(do.call(library, c(lib_opts)))
    }
    return(message(pkg, " ...loaded ", loaded))
}

core_pkgs <- c("devtools", "remotes", "BiocManager", "pak")
check_install(core_pkgs)
invisible(lapply(core_pkgs, load_lib))

bioc_params <- list(dependencies = core_deps, clean = TRUE, INSTALL_opts = "--no-multiarch",
                    ask = FALSE, update = FALSE)
utils::assignInNamespace("install", function(...) do.call(
    BiocManager::install, utils::modifyList(bioc_params, list(...))), ns = "BiocManager")

suppressMessages(pak::repo_add(.list = as.list(suppressMessages(BiocManager::repositories()))))

# Function: Returns the directory where brew is installed if not then NULL ----
getBrewDir <- function(module_name) tryCatch(
    system(sprintf("brew --prefix %s", module_name), intern = TRUE, ignore.stderr = TRUE),
    error   = function(e) NULL,
    warning = function(e) NULL
)

proj_path <- "/usr/share/proj"
sqli_path <- "/usr/lib/x86_64-linux-gnu"

if (Sys.info()[["sysname"]] == "Darwin") {
    proj_path_mac <- getBrewDir("proj")
    sqli_path_mac <- getBrewDir("sqlite")
    if (length(proj_path_mac) == 1L) {
        Sys.setenv(PROJ_LIB = file.path(proj_path_mac, "share", "proj"))
    }
    if (length(sqli_path_mac) == 1L) {
        Sys.setenv(SQLITE3_LIBS = file.path(sqli_path_mac, "lib"))
    }
} else {
    Sys.setenv(PROJ_LIB = proj_path, SQLITE3_LIBS = sqli_path)
}

# Custom packages and installer scripts -----------------------------------------------------------
cbioLn <- "/Volumes/CBioinformatics/Methylation/classifiers"
epicV2script <-
    "/Volumes/CBioinformatics/Methylation/Rscripts/install_epic_v2_classifier.R"
classPacks <- c(
    mgmtstp27 = "https://github.com/badozor/mgmtstp27/raw/master/archive/mgmtstp27_0.6-3.tar.gz",
    mnpqc = paste0(file.path(cbioLn, "mnpqc_0.1.0.tar.gz"))
)

# Cran Packages -----------------------------------------------------------------------------------
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
check_pak <- function(pkg) { if (not_installed(pkg)) pak::pkg_install(pkg, ask = FALSE) }

# FUNC: Installs package from github link
gh.inst <- function(git_repo, ...) {
    if (not_installed(basename(git_repo))) {
        tryCatch(
            pak::pak(repo = git_repo, dependencies = core_deps, ask = FALSE, ...),
            error = function(e) message("GitHub package could not install from repo:", git_repo)
        )
    }
    load_lib(basename(git_repo))
}

# FUNC: Installs package from Source link
source_inst <- function(pkg) {
    if (not_installed(names(pkg))) {
        params <- list(pkgs = pkg[[1]], repos = NULL, type = "source", dependencies = core_deps)
        tryCatch(
            do.call(install.packages, c(params)),
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
bio_install <- function(pkg) {
    if (not_installed(pkg)) {
        tryCatch(
            BiocManager::install(pkg),
            error = function(cond) pak::pkg_install(pkg, ask = FALSE)
        )
    }
    load_lib(pkg)
}

# FUNC: Installs and updates autoloading of the needs package
manage_needs <- function() {
    # Install Needs
    tryCatch(
        pak::pkg_install("needs", ask = FALSE),
        error = function(cond) {
            pak::pak("joshkatz/needs", ask = FALSE, dependencies = core_deps)
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
    
}

# FUNC: Applies each function install type (Github, CRAN, BioConductor)
load_install_pkgs <- function(pkgList, instFun) {
    invisible(lapply(X = 1:length(pkgList), function(X) instFun(pkgList[X])))
}

# FUNC: Ensures the correct minfi and manifest are installed for the classifier
install_minfi <- function() {
    if (isNamespaceLoaded("minfi")) try(unloadNamespace("minfi"), silent = TRUE)
    if (not_installed("IlluminaHumanMethylationEPICv2manifest"))
        pak::pak("mwsill/IlluminaHumanMethylationEPICv2manifest",
                 upgrade = FALSE, ask = FALSE, dependencies = core_deps)
    pak::pak("mwsill/minfi", upgrade = FALSE, ask = FALSE, dependencies = core_deps)
}

# FUNC: Checks if the package in the input list is already installed
checkNeeded <- function(pkgList) {
    neededPkgs <- pkgList[!pkgList %in% rownames(installed.packages())]
    if (length(neededPkgs) > 0) {
        pak::pkg_install(neededPkgs, ask = FALSE)
    }
}

# FUNC: Ensures correct version of Plotly to avoid issues loading ggplot2
check_plotly <- function(){
    plotly_version <- package_version("4.10.4")
    
    if (packageVersion("plotly") < plotly_version) {
        if ("package:plotly" %in% search()) {
            try(detach("package:plotly", unload = TRUE, character.only = TRUE), TRUE)
        }
        if ("plotly" %in% loadedNamespaces()) try(unloadNamespace("plotly"), TRUE)
        library("pak")
        pak::pkg_install("plotly", upgrade = TRUE, ask = FALSE, dependencies = core_deps)
    }
}

# FUNC: Loads the list of packages from CRAN, GitHub, and BioConductor needed
load_packages <- function() {
    check_install("zip")
    load_lib("zip")
    
    if (not_installed("needs")) manage_needs() else try(invisible(needs:::autoload(TRUE)), silent = TRUE)
    gh.inst("jakesherman/easypackages")
    check_install("tidyverse")
    
    checkNeeded(cranPkgs)
    
    if (not_installed("mdthemes")) {
        remotes::install_github("thomas-neitmann/mdthemes", upgrade = "never", ask = FALSE)
    }
    
    if (not_installed("minfi")) install_minfi()
    
    load_install_pkgs(gHubPkgs, gh.inst)
    
    if (not_installed("MethylAid")) {
        tryCatch(
            BiocManager::install("MethylAid", dependencies = core_deps, ask = FALSE, update = F),
            error = function(e) pak::pkg_install("MethylAid", ask = FALSE)
        )
    }
    check_plotly()
    checkNeeded(biocPkgs)
    
    suppressWarnings(load_install_pkgs(biocPkgs, bio_install))
    
    if (as.character(utils::packageVersion("minfi")) != "1.43.1") install_minfi()
    check_install("librarian")
    
    checkNeeded(extraPks)
    invisible(lapply(c('grid', 'pals', 'stringr'), load_lib))
}

# FUNC: Ensures the network drives are mounted
checkMounts <- function() {
    paths <- c("/Volumes/CBioinformatics/" = "CBioinformatics", "/Volumes/molecular/Molecular" = "Molecular")
    missed_dir <- names(paths)[!dir.exists(names(paths))]
    if (length(missed_dir) > 0)
        stop("You do not have the following drives mounted:\n", paste(missed_dir, collapse = "\n"))
}

# MAIN: Loads all packages and functions --------------------------------------
startLoadingAll <- function() {
    message("Your R library path(s):\n", paste(.libPaths(), collapse = "\n"))
    if (Sys.info()[["sysname"]] == "Darwin") { checkMounts() }
    load_packages()
    suppressWarnings(load_install_pkgs(classPacks, source_inst))
    
    mnp <- "mnp.v12epicv2"
    if (mnp %in% rownames(installed.packages())) {
        message(paste("Package", mnp, "is installed with version:", utils::packageVersion(mnp)))
    } else{
        message("Installing package from source:\n", epicV2script)
        source(epicV2script)
    }
}

startLoadingAll()
