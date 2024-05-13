#!/usr/bin/env Rscript
## ---------------------------
## Script name: all_installer.R
## Purpose: Functions to check if all required packages for the pipeline are installed
## Date Last Modified: January 11, 2024
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

options("install.packages.compile.from.source" = "No")
options("install.packages.check.source" = "no")

is_macos <- Sys.info()[['sysname']] == "Darwin"

# FUN: Check if brew installed ------------------------------------------------
install_homebrew <- function() {
    message("Homebrew is not installed. Installing Homebrew...")
    url <- "https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh"
    cmd <- sprintf('/bin/bash -c "$(curl -fsSL %s)"', url)
    message("Executing terminal command:\n", cmd)
    system(cmd, wait = TRUE)
}


# FUN: Install system package using homebrew ----------------------------------
brew_install <- function(pkg){
    message(pkg, " is not installed! Installing ", pkg, " via Homebrew...")
    system(paste("brew install", pkg), intern = TRUE, wait = TRUE)
}


# FUN: Checks if system compilers installed -----------------------------------
check_brew_pkgs <- function(){
    # Check if brew installed
    brew_exists <- file.exists("/usr/local/bin/brew")
    if (!brew_exists) {
        install_homebrew()
    }
    # Check if GCC installed
    gfortran_path <- try(Sys.which("gfortran")[[1]], T)
    if (gfortran_path == "") {
        message("gfortran is not installed. Please install gfortran.")
        brew_install("gcc")
    }
    # Check if LLVM installed
    llvm_installed <- file.exists("/usr/local/opt/llvm/bin/clang")
    if (!llvm_installed) {
        brew_install("llvm")
    }
    # Check if open-mpi installed
    mpi_installed <- file.exists("/usr/local/opt/open-mpi/bin")
    if (!mpi_installed) {
        brew_install("open-mpi")
    }

    clear_cmd <- "brew update && brew doctor && rm -rf $(brew --cache)"
    system(clear_cmd, wait = TRUE)
}


# FUN: Clears system environment flags ----------------------------------------
clear_enviro <- function(){
    Sys.setenv(CC = "")
    Sys.setenv(CFLAGS = "")
    Sys.setenv(CXX = "")
    Sys.setenv(CXXFLAGS = "")
    Sys.setenv(FC = "")
    Sys.setenv(FFLAGS = "")
    Sys.setenv(LDFLAGS = "")
    Sys.setenv(CPPFLAGS = "")
    Sys.setenv(SHLIB_CXXLD = "")
    Sys.setenv(SHLIB_LDFLAGS = "")
    Sys.setenv(OBJC = "")
}


# FUN: Sets system gfortran flags ---------------------------------------------
set_gfortran <- function(){
    gfortran_path <- Sys.which("gfortran")[[1]]

    tran_info <- system("gfortran -print-search-dirs", intern = TRUE)

    library_paths <- strsplit(tran_info[length(tran_info)], " ")[[1]][-1]
    flibs_value <- paste(sapply(library_paths, function(path) {
        paste("-L", path, sep = "")}), collapse = " ")

    Sys.setenv(FC = gfortran_path)
    Sys.setenv(FLIBS = flibs_value)
}


# FUN: Sets system openmpi flags ----------------------------------------------
set_openmpi <- function(){
    prte_path <- "/usr/local/opt/open-mpi/bin/prte"
    orte_path <- "/usr/local/opt/open-mpi/bin/orted"
    ln_cmd <- paste("ln -s", prte_path, orte_path)
    try(system(ln_cmd, wait = T), T)

    mpi_libs <- "-L/usr/local/opt/open-mpi/lib"
    mpi_flag <- paste(mpi_libs, Sys.getenv("LDFLAGS"))
    Sys.setenv(LDFLAGS = mpi_flag)

    cpp_libs <- "-I/usr/local/opt/open-mpi/include"
    cpp_flag <- paste(cpp_libs, Sys.getenv("CPPFLAGS"))
    Sys.setenv(CPPFLAGS = cpp_flag)
}


# FUN: Sets system compiler flags ---------------------------------------------
fix_compiler_flags <- function(){
    check_brew_pkgs()
    clear_enviro()

    Sys.setenv(CC = "/usr/local/opt/llvm/bin/clang")
    Sys.setenv(CXX = "/usr/local/opt/llvm/bin/clang++")
    Sys.setenv(CXX11 = "/usr/local/opt/llvm/bin/clang++")
    Sys.setenv(CXX14 = "/usr/local/opt/llvm/bin/clang++")
    Sys.setenv(CXX17 = "/usr/local/opt/llvm/bin/clang++")
    Sys.setenv(CXX1X = "/usr/local/opt/llvm/bin/clang++")
    Sys.setenv(OBJC = "/usr/local/opt/llvm/bin/clang")
    Sys.setenv(LDFLAGS = "-L/usr/local/opt/llvm/lib")
    Sys.setenv(CPPFLAGS = "-I/usr/local/opt/llvm/include")

    set_openmpi()
    set_gfortran()

    curr <- paste("/usr/local/opt/llvm/bin", Sys.getenv("PATH"), sep = ":")
    Sys.setenv(PATH = paste("/usr/local/opt/open-mpi/bin", curr, sep = ":"))
}

update_system_path <- function() {
    user_shell <- Sys.getenv("SHELL")
    
    shell_config_file <- if (grepl("zsh", user_shell)) {
        path.expand("~/.zshrc")
    } else if (grepl("bash", user_shell)) {
        if (file.exists(path.expand("~/.bash_profile"))) {
            path.expand("~/.bash_profile")
        } else {
            path.expand("~/.bashrc")
        }
    } else {
        stop("Your shell is not supported for automatic PATH updates by this script.")
    }
    
    path_command <- sprintf('echo "export PATH=\\"/usr/local/sbin:$PATH\\"" >> %s', shell_config_file)
    
    if (system(path_command, intern = FALSE) == 0) {
        message(sprintf("Successfully updated %s to include /usr/local/sbin in PATH.", shell_config_file))
    } else {
        warning("Failed to update the PATH in the shell configuration file.")
    }
}

if (is_macos) {
    local({
        path <- sub(":/opt/homebrew/bin", ":/usr/local/homebrew/bin", Sys.getenv("PATH"))
        Sys.setenv(PATH = path)
    })
    options(BioC_mirror = "https://packagemanager.rstudio.com/bioconductor")
    #options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2024-02-20"))
    options(warn = -1)
    options(repos = c(CRAN = 'https://cloud.r-project.org'))
    #options(pkgType = "binary")
    # Determine the SDK path using a system call in R
    fix_compiler_flags()
    update_system_path()
    #sdk_path <- system("xcrun --show-sdk-path", intern = TRUE)
    #Sys.setenv(CFLAGS = paste("-isysroot", sdk_path))
    #Sys.setenv(CXXFLAGS = paste("-isysroot", sdk_path))
}


# Get system architecture
arch <- Sys.info()[["machine"]]
# Check if architecture is 'arm64' or 'x86_64'
if (arch != "x86_64" & is_macos == T) {
    # Set JAVA_HOME environment variable
    #java_home <- system("/usr/libexec/java_home -v 11", intern = TRUE)
    system("brew tap homebrew/cask-versions && brew install --cask temurin17")
    java_home <- system("which java", intern = TRUE)
    Sys.setenv(JAVA_HOME = java_home)
    message("JAVA_HOME set to ", java_home)
    try(install.packages("rJava", type = "binary", dependencies = T, ask = F), T)
}

if (Sys.info()[['sysname']] == "Darwin" &
    !file.exists(file.path("~", ".Renviron"))) {
    system("touch ~/.Renviron")
    fileConn <- file("~/.Renviron")
    params <- c('PATH="/usr/local/gfortran/bin:${PATH}"')
    writeLines(params, fileConn)
    close(fileConn)
    #closeAllConnections()
}

# Setting US CRAN REPO
rlis = getOption("repos")
rlis["CRAN"] = "http://cran.us.r-project.org"
options(repos = rlis)


loadLibrary <- function(pkgName) {
    suppressPackageStartupMessages(library(
        pkgName,
        quietly = TRUE,
        logical.return = TRUE,
        warn.conflicts = FALSE,
        character.only = TRUE
    ))
}

checkRequire <- function(pkgName) {
    return(!requireNamespace(pkgName, quietly = TRUE))
}

if (checkRequire("devtools")) {
    install.packages(
        "devtools",
        dependencies = T,
        verbose = T,
        ask = F
    )
}

if (checkRequire("librarian")) {
    install.packages(
        "librarian",
        dependencies = T,
        verbose = T,
        ask = F
    )
}
loadLibrary("devtools")

# if (Sys.info()[['sysname']] == "Darwin") {
#     isOpen <- system("which openssl", intern = TRUE)
#     if (!exists("isOpen")) {
#         system("brew update")
#         system("brew install openssl")
#         system("ln -sf /usr/local/opt/openssl/lib/libcrypto.3.dylib /usr/local/lib/")
#         system("ln -sf /usr/local/opt/openssl/lib/libssl.3.dylib /usr/local/lib/")
#         system("ln -sf /usr/local/Cellar/openssl@3/3.0.5/bin/openssl /usr/local/bin/openssl")
#     } else{
#         if (isOpen != "/usr/bin/openssl") {
#             system("ln -sf /usr/local/opt/openssl/lib/libcrypto.3.dylib /usr/local/lib/")
#             system("ln -sf /usr/local/opt/openssl/lib/libssl.3.dylib /usr/local/lib/")
#             system("ln -sf /usr/local/Cellar/openssl@3/3.0.5/bin/openssl /usr/local/bin/openssl")

#         }
#     }
# }

# List Classifier Core Packages -------------------------------------------------------------------
corePkgs <- c("randomForest", "glmnet", "ggplot2", "gridExtra", "knitr", "pander", "gmp")

# List Prerequisite Packages ----------------------------------------------------------------------
preReqPkgs <- c(
    'MASS',
    'ade4',
    'dbplyr',
    'filelock',
    'png',
    'rjson',
    'generics',
    'tidyselect',
    'RCurl',
    'biomaRt',
    'XML',
    'bit64',
    'plyr',
    'GenomeInfoDbData',
    'base64',
    'GenomicAlignments',
    'readr',
    'dplyr',
    'tidyr',
    'R.utils',
    'BH',
    'affy',
    'GenomicFeatures',
    'GenomicRanges',
    'annotate',
    'nleqslv',
    'preprocessCore',
    'RSQLite',
    'DBI',
    'AnnotationDbi',
    'reshape2',
    'matrixStats',
    'FDb.InfiniumMethylation.hg19',
    'S4Vectors',
    'IRanges',
    'GenomeInfoDb',
    'SummarizedExperiment',
    'genefilter',
    'IlluminaHumanMethylationEPICanno.ilm10b2.hg19',
    'DNAcopy',
    'rtracklayer',
    'Biostrings',
    'bumphunter',
    'beanplot',
    'nor1mix',
    'limma',
    'mclust',
    'reshape',
    'quadprog',
    'data.table',
    'HDF5Array',
    'BiocParallel'
)

# List Remaining Packages -------------------------------------------------------------------------
pkgs <- c(
    "abind",
    "animation",
    "ape",
    "aplot",
    "arrow",
    "askpass",
    "assertr",
    "assertthat",
    "batchtools",
    "beepr",
    "bezier",
    "Biobase",
    "BiocCheck",
    "BiocManager",
    "BiocStyle",
    "BiocVersion",
    "biocViews",
    "bookdown",
    "boot",
    "brio",
    "broom",
    "bs4Dash",
    "BSgenome",
    "bslib",
    "cachem",
    "Cairo",
    "callr",
    "car",
    "carData",
    "caret",
    "cellranger",
    "changepoint",
    "checkmate",
    "chromote",
    "chron",
    "circlize",
    "class",
    "classInt",
    "cli",
    "clipr",
    "clustermq",
    "clusterProfiler",
    "cmprsk",
    "cn.mops",
    "codetools",
    "colorspace",
    "commonmark",
    "compiler",
    "config",
    "conquer",
    "conumee",
    "CopyNumberPlots",
    "cpp11",
    "crayon",
    "credentials",
    "curl",
    "datasets",
    "debugme",
    "desc",
    "devtools",
    "dichromat",
    "diffobj",
    "digest",
    "dlm",
    "docstring",
    "doParallel",
    "doSNOW",
    "downloader",
    "DT",
    "dtplyr",
    "easypackages",
    "enrichplot",
    "ensembldb",
    "evaluate",
    "exomeCopy",
    "extrafont",
    "extrafontdb",
    "fastmap",
    "fBasics",
    "fontawesome",
    "forecast",
    "foreign",
    "formattable",
    "Formula",
    "fracdiff",
    "fresh",
    "fs",
    "future",
    "gdata",
    "gdtools",
    "gert",
    "getopt",
    "ggeffects",
    "ggforce",
    "ggfortify",
    "ggfun",
    "ggnewscale",
    "ggplotify",
    "ggraph",
    "ggrepel",
    "ggtext",
    "ggthemes",
    "ggtree",
    "ggupset",
    "gh",
    "gitcreds",
    "GlobalOptions",
    "globals",
    "glue",
    "gmodels",
    "GO.db",
    "gplots",
    "graph",
    "graphics",
    "graphlayouts",
    "grDevices",
    "grid",
    "gridBase",
    "gridGraphics",
    "gridtext",
    "gsubfn",
    "gt",
    "gtools",
    "here",
    "hexbin",
    "highr",
    "Hmisc",
    "maptools",
    "hrbrthemes",
    "htmlTable",
    "htmltools",
    "htmlwidgets",
    "httr",
    "iterators",
    "igraph",
    "impute",
    "ini",
    "insight",
    "ipred",
    "jpeg",
    "jquerylib",
    "jsonlite",
    "kableExtra",
    "labelVector",
    "lars",
    "later",
    "lattice",
    "latticeExtra",
    "lazyeval",
    "lintr",
    "listenv",
    "lme4",
    "lmodel2",
    "lmtest",
    "lubridate",
    "magrittr",
    "mapdata",
    "mapproj",
    "maps",
    "markdown",
    "Matrix",
    "MatrixModels",
    "memoise",
    "methods",
    "MethylAid",
    "mlr",
    "modelr",
    "modeltools",
    "narray",
    "needs",
    "nlme",
    "nloptr",
    "nnet",
    "numDeriv",
    "openssl",
    "openxlsx",
    "optparse",
    "pals",
    "parallel",
    "parallelMap",
    "parameters",
    "ParamHelpers",
    "patchwork",
    "pathview",
    "pbkrtest",
    "performance",
    "pillar",
    "pingr",
    "pkgbuild",
    "pkgconfig",
    "pkgload",
    "plot3D",
    "plotly",
    "plotrix",
    "Polychrome",
    "polyclip",
    "polynom",
    "pracma",
    "praise",
    "prettyunits",
    "processx",
    "prodlim",
    "progressr",
    "promises",
    "ProtGenerics",
    "proto",
    "proxy",
    "ps",
    "purrr",
    "qdapTools",
    "qpdf",
    "qs",
    "quantreg",
    "qvalue",
    "R.cache",
    "raster",
    "RBGL",
    "Rcpp",
    "readxl",
    "recipes",
    "redcapAPI",
    "rematch",
    "rematch2",
    "remotes",
    "reprex",
    "reticulate",
    "rex",
    "Rgraphviz",
    "rio",
    "rJava",
    "rlang",
    "rmarkdown",
    "Rmpfr",
    "RMySQL",
    "robustbase",
    "roxygen2",
    "rpart",
    "RPMM",
    "rprojroot",
    "rstudioapi",
    "Rtsne",
    "rversions",
    "rvest",
    "sass",
    "scatterplot3d",
    "selectr",
    "sessioninfo",
    "shadowtext",
    "shiny",
    "sjmisc",
    "sourcetools",
    "spatial",
    "sqldf",
    "stats",
    "stats4",
    "stringdist",
    "stringi",
    "stringr",
    "strucchange",
    "svglite",
    "sys",
    "systemfonts",
    "targets",
    "tcltk",
    "tensorflow",
    "terra",
    "testit",
    "testthat",
    "tibble",
    "tidygraph",
    "tidytree",
    "tidyverse",
    "tiff",
    "timeDate",
    "timeSeries",
    "tinytex",
    "tools",
    "torch",
    "tseries",
    "units",
    "utf8",
    "utils",
    "uuid",
    "webshot",
    "webshot2",
    "websocket",
    "withr",
    "xfun",
    "xml2",
    "xopen",
    "xtable",
    "yaml",
    "zip",
    "zoo"
)


biocPkgs <-
    c(
        "lumi",
        "methylumi",
        "conumee",
        "minfi",
        "IlluminaHumanMethylation450kmanifest",
        "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
    )

if (checkRequire("BiocManager")) {
    install.packages("BiocManager", dependencies = T, verbose = T, ask = F)
}

loadLibrary("devtools")
loadLibrary("librarian")
loadLibrary("BiocManager")
if (checkRequire("Biobase")) {
    BiocManager::install("Biobase", update = F, ask = F, type = "binary")
}
loadLibrary("Biobase")

if (checkRequire("mapview")) {
     tryCatch(
        expr = remotes::install_github("r-spatial/mapview", dependencies = T, upgrade = "never"),
        error = install.packages("mapview", dependencies = T, verbose = T, ask = F, type = "binary")
    )
}

# Load/install missing pacakges without asking ----------------------------------------------------
supM <- function(pk) {return(suppressPackageStartupMessages(suppressWarnings(pk)))}

message("Loading packages:", paste0(capture.output(corePkgs), collapse = "\n"))

librarian::shelf(
    corePkgs,
    ask = F,
    update_all = F,
    quiet = FALSE
)

message("Loading packages:", paste0(capture.output(preReqPkgs), collapse = "\n"))

librarian::shelf(
    preReqPkgs,
    ask = F,
    update_all = F,
    quiet = FALSE
)

message("Loading BioConductor Packages and IlluminaHumanMethylation Manifest...")

if (checkRequire("IlluminaHumanMethylationEPICmanifest")) {
    devtools::install_github(repo = "mwsill/IlluminaHumanMethylationEPICmanifest",
                             dependencies = T,
                             upgrade = "never")
}

if (checkRequire("FDb.InfiniumMethylation.hg19")) {
    BiocManager::install(
        "FDb.InfiniumMethylation.hg19",
        update = F,
        ask = F,
        dependencies = T
    )
}

supM(librarian::shelf(
    biocPkgs,
    ask = F,
    update_all = F,
    quiet = F,
    dependencies = T
))

if (checkRequire("mgmtstp27")) {
    gitLink <-
        "https://github.com/badozor/mgmtstp27/raw/master/archive/mgmtstp27_0.6-3.tar.gz"
    install.packages(
        gitLink,
        repos = NULL,
        dependencies = T,
        verbose = T,
        type = "source",
        ask = F
    )
}

fixProf <- function() {
    txt1 <- "^[:blank:]*autoload\\(\"needs\", \"needs\"\\)"
    txt2 <- "\n\nautoload(\"needs\", \"needs\")\n\n"
    siteProf <- if (is.na(Sys.getenv("R_PROFILE", unset = NA))) {
        file.path(Sys.getenv("R_HOME"), "etc", "Rprofile.site")
    } else {
        Sys.getenv("R_PROFILE")
    }
    if (!file.exists(siteProf)) {
        try(file.create(siteProf), silent = T)
    }
    cxn <- try(file(siteProf), silent = T)
    lines <- try(base::readLines(cxn), silent = T)
    if (!any(grepl(txt1, lines))) {
        write(txt2, file = siteProf, append = T)
    }
    close(cxn)
    closeAllConnections()
}

fixNeeds <- function () {
    sysfile <- system.file("extdata", "promptUser", package = "needs")
    try(write(0, file = sysfile), silent = T)
    options(needs.promptUser = FALSE)
    invisible(needs:::autoload(TRUE))
    closeAllConnections()
}

checkNeeds <- function() {
    tryCatch(
        expr = {
            if (!("needs" %in% rownames(installed.packages()))) {
                install.packages(
                    "needs",
                    dependencies = T,
                    verbose = T,
                    Ncpus = 2
                )
                fixNeeds()
                try(fixProf(), silent = T)
            } else{
                fixNeeds()
                try(fixProf(), silent = T)
            }
        },
        error = function(cond) {
            devtools::install_github(
                "joshkatz/needs",
                ref = "development",
                dependencies = T,
                verbose = T,
                upgrade = "never"
            )
            fixNeeds()
            try(fixProf(), silent = T)
        },
        warning = function(cond) {
            devtools::install_github(
                "joshkatz/needs",
                ref = "development",
                dependencies = T,
                verbose = T,
                upgrade = "never"
            )
            fixNeeds()
            try(fixProf(), silent = T)
        }
    )
}

if (Sys.info()[['sysname']] == "Darwin") {
    checkNeeds()
    closeAllConnections()
    isGdal <- paste(system("echo `gdalinfo --version`", intern = T))
    if (isGdal == "") {
        system("brew install pkg-config")
        system("brew install gdal")
    }
    isProj <- system("which proj", intern = T)
    if (length(isProj) == 0) {
        system("brew install pkg-config")
        system("brew install proj")
    }
} else{
    if (!("needs" %in% rownames(installed.packages()))) {
        install.packages(
            "needs",
            dependencies = T,
            verbose = T,
            ask = F
        )
    }
    options(needs.promptUser = FALSE)
}

if (!(require("Rcpp"))) {
    remotes::install_github("RcppCore/Rcpp")
}
spat_config <-
    '--with-proj-lib=/usr/local/lib/ --with-proj-include=/usr/local/include/'
options(configure.args = c("sf" = spat_config))

if (checkRequire("sf")) {
    tryCatch(
        install.packages(
            c("sf"),
            type = "source",
            dependencies = T,
            verbose = T,
            Ncpus = 4
        ),
        error = function(e) {
            remotes::install_github(
                "r-spatial/sf",
                configure.args = "--with-proj-lib=/usr/local/lib/",
                dependencies = T,
                upgrade = "never"
            )
        }
    )
}

invisible(gc())

if (Sys.info()[['sysname']] == "Darwin") {
    system("export LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib")
}

options("install.packages.compile.from.source" = "No")
options("install.packages.check.source" = "no")
loadLibrary("BiocManager")
loadLibrary("Biobase")

terraDep <- c('tinytest', 'ncdf4', 'leaflet')
suppressWarnings(librarian::shelf(
    terraDep,
    ask = F,
    update_all = F,
    quiet = FALSE
))

if (checkRequire("terra")) {
    install.packages(
        'terra',
        repos = 'https://rspatial.r-universe.dev',
        dependencies = T,
        verbose = T,
        Ncpus = 4
    )
}

if (checkRequire("FField")) {
    gitLink <-
        "https://cran.r-project.org/src/contrib/Archive/FField/FField_0.1.0.tar.gz"
    install.packages(
        gitLink,
        repos = NULL,
        dependencies = T,
        verbose = T,
        type = "source",
        ask = F
    )
}

if (checkRequire("GenVisR")) {
    devtools::install_github("griffithlab/GenVisR", dependencies = T, upgrade = "never")
}

suppressWarnings(librarian::shelf(pkgs, ask = F, update_all = F, quiet = FALSE))
invisible(gc())

try(install.packages('foghorn', dependencies = T, type = 'binary'), silent = T)

cbioLn <- switch(Sys.info()[['sysname']],
                 "Darwin" = "/Volumes/CBioinformatics/Methylation/classifiers",
                 "Linux" = "~/molecpathlab/production/Methylation/classifiers")

try(librarian::shelf(c("mlr", "wateRmelon", "RPMM", "impute"),
                     ask = F, update_all = F, quiet = F),
    silent = T)

if(checkRequire("wateRmelon")){
    BiocManager::install("wateRmelon", dependencies = T, type = "source", update = F)
}

if (checkRequire("UniD")) {
    unidPath <- file.path(cbioLn, "UniD")
    try(install.packages(unidPath, type = "source", repos = NULL), silent = T)
}

closeAllConnections()

githubMain <- "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R"
devtools::source_url(file.path(githubMain, "LoadInstallPackages.R"))
