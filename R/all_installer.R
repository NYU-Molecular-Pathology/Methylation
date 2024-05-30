#!/usr/bin/env Rscript
## ---------------------------
## Script name: all_installer.R
## Purpose: Functions to check if all required packages for the pipeline are installed
## Date Created: August 9, 2022
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

options("install.packages.compile.from.source" = "Yes")
options("install.packages.check.source" = "no")
options(Ncpus = 4)

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
    # Check if GDAL installed
    isGdal <- paste(system("echo `gdalinfo --version`", intern = T))
    if (isGdal == "") {
        brew_install("pkg-config")
        brew_install("gdal")
    }
    # Check if proj installed
    isProj <- system("which proj", intern = T)
    if (length(isProj) == 0) {
        brew_install("pkg-config")
        brew_install("proj")
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
set_gfortran <- function() {
    gfortran_path <- Sys.which("gfortran")[[1]]
    Sys.setenv(FC = gfortran_path)
    Sys.setenv(FLIBS = "-L/usr/local/lib/gcc/current")

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
        message(sprintf("Updated %s to include /usr/local/sbin in PATH.", shell_config_file))
    } else {
        warning("Failed to update the PATH in the shell configuration file.")
    }
    Sys.setenv(PATH = paste("/usr/local/bin", Sys.getenv("PATH"), sep = ":"))
    Sys.setenv(LIBRARY_PATH = paste("/usr/local/lib", Sys.getenv("LIBRARY_PATH"), sep = ":"))
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
    Sys.setenv(LDFLAGS = "-L/usr/local/opt/llvm/lib -L/usr/local/opt/llvm/lib/c++ -Wl,-rpath,/usr/local/opt/llvm/lib/c++")
    Sys.setenv(CPPFLAGS = "-I/usr/local/opt/llvm/include")
    Sys.setenv(PKG_CFLAGS = "-I/usr/local/include -I/usr/local/opt/llvm/include -I/usr/local/opt/apache-arrow/include")
    Sys.setenv(PKG_LIBS = "-L/usr/local/lib -L/usr/local/opt/llvm/lib -L/usr/local/opt/apache-arrow/lib -larrow")
    Sys.setenv(LD_LIBRARY_PATH = "/usr/local/lib")
    Sys.setenv(R_LD_LIBRARY_PATH = "/usr/local/lib:/usr/local/opt/llvm/lib/c++")

    set_openmpi()
    set_gfortran()
    system("brew cleanup")
    update_system_path()

    curr <- paste("/usr/local/opt/llvm/bin", Sys.getenv("PATH"), sep = ":")
    Sys.setenv(PATH = paste("/usr/local/opt/open-mpi/bin", curr, sep = ":"))
}


set_compiler_paths <- function() {
    sdk_path <- system("xcrun --sdk macosx --show-sdk-path", intern = TRUE)
    cc_current <- Sys.getenv("CC")
    cc_flags <- paste("clang -isysroot", sdk_path)
    cxx_curr <- Sys.getenv("CXX")
    cxx_new <- paste("clang++ -isysroot", sdk_path)

    Sys.setenv(CC = paste(cc_current, cc_flags))
    Sys.setenv(CXX = paste(cxx_new, cxx_curr))
    Sys.setenv(CXXFLAGS = paste("-isysroot", sdk_path))
}

update_makevars <- function() {
    makevars_path <- file.path(Sys.getenv("HOME"), ".R", "Makevars")
    compiler_settings <- c(
        "CC = /usr/local/opt/llvm/bin/clang",
        "CXX = /usr/local/opt/llvm/bin/clang++",
        "CXX11 = /usr/local/opt/llvm/bin/clang++",
        "CXX14 = /usr/local/opt/llvm/bin/clang++",
        "CXX17 = /usr/local/opt/llvm/bin/clang++",
        "CXX1X = /usr/local/opt/llvm/bin/clang++",
        "OBJC = /usr/local/opt/llvm/bin/clang",
        "LDFLAGS= -L/usr/local/opt/llvm/lib -L/usr/local/opt/llvm/lib/c++ -Wl,-rpath,/usr/local/opt/llvm/lib/c++",
        "CPPFLAGS = -I/usr/local/opt/llvm/include",
        "PKG_CFLAGS = -I/usr/local/include -I/usr/local/opt/llvm/include -I/usr/local/opt/apache-arrow/include",
        "PKG_LIBS = -L/usr/local/lib -L/usr/local/opt/llvm/lib -L/usr/local/opt/apache-arrow/lib -larrow",
        "FC = /usr/local/gfortran/bin/gfortran",
        "FLIBS = -L/usr/local/lib/gcc/current",
        "LD_LIBRARY_PATH = /usr/local/lib",
        "R_LD_LIBRARY_PATH = /usr/local/lib:/usr/local/opt/llvm/lib/c++"
    )
    dir.create(dirname(makevars_path), showWarnings = F, recursive = T)
    writeLines(compiler_settings, makevars_path)
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
    #closeAllConnections()
}

if (is_macos) {
    options(BioC_mirror = "https://packagemanager.rstudio.com/bioconductor")
    #options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2024-02-20"))
    options(warn = -1)
    options(repos = c(CRAN = 'https://cloud.r-project.org'))
    # Setting US CRAN REPO
    rlis = getOption("repos")
    rlis["CRAN"] = "http://cran.us.r-project.org"
    options(repos = rlis)
    fix_compiler_flags()
    update_makevars()
    system("export LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib")
}


# Check if architecture is 'arm64' or 'x86_64' -------------------------------------
arch <- Sys.info()[["machine"]]

if (arch != "x86_64" & is_macos == T) {
    # Set JAVA_HOME environment variable
    #java_home <- system("/usr/libexec/java_home -v 11", intern = TRUE)
    system("brew tap homebrew/cask-versions && brew install --cask temurin17")
    java_home <- system("which java", intern = TRUE)
    Sys.setenv(JAVA_HOME = java_home)
    message("JAVA_HOME set to ", java_home)
    try(install.packages("rJava", type = "binary", dependencies = T, ask = F), T)
}

loadLibrary <- function(pkgName) {
    suppressPackageStartupMessages(library(
        pkgName,
        quietly = TRUE,
        logical.return = TRUE,
        warn.conflicts = FALSE,
        character.only = TRUE
    ))
}

checkPkg <- function(pkgName) {
    return(!requireNamespace(pkgName, quietly = TRUE))
}

if (checkPkg("devtools")) {
    install.packages("devtools", dependencies = T, verbose = T, ask = F)
}
if (checkPkg("librarian")) {
    install.packages("librarian", dependencies = T, verbose = T, ask = F)
}
if (checkPkg("BiocManager")) {
    install.packages("BiocManager", dependencies = T, verbose = T, ask = F)
}
if (checkPkg("Biobase")) {
    BiocManager::install("Biobase", update = F, ask = F, type = "binary")
}
if (loadLibrary("Biobase") == F) {
    BiocManager::install("Biobase", update = F, ask = F, type = "source")
}

stopifnot(loadLibrary("devtools"))
stopifnot(loadLibrary("librarian"))
stopifnot(loadLibrary("BiocManager"))
stopifnot(loadLibrary("Biobase"))

# FUN: Downloads Github repo locally then unzips it ---------------------------
download_pkg_unzip <- function(git_repo, zip_name = "main.zip") {
    repo_url <- file.path("https://github.com", git_repo, "archive/refs/heads")
    url_path <- file.path(repo_url, zip_name)
    local_dir <- file.path(fs::path_home(), "github_pkgs")
    out_file <- file.path(local_dir, paste0(basename(git_repo), "-", zip_name))
    utils::download.file(url = url_path, method = "curl", destfile = out_file)
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


# FUN: TryCatch for installing with devtools and install.packages -------------
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


fixNeeds <- function() {
    sysfile <- system.file("extdata", "promptUser", package = "needs")
    try(write(0, file = sysfile), silent = T)
    try(close(sysfile), silent = T)
    options(needs.promptUser = FALSE)
    invisible(needs:::autoload(TRUE))
    #closeAllConnections()
}

checkNeeds <- function() {
    tryCatch(
        expr = {
            if (!("needs" %in% rownames(installed.packages()))) {
                install.packages("needs", dependencies = T, verbose = T, Ncpus = 2)
                fixNeeds()
                try(fixProf(), silent = T)
            } else{
                fixNeeds()
                try(fixProf(), silent = T)
            }
        },
        error = function(cond) {
            try_github_inst("joshkatz/needs")
            fixNeeds()
            try(fixProf(), silent = T)
        },
        warning = function(cond) {
            try_github_inst("joshkatz/needs")
            fixNeeds()
            try(fixProf(), silent = T)
        }
    )
}


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
    "umap",
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


biocPkgs <- c(
    "lumi",
    "methylumi",
    "conumee",
    "minfi",
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
)

# Load/install missing pacakges without asking ----------------------------------------------------
supM <- function(pk) {return(suppressPackageStartupMessages(suppressWarnings(pk)))}

if (checkPkg("mapview")) {
    tryCatch(
        expr = remotes::install_github("r-spatial/mapview", dependencies = T, upgrade = "never"),
        error = install.packages("mapview", dependencies = T, verbose = T, ask = F, type = "binary")
    )
}

load_install <- function(pkg_list) {
    message("Loading packages:\n", paste0(capture.output(pkg_list), collapse = "\n"))
    librarian::shelf(pkg_list, ask = F, update_all = F, quiet = F, dependencies = T)
    missing_pkgs <- pkg_list[!(pkg_list %in% rownames(installed.packages()))]
    if (length(missing_pkgs) > 0) {
        for (pkg_missed in missing_pkgs) {
            tryCatch(
                expr = {
                    BiocManager::install(pkg_missed, update = F, ask = F, dependencies = T, type = "source")
                },
                error = {
                    install.packages(pkg_missed, dependencies = T, verbose = T, type = "source", ask = F)
                }
            )
        }
    }
}


load_install(corePkgs)
if (!requireNamespace("ff", quietly = TRUE)) {
    install.packages("ff", type = "binary", ask = F, dependencies = T)
}
load_install(preReqPkgs)

message("Loading BioConductor Packages and IlluminaHumanMethylation Manifest...")
if (checkPkg("IlluminaHumanMethylationEPICmanifest")) {
    tryCatch(
        expr = try_github_inst("mwsill/IlluminaHumanMethylationEPICmanifest"),
        error = function(e){
            message(e)
            message("You need to set your Git token to install Github packages")
            message("usethis::create_github_token()")
            message("Sys.setenv(GITHUB_PAT = 'YOUR_TOKEN')")
        }
    )
}

if (checkPkg("FDb.InfiniumMethylation.hg19")) {
    BiocManager::install(
        "FDb.InfiniumMethylation.hg19", update = F, ask = F, dependencies = T
    )
}

load_install(biocPkgs)

if (checkPkg("mgmtstp27")) {
    gitLink <-
        "https://github.com/badozor/mgmtstp27/raw/master/archive/mgmtstp27_0.6-3.tar.gz"
    install.packages(
        gitLink, repos = NULL, dependencies = T, verbose = T, type = "source", ask = F
    )
}
if (is_macos) {
    checkNeeds()
    #closeAllConnections()
} else{
    if (checkPkg("needs")) {
        install.packages("needs", dependencies = T, verbose = T, ask = F)
    }
}

try(options(needs.promptUser = FALSE), T)

if (checkPkg("Rcpp")) {
    try_github_inst("RcppCore/Rcpp")
}

spat_config <-
    '--with-proj-lib=/usr/local/lib/ --with-proj-include=/usr/local/include/'
options(configure.args = c("sf" = spat_config))

if (checkPkg("sf")) {
    tryCatch(
        install.packages("sf", type = "source", dependencies = T, verbose = T, Ncpus = 4),
        error = function(e) {
            remotes::install_github("r-spatial/sf", configure.args = spat_config,
                                    dependencies = T, upgrade = "never")})
}

invisible(gc())

options("install.packages.compile.from.source" = "No")
options("install.packages.check.source" = "no")

loadLibrary("BiocManager")
loadLibrary("Biobase")

terraDep <- c('tinytest', 'ncdf4', 'leaflet')
suppressWarnings(librarian::shelf(terraDep, ask = F, update_all = F, quiet = F))

if (checkPkg("terra")) {
    install.packages(
        'terra', repos = 'https://rspatial.r-universe.dev', ask = F,
        dependencies = T, verbose = T, Ncpus = 4)
}

if (checkPkg("FField")) {
    gitLink <-
        "https://cran.r-project.org/src/contrib/Archive/FField/FField_0.1.0.tar.gz"
    install.packages(
        gitLink, repos = NULL, dependencies = T, verbose = T, type = "source", ask = F
    )
}

if (checkPkg("GenVisR")) {
    try_github_inst("griffithlab/GenVisR")
}

load_install(pkgs)

invisible(gc())

try(install.packages('foghorn', dependencies = T, type = 'binary'), silent = T)

try(librarian::shelf(
    c("mlr", "wateRmelon", "RPMM", "impute"),
    ask = F, update_all = F, quiet = F), silent = T)

if (checkPkg("wateRmelon")) {
    BiocManager::install("wateRmelon", dependencies = T, type = "source", update = F)
}

class_pkg_path <- ifelse(
    is_macos,
    "/Volumes/CBioinformatics/Methylation/classifiers",
    "/gpfs/data/molecpathlab/production/Methylation/classifiers"
)

if (checkPkg("UniD")) {
    unidPath <- file.path(class_pkg_path, "UniD")
    try(install.packages(unidPath, type = "source", repos = NULL), silent = T)
}

githubMain <- "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R"
devtools::source_url(file.path(githubMain, "LoadInstallPackages.R"))
