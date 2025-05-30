#!/usr/bin/env Rscript
## Script name: all_installer.R
## Purpose: Functions to check if required packages are installed
## Date Created: August 9, 2022
## Version: 1.1.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025

github_main <-
    "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main"
arch <- Sys.info()[["machine"]]
is_macos <- Sys.info()[['sysname']] == "Darwin"

if (getRversion() <= "4.2.2") {
    stop("Your R version is ", R.version.string, ". Update to 4.4.0 or later.")
}

snapshot_date <- "2025-05-01"
bioc_version   <- "3.19"
options(
    repos = c(
        CRAN = sprintf("https://cran.microsoft.com/snapshot/%s/", snapshot_date),
        BioC = sprintf("https://packagemanager.rstudio.com/bioconductor/%s@%s",
                       bioc_version, snapshot_date)
        )
)

supM <- function(pk) {
    return(suppressPackageStartupMessages(suppressWarnings(pk)))
}

mac_tools <- c(
    'xz', 'tiff', 'libpng', 'openssl', 'jpeg', 'pcre2', 'cairo', 'texinfo',
    'libdeflate', 'zstd','apr', 'apr-util', 'pkgconfig', 'autoconf', 'm4',
    'automake', 'bdb', 'boost', 'libb2', 'lz4','freetype', 'fontconfig',
    'pixman', 'cgl', 'clp', 'osi', 'coinutils', 'dylp', 'emacs', 'expat',
    'ffi', 'fftw',  'flint', 'mpfr', 'gmp', 'freexl', 'fribidi', 'gdal',
    'proj', 'openjpeg', 'libwebp', 'libgeotiff', 'sqlite3', 'hdf4', 'hdf5',
    'szip', 'netcdf', 'geos', 'unixodbc', 'gettext', 'glib', 'pcre', 'glpk',
    'gsl', 'harfbuzz', 'icu', 'hwloc', 'isl', 'libarchive', 'libgit2', 'libssh2',
    'libsbml', 'libsndfile', 'mpc', 'nlopt', 'pango', 'poppler', 'poppler-data',
    'protobuf', 'qpdf', 'QuantLib', 'readline5', 'rsync', 'serf', 'subversion',
    'utf8proc', 'symphony', 'sys-stubs', 'texinfo', 'tidy', 'udunits', 'zeromq'
)

# Function to check for installed macOS binaries and install missing ones
check_and_install_mac_libs <- function(pkgs) {

    base_dir <- switch(Sys.info()[["machine"]],
                       "arm64" = "/opt/R/arm64",
                       "x86_64" = "/opt/R/x86_64",
                       "/usr/local")
    receipt_dir <- file.path(base_dir, "pkg")

    if (!dir.exists(receipt_dir)) {
        warning("Receipt directory not found: ", receipt_dir,
                "\nAssuming none of the binaries are installed.")
        missing_pkgs <- pkgs
    } else {
        installed_files <- dir(receipt_dir, all.files = TRUE)
        missing_pkgs <- pkgs[!sapply(pkgs, function(pkg) {
            any(grepl(paste0("^", pkg, "-.*\\.list$"), installed_files))
        })]
    }

    if (length(missing_pkgs) > 0) {
        message("Installing missing macOS binaries...")
        source("https://mac.R-project.org/bin/install.R")
        install.libs("all")
    }
}


fix_brew_path <- function() {
    # Define the target directories to add to the PATH
    target_dirs <- c("/opt/homebrew/bin", "/usr/local/bin")

    # Get the current PATH in the R session
    current_path <- strsplit(Sys.getenv("PATH"), ":")[[1]]

    # Add the target directories to the PATH if they are not already present
    new_path <- unique(c(target_dirs, current_path))
    Sys.setenv(PATH = paste(new_path, collapse = ":"))
    message("PATH updated for the current session.")

    # Persist the change in ~/.Renviron for future sessions
    renviron_file <- file.path(Sys.getenv("HOME"), ".Renviron")

    # Ensure the target directories are reflected in ~/.Renviron
    target_entry <- paste0('PATH="', paste(new_path, collapse = ":"), '"')
    if (file.exists(renviron_file)) {
        renviron_content <- readLines(renviron_file, warn = FALSE)
    } else {
        renviron_content <- character(0)
    }
    if (!any(grepl("^PATH=", renviron_content))) {
        # Append the new PATH if no PATH is defined
        writeLines(c(renviron_content, target_entry), renviron_file)
        message("PATH added to ~/.Renviron.")
    } else {
        # Update the PATH entry if it exists
        renviron_content <- sub("^PATH=.*", target_entry, renviron_content)
        writeLines(renviron_content, renviron_file)
        message("PATH updated in ~/.Renviron.")
    }
}


command_exists <- function(cmd) {
    path <- Sys.which(cmd)
    if (nzchar(path)) {
        return(TRUE)
    } else {
        warning(sprintf("Command '%s' not found in system PATH.", cmd), call. = FALSE)
        return(FALSE)
    }
}

# Ensures brew was added to the user shell profile
setup_brew <- function(brew_path) {
    zprof <- file.path(Sys.getenv("HOME"), ".zprofile")
    bashprof <- file.path(Sys.getenv("HOME"), ".bash_profile")
    commands <- c(
        sprintf("echo >> %s", shQuote(zprof)),
        sprintf("echo 'eval \"$(%s shellenv)\"' >> %s", brew_path, shQuote(zprof)),
        sprintf("echo >> %s", shQuote(bashprof)),
        sprintf("echo 'eval \"$(%s shellenv)\"' >> %s", brew_path, shQuote(bashprof)),
        sprintf("eval \"$(%s shellenv)\"", brew_path)
    )

    for (cmd in commands) system(cmd)
}

# Install Homebrew and packages if necessary
ensure_homebrew <- function() {
    pkgs <- c("gcc", "llvm", "lld", "open-mpi", "pkgconf", "gdal", "proj",
              "apache-arrow")
    arm_brew <- "/opt/homebrew/bin/brew"
    x64_brew <- "/usr/local/bin/brew"
    brew_installed <- file.exists(arm_brew) | file.exists(x64_brew)
    if (!brew_installed) {
        message("Installing Homebrew...")
        system('/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"', wait = TRUE)
        brew_path  <- ifelse(file.exists(arm_brew), arm_brew, x64_brew)
        setup_brew(brew_path)
    }

    fix_brew_path()

    installed_pkgs <- system2("brew", c("list", "--formula"),
                              stdout = T, stderr = NULL)

    for (pkg in pkgs) {
        if (!(pkg %in% installed_pkgs)) {
            message("Installing ", pkg, " via Homebrew...")
            system2("brew", c("install", pkg), wait = TRUE)
        }
    }
}

get_prefix <- function(pkg = "") {
    system2("brew", c("--prefix", pkg), stdout = TRUE, stderr = FALSE)
}

# Set environment variables dynamically
set_env_vars <- function() {
    brew_prefix = get_prefix()
    llvm_path = get_prefix("llvm")
    mpi_path = get_prefix("open-mpi")
    arrow_path = get_prefix("apache-arrow")

    # Dynamically find any versioned libarrow.*.dylib file
    libarrow_files <- Sys.glob(file.path(arrow_path, "lib", "libarrow.*.dylib"))
    if (length(libarrow_files) == 0) stop("No versioned libarrow.*.dylib file found in apache-arrow lib path.")

    # Select the latest version of libarrow.*.dylib (lexicographically last)
    libarrow_versioned <- sort(libarrow_files, decreasing = TRUE)[1]


    Sys.setenv(
        CC = file.path(llvm_path, "bin/clang"),
        CXX = file.path(llvm_path, "bin/clang++"),
        OBJC = file.path(llvm_path, "bin/clang"),
        LDFLAGS = paste(
            paste0("-L", file.path(llvm_path, "lib")),
            paste0("-L", file.path(llvm_path, "lib/c++")),
            paste0("-Wl,-rpath,", file.path(llvm_path, "lib/c++")),
            paste0("-L", file.path(llvm_path, "lib/unwind"), "-lunwind")
        ),
        CPPFLAGS = paste0("-I", file.path(llvm_path, "include")),
        PKG_CFLAGS = paste(
            paste0("-I", file.path(brew_prefix, "include")),
            if (nzchar(arrow_path))
                paste0("-I", file.path(arrow_path, "include"))
            else
                NULL
        ),
        PKG_LIBS = paste(
            paste0("-L", file.path(brew_prefix, "lib")),
            paste0("-L", file.path(llvm_path, "lib")),
            if (nzchar(arrow_path)) {
                paste0("-L", file.path(arrow_path, "lib"), " -larrow")
            }else { NULL }
        ),
        LD_LIBRARY_PATH = file.path(brew_prefix, "lib"),
        R_LD_LIBRARY_PATH = paste(
            file.path(brew_prefix, "lib"),
            file.path(llvm_path, "lib/c++"),
            sep = ":"
        ),
        DYLD_LIBRARY_PATH = paste(
            file.path(arrow_path, "lib"),
            libarrow_versioned,
            sep = ":"
        ),
        PATH = paste(file.path(mpi_path, "bin"), Sys.getenv("PATH"), sep = ":")
    )
    if (command_exists("gfortran")) Sys.setenv(FC = Sys.which("gfortran"))
}

# Ensure any required libraries are symlinked
get_brew_prefix <- function(pkg) {
    prefix <- system2("brew", c("--prefix", pkg), stdout = TRUE, stderr = FALSE)
    if (length(prefix) == 0 || !dir.exists(prefix)) {
        stop("Could not find Homebrew prefix for package: ", pkg)
    }
    return(prefix)
}

# Check for existing symlink and create if missing
ensure_symlink <- function(source, target) {
    if (!file.exists(source)) {
        stop("Source library not found: ", source)
    }
    if (!file.exists(target)) {
        message("Creating symlink: ", basename(target), " -> ", basename(source))
        owd <- setwd(dirname(source))
        on.exit(setwd(owd), add = TRUE)
        file.symlink(basename(source), basename(target))
    } else {
        message("Symlink or file already exists: ", target)
    }
}

# Locate current libarrow dylib and symlink it to the expected name
symlink_arrow <- function(target_name = "libarrow.1801.dylib") {
    arrow_prefix <- get_brew_prefix("apache-arrow")
    lib_dir      <- file.path(arrow_prefix, "lib")
    candidates   <- list.files(lib_dir,
                               pattern = "^libarrow\\.[0-9]+\\.dylib$",
                               full.names = TRUE)
    if (length(candidates) == 0) {
        stop("No libarrow.*.dylib files found in ", lib_dir)
    }
    versions   <- as.integer(gsub("^libarrow\\.([0-9]+)\\.dylib$", "\\1",
                                  basename(candidates)))
    # find the latest
    idx_max <- which.max(versions)
    source_lib <- candidates[idx_max]
    # all older versions get a symlink pointing at the latest
    for (i in seq_along(candidates)) {
        if (i == idx_max) next
        ver <- versions[i]
        target_lib <- file.path(lib_dir, sprintf("libarrow.%d.dylib", ver))
        ensure_symlink(source_lib, target_lib)
    }
}

# Checks if the symlink exists between a library and target path
check_symlink <- function(library_path, target_path = "libarrow.1700.dylib") {
    # Check if the library exists at the expected location
    if (!file.exists(library_path)) {
        stop("Library not found at\n", library_path)
    }

    # Check if the symlink exists and create it if it does not
    if (!file.exists(target_path)) {
        message("Creating symlink for ", basename(library_path), "...")
        system(paste("cd", dirname(library_path)))
        system(paste("ln -s", library_path, target_path))
    }
}


symlink_llvm <- function() {
    llvm_path = get_prefix("llvm")
    unwind_libpath <- file.path(llvm_path, "lib", "unwind", "libunwind.1.dylib")
    missing_path <- file.path(llvm_path, "lib", "libunwind.1.dylib")
    check_symlink(unwind_libpath, missing_path)

    arrow_path = get_prefix("apache-arrow")
    arrow_libs <- dir(file.path(arrow_path, "lib"), all.files = TRUE, full.names = TRUE)
    arrow_dylib <- arrow_libs[grepl("^libarrow\\.[0-9]{4}\\.dylib$", basename(arrow_libs))]
    miss_arrow <- file.path(arrow_path, "lib", "libarrow.1700.dylib")
    #check_symlink(arrow_dylib, miss_arrow)
}


# Main setup
setup_compilers <- function() {
    ensure_homebrew()
    comp_flags <- c("CC", "CXX", "OBJC", "LDFLAGS", "CPPFLAGS", "PKG_CFLAGS",
                    "PKG_LIBS", "LD_LIBRARY_PATH", "R_LD_LIBRARY_PATH")
    Sys.unsetenv(comp_flags)
    set_env_vars()
    options(Ncpus = 6)
    symlink_llvm()
}


loadLibrary <- function(pkgName) {
    suppressWarnings(suppressPackageStartupMessages(library(
        pkgName, verbose = FALSE, quietly = TRUE,
        logical.return = TRUE, warn.conflicts = FALSE, character.only = TRUE
    )))
}

# FUNC: Checks if a package is installed ---------------
not_installed <- function(pkgName) {
    return(!pkgName %in% rownames(installed.packages()))
}

# FUNC: Installs a package from binary file -----------------------------------
binary_install <- function(pkg) {
    if (not_installed(pkg)) {
        tryCatch(
            install.packages(pkg, dependencies = T, ask = F, type = "binary"),
            error = function(e) {
                message(e)
                message("trying to install as source")
                install.packages(pkg, dependencies = T, ask = F)
            }
        )
    }
}


# FUNC: Downloads Github repo locally then unzips it --------------------------
download_pkg_unzip <- function(git_repo, zip_name = "main.zip") {
    repo_url <- file.path("https://github.com", git_repo, "archive/refs/heads")
    url_path <- file.path(repo_url, zip_name)
    local_dir <- file.path(fs::path_home(), "github_pkgs")
    out_file <- file.path(local_dir, paste0(basename(git_repo), "-", zip_name))
    utils::download.file(url = url_path, method = "curl", destfile = out_file)
    utils::unzip(out_file, exdir = local_dir)
}

# FUNC: Downloads Github repo locally then installs ---------------------------
local_github_pkg_install <- function(git_repo) {
    repo_url <- file.path("https://github.com", git_repo, "archive/refs/heads")
    local_dir <- file.path(fs::path_home(), "github_pkgs")

    if (!dir.exists(local_dir)) dir.create(local_dir)
    tryCatch(
        expr = {
            download_pkg_unzip(git_repo, zip_name = "main.zip")
        },
        error = function(e) {
            download_pkg_unzip(git_repo, zip_name = "master.zip")
        }
    )

    unzipped_dir <- list.dirs(local_dir, full.names = TRUE, recursive = FALSE)
    pkg_dir <- unzipped_dir[grepl(basename(git_repo), unzipped_dir)]

    install.packages(pkg_dir, repo = NULL, type = "source", dependencies = T)
}


# FUNC: Installs package from a Github repository -----------------------------
install_repo <- function(git_repo) {
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


# FUNC: TryCatch for installing with devtools and install.packages ------------
try_github_inst <- function(git_repo) {
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


install_opts <- list(
    dependencies = c("Depends", "Imports", "LinkingTo"),
    ask = FALSE,
    update = "never",
    quiet = TRUE,
    repos = 'http://cran.us.r-project.org',
    Ncpus = 6
)


# FUNC: Quietly loads package library without messages ------------------------
quiet_load <- function(pkg_name) {
    if (isNamespaceLoaded(pkg_name) == F) {
        libLoad <- suppressWarnings(suppressPackageStartupMessages(
            library(pkg_name, character.only = T, logical.return = T, quietly = T)
        ))
    } else{
        libLoad <- TRUE
    }
    message(pkg_name, " loaded... ", libLoad)
    pkg_vec <- c(pkg_n = libLoad)
    names(pkg_vec) <- pkg_name
    return(pkg_vec)
}

# FUNC: Checks required package if not installs binary ------------------------
require_pkg <- function(pkg, pkg_type = "source") {
    install_opts$type <- pkg_type
    if (!pkg %in% rownames(installed.packages())) {
        do.call(install.packages, c(list(pkgs = pkg), install_opts))
    }
    quiet_load(pkg)
}


# Check devtools and install package dependencies -----------------------------
require_pkg("devtools", "binary")
require_pkg("utils", "binary")
require_pkg("BiocManager", "binary")

biocRepos <- suppressMessages(BiocManager::repositories())
avail_bioc_packs <- suppressMessages(BiocManager::available())
rbase_pkgs <- rownames(installed.packages(priority = "base"))
pkg_info <- utils::available.packages(repos = biocRepos)

if (not_installed("pak")) {
    install.packages("pak", dependencies = T, ask = F)
}
stopifnot(loadLibrary("pak"))

# FUNC: Returns all packages that are not installed ----------------------------
check_needed <- function(pkgs) {
    neededPkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
    return(neededPkgs)
}

# FUNC: Returns all package dependencies that are not installed ---------------
get_pkg_deps <- function(pkgs) {
    deps_list <- unique(pak::pkg_deps(pkgs)$package)
    all_pkgs <- setdiff(deps_list, pkgs)
    return(check_needed(all_pkgs))
}

# FUNC: Gets and installs any package dependencies ----------------------------
install_deps <- function(pkg) {
    any_deps <- get_pkg_deps(pkg)
    if (length(any_deps) > 0) {
        pak::pkg_install(any_deps, ask = F)
    }
}

# FUNC: Attempts source and binary package installation -----------------------
install_pkgs <- function(needed_pk, pkg_type = "source") {
    install_opts$type <- pkg_type
    tryCatch(
        do.call(install.packages, c(list(pkgs = needed_pk), install_opts)),
        error = function(e) {
            message("\nInitial installation failed! Trying binary install.\n")
            install_opts$type <- "binary"
            do.call(install.packages, c(list(pkgs = needed_pk), install_opts))
        }
    )
}

# FUNC: Attempts source and binary install of Bioconductor --------------------
install_bio_pkg <- function(new_pkg) {
    params <- list(dependencies = c("Depends", "Imports", "LinkingTo"),
                   ask = FALSE, update = TRUE)
    pkg_deps <- get_pkg_deps(new_pkg)
    if (length(pkg_deps) > 0) {
        new_pkg <- c(pkg_deps, new_pkg)
    }
    for (pkg in new_pkg) {
        message("Installing BioCpackage dependency:")
        tryCatch(
            do.call(BiocManager::install, c(list(pkgs = pkg), params)),
            error = function(e) {
                params$type <- "binary"
                do.call(BiocManager::install, c(list(pkgs = pkg), params))
            }
        )
    }
}

# FUNC: Installs a package from a tarball url ---------------------------------
install_url <- function(pkg_url){
    install.packages(
        pkg_url,
        repos = NULL,
        type = "source",
        ask = F,
        dependencies = T
    )
}


# FUNC: Installs any package dependencies and then the package ----------------
try_install <- function(new_pkg) {
    message("Trying to install required package: ", new_pkg)
    pkg_deps <- get_pkg_deps(new_pkg)
    if (length(pkg_deps) > 0) {
        new_pkg <- c(pkg_deps, new_pkg)
    }
    for (pkg in pkg_deps) {
        if (pkg %in% avail_bioc_packs) {
            install_bio_pkg(pkg)
        } else {
            install_pkgs(pkg)
        }
    }

}


# FUNC: Loads and installs necessary CRAN packages ----------------------------
check_pkg_install <- function(pkgs) {

    pkgs_needed <- check_needed(pkgs)

    if (length(pkgs_needed) > 0) {
        message("The following missing packages will be installed:\n",
                paste(pkgs_needed, collapse = "\n"))
        for (new_pkg in pkgs_needed) {
            tryCatch(
                pak::pkg_install(new_pkg, ask = F),
                error = function(cond){
                    try_install(new_pkg)
                }
            )
        }
    }

    not_loaded <- setdiff(pkgs, loadedNamespaces())

    if (length(not_loaded) == 0) {
        message("All packages loaded successfully!")
        return(NULL)
    }

    load_success <- sapply(not_loaded, quiet_load)

    if (any(load_success == F)) {
        failed_pkgs <- pkgs[load_success == F]
        message("\n>>The following package(s) failed to install:")
        message(paste(failed_pkgs, collapse = "\n"))
        return(failed_pkgs)
    } else{
        message("All packages loaded successfully!")
        return(NULL)
    }
}


manual_bioc <- function(bio_pkg) {
    if (not_installed(bio_pkg)) {
        install.packages(bio_pkg,
                         repos = "http://bioconductor.org/packages/release/bioc",
                         type = "binary",
                         ask = F)
    }
    if (!loadLibrary(bio_pkg)) {
        BiocManager::install(bio_pkg, update = F, ask = F, dependencies = T)
    }
    if (!loadLibrary(bio_pkg)) {
        try_github_inst(file.path("Bioconductor", bio_pkg))
    }
    loadLibrary(bio_pkg)
}


# Start execution -------------------------------------------------------------

if (is_macos) {
    # Execute the function to check and install macOS binaries
    check_and_install_mac_libs(mac_tools)
    setup_compilers()
    symlink_arrow()
}

if (not_installed("rJava")) {
    try(install.packages("rJava", type = "binary",
                         dependencies = T, ask = F), T)
}

java_inst <- system("java -version", ignore.stdout = T, ignore.stderr = T) == 0
if (!java_inst) {
    temurin <- "brew install --cask temurin"
    system(temurin, ignore.stdout = TRUE, ignore.stderr = TRUE)
}
java_home <- system("/usr/libexec/java_home", intern = TRUE)
Sys.setenv(JAVA_HOME = java_home)

if (not_installed("arrow")) {
    install.packages('arrow', type = "binary", ask = F, dependencies = T,
                     repos = c('https://apache.r-universe.dev',
                               'https://cloud.r-project.org'))
}

bin_pkgs <- c(
    # Core utilities and dependencies
    "curl", "jsonlite", "mime", "openssl", "R6",
    # Testing and coverage
    "covr", "testthat",
    # Web and image handling
    "httpuv", "jpeg", "png", "xml2",
    # Reporting and markdown (ensuring binary installation)
    "knitr", "rmarkdown", "readr"
)

for (pkg in bin_pkgs) {
    binary_install(pkg)
}

if (not_installed("httr")) {
    install.packages("httr", ask = F, dependencies = T)
}

binary_install("librarian")
binary_install("BiocManager")

if (not_installed("BiocGenerics")) {
    BiocManager::install("BiocGenerics", update = F, ask = F)
}

if (not_installed("Biobase")) {
    BiocManager::install("Biobase", update = F, ask = F, type = "source")
}

stopifnot(loadLibrary("devtools"))
stopifnot(loadLibrary("librarian"))
stopifnot(loadLibrary("BiocManager"))
stopifnot(loadLibrary("Biobase"))


# List Classifier Core Packages -----------------------------------------------
corePkgs <- c(
    'XML',
    "randomForest",
    "glmnet",
    "ggplot2",
    "gridExtra",
    "knitr",
    "pander",
    "gmp"
)

# List Prerequisite Packages --------------------------------------------------
preReqPkgs <- c(
    # Foundational packages and database interfaces
    "MASS", "BH", "filelock", "rjson", "RSQLite", "DBI", "AnnotationDbi",

    # Data manipulation and tidy selections
    "tidyselect", "dbplyr", "dplyr", "tidyr", "plyr",

    # Genomic and bioinformatics packages
    "biomaRt", "GenomicAlignments", "GenomicFeatures", "GenomicRanges",
    "annotate", "S4Vectors", "IRanges", "SummarizedExperiment",
    "genefilter", "DNAcopy", "rtracklayer",

    # Preprocessing and computational utilities
    "preprocessCore", "matrixStats", "nleqslv", "quadprog", "reshape2", "reshape",

    # Specialized bioinformatics and methylation analysis
    "FDb.InfiniumMethylation.hg19", "HDF5Array", "BiocParallel",

    # Additional statistical and plotting support
    "ade4", "beanplot", "nor1mix", "mclust"
)


# Group 1: Core Utilities and System Tools
core_util <- c(
    "rlang", "Rcpp", "cli", "crayon", "stringi", "stringr", "magrittr",
    "digest", "memoise", "abind", "askpass", "assertr", "assertthat", "brio",
    "cachem", "callr", "codetools", "commonmark", "clipr", "config", "cpp11",
    "credentials", "debugme", "desc", "diffobj", "dlm", "docstring",
    "downloader", "easypackages", "fBasics", "foreign", "formattable",
    "Formula", "fracdiff", "fresh", "gdata", "gdtools", "gert", "getopt",
    "coro", "gitcreds", "GlobalOptions", "globals", "gmodels", "GO.db", "httr",
    "iterators", "ini", "lazyeval", "openssl", "optparse", "R.cache", "sys"
)

# Group 2: Bioconductor and Biological Data Tools
bioc_tools <- c(
    "BiocVersion", "BiocManager", "BiocStyle", "biocViews", "BSgenome",
    "clusterProfiler", "cn.mops", "CopyNumberPlots", "ensembldb", "cellranger",
    "pathview", "Biobase", "MethylAid", "enrichplot"
)

# Group 3: Documentation, Reporting, Interface, and Dashboard Tools
doc_report <- c(
    "markdown", "evaluate", "tinytex", "rmarkdown", "roxygen2", "reprex",
    "bookdown", "boot", "broom", "bs4Dash", "bslib", "sass", "chromote",
    "jquerylib", "selectr", "htmlTable", "htmltools", "htmlwidgets", "webshot",
    "xtable", "chron", "later"
)

# Group 4: Statistical Modeling and Machine Learning
stat_model <- c(
    "checkmate", "carData", "car", "caret", "changepoint", "cmprsk",
    "classInt", "conquer", "dichromat", "ggfortify", "ggrepel", "lars", "lme4",
    "lmodel2", "lmtest", "mlr", "modelr", "modeltools", "nnet", "nloptr",
    "numDeriv", "pbkrtest", "performance", "prodlim", "quantreg", "qvalue",
    "robustbase", "RPMM", "rpart", "Rtsne", "strucchange"
)

# Group 5: Parallel, Batch, and Workflow Management
parallel_work <- c(
    "parallel", "batchtools", "clustermq", "doParallel", "doSNOW",
    "parallelMap", "parameters", "ParamHelpers", "progressr", "promises",
    "future", "targets"
)

# Group 6: Visualization, Graphics, and Plotting Extensions
viz_graphics <- c(
    "systemfonts", "Cairo", "svglite", "extrafont", "extrafontdb",
    "fontawesome", "animation", "ape", "aplot", "bezier", "circlize",
    "ggeffects", "ggforce", "ggfun", "ggnewscale", "ggplotify", "ggraph",
    "ggtext", "ggthemes", "ggtree", "ggupset", "graphlayouts", "gridBase",
    "gridGraphics", "gridtext", "gsubfn", "gt", "hexbin", "hrbrthemes",
    "plot3D", "plotrix", "polyclip", "polynom", "pracma", "scatterplot3d",
    "shadowtext", "latticeExtra", "mapdata", "mapproj", "maps", "colorspace",
    "gplots", "pals", "patchwork", "plotly", "Polychrome"
)

# Group 7: Data Manipulation, Import/Export, and I/O
data_io <- c(
    "dtplyr", "openxlsx", "recipes", "redcapAPI", "rio", "raster", "zip",
    "zoo", "qpdf", "qs", "pillar", "RMySQL", "MatrixModels", "jpeg", "readxl",
    "sqldf"
)

# Group 8: Development, Package Building, and Interoperability
dev_interop <- c(
    "pkgbuild", "pkgconfig", "pkgload", "reticulate", "rJava", "Rgraphviz",
    "rprojroot", "rstudioapi", "rex", "rematch", "rematch2", "lintr",
    "listenv", "processx", "ProtGenerics", "proto", "ps", "withr", "xfun",
    "devtools", "gh", "remotes", "sessioninfo", "sourcetools", "testthat"
)

# Group 9: Miscellaneous and Specialized Tools
misc_tools <- c(
    "beepr", "impute", "insight", "ipred", "labelVector", "narray", "pingr",
    "praise", "prettyunits", "proxy", "purrr", "qdapTools", "RBGL", "Rmpfr",
    "rversions", "rvest", "sjmisc", "tensorflow", "testit", "tidygraph",
    "tidytree", "tidyverse", "tiff", "timeDate", "timeSeries", "torch",
    "tseries", "umap", "units", "utf8", "uuid", "webshot2", "websocket",
    "xopen", "yaml", "stringdist"
)

biocPkgs <- c(
    "conumee",
    "lumi",
    "methylumi"
)


if (not_installed("mapview")) {
    tryCatch(
        devtools::install_github("r-spatial/mapview", dependencies = T,
                                 upgrade = "never", auth_token = NULL),
        error = function(e){
            install.packages("mapview", dependencies = T, verbose = T,
                             ask = F, type = "binary")
        }
    )
}

any_failed0 <- check_pkg_install(corePkgs)

if (not_installed("urca")) pak::pkg_install("urca", ask = F)
library("urca")

rhd_pkgs <- c("rhdf5", "Rhtslib", "Rhdf5lib", "HDF5Array", "rhdf5filters")
any_failed_rhd <- check_pkg_install(rhd_pkgs)

if (not_installed("GenomeInfoDbData")) {
    tryCatch(
        pak::pkg_install("GenomeInfoDbData", ask = F),
        error = function(cond){
            try_github_inst("Bioconductor/GenomeInfoDbData")
        }
    )
}
library("GenomeInfoDbData")

if (not_installed("GenomeInfoDb")) {
    tryCatch(
        pak::pkg_install("GenomeInfoDb", ask = F),
        error = function(cond){
            try_github_inst("Bioconductor/GenomeInfoDb")
        }
    )
}
library("GenomeInfoDb")

if (not_installed("ff")) {
    install.packages("ff", type = "binary", ask = F, dependencies = T)
}

if (not_installed("Hmisc")) {
    tryCatch(
        pak::pkg_install("Hmisc", ask = F),
        error = function(cond){
            bio_url <- "https://cran.r-project.org/src/contrib/Hmisc_5.1-3.tar.gz"
            install.packages(bio_url, repos = NULL, type = "source", ask = F, dependencies = T)
        }
    )
}

if (not_installed("karyoploteR")) {
    try(pak::pkg_install("karyoploteR", ask = F), silent = TRUE)
}

if (not_installed("karyoploteR")) {
    bio_url <- "https://www.bioconductor.org/packages/release/bioc/bin/macosx"
    karyo_tgz <- "karyoploteR_1.30.0.tgz"
    if (arch != "x86_64") {
        tgz_url <- file.path(bio_url, "big-sur-arm64/contrib/4.4", karyo_tgz)
    } else{
        tgz_url <- file.path(bio_url, "big-sur-x86_64/contrib/4.4", karyo_tgz)
    }
    install_url(tgz_url)
}

if (not_installed("Rsamtools")) {
    tryCatch(
        manual_bioc("Rsamtools"),
        error = function(cond) {
            pak::pkg_install("Rsamtools", ask = F)
        }
    )
}


bio_url <-
    "https://bioconductor.org/packages/release/data/annotation/src/contrib"

if (not_installed("FDb.InfiniumMethylation.hg19")) {
    tryCatch(
        pak::pkg_install("FDb.InfiniumMethylation.hg19", ask = F),
        error = function(cond) {
            install_url(file.path(bio_url, "org.Hs.eg.db_3.19.1.tar.gz"))
            install_url(file.path(bio_url, "TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2.tar.gz"))
            install_url(file.path(bio_url, "FDb.InfiniumMethylation.hg19_2.2.0.tar.gz"))
        }
    )
}

if (not_installed("IlluminaHumanMethylationEPICv2manifest")) {
    try_github_inst("mwsill/IlluminaHumanMethylationEPICv2manifest")
}

if (not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")) {
    install_url(file.path(bio_url, "IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1.tar.gz"))
}

if (not_installed("IlluminaHumanMethylation450kmanifest")) {
    install_url(file.path(bio_url, "IlluminaHumanMethylation450kmanifest_0.4.0.tar.gz"))
}

if (not_installed("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")) {
    install_url(file.path(bio_url, "IlluminaHumanMethylationEPICanno.ilm10b2.hg19_0.6.0.tar.gz"))
}

if (not_installed("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")) {
    install_url(file.path(bio_url, "IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0.tar.gz"))
}

if (not_installed("minfi")) {
    if (not_installed("sparseMatrixStats")) {
        pak::pkg_install("sparseMatrixStats", ask = F)
    }
    try_github_inst("mwsill/minfi")
}


if (not_installed("fields")) {
    install_deps("fields")
    tryCatch(
        pak::pkg_install("fields", ask = F),
        error = function(cond) {
            manual_bioc("fields")
        }
    )
}

if (not_installed("RnBeads")) {
    if (not_installed("gplots")) {
        install_deps("gplots")
        pak::pkg_install("gplots", ask = F)
    }
    rn_deps <- c(
        "BiocGenerics",
        "S4Vectors",
        "GenomicRanges",
        "MASS",
        "cluster",
        "fields",
        "ggplot2",
        "grid",
        "gridExtra",
        "limma",
        "matrixStats",
        "methods",
        "illuminaio",
        "methylumi",
        "plyr",
        "IRanges"
    )
    check_pkg_install(rn_deps)
    install_deps("RnBeads")
    pak::pkg_install("RnBeads", ask = F)
}

if (not_installed("conumee2.0")) {
    try_github_inst("hovestadtlab/conumee2")
}

any_fail <- check_pkg_install(preReqPkgs)
any_fail2 <- check_pkg_install(biocPkgs)

if (not_installed("mgmtstp27")) {
    gitLink <-
        "https://github.com/badozor/mgmtstp27/raw/master/archive/mgmtstp27_0.6-3.tar.gz"
    install_url(gitLink)
}

if (not_installed("Rcpp")) {
    tryCatch(
        pak::pkg_install("Rcpp", ask = F),
        error = function(cond) {
            try_github_inst("RcppCore/Rcpp")
        }
    )
}

spat_config <-
    '--with-proj-lib=/usr/local/lib/ --with-proj-include=/usr/local/include/'
options(configure.args = c("sf" = spat_config))

if (not_installed("sf")) {
    tryCatch(
        install.packages("sf", type = "source", dependencies = T, ask = F),
        error = function(e) {
            remotes::install_github("r-spatial/sf", configure.args = spat_config,
                                    dependencies = T, upgrade = "never")
        }
    )
}

invisible(gc())

loadLibrary("BiocManager")
loadLibrary("Biobase")

terraDep <- c('tinytest', 'ncdf4', 'leaflet')
check_pkg_install(terraDep)

if (not_installed("terra")) {
    tryCatch(
        pak::pkg_install("terra", ask = F),
        error = function(cond){
            install.packages(
                'terra', repos = 'https://rspatial.r-universe.dev', ask = F,
                dependencies = T, verbose = T, Ncpus = 6)
        }
    )
}

if (not_installed("FField")) {
    tryCatch(
        pak::pkg_install("FField", ask = F),
        error = function(cond) {
            cran_link <-
                "https://cran.r-project.org/src/contrib/Archive/FField/FField_0.1.0.tar.gz"
            install_url(cran_link)
        }
    )
}

if (not_installed("vdiffr")) {
    pak::pkg_install("vdiffr", ask = F)
}

if (not_installed("UCSC.utils")) {
    pak::pkg_install("UCSC.utils", ask = F)
}

if (not_installed("GenVisR")) {
    if (not_installed("VariantAnnotation")) {
        tryCatch(
            manual_bioc("VariantAnnotation"),
            error = function(cond){
                pak::pkg_install("VariantAnnotation", ask = F)
            }
        )
    }
    tryCatch(
        try_github_inst("griffithlab/GenVisR"),
        error = function(cond){
            pak::pkg_install("GenVisR", ask = F)
        }
    )
}

if (not_installed("forecast")) pak::pkg_install("forecast", ask = F)
if (not_installed("quantreg")) pak::pkg_install("quantreg", ask = F)

Sys.setenv(TORCH_INSTALL = "1")
options(needs.auto = TRUE)

msg_pkg <- function(li_name){
    message("Installing package list: ", li_name, "...", "\n")
}

msg_pkg("core_util")
any_fail3 <- check_pkg_install(core_util)
msg_pkg("bioc_tools")
any_fail3 <- check_pkg_install(bioc_tools)
msg_pkg("doc_report")
any_fail3 <- check_pkg_install(doc_report)
msg_pkg("parallel_work")
any_fail3 <- check_pkg_install(parallel_work)
msg_pkg("viz_graphics")
any_fail3 <- check_pkg_install(viz_graphics)
msg_pkg("data_io")
any_fail3 <- check_pkg_install(data_io)
msg_pkg("dev_interop")
any_fail3 <- check_pkg_install(dev_interop)
msg_pkg("misc_tools")
any_fail3 <- check_pkg_install(misc_tools)

invisible(gc())

if (not_installed("foghorn")) {
    pak::pkg_install("foghorn", ask = F)
}

water_pkgs <- c("mlr", "wateRmelon", "RPMM", "impute")
any_fail4 <- check_pkg_install(water_pkgs)

devtools::source_url(file.path(github_main, "R", "LoadInstallPackages.R"))
