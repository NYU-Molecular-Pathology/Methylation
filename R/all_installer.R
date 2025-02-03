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

# Execute the function to check and install macOS binaries
check_and_install_mac_libs(mac_tools)


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


setup_brew <- function() {
  # Define the commands
  commands <- c(
    'echo >> "$HOME/.zprofile"',
    'echo \'eval "$(/opt/homebrew/bin/brew shellenv)"\' >> "$HOME/.zprofile"',
    'eval "$(/opt/homebrew/bin/brew shellenv)"'
  )

  # Execute the commands
  for (cmd in commands) {
    system(cmd)
  }
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
    setup_brew()
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

# Main setup
setup_compilers <- function() {
  ensure_homebrew()
  Sys.unsetenv(c("CC", "CXX", "OBJC", "LDFLAGS", "CPPFLAGS", "PKG_CFLAGS",
                 "PKG_LIBS", "LD_LIBRARY_PATH", "R_LD_LIBRARY_PATH"))
  set_env_vars()
  options(Ncpus = 6)
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
  return(!pkgName %in% rownames(installed.packages()))
}

binary_install <- function(pkg) {
  if (!pkg %in% rownames(installed.packages())) {
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
    error = function(e) {
      download_pkg_unzip(git_repo, zip_name = "master.zip")
    }
  )

  unzipped_dir <- list.dirs(local_dir, full.names = TRUE, recursive = FALSE)
  pkg_dir <- unzipped_dir[grepl(basename(git_repo), unzipped_dir)]

  install.packages(pkg_dir, repo = NULL, type = "source", dependencies = T)
}


# FUN: Installs package from a Github repository ------------------------------
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


# FUN: TryCatch for installing with devtools and install.packages -------------
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


# FUN: Quietly loads package library without messages -------------------------
quiet_load <- function(pkg_name) {
  if (isNamespaceLoaded(pkg_name) == F){
    libLoad <- suppressWarnings(suppressPackageStartupMessages(
      library(pkg_name, character.only = T, logical.return = T, quietly = T)
    ))
  }else{libLoad <- TRUE}
  message(pkg_name, " loaded... ", libLoad)
  pkg_vec <- c(pkg_n = libLoad)
  names(pkg_vec) <- pkg_name
  return(pkg_vec)
}

# FUN: Checks required package if not installs binary -------------------------
require_pkg <- function(pkg, pkg_type = "source") {
  install_opts$type <- pkg_type
  if (!requireNamespace(pkg, quietly = T)) {
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

if (!requireNamespace("pak", quietly = T)) {
  install.packages("pak", dependencies = T, ask = F)
}
stopifnot(loadLibrary("pak"))

# FUN: Returns all packages that are not installed ----------------------------
check_needed <- function(pkgs) {
    neededPkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
    return(neededPkgs)
}


# FUN: Returns all package dependencies that are not installed ----------------
get_pkg_deps <- function(pkgs) {
  deps_list <- unique(pak::pkg_deps(pkgs)$package)
  all_pkgs <- setdiff(deps_list, pkgs)
  return(check_needed(all_pkgs))
}

# FUN: Attempts source and binary package installation ------------------------
install_pkgs <- function(needed_pkgs, pkg_type = "source") {
  install_opts$type <- pkg_type
  tryCatch(
    do.call(install.packages, c(list(pkgs = needed_pkgs), install_opts)),
    error = function(e) {
      message("\nInitial installation failed! Trying binary install.\n")
      install_opts$type <- "binary"
      do.call(install.packages, c(list(pkgs = needed_pkgs), install_opts))
    }
  )
}

# FUN: Attempts source and binary install of Bioconductor ---------------------
install_bio_pkg <- function(pkg_deps) {
  params <- list(dependencies = c("Depends", "Imports", "LinkingTo"),
                 ask = FALSE, update = TRUE)
  for (pkg in pkg_deps) {
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


# FUN: Checks if package installed from BioConductor --------------------------
check_bio_install <- function(pkgs) {
  needed_pkgs <- check_needed(pkgs)

  if (length(needed_pkgs) > 0) {
    for (new_pkg in needed_pkgs) {
      pkg_deps <- get_pkg_deps(new_pkg)
      if (length(pkg_deps) > 0) {
        install_bio_pkg(pkg_deps)
      }
      install_bio_pkg(new_pkg)
    }
  }
  sapply(pkgs, quiet_load)
}


# FUN: Installs any package dependencies and then the package -----------------
try_install <- function(new_pkg) {
  message("Trying to install required package: ", new_pkg)
  pkg_deps <- get_pkg_deps(new_pkg)
  if (length(pkg_deps) > 0) {
    message("The following missing dependencies will be installed:\n",
            paste(pkg_deps, collapse = "\n"))
    for (pkg in pkg_deps) {
      if (pkg %in% avail_bioc_packs) {
        check_bio_install(pkg)
      } else {
        install_pkgs(pkg)
      }
    }
  }
  install_pkgs(new_pkg)
}


# FUN: Loads and installs necessary CRAN packages -----------------------------
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

  load_success <- sapply(pkgs, quiet_load)
  if (any(load_success == F)) {
    failed_pkgs <- pkgs[load_success == F]
    message("\n>>The following package(s) failed to install:")
    message(paste(failed_pkgs, collapse = "\n"))
    return(failed_pkgs)
  } else {
    return(NULL)
  }
}


load_install <- function(pkg_list) {
  message("Loading packages:\n",
          paste0(capture.output(pkg_list), collapse = "\n"))
  librarian::shelf(
    pkg_list,
    ask = F,
    update_all = F,
    quiet = F,
    dependencies = T
  )
  missing_pkgs <- pkg_list[!(pkg_list %in% rownames(installed.packages()))]
  if (length(missing_pkgs) > 0) {
    for (pkg_missed in missing_pkgs) {
      tryCatch(expr = {
        BiocManager::install(
          pkg_missed,
          update = F,
          ask = F,
          dependencies = T,
          type = "source"
        )
      }, error = {
        install.packages(
          pkg_missed,
          dependencies = T,
          verbose = T,
          type = "source",
          ask = F
        )
      })
    }
  }
}


manual_bioc <- function(bio_pkg) {
  if (checkPkg(bio_pkg)) {
    install.packages(bio_pkg,
                     repos = "http://bioconductor.org/packages/release/bioc",
                     type = "binary",
                     ask = F)
  }
  if (!loadLibrary(bio_pkg)) {
    try_github_inst(file.path("Bioconductor", bio_pkg))
  }
  if (!loadLibrary(bio_pkg)) {
    BiocManager::install(bio_pkg, update = F, ask = F, dependencies = T)
  }
  loadLibrary(bio_pkg)
}


# Start execution -------------------------------------------------------------

if (is_macos) {
  setup_compilers()
}

if (!requireNamespace("rJava", quietly = T)) {
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

if (checkPkg("arrow")) {
  install.packages('arrow', type = "binary", ask = F, dependencies = T,
                   repos = c('https://apache.r-universe.dev',
                             'https://cloud.r-project.org'))
}

bin_pkgs <- c(
  "curl",
  "jsonlite",
  "mime",
  "openssl",
  "R6",
  "covr",
  "httpuv",
  "jpeg",
  "knitr",
  "png",
  "readr",
  "rmarkdown",
  "testthat",
  "xml2"
)

for (pkg in bin_pkgs) {
  if (checkPkg(pkg)) {
    install.packages(pkg, ask = F, dependencies = T, type = 'binary')
  }
}

if (checkPkg("httr")) {
  install.packages("httr", ask = F, dependencies = T)
}

binary_install("librarian")

if (checkPkg("BiocManager")) {
  binary_install("BiocManager")
}

if (checkPkg("BiocGenerics")) {
  BiocManager::install("BiocGenerics", update = F, ask = F)
}

if (checkPkg("Biobase")) {
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
  'plyr',
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
  'SummarizedExperiment',
  'genefilter',
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

# List Remaining Packages -----------------------------------------------------
pkgs <- c(
  "abind",
  "animation",
  "ape",
  "aplot",
  "askpass",
  "assertr",
  "assertthat",
  "batchtools",
  "beepr",
  "bezier",
  "Biobase",
  #"BiocCheck",
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
  "coro",
  "colorspace",
  "commonmark",
  "compiler",
  "config",
  "conquer",
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
  "extrafont",
  "extrafontdb",
  "fastmap",
  "fBasics",
  "fontawesome",
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
  "conumee",
  "lumi",
  "methylumi"
)


if (checkPkg("mapview")) {
  tryCatch(
    remotes::install_github("r-spatial/mapview", dependencies = T,
                            upgrade = "never", auth_token = NULL),
    error = install.packages("mapview", dependencies = T, verbose = T,
                             ask = F, type = "binary")
  )
}

any_failed0 <- check_pkg_install(corePkgs)

if (checkPkg("urca")) {
  pak::pkg_install("urca", ask = F)
  #install.packages("urca", ask = F, dependencies = T, verbose = T)
}
library("urca")

rhd_pkgs <- c("rhdf5", "Rhtslib", "Rhdf5lib", "HDF5Array", "rhdf5filters")
any_failed_rhd <- check_pkg_install(rhd_pkgs)

if (checkPkg("GenomeInfoDbData")) {
  tryCatch(
    pak::pkg_install("GenomeInfoDbData", ask = F),
    error = function(cond){
      try_github_inst("Bioconductor/GenomeInfoDbData")
    }
  )
}
library("GenomeInfoDbData")

if (checkPkg("GenomeInfoDb")) {
  tryCatch(
    pak::pkg_install("GenomeInfoDb", ask = F),
    error = function(cond){
      try_github_inst("Bioconductor/GenomeInfoDb")
    }
  )
}
library("GenomeInfoDb")

if (checkPkg("ff")) {
  install.packages("ff", type = "binary", ask = F, dependencies = T)
}

if (!requireNamespace("Hmisc", quietly = T)) {
  tryCatch(
    pak::pkg_install("Hmisc", ask = F),
    error = function(cond){
      bio_url <- "https://cran.r-project.org/src/contrib/Hmisc_5.1-3.tar.gz"
      install.packages(bio_url, repos = NULL, type = "source", ask = F, dependencies = T)
    }
  )
}

if (checkPkg("karyoploteR")) {
  try(pak::pkg_install("karyoploteR", ask = F), silent = TRUE)
}

if (checkPkg("karyoploteR")) {
  bio_url <- "https://www.bioconductor.org/packages/release/bioc/bin/macosx"
  karyo_tgz <- "karyoploteR_1.30.0.tgz"
  if (arch != "x86_64") {
    tgz_url <- file.path(bio_url, "big-sur-arm64/contrib/4.4", karyo_tgz)
  } else{
    tgz_url <- file.path(bio_url, "big-sur-x86_64/contrib/4.4", karyo_tgz)
  }
  install.packages(tgz_url, repos = NULL, type = "source", ask = F, dependencies = T)
}

if (checkPkg("Rsamtools")) {
  tryCatch(
    manual_bioc("Rsamtools"),
    error = function(cond) {
      pak::pkg_install("Rsamtools", ask = F)
    }
  )
}

if (checkPkg("FDb.InfiniumMethylation.hg19")) {
  tryCatch(
    pak::pkg_install("FDb.InfiniumMethylation.hg19", ask = F),
    error = function(cond) {
      bio_url <-
        "https://bioconductor.org/packages/release/data/annotation/src/contrib"
      pkg_url1 <- file.path(bio_url, "org.Hs.eg.db_3.19.1.tar.gz")
      pkg_url2 <- file.path(bio_url,
                            "TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2.tar.gz")
      pkg_url3 <- file.path(bio_url,
                            "FDb.InfiniumMethylation.hg19_2.2.0.tar.gz")
      install.packages(
        pkg_url1,
        repos = NULL,
        type = "source",
        ask = F,
        dependencies = T
      )
      install.packages(
        pkg_url2,
        repos = NULL,
        type = "source",
        ask = F,
        dependencies = T
      )
      install.packages(
        pkg_url3,
        repos = NULL,
        type = "source",
        ask = F,
        dependencies = T
      )
    }
  )
}

if (checkPkg("IlluminaHumanMethylationEPICv2manifest")) {
  try_github_inst("mwsill/IlluminaHumanMethylationEPICv2manifest")
}

if (checkPkg("IlluminaHumanMethylation450kanno.ilmn12.hg19")) {
  bio_url <- "https://bioconductor.org/packages/release/data/annotation/src/contrib/IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1.tar.gz"
  install.packages(bio_url, repos = NULL, type = "source", ask = F, dependencies = T)
}

if (checkPkg("IlluminaHumanMethylation450kmanifest")) {
  bio_url <- "https://bioconductor.org/packages/release/data/annotation/src/contrib/IlluminaHumanMethylation450kmanifest_0.4.0.tar.gz"
  install.packages(bio_url, repos = NULL, type = "source", ask = F, dependencies = T)
}

if (checkPkg("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")) {
  bio_url <- "https://bioconductor.org/packages/release/data/annotation/src/contrib/IlluminaHumanMethylationEPICanno.ilm10b2.hg19_0.6.0.tar.gz"
  install.packages(
    bio_url,
    repos = NULL,
    type = "source",
    ask = F,
    dependencies = T
  )
}

if (checkPkg("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")) {
  bio_url <- "https://www.bioconductor.org/packages/release/data/annotation/src/contrib/IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0.tar.gz"
  install.packages(bio_url, repos = NULL, type = "source", ask = F, dependencies = T)
}

if (checkPkg("minfi")) {
  if (checkPkg("sparseMatrixStats")) {
    pak::pkg_install("sparseMatrixStats", ask = F)
  }
  try_github_inst("mwsill/minfi")
}

install_deps <- function(pkg) {
  any_deps <- get_pkg_deps(pkg)
  if (length(any_deps) > 0) {
    pak::pkg_install(any_deps, ask = F)
  }
}

if (checkPkg("fields")) {
  install_deps("fields")
  tryCatch(
    pak::pkg_install("fields", ask = F),
    error = function(cond) {
      manual_bioc("fields")
    }
  )
}

if (checkPkg("RnBeads")) {
  if (checkPkg("gplots")) {
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

if (checkPkg("conumee2.0")) {
  try_github_inst("hovestadtlab/conumee2")
}

any_fail <- check_pkg_install(preReqPkgs)

any_fail2 <- check_pkg_install(biocPkgs)

if (checkPkg("mgmtstp27")) {
  gitLink <-
    "https://github.com/badozor/mgmtstp27/raw/master/archive/mgmtstp27_0.6-3.tar.gz"
  install.packages(
    gitLink, repos = NULL, dependencies = T, verbose = T, type = "source", ask = F
  )
}

if (checkPkg("Rcpp")) {
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

if (checkPkg("sf")) {
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

if (checkPkg("terra")) {
  tryCatch(
    pak::pkg_install("terra", ask = F),
    error = function(cond){
      install.packages(
        'terra', repos = 'https://rspatial.r-universe.dev', ask = F,
        dependencies = T, verbose = T, Ncpus = 6)
    }
  )
}

if (checkPkg("FField")) {
  tryCatch(
    pak::pkg_install("FField", ask = F),
    error = function(cond) {
      cran_link <-
        "https://cran.r-project.org/src/contrib/Archive/FField/FField_0.1.0.tar.gz"
      install.packages(
        cran_link,
        repos = NULL,
        dependencies = T,
        verbose = T,
        type = "source",
        ask = F
      )
    }
  )
}

if (checkPkg("vdiffr")) {
  pak::pkg_install("vdiffr", ask = F)
}

if (checkPkg("UCSC.utils")) {
  pak::pkg_install("UCSC.utils", ask = F)
}

if (checkPkg("GenVisR")) {
  if (checkPkg("VariantAnnotation")) {
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

if (checkPkg("forecast")) {
  pak::pkg_install("forecast", ask = F)
}

if (checkPkg("quantreg")) {
  pak::pkg_install("quantreg", ask = F)
}

Sys.setenv(TORCH_INSTALL = "1")
options(needs.auto = TRUE)


# Ensure any required libraries are symlinked
#library_path <- "/usr/local/opt/apache-arrow/lib/libarrow.1801.dylib"

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

llvm_path = get_prefix("llvm")
unwind_libpath <- file.path(llvm_path, "lib", "unwind", "libunwind.1.dylib")
missing_path <- file.path(llvm_path, "lib", "libunwind.1.dylib")
check_symlink(unwind_libpath, missing_path)

arrow_path = get_prefix("apache-arrow")
arrow_libs <- dir(file.path(arrow_path, "lib"), all.files = TRUE, full.names = TRUE)
arrow_dylib <- arrow_libs[grepl("^libarrow\\.[0-9]{4}\\.dylib$", basename(arrow_libs))]
miss_arrow <- file.path(arrow_path, "lib", "libarrow.1700.dylib")
#check_symlink(arrow_dylib, miss_arrow)

any_fail3 <- check_pkg_install(pkgs)

invisible(gc())

if (checkPkg("foghorn")) {
  pak::pkg_install("foghorn", ask = F)
}

water_pkgs <- c("mlr", "wateRmelon", "RPMM", "impute")
any_fail4 <- check_pkg_install(water_pkgs)

devtools::source_url(file.path(github_main, "R", "LoadInstallPackages.R"))
