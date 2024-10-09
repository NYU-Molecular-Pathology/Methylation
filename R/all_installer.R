#!/usr/bin/env Rscript
## Script name: all_installer.R
## Purpose: Functions to check if required packages are installed
## Date Created: August 9, 2022
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024

gb <- globalenv(); assign("gb", gb)
options("install.packages.compile.from.source" = "Yes")
options("install.packages.check.source" = "no")
options(Ncpus = 4)

github_main <-
    "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main"

is_macos <- Sys.info()[['sysname']] == "Darwin"

if (getRversion() <= "4.2.2") {
    stop("Your R version is ", R.version.string, ". Update to 4.4.0 or later.")
}

supM <- function(pk) {
    return(suppressPackageStartupMessages(suppressWarnings(pk)))
}

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

module_exists <- function(module_name) {
    mod_vers <- paste(module_name, "--version")
    mod_check <- tryCatch(system(mod_vers, intern = TRUE, ignore.stderr = TRUE),
                          error = function(e) return(FALSE),
                          warning = function(w) return(FALSE)
    )
    if (mod_check == FALSE) {
        return(FALSE)
    }else {
        return(TRUE)
    }
}

getBrewDir <- function(module_name) {
    brew_cmd <- paste("brew --prefix", module_name)
    module_path <- tryCatch(system(brew_cmd, intern = TRUE, ignore.stderr = TRUE),
                            error = function(e) NULL,
                            warning = function(e) NULL)
    return(module_path)
}

# FUN: Checks if system compilers installed -----------------------------------
check_brew_pkgs <- function(){
    # Check if brew installed
    brew_installed <- module_exists("brew")
    if (!brew_installed) {
        install_homebrew()
    }
    # Check if GCC installed
    gfortran_path <- try(Sys.which("gfortran")[[1]], T)
    if (gfortran_path == "") {
        message("gfortran is not installed.")
        brew_install("gcc")
    }
    # Check if LLVM installed
    llvm_path <- getBrewDir("llvm")
    # Check if LLVM installed
    if (is.null(llvm_path)) {
        brew_install("llvm")
    }
    # Check if open-mpi installed
    mpi_path <- getBrewDir("open-mpi")
    if (is.null(mpi_path)) {
        brew_install("open-mpi")
    }
    # Check if GDAL installed
    gdal_installed <- module_exists("gdalinfo")
    if (!gdal_installed) {
        brew_install("pkg-config")
        brew_install("gdal")
    }
    # Check if proj installed
    proj_check <- tryCatch(
        system("which proj", intern = T, ignore.stderr = T),
        error = function(e) return(FALSE),
        warning = function(w) return(FALSE)
    )
    if (proj_check == F) {
        brew_install("pkg-config")
        brew_install("proj")
    }

    clear_cmd <- "brew update && brew doctor && rm -rf $(brew --cache)"
    #system(clear_cmd, wait = TRUE)
}


locate_mod <- function(module_name) {
    locate_cmd <- paste("locate", module_name)
    module_path <- system(locate_cmd, intern = TRUE)
    if (length(module_path) > 1) {
        match_paths <- grep(paste0(module_name, "$"), module_path, value = TRUE)
        module_path <- match_paths[1]
    }
    return(module_path)
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
}


# FUN: Sets system openmpi flags ----------------------------------------------
set_openmpi <- function() {
    prte_path <- locate_mod("/bin/prte")
    orte_path <- locate_mod("/bin/orted")
    if (length(prte_path) > 0 & length(orte_path) > 0) {
        ln_cmd <- paste("ln -s", prte_path, orte_path)
        try(system(ln_cmd, wait = T), T)
    }

    mpi_path <- getBrewDir("open-mpi")
    mpi_libs <- paste0("-L", mpi_path, "/lib")
    mpi_flag <- paste(mpi_libs, Sys.getenv("LDFLAGS"))
    Sys.setenv(LDFLAGS = mpi_flag)

    cpp_libs <-  paste0("-I", mpi_path,"include")
    cpp_flag <- paste(cpp_libs, Sys.getenv("CPPFLAGS"))
    Sys.setenv(CPPFLAGS = cpp_flag)
}


update_system_path <- function() {
    user_shell <- Sys.getenv("SHELL")

    shell_config <- if (grepl("zsh", user_shell)) {
        path.expand("~/.zshrc")
    } else if (grepl("bash", user_shell)) {
        if (file.exists(path.expand("~/.bash_profile"))) {
            path.expand("~/.bash_profile")
        } else {
            path.expand("~/.bashrc")
        }
    } else {
        stop("Your shell does not support automatic PATH updates!")
    }

    profile_contents <- readLines(shell_config)
    contains_sbin <- any(grepl("/usr/local/sbin", profile_contents))
    path_command <- sprintf('echo "export PATH=\\"/usr/local/sbin:$PATH\\"" >> %s', shell_config)

    if (!contains_sbin) {
        system(path_command, intern = FALSE)
        message(sprintf(
            "Updated %s to include /usr/local/sbin in PATH.", shell_config))
    } else {
        message("PATH already contains sbin in configuration file.")
    }
    Sys.setenv(PATH = paste("/usr/local/bin", Sys.getenv("PATH"), sep = ":"))
    Sys.setenv(LIBRARY_PATH = paste("/usr/local/lib", Sys.getenv("LIBRARY_PATH"), sep = ":"))
}

# FUN: Sets system compiler flags ---------------------------------------------
fix_compiler_flags <- function(){
    check_brew_pkgs()
    #system("brew upgrade && brew update", intern = T, ignore.stderr = T)
    #system("brew doctor", intern = T, ignore.stderr = T)
    #system("brew cleanup", intern = T, ignore.stderr = T)
    clear_enviro()

    # Get paths dynamically using the previously defined functions
    llvm_path <- getBrewDir("llvm")
    arrow_dir <- getBrewDir("apache-arrow")
    gtran_path <- locate_mod("/bin/gfortran")
    flib_dir <- locate_mod("/lib/gcc/current")

    # Setting environment variables using dynamic paths
    Sys.setenv(CC = file.path(llvm_path, "bin/clang"))
    Sys.setenv(CXX = file.path(llvm_path, "bin/clang++"))
    Sys.setenv(CXX11 = file.path(llvm_path, "bin/clang++"))
    Sys.setenv(CXX14 = file.path(llvm_path, "bin/clang++"))
    Sys.setenv(CXX17 = file.path(llvm_path, "bin/clang++"))
    Sys.setenv(CXX1X = file.path(llvm_path, "bin/clang++"))
    Sys.setenv(OBJC = file.path(llvm_path, "bin/clang"))

    Sys.setenv(LDFLAGS = paste0(
        "-L", file.path(llvm_path, "lib"),
        " -L", file.path(llvm_path, "lib/c++"),
        " -Wl,-rpath,", file.path(llvm_path, "lib/c++")
    ))

    Sys.setenv(CPPFLAGS = paste0("-I", file.path(llvm_path, "include")))

    Sys.setenv(PKG_CFLAGS = paste0(
        "-I/usr/local/include -I", file.path(llvm_path, "include"),
        " -I", file.path(arrow_dir, "include")
    ))

    Sys.setenv(PKG_LIBS = paste0(
        "-L/usr/local/lib -L", file.path(llvm_path, "lib"),
        " -L", file.path(arrow_dir, "lib"),
        " -larrow"
    ))

    Sys.setenv(LD_LIBRARY_PATH = "/usr/local/lib")

    Sys.setenv(R_LD_LIBRARY_PATH = paste0(
        "/usr/local/lib:", file.path(llvm_path, "lib/c++")
    ))


    set_openmpi()
    set_gfortran()
    #system("brew cleanup")
    #update_system_path()
    mpi_path <- getBrewDir("open-mpi")
    curr <- paste(file.path(llvm_path, "bin"), Sys.getenv("PATH"), sep = ":")
    Sys.setenv(PATH = paste(file.path(mpi_path, "bin"), curr, sep = ":"))
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
    # Get paths for llvm, apache-arrow, gcc, gfortran, and gcc current lib directory
    llvm_path <- getBrewDir("llvm")
    arrow_dir <- getBrewDir("apache-arrow")
    gcc_dir <- getBrewDir("gcc")
    gtran_path <- locate_mod("/bin/gfortran")
    flib_dir <- locate_mod("/lib/gcc/current")

    # Dynamically constructing the compiler settings vector using file paths
    compiler_settings <- c(
        paste0("CC = ", file.path(llvm_path, "bin/clang")),
        paste0("CXX = ", file.path(llvm_path, "bin/clang++")),
        paste0("CXX11 = ", file.path(llvm_path, "bin/clang++")),
        paste0("CXX14 = ", file.path(llvm_path, "bin/clang++")),
        paste0("CXX17 = ", file.path(llvm_path, "bin/clang++")),
        paste0("CXX1X = ", file.path(llvm_path, "bin/clang++")),
        paste0("OBJC = ", file.path(llvm_path, "bin/clang")),
        paste0("LDFLAGS = -L", file.path(llvm_path, "lib"),
               " -L", file.path(llvm_path, "lib/c++"),
               " -Wl,-rpath,", file.path(llvm_path, "lib/c++")),
        paste0("CPPFLAGS = -I", file.path(llvm_path, "include")),
        paste0("PKG_CFLAGS = -I/usr/local/include -I",
               file.path(llvm_path, "include"),
               " -I", file.path(arrow_dir, "include")),
        paste0("PKG_LIBS = -L/usr/local/lib -L", file.path(llvm_path, "lib"),
               " -L", file.path(arrow_dir, "lib"), " -larrow"),
        paste0("FC = ", gtran_path),
        paste0("FLIBS = -L", flib_dir),
        "LD_LIBRARY_PATH = /usr/local/lib",
        paste0("R_LD_LIBRARY_PATH = /usr/local/lib:",
               file.path(llvm_path, "lib/c++"))
    )
    if (file.exists(makevars_path)) file.remove(makevars_path)
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

binary_install <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        tryCatch(
            install.packages(pkg, dependencies = T, ask = F, type = "binary"),
            error = function(e){
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
        error = function(e){
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


install_opts <- list(
    dependencies = c("Depends", "Imports", "LinkingTo"),
    ask = FALSE,
    update = "never",
    quiet = TRUE,
    repos = 'http://cran.us.r-project.org',
    Ncpus = 4
)

# FUN: Quietly loads package library without messages -------------------------
quiet_load <- function(pkg_name) {
    libLoad <- suppressWarnings(suppressPackageStartupMessages(
        library(pkg_name, character.only = T, logical.return = T, quietly = T)
    ))
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


# FUN: Returns all packages that are not installed ----------------------------
check_needed <- function(pkgs){
    return(pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)])
}


# FUN: Returns all package dependencies that are not installed ----------------
get_pkg_deps <- function(pkgs) {
    deps_list <- tools::package_dependencies(
        pkgs, db = pkg_info, recursive = T,
        which = c("Depends", "Imports", "LinkingTo")
    )
    all_deps <- unique(unlist(deps_list))
    all_pkgs <- setdiff(all_deps, rbase_pkgs)
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
            error = function(e){
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
            if (pkg %in% avail_bioc_packs){
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
            if (new_pkg %in% avail_bioc_packs){
                check_bio_install(new_pkg)
            } else {
                install_pkgs(new_pkg)
            }
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
    #system("export LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib")
}

# Check if architecture is 'arm64' or 'x86_64' --------------------------------
arch <- Sys.info()[["machine"]]

if (arch != "x86_64" & is_macos == T) {
    # Set JAVA_HOME environment variable
    #java_home <- system("/usr/libexec/java_home -v 11", intern = TRUE)
    #system("brew tap homebrew/cask-versions && brew install --cask temurin17")
    java_home <- system("which java", intern = TRUE)
    Sys.setenv(JAVA_HOME = java_home)
    message("JAVA_HOME set to ", java_home)
    if (!requireNamespace("rJava", quietly = T)) {
        try(install.packages("rJava", type = "binary", dependencies = T, ask = F), T)
    }
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

if (checkPkg("devtools")) {
    install.packages("devtools", ask = F, type = "source", dependencies = T)
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
corePkgs <- c("randomForest",
              "glmnet",
              "ggplot2",
              "gridExtra",
              "knitr",
              "pander",
              "gmp")

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
    'XML',
    'bit64',
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
    "conumee",
    "lumi",
    "methylumi"
)


if (checkPkg("mapview")) {
    tryCatch(
        expr = remotes::install_github("r-spatial/mapview", dependencies = T, upgrade = "never"),
        error = install.packages("mapview", dependencies = T, verbose = T, ask = F, type = "binary")
    )
}


any_failed0 <- check_pkg_install(corePkgs)

if (checkPkg("urca")) {
    install.packages("urca", ask = F, dependencies = T, verbose = T)
}
library("urca")

manual_bioc("rhdf5")
manual_bioc("Rhtslib")
manual_bioc("Rhdf5lib")
manual_bioc("HDF5Array")
manual_bioc("rhdf5filters")

if (checkPkg("GenomeInfoDbData")) {
    gb$try_github_inst("Bioconductor/GenomeInfoDbData")
    requireNamespace("GenomeInfoDbData")
}

if (checkPkg("GenomeInfoDb")) {
    gb$try_github_inst("Bioconductor/GenomeInfoDb")
    library("GenomeInfoDb")
}

if (checkPkg("ff")) {
    install.packages("ff", type = "binary", ask = F, dependencies = T)
}

bio_url <- "https://cran.r-project.org/src/contrib/Hmisc_5.1-3.tar.gz"
if (!requireNamespace("Hmisc", quietly = T)) {
    install.packages(bio_url, repos = NULL, type = "source", ask = F, dependencies = T)
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

manual_bioc("Rsamtools")

if (checkPkg("FDb.InfiniumMethylation.hg19")) {
    bio_url <-
        "https://bioconductor.org/packages/release/data/annotation/src/contrib"
    pkg_url1 <- file.path(bio_url, "org.Hs.eg.db_3.19.1.tar.gz")
    pkg_url2 <- file.path(bio_url, "TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2.tar.gz")
    pkg_url3 <- file.path(bio_url, "FDb.InfiniumMethylation.hg19_2.2.0.tar.gz")
    install.packages(pkg_url1, repos = NULL, type = "source", ask = F, dependencies = T)
    install.packages(pkg_url2, repos = NULL, type = "source", ask = F, dependencies = T)
    install.packages(pkg_url3, repos = NULL, type = "source", ask = F, dependencies = T)
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

if (checkPkg("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")) {
    bio_url <- "https://www.bioconductor.org/packages/release/data/annotation/src/contrib/IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0.tar.gz"
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

if (checkPkg("minfi")) {
    try_github_inst("mwsill/minfi")
}

if (checkPkg("fields")) {
    manual_bioc("fields")
}

if (checkPkg("RnBeads")) {
    manual_bioc("RnBeads")
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
    cran_link <-
        "https://cran.r-project.org/src/contrib/Archive/FField/FField_0.1.0.tar.gz"
    install.packages(
        cran_link, repos = NULL, dependencies = T, verbose = T, type = "source", ask = F
    )
}

if (checkPkg("vdiffr")) {
    install.packages("vdiffr", dependencies = T, ask = F, type = "binary")
}
manual_bioc("UCSC.utils")

if (checkPkg("GenVisR")) {
    manual_bioc("VariantAnnotation")
    try_github_inst("griffithlab/GenVisR")
}

if (checkPkg("forecast")) {
    cran_url <- "https://cran.r-project.org/bin/macosx"
    forcast_pkg <- "forecast_8.23.0.tgz"
    if (arch != "x86_64") {
        pkg_bin <- file.path(cran_url, "big-sur-arm64/contrib/4.4", forcast_pkg)
    } else {
        pkg_bin <- file.path(cran_url, "big-sur-x86_64/contrib/4.4", forcast_pkg)
    }
    install.packages(
        pkg_bin,
        repos = NULL,
        type = "source",
        ask = F,
        dependencies = T
    )
}

if (checkPkg("quantreg")) {
    install.packages('quantreg', ask = F, type = 'binary', dependencies = T)
}

any_failed <- gb$check_pkg_install(pkgs)

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

devtools::source_url(file.path(github_main, "R", "LoadInstallPackages.R"))
