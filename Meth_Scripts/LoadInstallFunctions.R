#!/usr/bin/env Rscript
## Script name: LoadInstallFunctions.R
## Purpose: source of global scripts to install packages
## Date Created: June 20, 2024
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024

install_opts <- list(
    dependencies = c("Depends", "Imports", "LinkingTo"),
    ask = FALSE,
    update = "never",
    quiet = TRUE,
    repos = 'http://cran.us.r-project.org',
    Ncpus = 4
)

# Logs error messages from function calls -------------------------------------
try_log_errors <- function(FUN) {
    tryCatch(
        FUN,
        error = function(e) {
            emsg <- paste("Error:", e$message)
            log_file <- file.path(Sys.getenv("HOME"), "errorlog.txt")
            message("The following error has been logged in ", log_file)
            message(emsg)
            write(emsg, log_file, append = TRUE)
        }
    )
}

# FUN: Quietly loads package library without messages -------------------------
quiet_load <- function(pkg_name) {
    libLoad <- suppressWarnings(suppressPackageStartupMessages(
        library(pkg_name, character.only = T, logical.return = T, quietly = T)
    ))
    message(paste0('"', pkg_name, '" package loaded... ', libLoad))
}

# FUN: Installs and loads package if not loaded -------------------------------
require_pkg <- function(pkg, pkg_type = "source") {
    install_opts$type <- pkg_type
    if (!requireNamespace(pkg, quietly = T)) {
        do.call(install.packages, c(list(pkgs = pkg), install_opts))
    }
    quiet_load(pkg)
}

# FUN: Returns all packages that are not installed ----------------------------
check_needed <- function(pkgs) {
    return(pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)])
}

# Check devtools and install package dependencies -----------------------------
require_pkg("devtools", "binary")
require_pkg("utils", "binary")
require_pkg("BiocManager", "binary")
biocRepos <- suppressMessages(BiocManager::repositories())
avail_bioc_packs <- suppressMessages(BiocManager::available())
rbase_pkgs <- rownames(installed.packages(priority = "base"))
pkg_info <- utils::available.packages(repos = biocRepos)

# FUN: Returns all package dependencies not installed as a vector -------------
get_pkg_deps <- function(pkgs) {
    deps_list <- tools::package_dependencies(
        pkgs,
        db = pkg_info,
        which = c("Depends", "Imports", "LinkingTo"),
        recursive = T
    )
    all_deps <- unique(unlist(deps_list))
    all_pkgs <- setdiff(all_deps, rbase_pkgs)
    return(check_needed(all_pkgs))
}


# FUN: Installs packages from a source URL, .tar, or .zip file ----------------
source_install <- function(pkg_name, pkg_url) {
    if (!requireNamespace(pkg_name, quietly = T)) {
        params <- list(
            pkgs = pkg_url,
            repos = NULL,
            type = "source",
            dependencies = T
        )
        tryCatch(
            expr = {
                do.call(install.packages, c(params, Ncpus = 6))
            },
            error = function(cond) {
                do.call(install.packages, c(params, method = "libcurl"))
            },
            warning = function(cond) {
                do.call(install.packages, c(params, method = "auto"))
            }
        )
    }
    quiet_load(pkg_name)
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
brew_install <- function(pkg) {
    message(pkg, " is not installed! Installing ", pkg, " via Homebrew...")
    system(paste("brew install", pkg), intern = TRUE, wait = TRUE)
}


# FUN: Checks if system compilers installed -----------------------------------
check_brew_pkgs <- function() {
    # Check if brew installed
    brew_exists <- file.exists("/usr/local/bin/brew")
    if (!brew_exists) {
        install_homebrew()
    }

    try(system("rm -rf /usr/local/Caskroom/gfortran"), T)

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
    update_br <- "brew update && brew upgrade"
    clear_cmd <- "brew doctor && rm -rf $(brew --cache)"
    system(update_br, wait = TRUE)
    system(clear_cmd, wait = TRUE)
    system("brew cleanup")
}


# FUN: Clears system environment flags ----------------------------------------
clear_enviro <- function() {
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

    tran_info <- system("gfortran -print-search-dirs", intern = TRUE)

    library_paths <- strsplit(tran_info[length(tran_info)], " ")[[1]][-1]
    flibs_value <- paste(sapply(library_paths, function(path) {
        paste("-L", path, sep = "")}), collapse = " ")
    flibs_value <- gsub("=", " ", flibs_value)
    Sys.setenv(FC = gfortran_path)
    Sys.setenv(FLIBS = flibs_value)
}


# FUN: Sets system openmpi flags ----------------------------------------------
set_openmpi <- function() {
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
fix_compiler_flags <- function() {
    #check_brew_pkgs()
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
    #Sys.setenv(GITHUB_PAT = "your_token_here")
}


# FUN: Attempts binary install of Bioconductor --------------------------------
bioc_install <- function(pkg) {
    params <- list(dependencies = c("Depends", "Imports", "LinkingTo"),
                   ask = F, update = F, type = "binary")
    do.call(BiocManager::install, c(list(pkgs = pkg), params))
}


# FUN: Attempts source and binary package installation ------------------------
install_pkgs <- function(curr_pkg, pkg_type = "source") {
    install_opts$type <- pkg_type
    if (any(curr_pkg %in% avail_bioc_packs)) {
        bio_curr <- curr_pkg[curr_pkg %in% avail_bioc_packs]
        return(bioc_install(bio_curr))
    }
    message(paste("Installing package:", curr_pkg), "...")
    tryCatch(
        do.call(install.packages, c(list(pkgs = curr_pkg), install_opts)),
        error = function(e) {
            message("\nInitial installation failed! Trying binary install.\n")
            install_opts$type <- "binary"
            do.call(install.packages, c(list(pkgs = curr_pkg), install_opts))
        }
    )
}

# FUN: Loops vector of packages and installs with Bioconductor if in repo -----
loop_install <- function(pkgs) {
    for (pkg in pkgs) {
        if (pkg %in% avail_bioc_packs) {
            bioc_install(pkg)
        } else {
            install_pkgs(pkg)
        }
    }
}


# FUN: Attempts binary install of Bioconductor if not CRAN --------------------
install_bio_pkg <- function(pkg_deps) {

    for (pkg in pkg_deps) {
        message("Installing BioCpackage: ", pkg)
        new_pkgs <- get_pkg_deps(pkg)
        if (length(new_pkgs) > 0) {
            message("Installing additional dependencies for: ", pkg)
            message(paste(new_pkgs, collapse = "\n"))
            loop_install(new_pkgs)
        }
        loop_install(pkg)
    }
}


# FUN: Checks if package installed from BioConductor --------------------------
check_bio_install <- function(pkgs) {
    needed_pkgs <- check_needed(pkgs)

    if (length(needed_pkgs) > 0) {
        for (new_pkg in needed_pkgs) {
            message("Trying to install required package: ", new_pkg)
            pkg_deps <- get_pkg_deps(new_pkg)
            if (length(pkg_deps) > 0) {
                message(paste("More dependencies found for package:", new_pkg))
                message(paste(pkg_deps, collapse = "\n"))
                install_bio_pkg(pkg_deps)
            }
            install_bio_pkg(new_pkg)
        }
    }

    sapply(pkgs, quiet_load)
}


# FUN: Installs any package dependencies and then the package -----------------
try_install <- function(new_pkg) {
    if (length(check_needed(new_pkg)) == 0) {
        return(try(quiet_load(new_pkg), T))
    }

    message("Trying to install required package: ", new_pkg)
    pkg_deps <- get_pkg_deps(new_pkg)
    if (length(pkg_deps) > 0) {
        message(paste("More dependencies found for package:", new_pkg))
        message(paste(pkg_deps, collapse = "\n"))
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
            if (new_pkg %in% avail_bioc_packs) {
                check_bio_install(new_pkg)
            } else {
                try_install(new_pkg)
            }
        }
    }
    invisible(sapply(pkgs, quiet_load))
}


# FUN: Downloads and unzips pkg download --------------------------------------
download_pkg_unzip <- function(git_repo, zip_name = "main.zip") {
    repo_url <- file.path("https://github.com", git_repo, "archive/refs/heads")
    url_path <- file.path(repo_url, zip_name)
    message("Trying to download from url:\n", url_path)
    down_dir <- file.path(fs::path_home(), "github_pkgs")
    if (!dir.exists(down_dir)) dir.create(down_dir, recursive = T)
    out_file <- file.path(down_dir, zip_name)
    utils::download.file(url = url_path, method = "libcurl", destfile = out_file)
    if (!file.exists(out_file)) {
        message("File failed to download:\n", out_file)
    } else {
        message("File download to:\n", out_file)
        message("Unzipping ", zip_name, " to ", down_dir)
        utils::unzip(out_file, exdir = down_dir)
    }
}


# FUN: Downloads Github repo locally then installs ----------------------------
local_github_pkg_install <- function(git_repo) {

    local_dir <- file.path(fs::path_home(), "github_pkgs")
    message("Github install failed, trying local install from: ", local_dir)
    if (!dir.exists(local_dir)) dir.create(local_dir)
    tryCatch(
        expr = {
            download_pkg_unzip(git_repo, zip_name = "main.zip")
        },
        error = function(e) {
            download_pkg_unzip(git_repo, zip_name = "master.zip")
        }
    )

    unzip_dir <- list.dirs(local_dir, full.names = TRUE, recursive = FALSE)
    stopifnot(length(unzip_dir) > 0)

    message("Downloaded Github packages:\n", paste(unzip_dir, collapse = '\n'))

    pkg_dir <- unzip_dir[grepl(basename(git_repo), unzip_dir)]
    if (length(pkg_dir) == 0) {
        message('Unzipped folder for package "',
                basename(git_repo), '" not found in:\n',
                paste(unzip_dir, collapse = '\n'))
    }

    message("Trying to install Github package saved in:\n", pkg_dir)

    deps <- devtools::dev_package_deps(pkg_dir)
    if (length(deps) > 0) {
        message("Checking Github package dependencies...")
        git_pkg_deps <- rownames(deps)
        check_pkg_install(git_pkg_deps)
    }

    install.packages(pkg_dir, repo = NULL, type = "source", dependencies = T)
}


# FUN: Installs package from a Github repository ------------------------------
install_repo <- function(git_repo) {
    message("Trying github install of ", git_repo)
    tryCatch(
        expr = {
            devtools::install_github(
                git_repo,  upgrade = "always", type = "source",
                dependencies = T, auth_token = NULL,
                subdir = basename(git_repo))
        },
        error = function(e) {
            devtools::install_github(
                git_repo, upgrade = "always", type = "source",
                dependencies = T, auth_token = NULL)
        }
    )
}


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


# FUN: Checks if package version is correct then installs from Github ---------
check_pkg_vers <- function(pkg_name, git_repo, min_vers) {
    
    message("Trying to load package: ", pkg_name)
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
        try_github_inst(git_repo)
    }
    
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
        message("The package ", pkg_name, " failed to install from Github")
        try_github_inst(git_repo)
    }
    
    current_vers <- as.character(utils::packageVersion(pkg_name))
    is_current <- utils::compareVersion(current_vers, min_vers) >= 0
    
    if (!is_current) {
        try(unloadNamespace(pkg_name), T)
        try(remove.packages(pkg_name), T)
        try_github_inst(git_repo)
    }
    
    quiet_load(pkg_name)
}

# FUN: Checks if user home has the local sbin in PATH -------------------------
check_rc_path <- function() {
    line <- 'export PATH="/usr/local/sbin:$PATH"'

    zshrc_check <- system(
        paste("grep -Fxq", shQuote(line), "~/.zshrc"),
        ignore.stdout = TRUE,
        ignore.stderr = TRUE
    )

    if (zshrc_check != 0) {
        system(
            paste("echo", shQuote(line), ">> ~/.zshrc"),
            ignore.stdout = TRUE,
            ignore.stderr = TRUE
        )
    }

    bashrc_check <- system(
        paste("grep -Fxq", shQuote(line), "~/.bashrc"),
        ignore.stdout = TRUE,
        ignore.stderr = TRUE
    )

    if (bashrc_check != 0) {
        system(
            paste("echo", shQuote(line), ">> ~/.bashrc"),
            ignore.stdout = TRUE,
            ignore.stderr = TRUE
        )
    }

}


fix_compiler_flags()
check_rc_path()
# check_pkg_install()
