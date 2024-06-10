#!/usr/bin/env Rscript
## Script name: get_dependencies.R
## Purpose: Source of global scripts for installing packages and dependencies
## Date Created: June 10, 2024
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

# FUN: Quietly loads package library without messages -------------------------
quiet_load <- function(pkg_name) {
    libLoad <- suppressWarnings(suppressPackageStartupMessages(
        library(pkg_name, character.only = T, logical.return = T, quietly = T)
    ))
    message(pkg_name, " loaded... ", libLoad)
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


# FUN: Installs any package dependencies and then the package -----------------
try_install <- function(new_pkg) {
    message("Trying to install required package: ", new_pkg)
    pkg_deps <- get_pkg_deps(new_pkg)
    if (length(pkg_deps) >= 1) {
        for (pkg in pkg_deps) {
            install_pkgs(pkg)
        }
    }
    install_pkgs(new_pkg)
}


# FUN: Loads and installs necessary CRAN packages -----------------------------
check_pkg_install <- function(pkgs) {
    pkg_deps <- check_needed(pkgs)
    if (length(pkg_deps) > 0) {
        for (new_pkg in pkg_deps) {
            try_install(new_pkg)
        }
    }
    invisible(sapply(pkgs, quiet_load))
}

# Example use:
# base_pkgs <- c("MASS", "cluster", "ff", "fields", "uwot", "BiocManager")
#check_pkg_install(base_pkgs)
