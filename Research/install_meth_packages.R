#!/usr/bin/env Rscript
## Script name: install_meth_packages.R
## Purpose: Source of functions to install packages using the pak package or BiocManager
## Date Created: May 28, 2026
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2026

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  ask = FALSE,
  timeout = max(600, getOption("timeout"))
)

system_type <- Sys.info()[["sysname"]]

not_installed <- function(pkgName) {
  return(!pkgName %in% rownames(utils::installed.packages()))
}

if (system_type == "Darwin") {
  message("macOS detected. If needed, install Xcode Command Line Tools: xcode-select --install")
} else if (system_type == "Windows") {
  message("Windows detected. If needed, install Rtools: https://cran.r-project.org/bin/windows/Rtools/")
} else if (system_type == "Linux") {
  message("Linux detected. If package compilation fails, install required system libraries.")
}

if (not_installed("pak")) {
  utils::install.packages("pak", dependencies = TRUE, ask = FALSE)
}

if (not_installed("BiocManager")) {
  utils::install.packages("BiocManager", dependencies = TRUE, ask = FALSE)
}

install_cran_pkg <- function(pkgName) {
  if (not_installed(pkgName)) {
    message("Installing CRAN package: ", pkgName)

    tryCatch(
      {
        utils::install.packages(pkgName, dependencies = TRUE, ask = FALSE)
      },
      error = function(e) {
        message("Base installer failed for ", pkgName, ". Trying pak...")
        pak::pkg_install(pkgName, ask = FALSE)
      }
    )
  }
}

install_bioc_pkg <- function(pkgName) {
  if (not_installed(pkgName)) {
    message("Installing Bioconductor package: ", pkgName)

    tryCatch(
      {
        BiocManager::install(pkgName, dependencies = TRUE, ask = FALSE, update = FALSE)
      },
      error = function(e) {
        message("BiocManager failed for ", pkgName, ". Trying pak...")
        pak::pkg_install(pkgName, ask = FALSE)
      }
    )
  }
}

cran_packages <- c(
  "tidyverse",
  "ggplot2",
  "dplyr",
  "readr",
  "stringr",
  "tibble",
  "tidyr",
  "purrr",
  "data.table",
  "janitor",
  "here",
  "fs",
  "remotes"
)

bioc_packages <- c(
  "minfi",
  "limma",
  "sesame",
  "sesameData",
  "GenomicRanges",
  "SummarizedExperiment",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19"
)

for (pkgName in cran_packages) {
  install_cran_pkg(pkgName)
}

for (pkgName in bioc_packages) {
  install_bioc_pkg(pkgName)
}

message("Methylation R setup complete.")
