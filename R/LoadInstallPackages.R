#!/usr/bin/env Rscript
## Script name: LoadInstallPackages.R
## Purpose: Check if the required dependencies are installed and loaded
## Date Created: August 5, 2022
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024

if (!require("devtools")) {
  install.packages("devtools", dependencies = T, ask = F)
}

isMacOS <- Sys.info()[['sysname']] == "Darwin"

cbioLn <- "/Volumes/CBioinformatics/Methylation/classifiers"
epicV2script <- 
  "/Volumes/CBioinformatics/Methylation/Rscripts/install_epic_v2_classifier.R"

# Function: Returns the directory where brew is installed if not then NULL ----
getBrewDir <- function(module_name) {
  brew_cmd <- paste("brew --prefix", module_name)
  module_path <- tryCatch(
    system(brew_cmd, intern = TRUE, ignore.stderr = TRUE),
    error = function(e)
      NULL,
    warning = function(e)
      NULL
  )
  return(module_path)
}

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
  'devtools',
  'remotes',
  'sjmisc',
  'stringi',
  'digest',
  'RCurl',
  'rlang',
  'parallel',
  'grid',
  'gridExtra',
  'knitr',
  'kableExtra',
  'htmltools',
  'ggplot2',
  'plotly',
  'ggfortify',
  'ggrepel',
  'gplots',
  'fastmatch',
  'pals',
  'Polychrome',
  'qdapTools',
  'beepr',
  'xtable',
  'pander',
  'grDevices',
  'graphics',
  'stats',
  'utils',
  'magick',
  'ade4',
  "MASS",
  "R.utils",
  "optparse",
  "targets",
  "usethis",
  "webshot",
  "reshape",
  "reshape2",
  "data.table",
  "DT",
  "scales",
  "RColorBrewer",
  "readxl",
  "stringr",
  "tinytex",
  "future.apply",
  "Rcpp",
  "httpuv",
  "shiny"
)

# GitHub Packages -------------------------------------------------------------
gHubPkgs <- data.frame(
  rmarkdown = 'rstudio/rmarkdown',
  docstring = 'dasonk/docstring',
  rstudioapi = 'rstudio/rstudioapi',
  easypackages = 'jakesherman/easypackages',
  redcapAPI = 'nutterb/redcapAPI',
  crayon = "r-lib/crayon",
  redcap = "epicentre-msf/redcap",
  mdthemes = "thomas-neitmann/mdthemes"
)

# BioConductor Packages  ------------------------------------------------------
biocPkgs <- c(
  'HDF5Array', 'bumphunter','GEOquery', 'lumi', 'rngtools', 'Rtsne',
  'methylumi', 'randomForest', 'glmnet','IlluminaHumanMethylation450kmanifest',
  'IlluminaHumanMethylation450kanno.ilmn12.hg19',
  'IlluminaHumanMethylationEPICmanifest',
  'IlluminaHumanMethylationEPICanno.ilm10b2.hg19',
  'IlluminaHumanMethylationEPICanno.ilm10b4.hg19',
  'MethylAid', 'conumee','BiocParallel', "Biobase", "limma","MatrixGenerics",
  "GenomeInfoDb"
)

# Extra Libraries -------------------------------------------------------------
extraPks <- c(
  "doSNOW",
  "doParallel",
  "foreach",
  "compiler",
  "tidyverse",
  "plyr",
  #"mdthemes",
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
sup <- function(x) return(suppressWarnings(suppressPackageStartupMessages(x)))

sw <- function(pkgOb) try(return(suppressMessages(suppressWarnings(pkgOb))), T)

ld <- function(libName) {
  lib.opts <- list(package = libName, character.only = T,verbose = T, warn.conflicts = F, quietly = F)
  suppressPackageStartupMessages(do.call(library, c(lib.opts)))
  message(libName, " ...load successful")
}

rq <- function(pkgName) {
  !requireNamespace(pkgName, quietly = TRUE)
}

mkred <- function(strMsg) {
  return(crayon::white$bgRed$bold(strMsg))
}

mkblu <- function(strMsg) {
  return(crayon::white$bgBlue$bold(strMsg))
}

mkGrn <- function(strMsg) {
  return(crayon::white$bgGreen$bold(strMsg))
}

msgCheck <- function(pkg, warn = F) {
  if (warn == F) {
    message("Checking ", pkg, "...")
  } else{
    message("\nError caught for package:\n", pkg, "\n-----------")
  }
}


#' pk.inst will load a package library or install it if it does not exist
pk.inst <- function(pkg) {
  msgCheck(pkg)
  if (rq(pkg)) {
    pk.opt <- list(
      pkgs = pkg,
      dependencies = T,
      verbose = T,
      Ncpus = 6
    )
    tryCatch(
      expr = {
        do.call(install.packages, c(pk.opt, list(type = bothType)))
      },
      warning = function(cond) {
        message("\nWarning on package:\n",
                pkg,
                "\n-----------\n",
                cond)
        do.call(install.packages, c(pk.opt, list(type = typeSrc)))
      },
      error = function(cond) {
        msgCheck(cond, T)
        do.call(install.packages, c(pk.opt, list(type = "source")))
      },
      custom_error = function(cond) {
        msgCheck(cond, T)
        message("\nCustom_error on package:\n",
                pkg,
                "\n-----------\n",
                cond)
        easypackages::packages(pkg, prompt = F, Ncpus = 6)
      },
      finally = {
        ld(pkg)
      }
    )
  } else {
    ld(pkg)
  }
}

# FUN: Installs package from github link
gh.inst <- function(pkNam, ...) {
  msgCheck(basename(pkNam[[1]]))
  if (rq(names(pkNam))) {
    gh.opt = list(pkNam[[1]], dependencies = T, upgrade = "never", ...)
    tryCatch(
      expr = {
        params <- c(gh.opt, list(type = "source"))
        do.call(devtools::install_github, params)
      },
      error = function(cond) {
        params  <- c(gh.opt, list(type = typeSrc))
        do.call(devtools::install_github, params)
      },
      finally = {
        ld(names(pkNam))}
    )} else {ld(names(pkNam))}
}

# FUN: Installs package from Source link
source_inst <- function(fn) {
  msgCheck(names(fn))
  if (rq(names(fn))) {
    params <- list(pkgs = fn[[1]], repos = NULL, type = "source", dependencies = T)
    tryCatch(
      expr = {
        do.call(install.packages, c(params, Ncpus = 6))
      },
      error = function(cond) { msgCheck(cond,T)
        do.call(install.packages, c(params, method = "libcurl"))
      },
      warning = function(cond) { msgCheck(cond,T)
        do.call(install.packages, c(params, method = "auto"))
      },
      finally = {ld(names(fn))}
    )} else {ld(names(fn))}
}

# FUN: Installs package from Bioconductor
bc.inst <- function(pknm) {
  message("Checking ", pknm, "...")
  if (rq(pknm)) {
    bio.opt <- list(
      pkgs = pknm,
      update = F,
      ask = F,
      dependencies = T
    )
    tryCatch(
      expr = {
        do.call(BiocManager::install, c(bio.opt))
      },
      error = function(cond) {
        message("Error with Package install:\n", pknm)
        message(cond)
      },
      finally = {
        ld(pknm)
      }
    )
  } else {
    ld(pknm)
  }
}

# FUN: Installs and updates autoloading of the needs package
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
      install.packages("needs", dependencies = TRUE, verbose = TRUE,
                       Ncpus = 4, ask = FALSE, type = "binary"),
      error = function(cond) {
        devtools::install_github(
          "joshkatz/needs",
          ref = "development",
          dependencies = TRUE,
          verbose = TRUE,
          upgrade = "always"
        )
      }
    )
  }
  
  if (!("needs" %in% rownames(installed.packages()))) {
    installNeeds()
    configureNeeds()
    updateRProfile()
  }
  
}

# FUN: ensures the main package installers are loaded
loadMainPkgs <- function() {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", dependencies = TRUE, ask = F)
  }
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", dependencies = TRUE, ask = F)
  }
  ld("BiocManager")
  if (!requireNamespace("zip", quietly = TRUE)) {
    install.packages("zip", dependencies = TRUE, ask = F)
  }
  ld("zip")
  if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak", type = "binary", ask = F, dependencies = T)
  }
  
  manageNeeds()
  options(needs.promptUser = FALSE)
  options(promptUser = FALSE)
  options(device.ask.default = FALSE)
  options(install.packages.check.source = "yes")
  
  gh.inst(data.frame(easypackages = "jakesherman/easypackages"))
  
  if (!requireNamespace("tidyverse", quietly = T)) {
    pak::pkg_install("tidyverse", ask = F)
  }
}

# FUN: Applies each function install type (Github, CRAN, BioConductor) to Pkgs
installAll <- function(pkgList, instFun) {
  invisible(lapply(X = 1:length(pkgList), function(X) {
    instFun(pkgList[X])
  }))
}

# FUN: Ensures the correct minfi and manifest are installed for the classifier
install_minfi <- function() {
  devtools::install_github(
    "mwsill/minfi",
    upgrade = "never",
    force = T,
    dependencies = T
  )
  devtools::install_github(
    "mwsill/IlluminaHumanMethylationEPICv2manifest",
    upgrade = "always",
    force = T,
    dependencies = T
  )
}

# FUN: Ensures if the version of minfi installed is correct for the classifier
check_minfi <- function(){
  minfiVers <- as.character(utils::packageVersion("minfi"))
  if (minfiVers != "1.43.1") {
    if (isNamespaceLoaded("minfi")) {
      try(unloadNamespace("minfi"), silent = TRUE)
    }
    install_minfi()
  }
}

# FUN: Checks if the package in the input list is already installed
checkNeeded <- function(pkgList) {
  installedPkgs <- rownames(installed.packages())
  neededPkgs <- pkgList[!pkgList %in% installedPkgs]
  return(neededPkgs)
}

# FUN: Loads the list of packages from CRAN, GitHub, and BioConductor needed
loadPacks <- function() {
  loadMainPkgs()
  neededPkgs <- checkNeeded(cranPkgs)
  
  if (length(neededPkgs) > 0) {
    pak::pkg_install(neededPkgs, ask = F)
  }
  
  if (rq("minfi")) {
    install_minfi()
  }
  sup(installAll(gHubPkgs, gh.inst))
  
  neededPkgs <- checkNeeded(biocPkgs)
  if (length(neededPkgs) > 0) {
    sup(installAll(neededPkgs, bc.inst))
  }
  check_minfi()
  if (rq("MethylAid")) {
    BiocManager::install("MethylAid", dependencies = T, ask = F, update = F)
  }
  if (rq("librarian")) {
    install.packages("librarian", dependencies = T, ask = F)
  }
  neededPkgs <- checkNeeded(extraPks)
  if (length(neededPkgs) > 0) {
    pak::pkg_install(neededPkgs, ask = F)
  }
  library('grid')
  library("ggplot2")
  library("pals")
  library("scales")
  library("stringr")
  
}

# FUN: Ensures the network drives are mounted
checkMounts <- function() {
  paths <- c(
    "/Volumes/CBioinformatics/" = "CBioinformatics",
    "/Volumes/molecular/Molecular" = "Molecular"
  )
  
  message("You have the following drives mounted:")
  system("ls /Volumes")
  
  missingPaths <- names(paths)[!dir.exists(names(paths))]
  if (length(missingPaths) > 0) {
    message("You do not have the following drives mounted:")
    message(paste(missingPaths, collapse = "\n"))
  }
}

# FUN: Ensures the latest classifier version is installed
checkEpicV2 <- function(pkg, epicV2script) {
  if (pkg %in% rownames(installed.packages())) {
    message(paste("Package", pkg, "is installed with version:",
                  utils::packageVersion(pkg)))
  } else{
    message("Installing package from source:\n", epicV2script)
    source(epicV2script)
  }
}


# FUN: Checks if the package UniD and it's requirements are installed
check_uniD_pkg <- function() {
  uniDpkgs <- c("ade4", "methylumi", "mlr")
  librarian::shelf(uniDpkgs, ask = F, update_all = F, quiet = F)
  try(bc.inst("impute"), silent = T)
  try(bc.inst("wateRmelon"), silent = T)
  if (rq("UniD")) {
    try(install.packages(uniD_path, type = "source", dependencies = T,
                         repo = NULL), silent = T)
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
  installAll(classPacks, source_inst)
  if (isMacOS) {
    check_uniD_pkg()
    checkEpicV2("mnp.v12epicv2", epicV2script)
  }
}

startLoadingAll(isMacOS, typeSrc, epicV2script)
