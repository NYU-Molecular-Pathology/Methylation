#!/usr/bin/env Rscript
## Script name: SetRunParams
## Purpose: source of global variables and install needed packages
## Date Created: September 20, 2021
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025

options(stringsAsFactors = FALSE)
gb <- globalenv()
assign("gb", gb)
library("base")
setRunLnk <- "https://github.com/NYU-Molecular-Pathology/Methylation/edit/main/SetRunParams.R"


msgFunName <- function(pthLnk, funNam) {
    message("Executing function: ", crayon::black$bgYellow(funNam),
            " from RScript in:\n", pthLnk, "\n")
}


not_installed <- function(pkgName) {
    return(!pkgName %in% rownames(installed.packages()))
}

load_pkgs <- function(pkgs){
    suppressPackageStartupMessages(
        sapply(pkgs, library, character.only = T, logical.return = T, quietly = T)
    )
}

if (not_installed("pak")) {
    install.packages("pak", ask = F, type = "binary", dependencies = T)
}

methylQCpkgs <- c(
    "kableExtra", "magick", "webshot", "plyr", "ggplot2", "knitr", "reshape2",
    "data.table", "DT", "plotly", "MethylAid", "minfi", "scales",
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
    "IlluminaHumanMethylationEPICmanifest",
    "IlluminaHumanMethylationEPICv2manifest",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    "Biobase", "RColorBrewer", "limma", "ggfortify", "Rtsne",
    "qdapTools", "gplots", "readxl", "stringr", "ggrepel", "Polychrome",
    "tinytex", "gridExtra", "rmarkdown", "BiocParallel", "grid", "pals"
)

reqPkg <- c("ggplot2", "pals", "stringr", "scales", "grid")


# FUNC: Installs version-specific minfi
minfi_install <- function() {
    Sys.setenv(R_COMPILE_AND_INSTALL_PACKAGES = "always")
    devtools::install_github(
        "mwsill/minfi",
        upgrade = "always",
        force = T, dependencies = T, type = "source", auth_token = NULL
    )
    if (!utils::packageVersion("IlluminaHumanMethylationEPICv2manifest") == "0.1.0") {
        try(unloadNamespace("IlluminaHumanMethylationEPICv2manifest"), silent = TRUE)
        try(remove.packages("IlluminaHumanMethylationEPICv2manifest"), silent = TRUE)
        devtools::install_github(
            "mwsill/IlluminaHumanMethylationEPICv2manifest",
            upgrade = "always",
            force = T, dependencies = T, type = "source", auth_token = NULL
        )
    }
}


# FUNC: Ensures the required packages for the QC report are installed/loaded
check_QC_pkgs <- function() {
    msgFunName(setRunLnk, "check_QC_pkgs")

    needed_pkgs <- sapply(methylQCpkgs, not_installed)

    if (any(needed_pkgs)) {
        missing_pkgs <- methylQCpkgs[needed_pkgs]
        library("pak")
        pak::pkg_install(missing_pkgs, ask = F, upgrade = F)
    }

    not_loaded <- !methylQCpkgs %in% loadedNamespaces()
    if (any(not_loaded)) {
        load_pkgs(methylQCpkgs[not_loaded])
    }

    if (!utils::packageVersion("minfi") == "1.43.1") {
        try(unloadNamespace("minfi"), silent = TRUE)
        try(remove.packages("minfi"), silent = TRUE)
        minfi_install()
    }


    load_pkgs(reqPkg)
}


# Helper functions to get and set global variables
setVar <- function(valueName, val) {
    msgFunName(setRunLnk, "setVar")
    return(assign(valueName, val, envir = .GlobalEnv))
}


assignVar <- function(varStr, assignedVal) {
    msgFunName(setRunLnk, "assignVar")
    return(
        tryCatch(
            expr = {
                if (!is.null(base::get(varStr))) {
                    cat("")
                }
            },
            error = {
                setVar(varStr, assignedVal)
            },
            finally = message(varStr, " = ", assignedVal)
        )
    )
}

# Checks if variable is null and assigns value
ckNull <- function(nullVar, subVar, varName) {
    msgFunName(setRunLnk, "ckNull")
    if (is.null(nullVar)) {
        setVar(as.character(varName), as.character(subVar))
        return(paste0(subVar))
    } else {
        return(paste0(nullVar))
    }
}

# FUNC: Loads the QC package dependencies checks if any are not installed
loadClassifierPacks <- function() {
    msgFunName(setRunLnk, "loadClassifierPacks")
    ldPkg <- tryCatch(
        expr = {
            check_QC_pkgs()
        },
        error = function(cond) {
            message("error in loading QC package dependency:\n")
            message(cond)
        },
        warning = function(cond) {
            message("Warning in loading QC package dependency:\n")
            message(cond)
        }
    )
    return(ldPkg)
}

# Sets default variable paths/names
getDefaults <- function() {
    msgFunName(setRunLnk, "getDefaults")
    cbVol <- switch(Sys.info()[["sysname"]],
        "Darwin" = "/Volumes/CBioinformatics/Methylation",
        "Linux" = "~/molecpathlab/production/Methylation"
    )
    moVol <- "/Volumes/molecular"
    rsVol <- "/Volumes/snudem01labspace"
    if (is.null(gb$ApiToken)) {
        gb$ApiToken <- gb$token
    }
    defaultParams <- data.frame(
        mnp.pk.loc = paste0(file.path(cbVol, "classifiers/mnp.v12epicv2")),
        ApiToken = gb$ApiToken,
        methDir = paste0(file.path(cbVol, "Clinical_Runs")),
        clinDrv = paste0(file.path(moVol, "MOLECULAR LAB ONLY/NYU-METHYLATION")),
        rschOut = paste0(file.path(rsVol, "FINAL_PDF_Reports_Brain")),
        clinOut = paste0(file.path(moVol, "MOLECULAR/MethylationClassifier")),
        rsch.idat = paste0(file.path(rsVol, "idats")),
        clin.idat = paste0(file.path(moVol, "MOLECULAR/iScan")),
        QC_file = paste0("~/Methyl_QC.Rmd"),
        baseDir = paste0(file.path(cbVol, "Clinical_Runs")),
        stringsAsFactors = F
    )
    return(defaultParams)
}

# FUNC: Returns pre-defined global variable values used in the pipeline
getSetvars <- function() {
    msgFunName(setRunLnk, "getSetvars")

    assignedVars <- data.frame(
        mnp.pk.loc = gb$mnp.pk.loc,
        ApiToken = gb$ApiToken,
        methDir = gb$methDir,
        clinDrv = gb$clinDrv,
        rschOut = gb$rschOut,
        clinOut = gb$clinOut,
        rsch.idat = gb$rsch.idat,
        clin.idat = gb$clin.idat,
        QC_file = gb$QC_file,
        baseDir = gb$baseDir,
        stringsAsFactors = F
    )
    return(assignedVars)
}

# Sets the default parameters for a methylation Run
defineParams <- function(mnp.pk.loc = NULL,
                         ApiToken = NULL,
                         methDir = NULL,
                         clinDrv = NULL,
                         rschOut = NULL,
                         clinOut = NULL,
                         rsch.idat = NULL,
                         clin.idat = NULL,
                         QC_file = NULL,
                         isMC = T,
                         baseDir = NULL,
                         runID = NULL,
                         loadClassifier = T) {
    msgFunName(setRunLnk, "defineParams")
    defVars <- getDefaults()
    inVars <- list(
        mnp.pk.loc,
        ApiToken,
        methDir,
        clinDrv,
        rschOut,
        clinOut,
        rsch.idat,
        clin.idat,
        QC_file,
        baseDir
    )
    i <- 1:length(inVars)
    for (x in 1:length(defVars)) {
        if (!is.null(inVars[[x]])) {
            defVars[1, x] <- inVars[[x]]
        }
    }
    invisible(lapply(i, function(x) {
        if (!is.null(inVars[[x]])) {
            setVar(names(defVars[x]), inVars[[x]])
        }
    }))
    if (loadClassifier == T) {
        loadClassifierPacks()
    }
    i <- 1:length(defVars)
    invisible(sapply(i, FUN = function(i) {
        assignVar(names((defVars[i])), paste0(defVars[, i]))
    }))
    cbVol <- switch(Sys.info()[["sysname"]],
        "Darwin" = "/Volumes/CBioinformatics/",
        "Linux" = "~/molecpathlab/production/"
    )
    if (!isMC) {
        methDir <- rschOut
        assign("workDir", cbVol)
    }
}

# Changes the working directory using the system CD command
setDirectory <- function(foldr) {
    msgFunName(setRunLnk, "setDirectory")
    bsDir <- paste("cd", foldr)
    if (dir.exists(foldr)) {
        system(bsDir)
        setwd(foldr)
        assign("workDir", foldr)
    } else {
        mm2 <- crayon::white$bgRed("Location Not Found:", foldr)
        warning(mm2)
    }
}

# Validate Input Parameters ---------------------------------------------------
CheckInputArg <- function(varValue, gb, defVal = NULL) {
    varStr <- deparse(substitute(varValue))
    if (length(varValue) == 0 | identical(varValue, NULL) | identical(varValue, "NULL")) {
        gb[[varStr]] <- varValue <- defVal
    } else {
        varValue <- ifelse(is.na(varValue), NULL, varValue)
    }
    if (varStr == "token") {
        message("\n~~~~~~~~~~~~~~~~~~~~~Parameters input~~~~~~~~~~~~~~~~~~~~~")
    }
    message(varStr, ": ", ifelse(is.null(varValue), "NULL", varValue))
    return(assign(varStr, varValue, envir = gb))
}


# Verifies the work directory exists
CheckBaseFolderInput <- function(baseFolder) {
    if (!is.null(baseFolder) & !identical(baseFolder, "NULL")) {
        stopifnot(
            "Input directory does not exist!" = dir.exists(baseFolder) == TRUE
            )
    } else {
        baseFolder <- NULL
    }
    return(baseFolder)
}
