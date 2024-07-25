#!/usr/bin/env Rscript
## ---------------------------
## Script name: SetRunParams
## Purpose: source of global variables and install needed packages
## Date Created: September 20, 2021
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

options(stringsAsFactors = FALSE)
gb <- globalenv(); assign("gb", gb)
library("base")
setRunLnk = "https://github.com/NYU-Molecular-Pathology/Methylation/edit/main/SetRunParams.R"

msgFunName <- function(pthLnk, funNam){
    message("Executing function: ", crayon::black$bgYellow(funNam), " from RScript in:\n", pthLnk,"\n")
}

SpSm <- function(pkg){return(suppressPackageStartupMessages(pkg))}

fix_compiler_flags <- function(){
    # Check if brew installed
    if (is.na(Sys.which("brew")["brew"][[1]])) {
        message("Homebrew is not installed. Installing Homebrew...")
        cmd <-
            '/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"'
        system(cmd, wait = TRUE)
    }
    # Check if LLVM installed
    llvm_installed <- file.exists("/usr/local/opt/llvm/bin/clang")
    if (!llvm_installed) {
        message("LLVM is not installed. Installing LLVM via Homebrew...")
        system("brew install llvm", intern = TRUE, wait = TRUE)
    }
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
    Sys.setenv(PATH = paste("/usr/local/opt/llvm/bin", Sys.getenv("PATH"), sep = ":"))
    Sys.setenv(CC = "/usr/local/opt/llvm/bin/clang")
    Sys.setenv(CXX = "/usr/local/opt/llvm/bin/clang++")
    Sys.setenv(CXX11 = "/usr/local/opt/llvm/bin/clang++")
    Sys.setenv(CXX14 = "/usr/local/opt/llvm/bin/clang++")
    Sys.setenv(CXX17 = "/usr/local/opt/llvm/bin/clang++")
    Sys.setenv(CXX1X = "/usr/local/opt/llvm/bin/clang++")
    Sys.setenv(OBJC = "/usr/local/opt/llvm/bin/clang")
    Sys.setenv(LDFLAGS = "-L/usr/local/opt/llvm/lib")
    Sys.setenv(CPPFLAGS = "-I/usr/local/opt/llvm/include")
}


# FUN: Loads all the packages used in the RMD Methylation QC file
checkQCpkg <- function(){
    msgFunName(setRunLnk,"checkQCpkg")
    methylQCpacks <- c(
        "kableExtra","magick","webshot","plyr","ggplot2","knitr","reshape2",
        "data.table","DT","plotly", "MethylAid","minfi","scales",
        "IlluminaHumanMethylation450kmanifest",
        "IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
        "IlluminaHumanMethylationEPICmanifest",
        "IlluminaHumanMethylationEPICv2manifest",
        "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
        "Biobase", "RColorBrewer", "limma","ggfortify","Rtsne",
        "qdapTools","gplots","readxl","stringr","ggrepel","Polychrome",
        "tinytex","gridExtra","rmarkdown", "BiocParallel", "grid"
    )
    
    #fix_compiler_flags()
    message("Now Loading...\n", paste(methylQCpacks, collapse=" "),"\n")
    
    needed_pkgs <- pkgs[!sapply(methylQCpacks, requireNamespace, quietly = TRUE)]
    params <- list(
        dependencies = TRUE,
        ask = FALSE,
        update = "never",
        quiet = TRUE,
        repos = 'http://cran.us.r-project.org',
        Ncpus = 4
    )
    # if (length(needed_pkgs) > 0) {
    #     do.call(install.packages, c(list(pkgs = needed_pkgs), params))
    # }
    
    reqPkg <- list("ggplot2", "pals", "stringr", "scales", "grid")
    needed_pkgs <- pkgs[!sapply(reqPkg, requireNamespace, quietly = TRUE)]
    # if (length(needed_pkgs) > 0) {
    #     do.call(install.packages, c(list(pkgs = needed_pkgs), params))
    # }
    sapply(pkgs, library, character.only = T, logical.return = T, quietly = T)
    sapply(reqPkg, library, character.only = T, logical.return = T, quietly = T)
}

# Helper functions to get and set global variables
setVar <- function(valueName,val){
    msgFunName(setRunLnk,"setVar")
    return(assign(valueName, val, envir=.GlobalEnv))
}

assignVar <- function(varStr, assignedVal){
    msgFunName(setRunLnk,"assignVar")
    return(
        tryCatch(
            expr = {if(!is.null(base::get(varStr))){cat("")}},
            error = {setVar(varStr,assignedVal)},
            finally = message(varStr," = ",assignedVal)
            )
        )
}

# Checks if variable is null and assigns value
ckNull <- function(nullVar, subVar, varName){
    msgFunName(setRunLnk,"ckNull")
    if (is.null(nullVar)){
        setVar(as.character(varName),as.character(subVar))
        return(paste0(subVar))} else {return(paste0(nullVar))}
}

# FUN: Loads the main packages and dependencies checks if any are not installed
loadClassifierPacks <- function(){
    msgFunName(setRunLnk,"loadClassifierPacks")
    ldPkg <- tryCatch(
        expr={checkQCpkg()},
        error=function(cond){message("error in loading QC package dependency:\n")
                             message(cond)},
        warning=function(cond){message("Warning in loading QC package dependency:\n")
                               message(cond)}
    )
    return(ldPkg)
}

# Sets default variable paths/names
getDefaults <- function() {
    msgFunName(setRunLnk,"getDefaults")
    cbVol <- switch (Sys.info()[['sysname']],
                  "Darwin" = "/Volumes/CBioinformatics/Methylation",
                  "Linux" = "~/molecpathlab/production/Methylation"
    )
    moVol = "/Volumes/molecular"
    rsVol = "/Volumes/snudem01labspace"
    if(is.null(gb$ApiToken)){
       gb$ApiToken <- gb$token
    }
    defaultParams <- data.frame(
        mnp.pk.loc = paste0(file.path(cbVol, "classifiers/mnp.v12epicv2")),
        ApiToken = gb$ApiToken,
        methDir = paste0(file.path(cbVol, "Clinical_Runs")),
        clinDrv = paste0(file.path(moVol, "MOLECULAR LAB ONLY/NYU-METHYLATION")),
        rschOut = paste0(file.path(rsVol, "FINAL_PDF_Reports_Brain")),
        clinOut = paste0(file.path(moVol, "MOLECULAR/MethylationClassifier")),
        rsch.idat = paste0(file.path(rsVol,"idats")),
        clin.idat = paste0(file.path(moVol, "MOLECULAR/iScan")),
        QC_file = paste0("~/Methyl_QC.Rmd"),
        baseDir = paste0(file.path(cbVol, "Clinical_Runs")),
        stringsAsFactors=F
    )
    return(defaultParams)
}

# returns assigned global variable values
getSetvars <- function() {
    msgFunName(setRunLnk,"getSetvars")

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
    i=1:length(inVars)
    for (x in 1:length(defVars)) {
        if(!is.null(inVars[[x]])){
            defVars[1,x] <- inVars[[x]]
        }
    }
    invisible(lapply(i,function(x){if(!is.null(inVars[[x]])){setVar(names(defVars[x]), inVars[[x]])}}))
    if(loadClassifier==T){loadClassifierPacks()}
    i=1:length(defVars)
    invisible(sapply(i,FUN=function(i){assignVar(names((defVars[i])), paste0(defVars[,i]))}))
    cbVol <- switch (Sys.info()[['sysname']],
                  "Darwin" = "/Volumes/CBioinformatics/",
                  "Linux" = "~/molecpathlab/production/"
    )
    if(!isMC) {methDir=rschOut;assign("workDir", cbVol)}
}

# Changes the working directory using the system CD command
setDirectory <- function(foldr) {
    msgFunName(setRunLnk,"setDirectory")
    bsDir = paste("cd", foldr)
    mm2 = crayon::white$bgRed("Location Not Found:", foldr)
    if (dir.exists(foldr)) {
        system(bsDir)
        setwd(foldr)
        assign("workDir", foldr)} else{warning(mm2)}
}

# Check Input Parameters -----------------------------------------------------------------------
CheckInputArg <- function(varValue, gb, defVal = NULL) {varStr <- deparse(substitute(varValue))
    if (length(varValue) == 0 | identical(varValue, NULL) | identical(varValue, "NULL")) {gb[[varStr]] <- varValue <- defVal} else{ varValue <- ifelse(is.na(varValue), NULL, varValue)}
    if(varStr=="token"){message("\n~~~~~~~~~~~~~~~~~~~~~Parameters input~~~~~~~~~~~~~~~~~~~~~")}
    message(varStr, ": " , ifelse(is.null(varValue), "NULL", varValue)); return(assign(varStr, varValue, envir = gb))
}


# LoadAndUnloadPacks <- function(){
#     unloadNamespace("mnp.v11b4"); unloadNamespace("mnp.v12b6"); unloadNamespace("sarc.v12b6")
#     loadNamespace("mnp.v11b6"); require("mnp.v11b6"); library("mnp.v11b6")
# }

CheckBaseFolderInput <- function(baseFolder){
    if(!is.null(baseFolder) & !identical(baseFolder, "NULL")) {
        message("Checking if custom run directory is valid: ", baseFolder, "\n")
        stopifnot("Input directory does not exist! Create it with mkdir" = dir.exists(baseFolder) == T)
    } else{baseFolder <- NULL}
    #LoadAndUnloadPacks()
    return(baseFolder)
}

