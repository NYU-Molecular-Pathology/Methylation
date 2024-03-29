#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
gb <- globalenv(); assign("gb", gb)
library("base")
setRunLnk = "https://github.com/NYU-Molecular-Pathology/Methylation/edit/main/SetRunParams.R"

msgFunName <- function(pthLnk, funNam){
    message("Executing function: ", crayon::black$bgYellow(funNam), " from RScript in:\n", pthLnk,"\n")
}

SpSm <- function(pkg){return(suppressPackageStartupMessages(pkg))}

# FUN: Loads all the packages used in the RMD Methylation QC file
checkQCpkg <- function(){
    msgFunName(setRunLnk,"checkQCpkg")
    methylQCpacks <- c(
        "kableExtra","magick","webshot","plyr","ggplot2","knitr","reshape2",
        "data.table","DT","plotly", "MethylAid","minfi","scales",
        "IlluminaHumanMethylation450kmanifest",
        "IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
        "IlluminaHumanMethylationEPICmanifest",
        "IlluminaHumanMethylationEPICanno.ilm10b4.hg19","Biobase", "RColorBrewer",
        "limma","ggfortify","Rtsne",
        "qdapTools","gplots","readxl","stringr","ggrepel","Polychrome",
        "tinytex","gridExtra","rmarkdown", "BiocParallel", "grid"
    )
    message("Now Loading...\n", paste(methylQCpacks, collapse=" "),"\n")
    SpSm(easypackages::packages(methylQCpacks,prompt=F))
    reqPkg <- list("ggplot2","pals","stringr","scales","grid")
    invisible(lapply(reqPkg, FUN = function(X) {SpSm(do.call("require", list(X)))}))
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
        mnp.pk.loc = paste0(file.path(cbVol, "classifiers/mnp.v11b6")),
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


LoadAndUnloadPacks <- function(){
    unloadNamespace("mnp.v11b4"); unloadNamespace("mnp.v12b6"); unloadNamespace("sarc.v12b6")
    loadNamespace("mnp.v11b6"); require("mnp.v11b6"); library("mnp.v11b6")
}

CheckBaseFolderInput <- function(baseFolder){
    if(!is.null(baseFolder) & !identical(baseFolder, "NULL")) {
        message("Checking if custom run directory is valid: ", baseFolder, "\n")
        stopifnot("Input directory does not exist! Create it with mkdir" = dir.exists(baseFolder) == T)
    } else{baseFolder <- NULL}
    LoadAndUnloadPacks()
    return(baseFolder)
}

