#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
gb <- globalenv(); assign("gb", gb)
setRunLnk = "https://github.com/NYU-Molecular-Pathology/Methylation/edit/main/SetRunParams.R"

msgFunName <- function(pthLnk, funNam){
    message("Executing function: ", funNam, " from RScript in:\n", pthLnk,"\n")
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
        "IlluminaHumanMethylationEPICanno.ilm10b4.hg19","Biobase", "RColorBrewer","limma","ggfortify","Rtsne",
        "qdapTools","gplots","readxl","stringr","ggrepel","Polychrome",
        "tinytex","gridExtra","rmarkdown", "BiocParallel", "grid"
    )
    message("Now Loading:\n", methylQCpacks,"\n")
    suppressPackageStartupMessages(easypackages::packages(methylQCpacks,prompt=F))
    reqPkg <- list("ggplot2","pals","stringr","scales","grid")
    invisible(lapply(reqPkg, FUN = function(X) {suppressPackageStartupMessages(do.call("require", list(X)))}))
    message("Packages loading finished\n")
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
            expr = {if(!is.null(get(varStr))){cat("")}},
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
    
    cbVol = "/Volumes/CBioinformatics"
    moVol = "/Volumes/molecular"
    rsVol = "/Volumes/snudem01labspace"
    defaultParams <- data.frame(
        mnp.pk.loc = paste0(file.path(cbVol, "Methylation/in_house/mnp.v116/mnp.v11b6")),
        ApiToken = "",
        methDir = paste0(file.path(cbVol, "Methylation/Clinical_Runs")),
        #methDir = paste0(file.path(cbVol)),
        clinDrv = paste0(file.path(moVol, "MOLECULAR LAB ONLY/NYU-METHYLATION")),
        rschOut = paste0(file.path(rsVol, "FINAL_PDF_Reports_Brain")),
        #rschOut = paste0(file.path(cbVol)),
        clinOut = paste0(file.path(moVol, "MOLECULAR/MethylationClassifier")),
        rsch.idat = paste0(file.path(rsVol,"idats")),
        clin.idat = paste0(file.path(moVol, "MOLECULAR/iScan")),
        QC_file = paste0(system.file('Methyl_QC.Rmd', package = "mnp.v11b6")),
        baseDir = paste0(file.path(cbVol, "Methylation/Clinical_Runs")),
        #baseDir = paste0(file.path(cbVol)),
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
defineParams <- function(
    mnp.pk.loc = NULL, ApiToken = NULL, methDir = NULL, clinDrv = NULL, rschOut = NULL, clinOut = NULL,
    rsch.idat = NULL, clin.idat = NULL, QC_file = NULL, isMC = T, baseDir = NULL, runID = NULL, loadClassifier=T
    ){
    msgFunName(setRunLnk,"defineParams")
    
    defVars <- getDefaults()
    inVars <- list(mnp.pk.loc, ApiToken, methDir, clinDrv, rschOut, clinOut,rsch.idat, clin.idat, QC_file, baseDir)
    i=1:length(inVars)
    invisible(lapply(i,function(x){if(!is.null(inVars[[x]])){setVar(names(defVars[x]), inVars[[x]])}}))
#    message("\nLoading Packages Silently...\n")
    if(loadClassifier==T){loadClassifierPacks()}
    i=1:length(defVars)
    invisible(sapply(i,FUN=function(i){assignVar(names((defVars[i])), paste0(defVars[,i]))}))
    if(!isMC) {methDir=rschOut;assign("workDir","/Volumes/CBioinformatics/")}
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
