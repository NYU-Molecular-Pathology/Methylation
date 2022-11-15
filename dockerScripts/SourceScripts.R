#!/usr/bin/env Rscript
if(!require(compiler)){install.packages("compiler", dependencies=T)}
library(compiler)
enableJIT(3)

cbioLn <- switch (Sys.info()[['sysname']],
                  "Darwin" = "/Volumes/CBioinformatics/Methylation/classifiers",
                  "Linux" = "./private_pkgs")

uniDpath <- file.path(cbioLn, "UniD")

# Classifier Packages and Versions

# Custom Classifier packages:
classPacks <- c(
    sest = "https://github.com/jungch/sest/raw/master/sest.tar",
    mgmtstp27 = "https://github.com/badozor/mgmtstp27/raw/master/archive/mgmtstp27_0.6-3.tar.gz",
    mnpqc = paste0(file.path(cbioLn, "mnpqc_0.1.0.tar.gz"))
)

# Cran Packages ----
cranPkgs <-
    c(
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
        "Rcpp", "httpuv", "shiny"
    )

# GitHub Packages ----
gHubPkgs <-
    data.frame(
        rmarkdown = 'rstudio/rmarkdown',
        docstring = 'dasonk/docstring',
        rstudioapi = 'rstudio/rstudioapi',
        easypackages = 'jakesherman/easypackages',
        redcapAPI = 'nutterb/redcapAPI',
        crayon = "r-lib/crayon",
        redcap = "epicentre-msf/redcap"
    )

# BioConductor Packages ----
biocPkgs <-
    c(
        'HDF5Array', 'rngtools', 'bumphunter','GEOquery', 'minfi', 'lumi', 'methylumi',
        'randomForest', 'glmnet','IlluminaHumanMethylation450kmanifest',
        'IlluminaHumanMethylation450kanno.ilmn12.hg19', 'Rtsne',
        'IlluminaHumanMethylationEPICmanifest',
        'IlluminaHumanMethylationEPICanno.ilm10b2.hg19',
        'IlluminaHumanMethylationEPICanno.ilm10b4.hg19',
        'MethylAid', 'conumee','BiocParallel', "Biobase","limma","MatrixGenerics", "GenomeInfoDb"
    )

# Extra Libraries ----
cpuPacks <- c("parallel","doSNOW","doParallel", "foreach","compiler")
easyPkgs <- c('tidyverse','sjmisc','stringi','digest','RCurl','gridExtra','needs')

# Helper Functions ----
sup <- function(x){return(suppressWarnings(suppressPackageStartupMessages(x)))}
sw <- function(pkgOb){try(return(suppressMessages(suppressWarnings(pkgOb))),silent=T)}
ld <- function(libName) {
    lib.opts <- list(
        package = libName, character.only = T,verbose = T, warn.conflicts = F, quietly = F)
    suppressPackageStartupMessages(do.call(library, c(lib.opts)))
    message(libName, " ...load successful")
}
up <- function(){update.packages(repos='http://cran.rstudio.com/', type = "source", ask=F, checkBuilt=T)}
rq <- function(pkgName){ifelse(pkgName %in% row.names(installed.packages()), F, T)}
mkred <- function(strMsg) {return(crayon::white$bgRed$bold(strMsg))}
mkblu <- function(strMsg) {return(crayon::white$bgBlue$bold(strMsg))}
mkGrn <- function(strMsg) {return(crayon::white$bgGreen$bold(strMsg))}

msgCheck <- function(pkg, warn = F) {
    if (warn == F) {message("Checking ", pkg, "...")} else{
        message("\nError caught for package:\n", pkg, "\n-----------")
    }
}

#' pk.inst will load a package library or install it if it does not exist
#' @param pkg a character string package name
#' @return Message if loading the library was successful
pk.inst <- function(pkg){
    msgCheck(pkg)
    if(rq(pkg)) {
        pk.opt <- list(pkgs=pkg,dependencies=T,verbose=T,Ncpus = 6)
        tryCatch(
            expr = {
                do.call(install.packages, c(pk.opt, list(type = "both")))
            },
            warning = function(cond) {
                message("\nWarning on package:\n", pkg, "\n-----------\n",cond)
                do.call(install.packages, c(pk.opt, list(type = "binary")))
            },
            error = function(cond) {
                msgCheck(cond,T)
                do.call(install.packages, c(pk.opt, list(type = "source")))
            },
            custom_error = function(cond) {
                msgCheck(cond,T)
                message("\nCustom_error on package:\n", pkg, "\n-----------\n",cond)
                easypackages::packages(pkg, prompt = F, Ncpus = 6)
            },
            finally={ld(pkg)}
        )
    } else {ld(pkg)}
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
                params  <- c(gh.opt, list(type = "binary"))
                do.call(devtools::install_github, params)
            },
            finally={ld(names(pkNam))}
        )} else {ld(names(pkNam))}
}

# FUN: Installs package from Source link
srcInst <- function(fn){
    msgCheck(names(fn))
    if (rq(names(fn))) {
        params <- list(pkgs=fn[[1]], repos = NULL, type = "source", dependencies=T)
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
bc.inst <- function(pknm){
    message("Checking ", pknm, "...")
    if(rq(pknm)) {
        bio.opt <- list(pkgs=pknm, update=F, ask=F)
        tryCatch(
            expr = {do.call(BiocManager::install, c(bio.opt))},
            error = function(cond) {message("Package already loaded")},
            finally = {ld(pknm)}
        )
    } else {ld(pknm)}
}

fixProf <- function(){
    txt1 <- "^[:blank:]*autoload\\(\"needs\", \"needs\"\\)"
    txt2 <- "\n\nautoload(\"needs\", \"needs\")\n\n"
    siteProf <- if(is.na(Sys.getenv("R_PROFILE", unset=NA))) {
        file.path(Sys.getenv("R_HOME"),"etc","Rprofile.site")
    } else {Sys.getenv("R_PROFILE")}
    if (!file.exists(siteProf)) {file.create(siteProf)}
    cxn <- file(siteProf); lines <- readLines(cxn)
    if (!any(grepl(txt1,lines))) {write(txt2, file=siteProf, append=T)}
    close(cxn)
}

fixNeeds <- function (){
    sysfile <- system.file("extdata", "promptUser", package="needs")
    write(0, file=sysfile); options(needs.promptUser=FALSE); invisible(needs:::autoload(TRUE))
}

checkNeeds <- function(){
    tryCatch(
        expr={
            if(!("needs" %in% rownames(installed.packages()))){
                install.packages("needs",dependencies=T,verbose=T, Ncpus = 6)
                fixNeeds();fixProf()
            }else{
                fixNeeds();fixProf()
            }
        },
        error=function(cond){
            devtools::install_github("joshkatz/needs", ref = "development",
                                     dependencies=T,verbose=T,upgrade="always")
            fixNeeds();fixProf()
        },
        warning=function(cond){
            devtools::install_github("joshkatz/needs", ref = "development",
                                     dependencies=T,verbose=T,upgrade="always")
            fixNeeds();fixProf()
        }
    )
}

setOptions <- function(){
    options(needs.promptUser=F)
    #Sys.setenv(R_ENABLE_JIT=T)
    options("needs.promptUser"=F); options("promptUser"=F)
    options("device.ask.default" = FALSE)
    options("install.packages.compile.from.source"="Yes")
    options("install.packages.check.source"="yes")
    #compiler::enableJIT(3); 
    compiler::compilePKGS(enable=T);
    compiler::setCompilerOptions(suppressAll=T, optimize=3)
    fixProf()
}

cpMakeV <- function(editFile=F){
    #mkFi <- "/Volumes/CBioinformatics/jonathan/Rprojects/Methylation-scripts/Makevars"
    #message("Copying Makevars to local dotfiles: ", mkFi)
    #message("Copying Makevars to local dotfiles")
    if(!dir.exists("~/.R/")){
        message("No Makevars file in ~/.R/")
        #dir.create("~/.R/")
    }else{
        system("rm -R ~/.R/")
    }
    #if(editFile){file.edit(mkFi)}
    #if(file.exists(mkFi)){file.copy(mkFi,"~/.R/", overwrite=T )}
}

# FUN: Load Classifier Versions
classifierInstall <- function(pathtoFile=NULL, instNew=T, rmpkg=F) {
    mnp.pk.loc <- file.path(cbioLn, pathtoFile)
    if(rmpkg){remove.packages(basename(mnp.pk.loc), lib=.libPaths()[[1]][1])}
    if(instNew){
        message("Installing package: ", basename(mnp.pk.loc))
        message("Classifier install may take several minutes!  Note the stage '** inst' may take the longest.")
        install.packages(mnp.pk.loc, repos=NULL, type="source", force=T, Ncpus = 6)
        classifierInstall(pathtoFile,F,F)
    } else {ld(basename(mnp.pk.loc))}
}

checkBioC <- function(){
    if(sw(suppressPackageStartupMessages((!require("BiocManager"))))){
        install.packages("BiocManager", Ncpus = 4)
        BiocManager::install(version="3.10", update=T, ask=F, type="source")
    } else{ld("BiocManager")}
    if(rq("zip")){install.packages("zip", dependencies=T, type="binary")
    } else{ld("zip")}
}

loadMainPkgs <- function(){
    setOptions()
    if(rq("devtools")){install.packages("devtools",dependencies=T,verbose=T)}else{ld("devtools")}
    checkNeeds()
    gh.inst(data.frame(easypackages="jakesherman/easypackages"))
    pk.inst("tidyverse")
    suppressPackageStartupMessages(easypackages::packages(cpuPacks, prompt=F))
    checkBioC()
}

readyPkgs <- function(ezLibs){
    ready <- ezLibs %in% sessioninfo::package_info()[["package"]]
    if(any(!ready)){lapply(ezLibs[!ready], require, character.only=TRUE)}
}

installAll <- function(pkgList, instFun){
    invisible(lapply(X=1:length(pkgList), function(X){instFun(pkgList[X])}))
}

# Loads default packages or custom if input provided
#' @return package installs or loads
#' @export
loadPacks <- function(pkgs=cranPkgs, ezLibs=easyPkgs, ghPk=gHubPkgs, bcPks=biocPkgs) {
    tryCatch(
        expr = {
            loadMainPkgs()
            sup(installAll(pkgs, pk.inst))
            if(!require("BiocManager")){install.packages("BiocManager", dependencies=T, quiet=F)}
            sup(installAll(ghPk, gh.inst))
            sup(installAll(bcPks, bc.inst))
            readyPkgs(ezLibs)

            if(!require("MethylAid")){BiocManager::install("MethylAid",update=F, ask=F)}
            if(!require("librarian")){
                install.packages("librarian", dependencies=T, verbose=T, Ncpus = 4, quiet=F)
                }
            pkgs <- c(
                "knitr",
                "kableExtra",
                "magick",
                "webshot",
                "plyr",
                "ggplot2",
                "knitr",
                "reshape2",
                "data.table",
                "DT",
                "plotly",
                "MethylAid",
                "minfi",
                "scales",
                #"IlluminaHumanMethylation450kmanifest",
                #"IlluminaHumanMethylationEPICmanifest",
                #"IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
                "Biobase",
                "RColorBrewer",
                "limma",
                "ggfortify",
                "Rtsne",
                "qdapTools",
                "gplots",
                "readxl",
                "stringr",
                "ggrepel",
                "Polychrome",
                "tinytex",
                "gridExtra",
                "rmarkdown",
                "tinytex",
                "BiocParallel",
                'grid',
                'grDevices',
                "mdthemes",
                "dplyr",
                "sqldf",
                "gridExtra",
                "reshape"
            )
            librarian::shelf(pkgs, ask=F, update_all = FALSE, quiet = FALSE)
            require('grid')
            require("ggplot2")
            require("pals")
            require("scales")
            require("stringr")
            require("scales")
        },
        error = function(cond){
            message("\n~~~You encountered the following error during install:\n", cond)
            message("If there is a compile error, try running gb$fixCompiles() and then try gb$loadPacks() again")
        }
    )
}


colorMsg <- function(){
    ms1 <- paste0(mkGrn("Updating in-house classifier to current version:"),"\n",
                  mkblu("classifierInstall(pathtoFile=NULL, instNew=F, rmpkg=F)"),"\n"
                  )
    ms2 <- paste0(mkred("Classifier package is not installed installing classifier"),"\n")
    ms3 <- paste0(mkGrn("Your classifier package is up-to-date and loading"),"\n")
    return(c(ms1,ms2,ms3))
}

setEnviron <- function(){
    options(Ncpus = 6)
    library("base")
    #try(cpMakeV(),silent=T)
    rlis = getOption("repos")
    rlis["CRAN"] = "http://cran.us.r-project.org"
    options(repos = rlis)
    if(Sys.info()[['sysname']]=="Darwin"){
        Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/MacOS/pandoc")
        system("export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc")        
        system("export RSTUDIO_WHICH_R=/usr/local/bin/R")
        Sys.setenv(RSTUDIO_WHICH_R = "/usr/local/bin/R")
        }
    options("install.packages.compile.from.source" = "never")
    options("install.packages.check.source"="no")
    if(!require("devtools")){install.packages("devtools", dependencies=T)}
}

# classifier mnpPath currVers
checkClassifier <- function(mnpClass) {
    ms <- colorMsg()
    isInstalled <- mnpClass[,1] %in% rownames(installed.packages())
    if (isInstalled==F) {
        cat(ms[2])
        cat(mnpClass[, 1], sep = "\n")
        classifierInstall(mnpClass[, 2], T, F)
    } else {
        currentVers <- as.character(utils::packageVersion(mnpClass[, 1]))
        latestVers <- as.character(mnpClass[, 3])
        message("Current Version Installed: ", currentVers, "\nNewest Package Version: ",latestVers,"\n")
        if (currentVers==latestVers) {
            cat(ms[3])
            classifierInstall(mnpClass[, 2], F, F)
        } else {
            cat(ms[1])
            classifierInstall(mnpClass[, 2], T, T)
        }
    }
}

# Load all Functions ---------------------
startLoadingAll <- function() {
    setEnviron()
    loadPacks()
    installAll(classPacks, srcInst)
}

startLoadingAll()

mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/"
script.list = c("SetRunParams.R", "MakeSampleSheet.R", "CopyInputs.R", "CopyOutput.R", "pipelineHelper.R", "CustomRuns.R")
scripts = paste0(mainHub, script.list)
lapply(scripts, function(i){message("Sourcing: ", i); devtools::source_url(i)})

mainPage = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/Report-Scripts/"
scripts <- c("ClassTables.R", "MLH1_Functions.R", "PipeLineU.R", "RedcapOutput.R", "TsneFunctions.R", "cnvggplotly.R")
scripts <- paste0(mainPage, scripts)
lapply(scripts, function(i) {devtools::source_url(i); return(message("Sourced: ", i))})
