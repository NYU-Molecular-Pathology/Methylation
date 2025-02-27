gb <- .GlobalEnv; assign("gb", gb)
formals(library)$quietly <- T
formals(library)$warn.conflicts <- F
formals(require)$warn.conflicts <- F
formals(install.packages)$dependencies <- T
formals(install.packages)$verbose <- T
formals(install.packages)$ask <- F 


library("grid")
if(!require("GenVisR")){BiocManager::install("GenVisR")}
library("minfi")
library("conumee")
#library("mnp.v11b6")
library("foreach")
library("plotly")
library("graphics")
if(!require("clusterProfiler")){BiocManager::install("clusterProfiler", ask=F)}
if(!require("org.Hs.eg.db")){BiocManager::install("org.Hs.eg.db", ask=F)}
if(!require("pathview")){BiocManager::install("pathview", ask=F)}
library("org.Hs.eg.db")
library("clusterProfiler")
library("minfi")
library("pathview")


cranPkgs <-
    c(
        "Rtsne",
        "parallel",
        "doParallel",
        "foreach",
        "pkgbuild",
        "ggrepel",
        "pals",
        "ggfortify",
        "plotly",
        "compiler",
        "circlize",
        "grid",
        "gridExtra",
        "gridtext",
        "feather",
        "sqldf",
        "reshape2",
        "ggplot2",
        "ggnewscale",
        "ggupset",
        "umap"
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
        "minfi",
        "conumee",
        "IlluminaHumanMethylationEPICmanifest",
        "limma",
        "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
        "ComplexHeatmap",
        "DNAcopy",
        "CopyNumberPlots"
    )

# Extra Libraries ----
cpuPacks <- c("parallel","doSNOW","doParallel", "foreach","compiler")
easyPkgs <- c('tidyverse','sjmisc','stringi','digest','RCurl','gridExtra','needs')

# Helper Functions ----
sw <- function(pkgOb){try(return(suppressMessages(suppressWarnings(pkgOb))),silent=T)}
ld <- function(libName) {
    lib.opts <- list(package = libName, character.only = T, verbose = F, warn.conflicts = F)
    sw(suppressPackageStartupMessages(do.call(library, c(lib.opts))))
    message(libName, " ...load successful")
}
up <- function(){update.packages(repos='http://cran.rstudio.com/',
                                 type = "source", ask=F, checkBuilt=T)}
rq <- function(pkgName){
    ifelse(pkgName %in% row.names(installed.packages()), F, T)
}

msgCheck <- function(pkg, warn = F) {
    if (warn == F) {
        message("Checking ", pkg, "...")
    } else{
        message("\nError caught:\n", cond, "\n-----------")
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
            warning = function(cond) { msgCheck(cond,T)
                do.call(install.packages, c(pk.opt, list(type = "source")))
            },
            error = function(cond) { msgCheck(cond,T)
                do.call(install.packages, c(pk.opt, list(type = "binary")))
            },
            custom_error = function(cond) { msgCheck(cond,T)
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
        params <- list(pkgs=fn[[1]], repos = NULL, type = "source")
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
        bio.opt <- list(pkgs=pknm, ask=F,update=F)
        tryCatch(
            expr = {sw(do.call(BiocManager::install, c(bio.opt)))},
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
    cxn <- file(siteProf); lines <- base::readLines(cxn)
    if (!any(grepl(txt1,lines))) {write(txt2, file=siteProf, append=T)}
    close(cxn)
#    closeAllConnections()
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
            }else{fixNeeds();fixProf()}
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
    options(needs.promptUser=F); Sys.setenv(R_ENABLE_JIT=T)
    options("needs.promptUser"=F); options("promptUser"=F)
    options("device.ask.default" = FALSE)
    options("install.packages.compile.from.source" = "No")
    options("install.packages.check.source" = "no")
    compiler::enableJIT(3); compiler::compilePKGS(enable=T);
    compiler::setCompilerOptions(suppressAll=T, optimize=3)
    fixProf()
}

cpMakeV <- function(editFile=F){
    mkFi <- "/Volumes/CBioinformatics/jonathan/Rprojects/Methylation-scripts/Makevars"
    message("Copying Makevars to local dotfiles: ", mkFi)
    if(editFile){file.edit(mkFi)}
    if(file.exists(mkFi)){file.copy(mkFi,"~/.R/", overwrite=T )}
}

checkBioC <- function(){
    if(sw(suppressPackageStartupMessages((!require("BiocManager"))))){
        install.packages("BiocManager", Ncpus = 6)
        BiocManager::install(version="3.10", update=F, ask=F, type="source")
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
    sw(suppressPackageStartupMessages(easypackages::packages(cpuPacks, prompt=F)))
    checkBioC()
}

readyPkgs <- function(ezLibs){
    ready <- ezLibs %in% sessioninfo::package_info()[["package"]]
    if(any(!ready)){lapply(ezLibs[!ready], require, character.only=TRUE)}
}

installAll <- function(pkgList, instFun){
    invisible(lapply(X=1:length(pkgList),function(X){instFun(pkgList[X])}))
}

# Loads default packages or custom if input provided
#' @return package installs or loads
#' @export
loadPacks <- function(pkgs=cranPkgs, ezLibs=easyPkgs, ghPk=gHubPkgs, bcPks=biocPkgs) {
    tryCatch(
        expr = {
            loadMainPkgs()
            installAll(pkgs, pk.inst)
            installAll(ghPk, gh.inst)
            installAll(bcPks, bc.inst)
            readyPkgs(ezLibs)
        },
        error = function(cond){
            message("\n~~~You encountered the following error during install:\n", cond,"\n")
            message("If there is a compile error, try running fixCompiles() and then try loadPacks() again")
        }
    )
    knitr::opts_chunk$set(echo = T)
    Sys.setenv(R_ENABLE_JIT = T)
    compiler::enableJIT(3); compiler::compilePKGS(enable = T); compiler::setCompilerOptions(suppressAll = T, optimize = 3)
    
}

setDirectory <- function(foldr) {
    bsDir = paste("cd", foldr); mm2 = crayon::white$bgRed("Location Not Found:", foldr)
    if (dir.exists(foldr)) {system(bsDir);setwd(foldr);assign("runDir", foldr)} else{warning(mm2)}
}

loadPacks()
