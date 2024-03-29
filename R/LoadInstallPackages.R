#!/usr/bin/env Rscript
## ---------------------------
## Script name: LoadInstallPackages.R
## Purpose: Functions that check if the required classifier packages and dependencies are installed and loaded
## Date Created: August 5, 2022
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

formals(library)$quietly <- T
formals(library)$warn.conflicts <- F
formals(require)$warn.conflicts <- F
if(!require("devtools")){install.packages("devtools", dependencies=T, ask=F)}
message("Your library path(s):")
message(paste(.libPaths()))
yourHome <- fs::path_home()

cbioLn <- switch (
    Sys.info()[['sysname']],
    "Darwin" = "/Volumes/CBioinformatics/Methylation/classifiers",
    "Linux" = "/Volumes/CBioinformatics/Methylation/classifiers"
    #"Linux" = file.path(yourHome,"molecpathlab/production/Methylation/classifiers") # BigPurple
)

if (Sys.info()[['sysname']]=="Darwin") {
    Sys.setenv(PROJ_LIBS = "/opt/homebrew/opt/proj/lib")
    Sys.setenv(SQLITE3_LIBS = "/opt/homebrew/opt/sqlite/lib")
} else{
    options(BioC_mirror = "https://packagemanager.rstudio.com/bioconductor")
    options(repos = c(CRAN = "https://packagemanager.rstudio.com/cran/__linux__/jammy/latest"))
}

if (Sys.info()[['sysname']]=="Darwin") {
    typeSrc <- "binary"
} else{
    typeSrc <- "source"
}

uniDpath <- file.path(cbioLn, "UniD")

# Classifier Packages and Versions
mnpV4 <- data.frame(mnpVers = "mnp.v11b4", mnpPath = "mnp.v11b4", mnpNumb = "0.1.126")
mnpV6 <- data.frame(mnpVers = "mnp.v11b6", mnpPath = "mnp.v11b6", mnpNumb = "0.1.140")
mnpV12 <- data.frame(mnpVers = "mnp.v12b6", mnpPath = "mnp.v12b6", mnpNumb = "0.1.132")
srcV12 <- data.frame(mnpVers = "sarc.v12b6", mnpPath = "sarc.v12b6", mnpNumb = "0.1.129")

# Custom Classifier packages:
if (Sys.info()[['sysname']]=="Darwin") {
    classPacks <- c(
    sest = "https://github.com/jungch/sest/raw/master/sest.tar",
    mgmtstp27 = "https://github.com/badozor/mgmtstp27/raw/master/archive/mgmtstp27_0.6-3.tar.gz",
    mnpqc = paste0(file.path(cbioLn, "mnpqc_0.1.0.tar.gz"))
)
} else{
    classPacks <- c(
    sest = "https://github.com/jungch/sest/raw/master/sest.tar",
    mgmtstp27 = "https://github.com/badozor/mgmtstp27/raw/master/archive/mgmtstp27_0.6-3.tar.gz"
    )
}

# Cran Packages ----
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

# GitHub Packages ----
gHubPkgs <- data.frame(
    rmarkdown = 'rstudio/rmarkdown',
    docstring = 'dasonk/docstring',
    rstudioapi = 'rstudio/rstudioapi',
    easypackages = 'jakesherman/easypackages',
    redcapAPI = 'nutterb/redcapAPI',
    crayon = "r-lib/crayon",
    redcap = "epicentre-msf/redcap"
)

# BioConductor Packages ----
biocPkgs <- c(
    'HDF5Array', 'bumphunter','GEOquery', 'minfi', 'lumi', 'rngtools',
    'methylumi', 'randomForest', 'glmnet','IlluminaHumanMethylation450kmanifest',
    'IlluminaHumanMethylation450kanno.ilmn12.hg19', 'IlluminaHumanMethylationEPICmanifest', 'Rtsne',
    'IlluminaHumanMethylationEPICanno.ilm10b2.hg19', 'IlluminaHumanMethylationEPICanno.ilm10b4.hg19',
    'MethylAid', 'conumee','BiocParallel', "Biobase", "limma","MatrixGenerics", "GenomeInfoDb"
)

# Extra Libraries ----
cpuPacks <- c("parallel","doSNOW","doParallel", "foreach","compiler")
easyPkgs <- c('tidyverse','sjmisc','stringi','digest','RCurl','gridExtra')

# Helper Functions ----
sup <- function(x){return(suppressWarnings(suppressPackageStartupMessages(x)))}
sw <- function(pkgOb){try(return(suppressMessages(suppressWarnings(pkgOb))),silent=T)}

ld <- function(libName) {
    lib.opts <- list(package = libName, character.only = T,verbose = T, warn.conflicts = F, quietly = F)
    suppressPackageStartupMessages(do.call(library, c(lib.opts)))
    message(libName, " ...load successful")
}

up <- function(){update.packages(repos='http://cran.rstudio.com/', type = "source", ask=F, checkBuilt=T)}
rq <- function(pkgName) {!requireNamespace(pkgName, quietly = TRUE)}
mkred <- function(strMsg) {return(crayon::white$bgRed$bold(strMsg))}
mkblu <- function(strMsg) {return(crayon::white$bgBlue$bold(strMsg))}
mkGrn <- function(strMsg) {return(crayon::white$bgGreen$bold(strMsg))}

msgCheck <- function(pkg, warn = F) {
    if (warn == F) {message("Checking ", pkg, "...")} else{
        message("\nError caught for package:\n", pkg, "\n-----------")
    }
}

if (Sys.info()[['sysname']]=="Darwin") {
    bothType <- "both"
}else{
    bothType <- "source"
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
                do.call(install.packages, c(pk.opt, list(type = bothType)))
            },
            warning = function(cond) {
                message("\nWarning on package:\n", pkg, "\n-----------\n",cond)
                do.call(install.packages, c(pk.opt, list(type = typeSrc)))
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
                params  <- c(gh.opt, list(type = typeSrc))
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
    profPath <- Sys.getenv("R_PROFILE", unset = NA)
    siteProf <- if (is.na(profPath)) {
        file.path(Sys.getenv("R_HOME"), "etc", "Rprofile.site")
    } else {
        Sys.getenv("R_PROFILE")
    }
    if (siteProf=="") {
        siteProf <- file.path(Sys.getenv("HOME"), "Rprofile.site")
    } else{
        if (!file.exists(siteProf)) {
            try(file.create(siteProf), silent = T)
        }
    }
    if (!file.exists(siteProf)) {
        siteProf <- file.path(Sys.getenv("HOME"), "Rprofile.site")
        if (!file.exists(siteProf)) {
            try(file.create(siteProf), silent = T)
        }
    }
    if (file.exists(siteProf)) {
        cxn <- file(siteProf)
        lines <- base::readLines(cxn)
        if (!any(grepl(txt1, lines))) {
            try(write(txt2, file = siteProf, append = T), silent=T)
        }
        close(cxn)
    }
}

fixNeeds <- function(){
  sysfile <- system.file("extdata", "promptUser", package = "needs")
  options(needs.promptUser = FALSE)
  invisible(needs:::autoload(TRUE))
  if (file.exists(sysfile)) {
    try(write(0, file = sysfile), silent = T)
  }
}


checkNeeds <- function(){
    tryCatch(
        expr={
            if(!("needs" %in% rownames(installed.packages()))){
                install.packages("needs",dependencies=T,verbose=T, Ncpus = 6)
                invisible(needs:::autoload(TRUE))
                fixNeeds()
                try(fixProf(),T)
            }else{
                fixNeeds()
                try(fixProf(),T)
            }
        },
        error=function(cond){
            devtools::install_github("joshkatz/needs", ref = "development",
                                     dependencies=T,verbose=T,upgrade="always")
            invisible(needs:::autoload(TRUE))
            fixNeeds()
            try(fixProf(),T)
        },
        warning=function(cond){
            devtools::install_github("joshkatz/needs", ref = "development",
                                     dependencies=T,verbose=T,upgrade="always")
            invisible(needs:::autoload(TRUE))
            fixNeeds()
            try(fixProf(),T)
        }
    )
}

setOptions <- function(){
    options(needs.promptUser=F)
    Sys.setenv(R_ENABLE_JIT=T)
    options("needs.promptUser"=F); options("promptUser"=F)
    options("device.ask.default" = FALSE)
    options("install.packages.compile.from.source"="Yes")
    options("install.packages.check.source"="yes")
    compiler::enableJIT(3);
    compiler::compilePKGS(enable=T);
    compiler::setCompilerOptions(suppressAll=T, optimize=3)
    fixProf()
}

cpMakeV <- function(editFile=F){
    #mkFi <- "/Volumes/CBioinformatics/jonathan/Rprojects/Methylation-scripts/Makevars"
    #message("Copying Makevars to local dotfiles: ", mkFi)
    #message("Copying Makevars to local dotfiles")
    if(!dir.exists(file.path(yourHome,".R/"))){
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
    if(sw(suppressPackageStartupMessages((!require("BiocManager", warn.conflicts = F))))){
        install.packages("BiocManager", Ncpus = 4)
        BiocManager::install(version="3.10", update=T, ask=F, type="source")
    } else{ld("BiocManager")}
    
    if (rq("zip")){
        install.packages("zip", dependencies=T, type=typeSrc)
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
            if(!require("BiocManager", warn.conflicts = F)){install.packages("BiocManager", dependencies=T, quiet=F)}
            sup(installAll(ghPk, gh.inst))
            sup(installAll(bcPks, bc.inst))
            readyPkgs(ezLibs)

            if(!require("MethylAid", warn.conflicts = F)){BiocManager::install("MethylAid",update=F, ask=F)}
            if(!require("librarian", warn.conflicts = F)){
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
                #"gridExtra",
                #"rmarkdown",
                "tinytex",
                "BiocParallel",
                'grid',
                'grDevices',
                "mdthemes",
                "dplyr",
                "sqldf",
                "gridExtra",
                "reshape",
                "mlr", "wateRmelon", "RPMM", "impute"
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


runAllBrew <- function(){
    # export RSTUDIO_WHICH_R=/usr/local/bin/R
    system("brew install aspell gdal autoconf automake bash cairo cmake coreutils dos2unix exa fd fontconfig freetype fribidi fzf gawk gcc gdbm gdk-pixbuf gettext giflib git glib gmp gnu-sed gnu-tar gnutls graphite2 harfbuzz htop icu4c ilmbase imagemagick isl jpeg libcroco libde265 libevent libffi libgcrypt libgpg-error libheif libidn2 libmpc libomp libpng librsvg libtasn1 libtermkey libtiff libtool libunistring libuv libvterm libyaml little-cms2 llvm lua lua@5.1 luajit luarocks lzo mpfr msgpack multimarkdown ncurses neovim nettle nspr nss oniguruma openconnect openexr openjpeg openssl@1.1 p11-kit pandoc pandoc-citeproc pandoc-crossref pango pcre pcre2 pdfgrep perl pixman pkg-config poppler qt readline rename ripgrep ruby shared-mime-info shellcheck sqlite stoken swig texinfo tmux tree unbound unibilium vim webp wget x265")
    system('export CFLAGS+=-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
  export CCFLAGS+=-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
  export CXXFLAGS+=-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
  export CPPFLAGS+=-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk')
    system("brew install golang libstfl")
    if(dir.exists("/Library/Developer/CommandLineTools")){
        if(length(list.dirs("/Library/Developer/CommandLineTools", recursive = F))>0){
            system("ln -s /Applications/Xcode.app/Contents/Developer /Library/Developer/CommandLineTools")}}
    system("brew unlink gcc && brew link gcc --force")
    system('brew link --overwrite gcc')
    system('brew doctor')
}


fixCompiles <- function(brewExtra=F){
    #system("/bin/bash -c $(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)")
    #cpMakeV()
    system("brew install llvm gcc")
    system("defaults write org.R-project.R NSAppSleepDisabled -bool YES")
    system('echo "export PATH="/usr/local/sbin:$PATH"" >> ~/.zshrc >')
    system('echo "export PATH="/usr/local/sbin:$PATH"" >> ~/.bashrc >')
    system("export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES")
    system('echo "export OBJC_DISABLE_INITIALIZE_FORK_SAFETY="YES"" >> ~/.zshrc')
    system('echo "export OBJC_DISABLE_INITIALIZE_FORK_SAFETY="YES"" >> ~/.bashrc')
    system("xcode-select --install")
    system("xcodebuild -runFirstLaunch")
    system("brew install dpkg openssl libffi go python")
    system('export PATH="/usr/local/opt/python/libexec/bin:$PATH"')
    system('export GOROOT=/usr/local/opt/go/libexec')
    system("export PATH=$PATH:$GOROOT/bin")
    system("brew install cmake gsl gmp libopenmpt open-mpi")
    if (brewExtra == T) {
      runAllBrew()
    }
}

colorMsg <- function(){
    ms1 <-
      paste(
        mkGrn("Updating in-house classifier to current version:"),
        mkblu("classifierInstall(pathtoFile=NULL, instNew=F, rmpkg=F)"),
        sep = "\n"
      )
    ms2 <- paste0(mkred("Classifier package is not installed installing classifier"),"\n")
    ms3 <- paste0(mkGrn("Your classifier package is up-to-date and loading"),"\n")
    return(c(ms1,ms2,ms3))
}

startmsg <- function(){
    cbio = "/Volumes/CBioinformatics/"; zdriv = "/Volumes/molecular/Molecular"
    wmm = "You do not have this path mounted:\n"
    wmm2 = "If any are missing, try running:\ndiskutil mountdisk "
    message("You have the following drives mounted:")
    system("ls /Volumes")
    if (!dir.exists(cbio)) {
        message(paste(wmm, cbio))
        message(paste(wmm2, cbio))
    }
    if (!dir.exists(zdriv)) {
        message(paste(wmm, zdriv))
        message(paste(wmm2, zdriv))
    }
}

setEnviron <- function(){
    options(Ncpus = 6)
    startmsg()
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
        try(options(gsubfn.engine = "R"), T)
        }
    options("install.packages.compile.from.source" = "never")
    options("install.packages.check.source"="no")
    if(!require("devtools", warn.conflicts = F)){install.packages("devtools", dependencies=T)}
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


checkEpicV2 <- function(pkg = "mnp.v12epicv2"){
    if (pkg %in% rownames(installed.packages())) {
      message(paste("Package", pkg, "is installed with version:", utils::packageVersion(pkg)))
    } else{
        epicV2script <- "/Volumes/CBioinformatics/Methylation/Rscripts/install_epic_v2_classifier.R"
        message("Installing package from source:\n", epicV2script)
        source(epicV2script)
    }
}


# Load all Functions ---------------------
startLoadingAll <- function() {
    setEnviron()
    loadPacks()

    if(Sys.info()[['sysname']]=="Darwin") {
        typeSrc <- "binary"
    } else{
        typeSrc <- "source"
        classPath = "molecpathlab/production/Methylation/classifiers/"
        #mgmtstp27 = file.path(yourHome, classPath,"mgmtstp27_0.6-4.tar.gz")
        #if(!require("mgmtstp27")){install.packages(mgmtstp27, repos=NULL, type="source", dependencies=T)}
    }
    
    installAll(classPacks, srcInst)
    
    if(Sys.info()[['sysname']] == "Darwin") {
        checkClassifier(mnpV4)
        checkClassifier(mnpV6)
        checkClassifier(mnpV12)
        checkClassifier(srcV12)
        #checkEpicV2()
    }
    
    uniDpkgs <- c("lumi","ade4","methylumi","mlr")
    librarian::shelf(uniDpkgs, ask=F, update_all = F, quiet = F)
    
    if(!require("impute")){
        try(BiocManager::install("impute", update=F, ask=F, dependencies=T, type=typeSrc), silent=T)
    }

    if(!require("wateRmelon")){
        try(BiocManager::install("wateRmelon", update=F, ask=F, dependencies=T, type=typeSrc), silent=T)
    }
    if(Sys.info()[['sysname']]=="Darwin") {
    if (!requireNamespace("UniD", quietly = T)) {
        try(install.packages(uniDpath, type = "source", dependencies = T, repo = NULL), silent = T)
    }
        }

}

startLoadingAll()
