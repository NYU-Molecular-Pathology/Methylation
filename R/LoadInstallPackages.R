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
if (!require("devtools")) {
    install.packages("devtools", dependencies = T, ask = F)
}
message("Your library path(s):")
message(paste(.libPaths()))
yourHome <- fs::path_home()

isMacOS <- Sys.info()[['sysname']] == "Darwin"

if (isMacOS) {
    cbioLn <-  "/Volumes/CBioinformatics/Methylation/classifiers"
}else {
    cbioLn <- "/Volumes/CBioinformatics/Methylation/classifiers"
}

getBrewDir <- function(module_name) {
    brew_cmd <- paste("brew --prefix", module_name)
    module_path <- tryCatch(system(brew_cmd, intern = TRUE, ignore.stderr = TRUE),
                            error = function(e) NULL,
                            warning = function(e) NULL)
    return(module_path)
}


proj_path <- getBrewDir("proj")
sqli_path <- getBrewDir("sqlite")

if (isMacOS) {
    Sys.setenv(PROJ_LIBS = file.path(proj_path,"lib"))
    Sys.setenv(SQLITE3_LIBS = file.path(sqli_path,"lib"))
} else{
    options(BioC_mirror = "https://packagemanager.rstudio.com/bioconductor")
    options(repos = c(CRAN = "https://packagemanager.rstudio.com/cran/__linux__/jammy/latest"))
}

if (isMacOS) {
    typeSrc <- "binary"
} else{
    typeSrc <- "source"
}

uniDpath <- file.path(cbioLn, "UniD")

# Custom Classifier packages:
if (isMacOS) {
    classPacks <- c(
       # sest = "https://github.com/jungch/sest/raw/master/sest.tar",
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
    'HDF5Array', 'bumphunter','GEOquery', 'lumi', 'rngtools',
    'methylumi', 'randomForest', 'glmnet','IlluminaHumanMethylation450kmanifest',
    'IlluminaHumanMethylation450kanno.ilmn12.hg19', 'IlluminaHumanMethylationEPICmanifest', 'Rtsne',
    'IlluminaHumanMethylationEPICanno.ilm10b2.hg19', 'IlluminaHumanMethylationEPICanno.ilm10b4.hg19',
    'MethylAid', 'conumee','BiocParallel', "Biobase", "limma","MatrixGenerics", "GenomeInfoDb"
)

# Extra Libraries ----
cpuPacks <- c("parallel","doSNOW","doParallel", "foreach","compiler")
easyPkgs <- c('tidyverse','sjmisc','stringi','digest','RCurl','gridExtra')

# Helper Functions ----
sup <- function(x) {return(suppressWarnings(suppressPackageStartupMessages(x)))}
sw <- function(pkgOb) {
    try(return(suppressMessages(suppressWarnings(pkgOb))), silent = T)
}

ld <- function(libName) {
    lib.opts <- list(package = libName, character.only = T,verbose = T, warn.conflicts = F, quietly = F)
    suppressPackageStartupMessages(do.call(library, c(lib.opts)))
    message(libName, " ...load successful")
}

up <- function() {
    update.packages(
        repos = 'http://cran.rstudio.com/',
        type = "source",
        ask = F,
        checkBuilt = T
    )
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


if (isMacOS) {
    bothType <- "both"
} else{
    bothType <- "source"
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
srcInst <- function(fn) {
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
        bio.opt <- list(pkgs = pknm, update = F, ask = F, dependencies = T)
        tryCatch(
            expr = {do.call(BiocManager::install, c(bio.opt))},
            error = function(cond) {
                message("Error with Package install:\n", pknm)
                message(cond)
            },
            finally = {ld(pknm)}
        )
    } else {ld(pknm)}
}


fixProf <- function() {
    txt1 <- "^[:blank:]*autoload\\(\"needs\", \"needs\"\\)"
    txt2 <- "\n\nautoload(\"needs\", \"needs\")\n\n"
    profPath <- Sys.getenv("R_PROFILE", unset = NA)
    siteProf <- if (is.na(profPath)) {
        file.path(Sys.getenv("R_HOME"), "etc", "Rprofile.site")
    } else {
        Sys.getenv("R_PROFILE")
    }
    if (siteProf == "") {

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
            try(write(txt2, file = siteProf, append = T), silent = T)
        }
        close(cxn)
    }
}

fixNeeds <- function() {
    sysfile <- system.file("extdata", "promptUser", package = "needs")
    options(needs.promptUser = FALSE)
    invisible(needs:::autoload(TRUE))
    if (file.exists(sysfile)) {
        try(write(0, file = sysfile), silent = T)
    }
}


checkNeeds <- function() {
    tryCatch(
        expr = {
            if (!("needs" %in% rownames(installed.packages()))) {
                install.packages(
                    "needs",
                    dependencies = T,
                    verbose = T,
                    Ncpus = 6
                )
                invisible(needs:::autoload(TRUE))
                fixNeeds()
                try(fixProf(), T)
            } else{
                fixNeeds()
                try(fixProf(), T)
            }
        },
        error = function(cond) {
            devtools::install_github(
                "joshkatz/needs",
                ref = "development",
                dependencies = T,
                verbose = T,
                upgrade = "always"
            )
            invisible(needs:::autoload(TRUE))
            fixNeeds()
            try(fixProf(), T)
        },
        warning = function(cond) {
            devtools::install_github(
                "joshkatz/needs",
                ref = "development",
                dependencies = T,
                verbose = T,
                upgrade = "always"
            )
            invisible(needs:::autoload(TRUE))
            fixNeeds()
            try(fixProf(), T)
        }
    )
}

setOptions <- function() {
    options(needs.promptUser = F)
    options("needs.promptUser" = F)
    options("promptUser" = F)
    options("device.ask.default" = FALSE)
    options("install.packages.compile.from.source" = "Yes")
    options("install.packages.check.source" = "yes")
    fixProf()
}

cpMakeV <- function(editFile = F) {
    if (!dir.exists(file.path(yourHome, ".R/"))) {
        message("No Makevars file in ~/.R/")
        #dir.create("~/.R/")
    } else{
        system("rm -R ~/.R/")
    }
}


# FUN: Install or Load the Required Classifier Version(s) ---------------------
classifierInstall <- function(pathtoFile = NULL,
                              instNew = T,
                              rmpkg = F) {
    mnp.pk.loc <- file.path(cbioLn, pathtoFile)
    if (rmpkg) {
        remove.packages(basename(mnp.pk.loc), lib = .libPaths()[[1]][1])
    }
    if (instNew) {
        message("Installing package: ", basename(mnp.pk.loc))
        message("Classifier install may take several minutes!")
        message("Note the stage '** inst' may take the longest.")
        install.packages(
            mnp.pk.loc,
            repos = NULL,
            type = "source",
            force = T,
            Ncpus = 6
        )
        classifierInstall(pathtoFile, F, F)
    } else {
        ld(basename(mnp.pk.loc))
    }
}

checkBioC <- function() {
    if (sw(suppressPackageStartupMessages((!require(
        "BiocManager", warn.conflicts = F
    ))))) {
        install.packages("BiocManager", Ncpus = 4)
    } else{
        ld("BiocManager")
    }

    if (rq("zip")) {
        install.packages("zip", dependencies = T, type = typeSrc)
    } else{
        ld("zip")
    }
}

loadMainPkgs <- function() {
    setOptions()
    if (rq("devtools")) {
        install.packages("devtools",
                         dependencies = T,
                         verbose = T)
    } else{
        ld("devtools")
    }
    checkNeeds()
    gh.inst(data.frame(easypackages = "jakesherman/easypackages"))
    pk.inst("tidyverse")
    suppressPackageStartupMessages(easypackages::packages(cpuPacks, prompt = F))
    checkBioC()
}

readyPkgs <- function(ezLibs) {
    ready <- ezLibs %in% sessioninfo::package_info()[["package"]]
    if (any(!ready)) {lapply(ezLibs[!ready], require, character.only = TRUE)}
}

installAll <- function(pkgList, instFun) {
    invisible(lapply(X = 1:length(pkgList), function(X) {
        instFun(pkgList[X])
    }))
}

# Loads default packages or custom if input provided
loadPacks <- function(pkgs = cranPkgs,
                      ezLibs = easyPkgs,
                      ghPk = gHubPkgs,
                      bcPks = biocPkgs) {

    tryCatch(
        expr = {
            loadMainPkgs()
            sup(installAll(pkgs, pk.inst))
            if (!require("minfi")) {
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
            if (!require("BiocManager", warn.conflicts = F)) {
                install.packages("BiocManager", dependencies = T, quiet = F)
            }
            sup(installAll(ghPk, gh.inst))
            sup(installAll(bcPks, bc.inst))
            readyPkgs(ezLibs)
            minfiVers <- as.character(utils::packageVersion("minfi"))
            if (minfiVers != "1.43.1") {
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
            if (!require("MethylAid", warn.conflicts = F)) {
                BiocManager::install("MethylAid", update = F, ask = F)
            }
            if (!require("librarian", warn.conflicts = F)) {
                install.packages(
                    "librarian",
                    dependencies = T,
                    verbose = T,
                    Ncpus = 4,
                    quiet = F
                )
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
            librarian::shelf(pkgs, ask = F,  update_all = FALSE, quiet = FALSE)
            require('grid')
            require("ggplot2")
            require("pals")
            require("scales")
            require("stringr")
            require("scales")
        },
        error = function(cond) {
            message("\n~~~You encountered the following error during install:\n", cond)
        }
    )
}

colorMsg <- function() {
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

startmsg <- function() {
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

setEnviron <- function() {
    options(Ncpus = 6)
    startmsg()
    library("base")
    #try(cpMakeV(),silent=T)
    rlis = getOption("repos")
    rlis["CRAN"] = "http://cran.us.r-project.org"
    options(repos = rlis)
    if (isMacOS) {
        Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/MacOS/pandoc")
        system("export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc")
        system("export RSTUDIO_WHICH_R=/usr/local/bin/R")
        Sys.setenv(RSTUDIO_WHICH_R = "/usr/local/bin/R")
        try(options(gsubfn.engine = "R"), T)
    }
    options("install.packages.compile.from.source" = "never")
    options("install.packages.check.source"="no")
    if (!require("devtools", warn.conflicts = F)) {
        install.packages("devtools", dependencies = T)
    }
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


checkEpicV2 <- function(pkg = "mnp.v12epicv2") {
    if (pkg %in% rownames(installed.packages())) {
        message(paste("Package", pkg, "is installed with version:", utils::packageVersion(pkg)))
    } else{
        epicV2script <- "/Volumes/CBioinformatics/Methylation/Rscripts/install_epic_v2_classifier.R"
        message("Installing package from source:\n", epicV2script)
        source(epicV2script)
    }
}


# Load all Functions ---------------------
check_uniD_pkg <- function(typeSrc) {

    uniDpkgs <- c("lumi", "ade4", "methylumi", "mlr")

    librarian::shelf(uniDpkgs, ask = F, update_all = F, quiet = F)
    try(bc.inst("impute"), silent = T)
    try(bc.inst("wateRmelon"), silent = T)

    if (rq("UniD")) {
        try(install.packages(uniDpath, type = "source", dependencies = T, repo = NULL), silent = T)
    }
}


startLoadingAll <- function() {
    setEnviron()
    loadPacks()

    if (isMacOS) {
        typeSrc <- "binary"
    } else{
        typeSrc <- "source"
        classPath = "molecpathlab/production/Methylation/classifiers/"
    }

    installAll(classPacks, srcInst)

    if (isMacOS) {
        check_uniD_pkg(typeSrc)
        checkEpicV2()
    }

}

startLoadingAll()
