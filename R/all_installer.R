options("install.packages.compile.from.source" = "No")
options("install.packages.check.source" = "no")

if (identical(Sys.info()[["machine"]], "x86_64")) {
    local({
        path <- sub(":/opt/homebrew/bin", ":/usr/local/homebrew/bin", Sys.getenv("PATH"))
        Sys.setenv(PATH = path)
    })
}

if (!dir.exists(file.path("~", ".R"))) {
    message("No Makevars file in ~/.R")
    system("mkdir -p ~/.R")
    system("touch ~/.R/MakeVars")
    fileConn <- file("~/.R/MakeVars")
    message("Creating MakeVars in ~/.R/MakeVars")
    params <- c(
        'FLIBS=-L/usr/local/gfortran/lib',
        'CC = gcc',
        'CXX = g++',
        'CXX98 = g++',
        'CXX11 = g++',
        'CXX14 = g++',
        'CXX17 = g++',
        'CXX20 = g++',
        'CXXCPP = g++',
        'FC = gfortran',
        'F77 = gfortran',
        'OBJC = gcc',
        'OBJCXX = g++',
        'CXX1X=/usr/local/gfortran/bin/g++',
        'CXX98=/usr/local/gfortran/bin/g++',
        'CXX11=/usr/local/gfortran/bin/g++',
        'CXX14=/usr/local/gfortran/bin/g++',
        'CXX17=/usr/local/gfortran/bin/g++',
        'LLVM_LOC = /usr/local/opt/llvm',
        'CC=/usr/local/gfortran/bin/gcc -fopenmp',
        'CXX=/usr/local/gfortran/bin/g++ -fopenmp',
        'CFLAGS=-g -O3 -Wall -pedantic -std=gnu99 -mtune=native -pipe',
        'CXXFLAGS=-g -O3 -Wall -pedantic -std=c++11 -mtune=native -pipe',
        'LDFLAGS=-L/usr/local/opt/gettext/lib -L$(LLVM_LOC)/lib -Wl,-rpath,$(LLVM_LOC)/lib',
        'LDFLAGS=-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib',
        'CPPFLAGS=-I/usr/local/opt/gettext/include -I$(LLVM_LOC)/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include'
    )
    writeLines(params, fileConn)
    close(fileConn)
    cmd = 'devtools::source_url("https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/all_installer.R")'
    #rstudioapi::restartSession(cmd)
}

if (!file.exists(file.path("~", ".Renviron"))) {
    system("touch ~/.Renviron")
    fileConn <- file("~/.Renviron")
    params <- c('PATH="/usr/local/gfortran/bin:${PATH}"')
    writeLines(params, fileConn)
    close(fileConn)
    cmd = 'devtools::source_url("https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/all_installer.R")'
    #rstudioapi::restartSession(cmd)
}

# Setting US CRAN REPO
rlis = getOption("repos")
rlis["CRAN"] = "http://cran.us.r-project.org"
options(repos = rlis)


loadLibrary <- function(pkgName) {
    suppressPackageStartupMessages(library(
        pkgName,
        quietly = TRUE,
        logical.return = TRUE,
        warn.conflicts = FALSE,
        character.only = TRUE
    ))
}

checkRequire <- function(pkgName){
    return(
        suppressWarnings(!require(pkgName,character.only = T,warn.conflicts = FALSE))
    )
}

if(checkRequire("devtools")){install.packages("devtools", dependencies=T, verbose=T, quiet=T, ask=F)}
if(checkRequire("librarian")){install.packages("librarian", dependencies=T,verbose=T,ask=F, quiet = FALSE)}
loadLibrary("devtools")

#if(checkRequire("rstudioapi")){devtools::install_github("rstudio/rstudioapi", dependencies = T, upgrade = "never")}
# if(rstudioapi::buildToolsCheck()==FALSE){
#     system("xcode-select --install")
#     system("xcode-select -s /Library/Developer/CommandLineTools")
# }

isOpen <- system("which openssl", intern = TRUE)
if(!exists("isOpen") | isOpen != "/usr/bin/openssl"){
    system("brew update")
    system("brew install openssl")
    system("ln -sf /usr/local/opt/openssl/lib/libcrypto.3.dylib /usr/local/lib/")
    system("ln -sf /usr/local/opt/openssl/lib/libssl.3.dylib /usr/local/lib/")
    system("ln -sf /usr/local/Cellar/openssl@3/3.0.5/bin/openssl /usr/local/bin/openssl")
}

corePkgs <- c("randomForest","glmnet","ggplot2","gridExtra","knitr","pander")

# List of all Packages
pkgs <-
    c(
        "abind",
        "affy",
        "affyio",
        "animation",
        "annotate",
        "AnnotationDbi",
        "AnnotationFilter",
        "AnnotationForge",
        "AnnotationHub",
        "AnnotationHubData",
        "ape",
        "aplot",
        "arrow",
        "askpass",
        "assertr",
        "assertthat",
        "audio",
        "aws.s3",
        "aws.signature",
        "backports",
        "bamsignals",
        "base",
        "base64",
        "base64enc",
        "base64url",
        "batchtools",
        "bayestestR",
        "BBmisc",
        "beanplot",
        "beepr",
        "bezier",
        "BH",
        "Biobase",
        "BiocCheck",
        "BiocFileCache",
        "BiocIO",
        "BiocManager",
        "BiocParallel",
        "BiocStyle",
        "BiocVersion",
        "biocViews",
        "biomaRt",
        "Biostrings",
        "biovizBase",
        "bit",
        "bit64",
        "bitops",
        "blob",
        "bookdown",
        "boot",
        "brio",
        "broom",
        "bs4Dash",
        "BSgenome",
        "BSgenome.Hsapiens.NCBI.GRCh38",
        "bslib",
        "bumphunter",
        "cachem",
        "Cairo",
        "callr",
        "car",
        "carData",
        "caret",
        "caTools",
        "cellranger",
        "changepoint",
        "checkmate",
        "chromote",
        "chron",
        "circlize",
        "class",
        "classInt",
        "cli",
        "clipr",
        "clue",
        "cluster",
        "clustermq",
        "clusterProfiler",
        "cmprsk",
        "cn.mops",
        "codetools",
        "colorspace",
        "commonmark",
        "compiler",
        "ComplexHeatmap",
        "config",
        "conquer",
        "conumee",
        "CopyNumberPlots",
        "coro",
        "corpcor",
        "covr",
        "cowplot",
        "cpp11",
        "crayon",
        "credentials",
        "cronR",
        "crosstalk",
        "cubature",
        "curl",
        "cyclocomp",
        "data.table",
        "datasets",
        "datawizard",
        "DBI",
        "dbplyr",
        "debugme",
        "DelayedArray",
        "DelayedMatrixStats",
        "dendextend",
        "DEoptimR",
        "desc",
        "deSolve",
        "devtools",
        "dichromat",
        "diffobj",
        "digest",
        "dlm",
        "DNAcopy",
        "DO.db",
        "docstring",
        "doParallel",
        "doRNG",
        "DOSE",
        "doSNOW",
        "downloader",
        "dplyr",
        "DT",
        "dtplyr",
        "dygraphs",
        "e1071",
        "easypackages",
        "effectsize",
        "egg",
        "ellipse",
        "ellipsis",
        "emmeans",
        "enrichplot",
        "ensembldb",
        "estimability",
        "evaluate",
        "exomeCopy",
        "ExperimentHub",
        "ExperimentHubData",
        "extrafont",
        "extrafontdb",
        "fansi",
        "fastICA",
        "fastmap",
        "fastmatch",
        "fBasics",
        "FDb.InfiniumMethylation.hg19",
        "feather",
        "fGarch",
        "fgsea",
        "filelock",
        "flexmix",
        "FNN",
        "fontawesome",
        "forcats",
        "forecast",
        "foreign",
        "formatR",
        "formattable",
        "Formula",
        "fracdiff",
        "fresh",
        "fs",
        "fst",
        "furrr",
        "futile.logger",
        "futile.options",
        "future",
        "future.apply",
        "future.callr",
        "gamm4",
        "gargle",
        "gbm",
        "gdata",
        "gdtools",
        "genefilter",
        "generics",
        "GenomeInfoDb",
        "GenomeInfoDbData",
        "GenomicAlignments",
        "GenomicFeatures",
        "GenomicRanges",
        "GEOquery",
        "gert",
        "getopt",
        "GetoptLong",
        "ggeffects",
        "ggforce",
        "ggfortify",
        "ggfun",
        "ggnewscale",
        "ggplot2",
        "ggplotify",
        "ggraph",
        "ggrepel",
        "ggtext",
        "ggthemes",
        "ggtree",
        "ggupset",
        "gh",
        "gitcreds",
        "glmnet",
        "GlobalOptions",
        "globals",
        "glue",
        "gmodels",
        "gmp",
        "GO.db",
        "googledrive",
        "googlesheets4",
        "GOSemSim",
        "gower",
        "gplots",
        "graph",
        "graphics",
        "graphlayouts",
        "grDevices",
        "grid",
        "gridBase",
        "gridExtra",
        "gridGraphics",
        "gridtext",
        "gss",
        "gsubfn",
        "gt",
        "gtools",
        "haven",
        "HDF5Array",
        "here",
        "hexbin",
        "highr",
        "Hmisc",
        "hms",
        "hrbrthemes",
        "htmlTable",
        "htmltools",
        "htmlwidgets",
        "httpuv",
        "httr",
        "igraph",
        "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        "IlluminaHumanMethylation450kmanifest",
        "IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
        "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
        "illuminaio",
        "impute",
        "ini",
        "insight",
        "interactiveDisplayBase",
        "ipred",
        "IRanges",
        "irlba",
        "jpeg",
        "jquerylib",
        "jsonlite",
        "kableExtra",
        "karyoploteR",
        "KEGGgraph",
        "KEGGREST",
        "keras",
        "kernlab",
        "KernSmooth",
        "KFAS",
        "knitr",
        "ks",
        "labelVector",
        "lambda.r",
        "lars",
        "later",
        "lattice",
        "latticeExtra",
        "lava",
        "lazyeval",
        "leaps",
        "lfda",
        "lifecycle",
        "limma",
        "lintr",
        "listenv",
        "lme4",
        "lmodel2",
        "lmtest",
        "locfit",
        "lubridate",
        "lumi",
        "magick",
        "magrittr",
        "mapdata",
        "mapproj",
        "maps",
        "maptools",
        "markdown",
        "MASS",
        "Matrix",
        "MatrixGenerics",
        "MatrixModels",
        "matrixStats",
        "mclust",
        "mdthemes",
        "memoise",
        "methods",
        "MethylAid",
        "methylumi",
        "mgcv",
        "mime",
        "minfi",
        "miniUI",
        "minqa",
        "misc3d",
        "mlr",
        "mockery",
        "ModelMetrics",
        "modelr",
        "modeltools",
        "MSwM",
        "multcomp",
        "multcompView",
        "multicool",
        "multtest",
        "mvtnorm",
        "narray",
        "needs",
        "nleqslv",
        "nlme",
        "nloptr",
        "nnet",
        "nor1mix",
        "numDeriv",
        "openssl",
        "openxlsx",
        "optparse",
        "org.Hs.eg.db",
        "OrganismDbi",
        "packrat",
        "pals",
        "pander",
        "parallel",
        "parallelly",
        "parallelMap",
        "parameters",
        "ParamHelpers",
        "patchwork",
        "pathview",
        "pbkrtest",
        "pcaPP",
        "pdftools",
        "performance",
        "pillar",
        "pingr",
        "pkgbuild",
        "pkgconfig",
        "pkgload",
        "plogr",
        "plot3D",
        "plotly",
        "plotrix",
        "plyr",
        "png",
        "polspline",
        "Polychrome",
        "polyclip",
        "polynom",
        "pracma",
        "praise",
        "preprocessCore",
        "prettyunits",
        "pROC",
        "processx",
        "prodlim",
        "progress",
        "progressr",
        "promises",
        "ProtGenerics",
        "proto",
        "proxy",
        "ps",
        "purrr",
        "qdapTools",
        "qpdf",
        "qs",
        "quadprog",
        "quantmod",
        "quantreg",
        "qvalue",
        "R.cache",
        "R.methodsS3",
        "R.oo",
        "R.utils",
        "R6",
        "randomForest",
        "RApiSerialize",
        "rappdirs",
        "rARPACK",
        "raster",
        "RBGL",
        "rcmdcheck",
        "Rcpp",
        "RcppArmadillo",
        "RcppParallel",
        "RCurl",
        "readr",
        "readxl",
        "recipes",
        "redcapAPI",
        "regioneR",
        "rematch",
        "rematch2",
        "remotes",
        "reprex",
        "reshape",
        "reshape2",
        "restfulr",
        "reticulate",
        "rex",
        "Rgraphviz",
        "RH2",
        "rhdf5",
        "rhdf5filters",
        "Rhdf5lib",
        "Rhtslib",
        "rio",
        "rJava",
        "RJDBC",
        "rjson",
        "rlang",
        "rmarkdown",
        "Rmpfr",
        "rms",
        "RMySQL",
        "rngtools",
        "robustbase",
        "ROC",
        "ROCR",
        "roxygen2",
        "rpart",
        "RPMM",
        "rprojroot",
        "Rsamtools",
        "rsconnect",
        "rsm",
        "RSpectra",
        "RSQLite",
        "rstudioapi",
        "rtracklayer",
        "Rtsne",
        "Rttf2pt1",
        "RUnit",
        "rversions",
        "rvest",
        "s2",
        "sass",
        "scatterpie",
        "scatterplot3d",
        "scrime",
        "selectr",
        "sesame",
        "sesameData",
        "sessioninfo",
        "shadowtext",
        "shiny",
        "siggenes",
        "sjlabelled",
        "sjmisc",
        "snow",
        "sourcetools",
        "sparkline",
        "spatial",
        "splines",
        "sqldf",
        "SQUAREM",
        "stabledist",
        "stats",
        "stats4",
        "stringdist",
        "stringfish",
        "stringi",
        "stringr",
        "strucchange",
        "styler",
        "SummarizedExperiment",
        "survival",
        "svglite",
        "svUnit",
        "sys",
        "systemfonts",
        "targets",
        "tcltk",
        "tensorflow",
        "terra",
        "testit",
        "testthat",
        "tfautograph",
        "tfruns",
        "TH.data",
        "tibble",
        "tidygraph",
        "tidyr",
        "tidyselect",
        "tidytree",
        "tidyverse",
        "tiff",
        "timeDate",
        "timeSeries",
        "tinytex",
        "tools",
        "torch",
        "treeio",
        "truncnorm",
        "tseries",
        "TxDb.Hsapiens.UCSC.hg19.knownGene",
        "tzdb",
        "units",
        "urca",
        "utf8",
        "utils",
        "uuid",
        "waiter",
        "waldo",
        "wateRmelon",
        "webshot",
        "webshot2",
        "websocket",
        "wheatmap",
        "whisker",
        "withr",
        "wk",
        "xfun",
        "XML",
        "xml2",
        "xmlparsedata",
        "xopen",
        "xtable",
        "xts",
        "XVector",
        "yaml",
        "zip",
        "zlibbioc",
        "zoo"
    )

preReqPkgs <- c(
    'MASS',
    'ade4',
    'dbplyr',
    'filelock',
    'formatR',
    'bitops',
    'progress',
    'BiocFileCache',
    'bit',
    'png',
    'Rhtslib',
    'rjson',
    'rngtools',
    'hms',
    'vroom',
    'tzdb',
    'generics',
    'tidyselect',
    'R.oo',
    'R.methodsS3',
    'lambda.r',
    'futile.options',
    'affyio',
    'zlibbioc',
    'RCurl',
    'XVector',
    'BiocIO',
    'biomaRt',
    'XML',
    'bit64',
    'blob',
    'plogr',
    'KEGGREST',
    'plyr',
    'TxDb.Hsapiens.UCSC.hg19.knownGene',
    'org.Hs.eg.db',
    'GenomeInfoDbData',
    'MatrixGenerics',
    'base64',
    'Rsamtools',
    'GenomicAlignments',
    'restfulr',
    'locfit',
    'doRNG',
    'multtest',
    'scrime',
    'sparseMatrixStats',
    'readr',
    'dplyr',
    'tidyr',
    'R.utils',
    'rhdf5',
    'rhdf5filters',
    'Rhdf5lib',
    'futile.logger',
    'snow',
    'BH',
    'affy',
    'GenomicFeatures',
    'GenomicRanges',
    'annotate',
    'nleqslv',
    'preprocessCore',
    'RSQLite',
    'DBI',
    'AnnotationDbi',
    'reshape2',
    'matrixStats',
    'FDb.InfiniumMethylation.hg19',
    'S4Vectors',
    'IRanges',
    'GenomeInfoDb',
    'SummarizedExperiment',
    'genefilter',
    'illuminaio',
    'IlluminaHumanMethylationEPICanno.ilm10b2.hg19',
    'DNAcopy',
    'rtracklayer',
    'Biostrings',
    'bumphunter',
    'beanplot',
    'nor1mix',
    'siggenes',
    'limma',
    'DelayedMatrixStats',
    'mclust',
    'reshape',
    'quadprog',
    'data.table',
    'GEOquery',
    'DelayedArray',
    'HDF5Array',
    'BiocParallel'
)


biocPkgs <-
    c(
        "lumi",
        "methylumi",
        "conumee",
        "minfi",
        "IlluminaHumanMethylation450kmanifest",
        "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
    )

message("Version of Clang:"); system("gcc --version")

if (checkRequire("BiocManager")) {
    install.packages(
        "BiocManager",
        dependencies = T,
        verbose = T,
        ask = F
    )
}
loadLibrary("devtools")
loadLibrary("librarian")
loadLibrary("BiocManager")
if (checkRequire("Biobase")) {BiocManager::install("Biobase", update = F, ask = F)}
loadLibrary("Biobase")

if (checkRequire("mapview")) {
remotes::install_github("r-spatial/mapview")
    }

#if(checkRequire("S4Vectors")){BiocManager::install("S4Vectors", update = F, ask = F)}
#if(checkRequire("MASS")){install.packages("MASS",type = "binary", update = F, ask = F)}
# Load/install missing pacakges without asking

suppressWarnings(librarian::shelf(corePkgs, ask = F, update_all = F, quiet = FALSE))
suppressWarnings(librarian::shelf(preReqPkgs, ask = F, update_all = F, quiet = FALSE))
invisible(gc())
message("Loading BioConductor Packages and IlluminaHumanMethylation Manifest...")
suppressWarnings(librarian::shelf(biocPkgs, ask = F, update_all = F, quiet = FALSE))
if(checkRequire("IlluminaHumanMethylationEPICmanifest")){
    devtools::install_github(repo = "mwsill/IlluminaHumanMethylationEPICmanifest",
                             dependencies = T,
                             upgrade = "never")
}

if(checkRequire("mgmtstp27")){
    gitLink <- "https://github.com/badozor/mgmtstp27/raw/master/archive/mgmtstp27_0.6-3.tar.gz"
    install.packages(gitLink, repos = NULL, dependencies = T, verbose = T, type = "source", ask = F)
}
invisible(gc())

fixProf <- function() {
    txt1 <- "^[:blank:]*autoload\\(\"needs\", \"needs\"\\)"
    txt2 <- "\n\nautoload(\"needs\", \"needs\")\n\n"
    siteProf <- if (is.na(Sys.getenv("R_PROFILE", unset = NA))) {
        file.path(Sys.getenv("R_HOME"), "etc", "Rprofile.site")
    } else {
        Sys.getenv("R_PROFILE")
    }
    if (!file.exists(siteProf)) {
        try(file.create(siteProf), silent=T)
    }
    cxn <- try(file(siteProf), silent=T)
    lines <- try(base::readLines(cxn),silent=T)
    if (!any(grepl(txt1, lines))) {
        write(txt2, file = siteProf, append = T)
    }
    close(cxn)
}

fixNeeds <- function (){
    sysfile <- system.file("extdata", "promptUser", package="needs")
    try(write(0, file=sysfile),silent=T)
    options(needs.promptUser=FALSE)
    invisible(needs:::autoload(TRUE))
}

checkNeeds <- function(){
    tryCatch(
        expr={
            if(!("needs" %in% rownames(installed.packages()))){
                install.packages("needs",dependencies=T,verbose=T, Ncpus = 2)
                fixNeeds()
                try(fixProf(), silent=T)
            }else{
                fixNeeds()
                try(fixProf(), silent=T)
            }
        },
        error=function(cond){
            devtools::install_github("joshkatz/needs", ref = "development",
                                     dependencies=T,verbose=T,upgrade="always")
            fixNeeds()
            try(fixProf(), silent=T)
        },
        warning=function(cond){
            devtools::install_github(
                "joshkatz/needs", ref = "development",
                                     dependencies=T,verbose=T,upgrade="always")
            fixNeeds()
            try(fixProf(), silent=T)
        }
    )
}

checkNeeds()

spat_config <- '--with-proj-lib=/usr/local/lib/ --with-proj-include=/usr/local/include/'
options(configure.args = c("sf" = spat_config, "rgdal" = spat_config))
if(checkRequire("sf")){
    install.packages(c("sf"), type = "source", dependencies=T, verbose=T, Ncpus = 4)
}
if(checkRequire("rgdal")){
    install.packages(c("rgdal"), type = "source", dependencies=T, verbose=T, Ncpus = 4)
}
invisible(gc())
system("export LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib")
options("install.packages.compile.from.source" = "No")
options("install.packages.check.source" = "no")
loadLibrary("BiocManager")
loadLibrary("Biobase")

terraDep <- c('tinytest', 'ncdf4', 'leaflet')
suppressWarnings(librarian::shelf(terraDep, ask = F, update_all = F, quiet = FALSE))
if(checkRequire("terra")) {
    install.packages(
        'terra',
        repos = 'https://rspatial.r-universe.dev',
        dependencies = T,
        verbose = T,
        Ncpus = 4
    )
}

if(checkRequire("FField")){
    gitLink <- "https://cran.r-project.org/src/contrib/Archive/FField/FField_0.1.0.tar.gz"
    install.packages(gitLink, repos = NULL, dependencies = T, verbose = T, type = "source", ask = F)
}

if (checkRequire("GenVisR")) {
    devtools::install_github("griffithlab/GenVisR", dependencies = T, upgrade = "never")
}
invisible(gc())
suppressWarnings(librarian::shelf(pkgs, ask = F, update_all = F, quiet = FALSE))
invisible(gc())
if (checkRequire("UniD")) {
try(install.packages("/Volumes/CBioinformatics/Methylation/UniD/", type="source", repos=NULL), silent=T)
}
