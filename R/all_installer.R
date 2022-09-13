options("install.packages.compile.from.source" = "No")
options("install.packages.check.source" = "no")

if(Sys.info()[['sysname']]=="Darwin"){
    local({
        path <- sub(":/opt/homebrew/bin", ":/usr/local/homebrew/bin", Sys.getenv("PATH"))
        Sys.setenv(PATH = path)
    })
}

if (Sys.info()[['sysname']]=="Darwin" & !dir.exists(file.path("~", ".R"))) {
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
closeAllConnections()
    cmd = 'devtools::source_url("https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/all_installer.R")'
    #rstudioapi::restartSession(cmd)
}

if (Sys.info()[['sysname']]=="Darwin" & !file.exists(file.path("~", ".Renviron"))) {
    system("touch ~/.Renviron")
    fileConn <- file("~/.Renviron")
    params <- c('PATH="/usr/local/gfortran/bin:${PATH}"')
    writeLines(params, fileConn)
    close(fileConn)
closeAllConnections()
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
if (Sys.info()[['sysname']]=="Darwin"){
    isOpen <- system("which openssl", intern = TRUE)
    if(!exists("isOpen")){
        system("brew update")
        system("brew install openssl")
        system("ln -sf /usr/local/opt/openssl/lib/libcrypto.3.dylib /usr/local/lib/")
        system("ln -sf /usr/local/opt/openssl/lib/libssl.3.dylib /usr/local/lib/")
        system("ln -sf /usr/local/Cellar/openssl@3/3.0.5/bin/openssl /usr/local/bin/openssl")
    }else{
        if(isOpen != "/usr/bin/openssl"){
            system("ln -sf /usr/local/opt/openssl/lib/libcrypto.3.dylib /usr/local/lib/")
            system("ln -sf /usr/local/opt/openssl/lib/libssl.3.dylib /usr/local/lib/")
            system("ln -sf /usr/local/Cellar/openssl@3/3.0.5/bin/openssl /usr/local/bin/openssl")

        }
    }
}

corePkgs <- c("randomForest","glmnet","ggplot2","gridExtra","knitr","pander")

# List of all Packages
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

pkgs <- c(
"abind", "animation", "AnnotationFilter", "AnnotationHub", "AnnotationHubData", "ape", "aplot", "arrow", "askpass", "assertr", "assertthat", "batchtools", "beepr", "bezier", "Biobase", "BiocCheck", "BiocManager", "BiocStyle", "BiocVersion", "biocViews", "bookdown", "boot", "brio", "broom", "bs4Dash", "BSgenome", "bslib", "cachem", "Cairo", "callr", "car", "carData", "caret", "cellranger", "changepoint", "checkmate", "chromote", "chron", "circlize", "class", "classInt", "cli", "clipr", "clustermq", "clusterProfiler", "cmprsk", "cn.mops", "codetools", "colorspace", "commonmark", "compiler", "ComplexHeatmap", "config", "conquer", "conumee", "CopyNumberPlots", "cpp11", "crayon", "credentials", "cronR", "curl", "datasets", "debugme", "DEoptimR", "desc", "devtools", "dichromat", "diffobj", "digest", "dlm", "docstring", "doParallel", "DOSE", "doSNOW", "downloader", "DT", "dtplyr", "easypackages", "enrichplot", "ensembldb", "evaluate", "exomeCopy", "ExperimentHub", "ExperimentHubData", "extrafont", "extrafontdb", "fansi", "fastmap", "fastmatch", "fBasics", "feather", "fontawesome", "forcats", "forecast", "foreign", "formattable", "Formula", "fracdiff", "fresh", "fs", "future", "future.apply", "future.callr", "gamm4", "gdata", "gdtools", "gert", "getopt", "GetoptLong", "ggeffects", "ggforce", "ggfortify", "ggfun", "ggnewscale", "ggplotify", "ggraph", "ggrepel", "ggtext", "ggthemes", "ggtree", "ggupset", "gh", "gitcreds", "GlobalOptions", "globals", "glue", "gmodels", "gmp", "GO.db", "gplots", "graph", "graphics", "graphlayouts", "grDevices", "grid", "gridBase", "gridGraphics", "gridtext", "gss", "gsubfn", "gt", "gtools", "haven", "here", "hexbin", "highr", "Hmisc", "hrbrthemes", "htmlTable", "htmltools", "htmlwidgets", "httpuv", "httr", "igraph", "impute", "ini", "insight", "interactiveDisplayBase", "ipred", "irlba", "jpeg", "jquerylib", "jsonlite", "kableExtra", "KEGGgraph", "keras", "kernlab", "KernSmooth", "KFAS", "ks", "labelVector", "lars", "later", "lattice", "latticeExtra", "lazyeval", "lintr", "listenv", "lme4", "lmodel2", "lmtest", "lubridate", "magick", "magrittr", "mapdata", "mapproj", "maps", "maptools", "markdown", "Matrix", "MatrixModels", "mdthemes", "memoise", "methods", "MethylAid", "mlr", "modelr", "modeltools", "narray", "needs", "nlme", "nloptr", "nnet", "numDeriv", "openssl", "openxlsx", "optparse", "OrganismDbi", "packrat", "pals", "parallel", "parallelMap", "parameters", "ParamHelpers", "patchwork", "pathview", "pbkrtest", "pdftools", "performance", "pillar", "pingr", "pkgbuild", "pkgconfig", "pkgload", "plot3D", "plotly", "plotrix", "Polychrome", "polyclip", "polynom", "pracma", "praise", "prettyunits", "pROC", "processx", "prodlim", "progressr", "promises", "ProtGenerics", "proto", "proxy", "ps", "purrr", "qdapTools", "qpdf", "qs", "quantreg", "qvalue", "R.cache", "raster", "RBGL", "Rcpp", "RcppArmadillo", "RcppParallel", "readxl", "recipes", "redcapAPI", "rematch", "rematch2", "remotes", "reprex", "reticulate", "rex", "Rgraphviz", "rio", "rJava", "rlang", "rmarkdown", "Rmpfr", "RMySQL", "robustbase", "ROC", "ROCR", "roxygen2", "rpart", "RPMM", "rprojroot", "RSpectra", "rstudioapi", "Rtsne", "rversions", "rvest", "s2", "sass", "scatterplot3d", "selectr", "sessioninfo", "shadowtext", "shiny", "sjmisc", "sourcetools", "spatial", "sqldf", "stats", "stats4", "stringdist", "stringi", "stringr", "strucchange", "svglite", "sys", "systemfonts", "targets", "tcltk", "tensorflow", "terra", "testit", "testthat", "TH.data", "tibble", "tidygraph", "tidytree", "tidyverse", "tiff", "timeDate", "timeSeries", "tinytex", "tools", "torch", "treeio", "tseries", "units", "utf8", "utils", "uuid", "webshot", "webshot2", "websocket", "withr", "xfun", "xml2", "xmlparsedata", "xopen", "xtable", "yaml", "zip", "zoo"
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

#message("Version of Clang:"); system("gcc --version")

if (checkRequire("BiocManager")) {
    install.packages("BiocManager", dependencies = T, verbose = T, ask = F)
}
loadLibrary("devtools")
loadLibrary("librarian")
loadLibrary("BiocManager")
if (checkRequire("Biobase")) {BiocManager::install("Biobase", update = F, ask = F)}
loadLibrary("Biobase")

if (checkRequire("mapview")) {remotes::install_github("r-spatial/mapview", dependencies = T, upgrade="never")}

#if(checkRequire("S4Vectors")){BiocManager::install("S4Vectors", update = F, ask = F)}
#if(checkRequire("MASS")){install.packages("MASS",type = "binary", update = F, ask = F)}
# Load/install missing pacakges without asking
supM <- function(pk){return(suppressPackageStartupMessages(suppressWarnings(pk)))}

message("Loading packages:", paste0(capture.output(corePkgs), collapse = "\n"))
message("...loading")
supM(librarian::shelf(corePkgs, ask = F, update_all = F, quiet = FALSE))

message("Loading packages:",paste0(capture.output(preReqPkgs), collapse = "\n"))
message("...loading")
supM(librarian::shelf(preReqPkgs, ask = F, update_all = F, quiet = FALSE))

message("Loading BioConductor Packages and IlluminaHumanMethylation Manifest...")
supM(librarian::shelf(biocPkgs, ask = F, update_all = F, quiet = FALSE))

if(checkRequire("IlluminaHumanMethylationEPICmanifest")){
    devtools::install_github(
        repo = "mwsill/IlluminaHumanMethylationEPICmanifest",
        dependencies = T, upgrade = "never"
    )
}

if(checkRequire("mgmtstp27")){
    gitLink <- "https://github.com/badozor/mgmtstp27/raw/master/archive/mgmtstp27_0.6-3.tar.gz"
    install.packages(gitLink, repos = NULL, dependencies = T, verbose = T, type = "source", ask = F)
}


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
    closeAllConnections()
}

fixNeeds <- function (){
    sysfile <- system.file("extdata", "promptUser", package="needs")
    try(write(0, file=sysfile),silent=T)
    options(needs.promptUser=FALSE)
    invisible(needs:::autoload(TRUE))
    closeAllConnections()
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
                                     dependencies=T,verbose=T,upgrade="never")
            fixNeeds()
            try(fixProf(), silent=T)
        },
        warning=function(cond){
            devtools::install_github(
                "joshkatz/needs", ref = "development",
                                     dependencies=T,verbose=T,upgrade="never")
            fixNeeds()
            try(fixProf(), silent=T)
        }
    )
}

checkNeeds()
closeAllConnections()
if(Sys.info()[['sysname']]=="Darwin"){
    isGdal <- paste(system("echo `gdalinfo --version`", intern = T))

    if(!exists("isGdal")){
        system("brew install pkg-config")
        system("brew install gdal")
    }

    isProj <- system("which proj", intern=T)

    if(!exists("isProj")){
        system("brew install pkg-config")
        system("brew install proj")
    }
}

spat_config <- '--with-proj-lib=/usr/local/lib/ --with-proj-include=/usr/local/include/'
options(configure.args = c("sf" = spat_config, "rgdal" = spat_config))
if(checkRequire("sf")){
    tryCatch(
        install.packages(c("sf"), type = "source", dependencies=T, verbose=T, Ncpus = 4),
        error=function(e){
            remotes::install_github("r-spatial/sf", configure.args = "--with-proj-lib=/usr/local/lib/", dependencies=T, upgrade="never")
        }
    )
}

if(checkRequire("rgdal")){
    install.packages("rgdal", configure.args = c("--with-proj-lib=/usr/local/lib/", "--with-proj-include=/usr/local/include/"), type = "source", dependencies=T, verbose=T, Ncpus = 4)

}
invisible(gc())
if(Sys.info()[['sysname']]=="Darwin"){
    system("export LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib")
}

options("install.packages.compile.from.source" = "No")
options("install.packages.check.source" = "no")
loadLibrary("BiocManager")
loadLibrary("Biobase")

terraDep <- c('tinytest', 'ncdf4', 'leaflet')
suppressWarnings(librarian::shelf(terraDep, ask = F, update_all = F, quiet = FALSE))
if(checkRequire("terra")) {
    install.packages('terra', repos = 'https://rspatial.r-universe.dev', dependencies = T, verbose = T, Ncpus = 4)
}

if(checkRequire("FField")){
    gitLink <- "https://cran.r-project.org/src/contrib/Archive/FField/FField_0.1.0.tar.gz"
    install.packages(gitLink, repos = NULL, dependencies = T, verbose = T, type = "source", ask = F)
}

if (checkRequire("GenVisR")) {
    devtools::install_github("griffithlab/GenVisR", dependencies = T, upgrade = "never")
}
suppressWarnings(librarian::shelf(pkgs, ask = F, update_all = F, quiet = FALSE))
invisible(gc())


cbioLn <- switch (Sys.info()[['sysname']],
                  "Darwin" = "/Volumes/CBioinformatics/Methylation/classifiers",
                  "Linux" = "~/molecpathlab/production/Methylation/classifiers"
)

if (checkRequire("UniD")) {
    try(install.packages(file.path(cbioLn, "UniD"), type="source", repos=NULL), silent=T)
}

closeAllConnections()
