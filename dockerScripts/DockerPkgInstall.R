
corePkgs <- c(
    "randomForest",
    "glmnet",
    "ggplot2",
    "gridExtra",
    "knitr",
    "pander",
    "gmp",
    "minfi",
    "locfit",
    "parallel",
    "iterators",
    "foreach",
    "Biostrings",
    "matrixStats",
    "GenomicRanges",
    "GenomeInfoDb",
    "stats4",
    "stringr",
    "tidyr",
    "readxl",
    "Biobase",
    "BiocGenerics",
    "stats",
    "graphics",
    "grDevices",
    "utils",
    "Rhdf5lib"
)

preReqPkgs <- c(
    'MASS',
    'ade4',
    'dbplyr',
    'filelock',
    'png',
    'rjson',
    'generics',
    'tidyselect',
    'RCurl',
    'XML',
    'bit64',
    'plyr',
    'base64',
    'GenomicAlignments',
    'readr',
    'dplyr',
    'tidyr',
    'R.utils',
    'BH',
    'affy',
    'GenomicFeatures',
    'GenomicRanges',
    'annotate',
    'preprocessCore',
    'RSQLite',
    'DBI',
    'AnnotationDbi',
    'reshape2',
    'matrixStats',
    'S4Vectors',
    'IRanges',
    'GenomeInfoDb',
    'SummarizedExperiment',
    'genefilter',
    'rtracklayer',
    'Biostrings',
    'bumphunter',
    'limma',
    'reshape',
    'quadprog',
    'data.table',
    'HDF5Array',
    'BiocParallel'
)

biocPkgs <- c(
    "lumi",
    "methylumi",
    "conumee",
    "minfi",
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
)
# Options Set ------------------------------------------------
options("install.packages.compile.from.source" = "No")
options("install.packages.check.source" = "no")
options(repos="https://packagemanager.rstudio.com/all/__linux__/focal/latest")
options(warn = -1)
# Load/install missing pacakges without asking ------------------------------------------------
supM <- function(pk){return(suppressPackageStartupMessages(suppressWarnings(pk)))}
loadLibrary <- function(pkgName) {supM(library(pkgName, quietly = T, logical.return = T, warn.conflicts = F, character.only = T))}
checkRequire <- function(pkgName){return(suppressWarnings(!require(pkgName, character.only = T, warn.conflicts = F)))}
in.pkg <- function(pkg, character.only=T){install.packages(pkg, dependencies=T, verbose=T, ask=F, quiet = F)}
# Check Main Packages ------------------------------------------------
if(checkRequire("devtools")){in.pkg("devtools")};loadLibrary("devtools")
if(checkRequire("remotes")){in.pkg("remotes")}
if(checkRequire("librarian")){in.pkg("librarian")}
if(checkRequire("pak")) {install.packages(lib = '/usr/local/lib/R/site-library/', 'pak', repos = 'https://r-lib.github.io/p/pak/dev/', dependencies = T)}
if(checkRequire("BiocManager")){in.pkg("BiocManager")};loadLibrary("BiocManager")
if(checkRequire("Biobase")) {BiocManager::install("Biobase", update = F, ask = F)}
loadLibrary("Biobase");loadLibrary("librarian");loadLibrary("BiocManager");loadLibrary("Biobase")

supM(librarian::shelf(c('tinytest', 'ncdf4', 'leaflet'), ask = F, update_all = F, quiet = FALSE))
if(checkRequire("terra")) {install.packages('terra', repos = 'https://rspatial.r-universe.dev', dependencies = T, verbose = T)}

if(checkRequire("FField")){
    gitLink <- "https://cran.r-project.org/src/contrib/Archive/FField/FField_0.1.0.tar.gz"
    install.packages(gitLink, repos = NULL, dependencies = T, verbose = T, type = "source", ask = F)
}

CatchPak <- function(pkgLiSub) {
    tryCatch(
        pak::pkg_install(pkgLiSub, ask = F, lib = '/usr/local/lib/R/site-library/'),
        error = function(e) {tryCatch(in.pkg(pkgLiSub), error = function(e) {
            BiocManager::install(pkgLiSub, update = F, ask = F)
            })}
        )
}


CheckPackages <- function(pkgList) {
    toDrop <- pkgList %in% rownames(installed.packages())
    pkgLiSub <- pkgList[!toDrop]
    if (length(pkgLiSub) > 0) {
        CatchPak(pkgLiSub)
    }
    supM(librarian::shelf(
        pkgList,
        ask = F,
        update_all = F,
        quiet = F
    ))
}


message("Librarian Installing corePkgs...")
CheckPackages(corePkgs)
message("Librarian Installing preReqPkgs...")
CheckPackages(preReqPkgs)
message("Librarian Installing biocPkgs...")
CheckPackages(biocPkgs)
invisible(gc())

if(checkRequire("illuminaio")){BiocManager::install("illuminaio", update = F, ask = F, dependencies = T)}
if(checkRequire("mapview")){remotes::install_github("r-spatial/mapview", dependencies = T, upgrade="never")}
if(checkRequire("IlluminaHumanMethylationEPICmanifest")){
    pak::pkg_install("mwsill/IlluminaHumanMethylationEPICmanifest", ask=F, dependencies = T)
}
if(checkRequire("mgmtstp27")){
    gitLink <- "https://github.com/badozor/mgmtstp27/raw/master/archive/mgmtstp27_0.6-3.tar.gz"
    install.packages(gitLink, repos = NULL, dependencies = T, verbose = T, type = "source", ask = F)
}
if (!("needs" %in% rownames(installed.packages()))) {install.packages("needs", dependencies = T, verbose = T, ask = F)}
options(needs.promptUser = FALSE)
spat_config <- '--with-proj-lib=/usr/local/lib/ --with-proj-include=/usr/local/include/'
options(configure.args = c("sf" = spat_config, "rgdal" = spat_config))
if(checkRequire("sf")){tryCatch(
    install.packages(c("sf"), type = "source", dependencies=T, verbose=T),
    error=function(e){remotes::install_github("r-spatial/sf", configure.args = "--with-proj-lib=/usr/local/lib/", dependencies=T, upgrade="never")})
}
invisible(gc())

