## Script name: DockerPkgInstall.R
## Purpose of script: Install Packages for methylation Pipeline
## Author: Jonathan Serrano
## Date Created: Dec-16-2022

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

ChromePkgs <- c("webshot2", "pagedown","curl","websocket", "curl", "jsonlite", "processx", "R6", "later", "promises", "magrittr", "rlang", "fastmap")

# Options Set ------------------------------------------------
options("install.packages.compile.from.source" = "No")
options("install.packages.check.source" = "no")
options(BioC_mirror = "https://packagemanager.rstudio.com/bioconductor")
options(repos = c(CRAN = "https://packagemanager.rstudio.com/cran/__linux__/focal/latest"))
options(warn = -1)

# Load/install missing pacakges without asking ------------------------------------------------
supM <- function(pk){return(suppressPackageStartupMessages(suppressWarnings(pk)))}
LoadLibrary <- function(pkgName) {supM(library(pkgName, quietly = T, logical.return = T, warn.conflicts = F, character.only = T))}
CheckReq <- function(pkgName){return(suppressWarnings(!require(pkgName, character.only = T, warn.conflicts = F)))}
InstPkg <- function(pkg, character.only=T){install.packages(pkg, dependencies=T, verbose=T, ask=F, quiet = F)}
# Check Main Packages ------------------------------------------------
if(CheckReq("devtools")){InstPkg("devtools")}
if(CheckReq("remotes")){InstPkg("remotes")}
if(CheckReq("librarian")){InstPkg("librarian")}
if(CheckReq("pak")) {install.packages(lib = '/usr/local/lib/R/site-library/', 'pak', repos = 'https://r-lib.github.io/p/pak/dev/', dependencies = T)}
if(CheckReq("BiocManager")){InstPkg("BiocManager")}
if(CheckReq("Biobase")) {BiocManager::install("Biobase", update = F, ask = F)}

LoadLibrary("devtools");LoadLibrary("librarian");LoadLibrary("BiocManager");LoadLibrary("Biobase")

supM(librarian::shelf(c('tinytest', 'ncdf4', 'leaflet'), ask = F, update_all = F, quiet = FALSE))
if(CheckReq("terra")) {install.packages('terra', repos = 'https://rspatial.r-universe.dev', dependencies = T, verbose = T)}
fieldLink <- "https://cran.r-project.org/src/contrib/Archive/FField/FField_0.1.0.tar.gz"
if(CheckReq("FField")){install.packages(fieldLink, repos = NULL, dependencies = T, verbose = T, type = "source", ask = F)}

CatchPak <- function(pkgLiSub) {
    tryCatch(pak::pkg_install(pkgLiSub, ask = F, lib = '/usr/local/lib/R/site-library/'),
             error = function(e) {tryCatch(InstPkg(pkgLiSub), error = function(e) {BiocManager::install(pkgLiSub, update = F, ask = F)})})
}


CheckPackages <- function(pkgList) {
    pkgLiSub <- pkgList[!(pkgList %in% rownames(installed.packages()))]
    if (length(pkgLiSub) > 0) {CatchPak(pkgLiSub)}
    supM(librarian::shelf(pkgList, ask = F, update_all = F, quiet = F))
}


message("Librarian Installing corePkgs...")
CheckPackages(corePkgs)
message("Librarian Installing preReqPkgs...")
CheckPackages(preReqPkgs)
message("Librarian Installing biocPkgs...")
CheckPackages(biocPkgs)
message("Librarian Installing ChromePkgs...")
CheckPackages(ChromePkgs)

if(CheckReq("illuminaio")){BiocManager::install("illuminaio", update = F, ask = F, dependencies = T)}
if(CheckReq("arrow")){install.packages("arrow", dependencies = T, verbose = T, ask = F)}
if(CheckReq("mapview")){remotes::install_github("r-spatial/mapview", dependencies = T, upgrade="never")}
if(CheckReq("IlluminaHumanMethylationEPICmanifest")){pak::pkg_install("mwsill/IlluminaHumanMethylationEPICmanifest", ask=F, dependencies = T)}
mgmLink <- "https://github.com/badozor/mgmtstp27/raw/master/archive/mgmtstp27_0.6-3.tar.gz"
if(CheckReq("mgmtstp27")){install.packages(mgmLink, repos = NULL, dependencies = T, verbose = T, type = "source", ask = F)}
if(CheckReq("needs")){install.packages("needs", dependencies = T, verbose = T, ask = F)}; options(needs.promptUser = FALSE)
if(CheckReq("chromote")){remotes::install_github("rstudio/chromote", dependencies = T, upgrade="never")}
spat_config <- '--with-proj-lib=/usr/local/lib/ --with-proj-include=/usr/local/include/'
options(configure.args = c("sf" = spat_config, "rgdal" = spat_config))
if(CheckReq("sf")){
    tryCatch(install.packages("sf", type = "source", dependencies=T, verbose=T), 
             error=function(e){
                 remotes::install_github("r-spatial/sf", configure.args = "--with-proj-lib=/usr/local/lib/", 
                                         dependencies=T, upgrade="never")})
}

gc()

