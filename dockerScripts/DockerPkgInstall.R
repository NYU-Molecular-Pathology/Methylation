options("install.packages.compile.from.source" = "No")
options("install.packages.check.source" = "no")
options(warn = -1)
# Load/install missing pacakges without asking
supM <- function(pk){return(suppressPackageStartupMessages(suppressWarnings(pk)))}

# Setting US CRAN REPO
rlis = getOption("repos")
rlis["CRAN"] = "http://cran.us.r-project.org"
options(repos = rlis)

loadLibrary <- function(pkgName) {
    supM(library(pkgName, quietly = T, logical.return = T, warn.conflicts = F, character.only = T))
}

checkRequire <- function(pkgName){
    return(suppressWarnings(!require(pkgName, character.only = T, warn.conflicts = F)))
}

in.pkg <- function(pkg, character.only=T){
    install.packages(pkg, dependencies=T, verbose=T, ask=F, quiet = F)
}

corePkgs <-
    c(
        #"randomForest",
        #"glmnet",
        "ggplot2",
        "gridExtra",
        "knitr",
        "pander",
        "gmp",
        #"Matrix",
        "minfi",
        #"bumphunter",
        "locfit",
        "parallel",
        "iterators",
        "foreach",
        "Biostrings",
        #"XVector",
        #"MatrixGenerics",
        "matrixStats",
        "GenomicRanges",
        "GenomeInfoDb",
        "stats4",
        #"remotes",
        "stringr",
        "tidyr",
        "readxl",
        "Biobase",
        "BiocGenerics",
        #"usethis",
        "stats",
        "graphics",
        "grDevices",
        "utils",
        #"datasets",
        #"methods",
        #"base",
        "Rhdf5lib"
    )

preReqPkgs <- c(
    'MASS',
    'ade4',
    'dbplyr',
    'filelock',
    #'formatR',
    #'bitops',
    #'progress',
    #'BiocFileCache',
    #'bit',
    'png',
    #'Rhtslib',
    'rjson',
    #'rngtools',
    #'hms',
    #'vroom',
    #'tzdb',
    'generics',
    'tidyselect',
    #'R.oo',
    #'R.methodsS3',
    #'lambda.r',
    #'futile.options',
    #'affyio',
    #'zlibbioc',
    'RCurl',
    #'XVector',
    #'BiocIO',
    #'biomaRt',
    'XML',
    'bit64',
    #'blob',
    #'plogr',
    #'KEGGREST',
    'plyr',
    #'TxDb.Hsapiens.UCSC.hg19.knownGene',
    #'org.Hs.eg.db',
    #'GenomeInfoDbData',
    #'MatrixGenerics',
    'base64',
    #'Rsamtools',
    'GenomicAlignments',
    #'restfulr',
    #'locfit',
    #'doRNG',
    #'multtest',
    #'scrime',
    #'sparseMatrixStats',
    'readr',
    'dplyr',
    'tidyr',
    'R.utils',
    #'rhdf5',
    #'rhdf5filters',
    #'Rhdf5lib',
    #'futile.logger',
    #'snow',
    'BH',
    'affy',
    'GenomicFeatures',
    'GenomicRanges',
    'annotate',
    #'nleqslv',
    'preprocessCore',
    'RSQLite',
    'DBI',
    'AnnotationDbi',
    'reshape2',
    'matrixStats',
    #'FDb.InfiniumMethylation.hg19',
    'S4Vectors',
    'IRanges',
    'GenomeInfoDb',
    'SummarizedExperiment',
    'genefilter',
    #'illuminaio',
    #'IlluminaHumanMethylationEPICanno.ilm10b2.hg19',
    #'DNAcopy',
    'rtracklayer',
    'Biostrings',
    'bumphunter',
    #'beanplot',
    #'nor1mix',
    #'siggenes',
    #'tvthemes',
    'limma',
    #'DelayedMatrixStats',
    #'mclust',
    'reshape',
    'quadprog',
    'data.table',
    #'GEOquery',
    #'DelayedArray',
    'HDF5Array',
    'BiocParallel'
)

pkgs1 <- c(
    "abind",
    "animation",
    "arrow",
    #"askpass",
    "beepr",
    "BiocStyle",
    "biocViews",
    "bookdown",
    "brio",
    "Cairo",
    "caret",
    "chromote",
    "chron",
    "cli",
    "clipr",
    "colorspace",
    "compiler",
    "CopyNumberPlots",
    "cpp11",
    "crayon",
    "curl",
    "dichromat",
    "diffobj",
    "digest",
    "doParallel",
    "DT"
    )

pkgs2 <- c(
    "easypackages",
    "evaluate",
    "extrafont",
    "extrafontdb",
    "fontawesome",
    "forecast",
    "formattable",
    "fs",
    "gdata",
    "gdtools",
    "getopt",
    "ggforce",
    "ggfortify",
    "ggnewscale",
    "ggplotify",
    "ggraph",
    "ggrepel",
    "ggtext",
    "ggthemes",
    "gh",
    "GlobalOptions",
    "graphlayouts",
    "grid",
    "gridBase",
    "gridGraphics",
    "gridtext")

pkgs3 <- c(
    "Hmisc",
    "htmlTable",
    "htmltools",
    "htmlwidgets",
    "httr",
    "igraph",
    "jpeg",
    "jquerylib",
    "jsonlite",
    "kableExtra",
    "magrittr",
    "markdown",
    "MethylAid",
    "needs",
    "openxlsx",
    "pals",
    "parameters",
    "plotly",
    "Polychrome",
    "prettyunits",
    "raster",
    "Rcpp",
    "redcapAPI",
    "rmarkdown",
    "Rtsne",
    "sjmisc",
    "sourcetools",
    "spatial",
    "sqldf",
    "stringi",
    "systemfonts",
    "targets",
    "terra",
    "testit",
    "tibble",
    "tidyverse",
    "tinytex",
    "utf8",
    "uuid",
    "webshot2",
    "websocket",
    "xtable",
    "yaml",
    "zip",
    "zoo"
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

if(checkRequire("devtools")){in.pkg("devtools")};loadLibrary("devtools")
if(checkRequire("remotes")){in.pkg("remotes")}
if(checkRequire("librarian")){in.pkg("librarian")}

if (checkRequire("BiocManager")){in.pkg("BiocManager")};loadLibrary("BiocManager")
if (checkRequire("Biobase")) {BiocManager::install("Biobase", update = F, ask = F)};loadLibrary("Biobase")

loadLibrary("librarian")

if(checkRequire("mapview")){remotes::install_github("r-spatial/mapview", dependencies = T, upgrade="never")}
message("Librarian Installing corePkgs...")
supM(librarian::shelf(corePkgs, ask = F, update_all = F, quiet = FALSE))
message("Librarian Installing preReqPkgs...")
supM(librarian::shelf(preReqPkgs, ask = F, update_all = F, quiet = FALSE))
message("Librarian Installing biocPkgs...")
supM(librarian::shelf(biocPkgs, ask = F, update_all = F, quiet = FALSE))

if(checkRequire("IlluminaHumanMethylationEPICmanifest")){
    devtools::install_github(repo = "mwsill/IlluminaHumanMethylationEPICmanifest",
                             dependencies = T, upgrade = "never")
}

if(checkRequire("mgmtstp27")){
    gitLink <- "https://github.com/badozor/mgmtstp27/raw/master/archive/mgmtstp27_0.6-3.tar.gz"
    install.packages(gitLink, repos = NULL, dependencies = T, verbose = T, type = "source", ask = F)
}

if (!("needs" %in% rownames(installed.packages()))) {
    install.packages("needs", dependencies = T, verbose = T, ask = F)
}
options(needs.promptUser = FALSE)

spat_config <- '--with-proj-lib=/usr/local/lib/ --with-proj-include=/usr/local/include/'
options(configure.args = c("sf" = spat_config, "rgdal" = spat_config))

if(checkRequire("sf")){
    tryCatch(install.packages(c("sf"), type = "source", dependencies=T, verbose=T),
             error=function(e){
                 remotes::install_github("r-spatial/sf", configure.args = "--with-proj-lib=/usr/local/lib/",
                                         dependencies=T, upgrade="never")})
}

invisible(gc())

options("install.packages.compile.from.source" = "No")
options("install.packages.check.source" = "no")
loadLibrary("BiocManager")
loadLibrary("Biobase")

terraDep <- c('tinytest', 'ncdf4', 'leaflet')
supM(librarian::shelf(terraDep, ask = F, update_all = F, quiet = FALSE))

if(checkRequire("terra")) {
    install.packages('terra', repos = 'https://rspatial.r-universe.dev', dependencies = T, verbose = T)
}

if(checkRequire("FField")){
    gitLink <- "https://cran.r-project.org/src/contrib/Archive/FField/FField_0.1.0.tar.gz"
    install.packages(gitLink, repos = NULL, dependencies = T, verbose = T, type = "source", ask = F)
}
message("Librarian Installing pkgs1...")
supM(librarian::shelf(pkgs1, ask = F, update_all = F, quiet = T))
message("Librarian Installing pkgs2...")
supM(librarian::shelf(pkgs2, ask = F, update_all = F, quiet = T))
message("Librarian Installing pkgs3...")
supM(librarian::shelf(pkgs3, ask = F, update_all = F, quiet = T))
invisible(gc())
