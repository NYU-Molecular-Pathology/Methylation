options("install.packages.compile.from.source" = "No")
options("install.packages.check.source" = "no")
options(warn = -1)
# Load/install missing pacakges without asking
supM <- function(pk){return(suppressPackageStartupMessages(suppressWarnings(pk)))}

# Setting US CRAN REPO
#rlis = getOption("repos")
#rlis["CRAN"] = "http://cran.us.r-project.org"
#options(repos = rlis)
options(BioC_mirror = "https://packagemanager.rstudio.com/bioconductor")
options(repos = c(CRAN = "https://packagemanager.rstudio.com/cran/__linux__/focal/latest"))

loadLibrary <- function(pkgName) {
    supM(library(pkgName, quietly = T, logical.return = T, warn.conflicts = F, character.only = T))
}

checkRequire <- function(pkgName){
    return(suppressWarnings(!require(pkgName, character.only = T, warn.conflicts = F)))
}

in.pkg <- function(pkg, character.only=T){
    install.packages(pkg, dependencies=T, verbose=T, ask=F, quiet = F)
}

CheckPackages <- function(pkgList) {
    toDrop <- pkgList %in% rownames(installed.packages())
    pkgLiSub <- pkgList[!toDrop]
    if (length(pkgLiSub) > 0) {
        tryCatch(
            pak::pkg_install(pkgLiSub, ask = F, lib = '/usr/local/lib/R/site-library/', dependencies=T),
            error = function(e) {
                tryCatch(
                    in.pkg(pkgLiSub),
                    error = function(e) {
                        BiocManager::install(pkgLiSub, update = F, ask = F, dependencies=T)
                    }
                )
            }
        )
    }
    supM(librarian::shelf(
        pkgList,
        ask = F,
        dependencies=T,
        update_all = F,
        quiet = FALSE
    ))
}

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
    "gridtext"
)

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
    #"prettyunits",
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
    "zoo",
    "mlr", 
    "wateRmelon", 
    "RPMM", 
    "impute"
)

loadLibrary("librarian")
loadLibrary("BiocManager")
loadLibrary("Biobase")
options(BioC_mirror = "https://packagemanager.rstudio.com/bioconductor")
options(repos = c(CRAN = "https://packagemanager.rstudio.com/cran/__linux__/focal/latest"))
if(checkRequire("arrow")) {CheckPackages('arrow')}

message("Librarian Installing pkgs1...")
CheckPackages(pkgs1)
message("Librarian Installing pkgs2...")
CheckPackages(pkgs2)
message("Librarian Installing pkgs3...")
CheckPackages(pkgs3)
pak::pkg_install("mwsill/IlluminaHumanMethylationEPICmanifest" , ask=F , dependencies = T)
invisible(gc())


