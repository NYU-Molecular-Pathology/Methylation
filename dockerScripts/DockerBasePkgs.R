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

CheckPackages <- function(pkgList) {
    toDrop <- pkgList %in% rownames(installed.packages())
    pkgLiSub <- pkgList[!toDrop]
    if (length(pkgLiSub) > 0) {
        tryCatch(
            pak::pkg_install(pkgLiSub, ask = F, lib = '/usr/local/lib/R/site-library/'),
            error = function(e) {
                tryCatch(
                    in.pkg(pkgLiSub),
                    error = function(e) {
                        BiocManager::install(pkgLiSub, update = F, ask = F)
                    }
                )
            }
        )
    }
    supM(librarian::shelf(
        pkgList,
        ask = F,
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
    "mlr", "wateRmelon", "RPMM", "impute"
)

loadLibrary("librarian")

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
CheckPackages(pkgs1)
message("Librarian Installing pkgs2...")
CheckPackages(pkgs2)
message("Librarian Installing pkgs3...")
CheckPackages(pkgs3)
invisible(gc())
