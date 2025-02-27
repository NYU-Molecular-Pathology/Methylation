# Initialize
library(parallel)
no_cores <- detectCores(logical = FALSE)
formals(library)$quietly <- TRUE
formals(library)$warn.conflicts <- FALSE
formals(require)$warn.conflicts <- FALSE

if(Sys.info()[['sysname']]!="Darwin"){
    options(BioC_mirror = "https://packagemanager.rstudio.com/bioconductor")
    options(repos = c(CRAN = "https://packagemanager.rstudio.com/all/__linux__/jammy/latest"))
}

# Dataframe to store failed package installations
failed_packages <- data.frame(Package = character(0), Source = character(0), stringsAsFactors = FALSE)

# Install and load packages from CRAN
install_load_cran <- function(pkg) {
    tryCatch({
        if (!require(pkg, character.only = TRUE)) {
            install.packages(pkg, dependencies = TRUE, ask = FALSE)
            library(pkg, character.only = TRUE)
            message(pkg, " loaded successfully")
        }
    }, error = function(e) {
        failed_packages <<- rbind(failed_packages, data.frame(Package = pkg, Source = "CRAN"))
        message(pkg, " failed to load or install")
    })
}

# Install and load packages from Bioconductor
install_load_bioc <- function(pkg) {
    tryCatch({
        if (!require(pkg, character.only = TRUE)) {
            BiocManager::install(pkg, ask = FALSE, dependencies = TRUE)
            library(pkg, character.only = TRUE)
        }
    }, error = function(e) {
        failed_packages <<-
            rbind(failed_packages,
                  data.frame(Package = pkg, Source = "Bioconductor"))
    })
}

# Install and load packages from GitHub
install_load_github <- function(pkg, repo) {
    tryCatch({
        if (!require(pkg, character.only = TRUE)) {
            devtools::install_github(repo, dependencies = TRUE, ask = FALSE)
            library(pkg, character.only = TRUE)
        }
    }, error = function(e) {
        failed_packages <<- rbind(failed_packages, data.frame(Package = pkg, Source = "GitHub"))
    })
}

# Bioconductor packages
bioc_packages <- c('HDF5Array', 'bumphunter', 'GEOquery', 'minfi', 'lumi', 'rngtools', 'methylumi', 'randomForest', 'glmnet','IlluminaHumanMethylation450kmanifest', 'IlluminaHumanMethylation450kanno.ilmn12.hg19', 'IlluminaHumanMethylationEPICmanifest', 'Rtsne', 'IlluminaHumanMethylationEPICanno.ilm10b2.hg19', 'IlluminaHumanMethylationEPICanno.ilm10b4.hg19', 'MethylAid', 'conumee', 'BiocParallel', 'Biobase', 'limma', 'MatrixGenerics', 'GenomeInfoDb')

# CRAN packages
cran_packages <- c('devtools','doSnow', 'remotes', 'sjmisc', 'stringi', 'digest', 'RCurl', 'rlang', 'parallel', 'grid', 'gridExtra', 'knitr', 'kableExtra', 'ggplot2', 'plotly', 'ggfortify', 'ggrepel', 'gplots', 'fastmatch', 'pals', 'Polychrome', 'qdapTools', 'beepr', 'xtable', 'pander', 'grDevices', 'graphics', 'stats', 'utils', 'magick', 'ade4', 'MASS', 'R.utils', 'optparse', 'targets', 'usethis', 'webshot', 'reshape', 'reshape2', 'data.table', 'DT', 'scales', 'RColorBrewer', 'readxl', 'stringr', 'tinytex', 'future.apply', 'Rcpp', 'httpuv', 'shiny')

# GitHub packages
github_repos <- c(
    rmarkdown = 'rstudio/rmarkdown',
    docstring = 'dasonk/docstring',
    rstudioapi = 'rstudio/rstudioapi',
    easypackages = 'jakesherman/easypackages',
    redcapAPI = 'nutterb/redcapAPI',
    crayon = 'r-lib/crayon',
    redcap = 'epicentre-msf/redcap'
)

message("Loading Packages...")

# Parallel installation
if (!require("BiocManager")) install.packages("BiocManager", dependencies = TRUE, ask = FALSE)
install_load_cran("foreach")
install_load_cran("doParallel")

library('doParallel')
registerDoParallel(cores = no_cores - 1)

# CRAN Packages
unlist(foreach(pkg = cran_packages) %dopar% {install_load_cran(pkg)})

# Bioconductor Packages
unlist(foreach(pkg = bioc_packages) %dopar% {install_load_bioc(pkg)})

# GitHub Packages
unlist(foreach(pkg = names(github_repos)) %dopar% {install_load_github(pkg, github_repos[[pkg]])})

# Unregister parallel backend and print failed packages if any
stopImplicitCluster()
if (nrow(failed_packages) > 0) {
    message("Failed to install/load the following packages:")
    message(paste0(capture.output(failed_packages), collapse="\n"))
    stop("Failed to install/load packages!")
}

# Update any packages as needed
# if(Sys.info()[['sysname']]!="Darwin"){
#     update.packages(ask = FALSE, repos = 'https://cran.rstudio.com/', lib='/usr/local/lib/R/site-library', dependencies=TRUE)
#     base_pkgs <- c('base', 'compiler', 'datasets', 'graphics', 'grDevices', 'methods', 'stats', 'utils')
#     install.packages(base_pkgs, repos = 'https://cran.rstudio.com/', ask = FALSE, type = 'source', dependencies=TRUE, lib='/usr/local/lib/R/site-library')
# }
                    
