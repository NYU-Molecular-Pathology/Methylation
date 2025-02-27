#!/usr/bin/env Rscript

set_environment <- function() {
    if (Sys.info()[['sysname']]=="Darwin") {
        options(Ncpus = 4)
        jdk_opt <- "/usr/local/opt/openjdk/bin/java"
        cmd <- "-version 2>&1 | awk '/version/ {print $3}' | tr -d '\"'"
        java_version <- system(paste(jdk_opt, cmd), intern = TRUE)
        jdk_path <- "/usr/local/Cellar/openjdk"
        jdk_exec <- "libexec/openjdk.jdk/Contents/Home"
        java_home_path <- file.path(jdk_path, java_version, jdk_exec)
        Sys.setenv(JAVA_HOME = java_home_path)
        arrow_path <- "/usr/local/Cellar/apache-arrow"
        arrow_version <- system(paste("ls", arrow_path, "|", "head -n 1"), intern = TRUE)
        pkg_config_path <- file.path(arrow_path, arrow_version, "lib/pkgconfig")
        Sys.setenv(PATH = paste("/usr/local/opt/gcc/bin", Sys.getenv("PATH"), sep = ":"))
        Sys.setenv(PKG_CONFIG_PATH = pkg_config_path)
    }
}


is_installed <- function(package_name) {
    tryCatch(
        expr = {
            return(length(find.package(package_name, quiet = TRUE)) > 0)
        },
        error = function(e) {
            return(FALSE)
        }
    )
}


install_pkgs <- function(pkg_list, ...) {
    inst_params <- list(dependencies = c("Depends", "Imports", "LinkingTo"),
                        verbose = TRUE, clean = TRUE, type = "both", quietly = TRUE,
                        ask = FALSE, lib = .libPaths()[1], update = F)
    final_params <- c(list(pkgs = pkg_list), list(...), inst_params)
    do.call("install.packages", final_params)
}


library_inst <- function(pkg_list) {
    inst_params <- list(pkgs = pkg_list, ask = F, update_all = F, quiet = F, type = "both",
                       dependencies = c("Depends", "Imports", "LinkingTo"))
    do.call("shelf", inst_params)
}


inst_load_pkg <- function(pkg){
    if (!is_installed(pkg)) {
        install_pkgs(pkg, repos = "https://cran.r-project.org")
    }
    library(pkg, character.only = T)
}


Setup_install <- function(){
    set_environment()
    inst_load_pkg("devtools")
    inst_load_pkg("BiocManager")
    inst_load_pkg("librarian")
}

your_packages <-
    c("glmnet",
      "ggplot2",
      "gridExtra",
      "knitr"
      )

Setup_install()
library_inst(your_packages)

