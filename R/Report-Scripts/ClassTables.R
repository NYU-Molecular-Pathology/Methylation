#!/usr/bin/env Rscript
## Script name: ClassTables.R
## Purpose: source of global functions used for methylation report.Rmd knitting
## Date Created: Aug 5, 2021
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH, 2025

gb <- globalenv() 
assign("gb", gb)
supM <- function(objTing){return(suppressMessages(suppressWarnings(objTing)))}
assign(x = "supM", value = supM, envir = .GlobalEnv)
options(install.packages.compile.from.source = "--no-multiarch")
options(Ncpus = 4)

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


pkgs <- c(
    "needs",
    "knitr",
    "jsonlite",
    "crayon",
    "RCurl",
    "ggplot2",
    "kableExtra",
    "tidyverse",
    "plotly",
    "pkgdown",
    "magrittr",
    "compiler",
    "minfi",
    "pander",
    "parallel",
    "doParallel",
    "randomForest",
    "glmnet",
    "gridExtra",
    "Rtsne",
    "dplyr"
)

optsLi <- list(
    error = FALSE,
    echo = FALSE,
    message = FALSE,
    warning = FALSE,
    fig.align = "left",
    self.contained = TRUE,
    comment = '',
    label_list = FALSE,
    highlight = FALSE,
    quiet = TRUE,
    verbose = FALSE,
    progress = TRUE
)

chunkOpts <- list(error = FALSE, echo = FALSE, message = FALSE, warning = FALSE, self.contained = TRUE, comment = '')

LoadReportPkgs <- function(pkgs, optsLi, chunkOpts){
    rlis = getOption("repos")
    rlis["CRAN"] = "http://cran.us.r-project.org"
    try(options(repos = rlis), silent=T)
    options("install.packages.compile.from.source" = "never")
    options("install.packages.check.source"="no")
    
    if(!is_installed("devtools")){
        install.packages("devtools", dependencies = c("Depends", "Imports", "LinkingTo"), quiet = T,
                         ask = FALSE, lib = .libPaths()[1], update = F)
        }
    if (!is_installed("BiocManager")) {
        install.packages("BiocManager", dependencies = c("Depends", "Imports", "LinkingTo"), quiet = T,
                         ask = FALSE, lib = .libPaths()[1], update = F)
        }
    if(!require("librarian")){
        install.packages("librarian", dependencies = c("Depends", "Imports", "LinkingTo"), quiet = T,
                         ask = FALSE, lib = .libPaths()[1], update = F)
        }
    
    library("BiocManager")
    librarian::shelf(pkgs, ask=F, verbose=F, warn.conflicts = F, quietly = T)
    require("needs", quietly = T, warn.conflicts=F)
    library(verbose=F, warn.conflicts = F, quietly = T, package = "knitr")
    knitr::opts_chunk$set(optsLi)
    knitr::opts_knit$set(optsLi)
    knitr::opts_current$set(optsLi)
    knitr::opts_knit$set(root.dir = getwd())
    options(width = 300, scipen = 5)
    library(verbose = F, warn.conflicts = F, quietly = T, package = "pander")
    pander::panderOptions('table.alignment.default', "left")
    require("compiler")
    invisible(supM(compiler::compilePKGS(enable = TRUE)))
    knitr::opts_chunk$set(chunkOpts)
    supM(library(verbose=F, warn.conflicts = F, quietly = T, package="ggplot2"))
}

LoadReportPkgs(pkgs, optsLi, chunkOpts)


FormatSuppInfo <- function(suppinfo){
    names(suppinfo) <- c(
        "ID",
        "Batch ID",
        "Tech",
        "B-number",
        "MP-TN#",
        "Sentrix ID",
        "Array type",
        "Material",
        "Gender"
    )
    suppinfo <- as.data.frame(suppinfo)
    colnames(suppinfo) <- ""
    return(suppinfo)
}


GetSuppInfo <- function(dat, RGset, msetDat) {
    suppinfo <- c(
        paste(dat$sampleID),
        paste(dat$run_id),
        paste(dat$tech),
        paste(dat$bnumber),
        paste(dat$mp_number),
        colnames(RGset),
        minfi::annotation(RGset)[[1]],
        msetDat$FFPE,
        msetDat$sex
    )
    suppinfo <- FormatSuppInfo(suppinfo)
    return(suppinfo)
}


GetUniDTables <- function(predU){
    be = c(booktabs = T, escape = F, linesep = "")
    btso = c("bordered")
    kgb <- c("striped",font_size = 14, bootstrap_options = btso, position = "left")
    uniTable1 <- predU[,2:9]
    uniTable2 <- predU[,10:14]
    colnames(uniTable2) <- c("Predicted Subtype", "C_prob", "M_prob", "P_prob", "missing.probe.subtype")
    unitab1 <- uniTable1 %>% mutate_all(as.character) %>% knitr::kable("html", be, align = 'clc') %>%
        kableExtra::kable_styling(kgb, full_width = F, position="left")
    unitab2 <- uniTable2 %>% mutate_all(as.character) %>% knitr::kable("html", be, align = 'clc') %>%
        kableExtra::kable_styling(kgb, full_width = F, position="left")
    return(list("tab1"=unitab1,"tab2"=unitab2))
}


GetClassProbTables <-  function(outList){
    xtraCss1 = "border-radius:0px;border-width:1px;border-style:solid;border-color:rgb(192,192,192);"
    xtraCss2 = "border-radius:0px;border-width:3px;border-style:solid;border-color:rgb(26,105,16);"
    xtraCss3 = "border-radius:0px;border-width:3px;border-style:solid;border-color:rgb(16,28,105);"
    txtc = "text-align:center;"
    btso = c("bordered")
    kgb <- c("striped",font_size = 14, bootstrap_options = btso, position = "left")
    kgh <- c(booktabs = T, escape = F, linesep = "")
    # all(names(outList) == c("super family", "family", "class", "subclass"))
    if (all(c("predicted", "maxscore") %in% names(outList))) {
        outList$maxscore <- as.character(outList$maxscore)
        return(outList %>%
                   knitr::kable("html", kgh, align='clc') %>%
                   kableExtra::kable_styling(kgb, full_width = F, position="left") %>%
                   kableExtra::column_spec(column = c(1, 2), extra_css = xtraCss1) %>%
                   kableExtra::column_spec(column = 2, background = "rgb(204,230,255)", extra_css = txtc) %>%
                   kableExtra::row_spec(row = 0, font_size = 16, background = "rgb(135,174,237)", color = "black")
               )
    } else{
        out_class_family=outList$out_class_family
        out <- outList$out
    }
    stopifnot(!is.null(outList) & !is.null(out_class_family))
    
    famTable <- out_class_family %>%
        knitr::kable("html", kgh, align = 'clc') %>%
        kableExtra::kable_styling(kgb, full_width = F, position="float_left") %>%
        kableExtra::column_spec(column=c(1,2),extra_css=xtraCss1) %>%
        kableExtra::column_spec(column = 2, background = "rgb(204,255,204)", extra_css = txtc) %>%
        kableExtra::row_spec(row = 0, font_size = 16, background = "rgb(127,217,126)", color = "black") %>%
        kableExtra::row_spec(row = 1, extra_css = xtraCss2)
    
    grpTable <- out %>%
        knitr::kable("html", kgh, align='clc') %>%
        kableExtra::kable_styling(kgb, full_width = F) %>%
        kableExtra::column_spec(column = c(1, 2), extra_css = xtraCss1) %>%
        kableExtra::column_spec(column = 2, background = "rgb(204,230,255)", extra_css = txtc) %>%
        kableExtra::row_spec(row = 0, font_size = 16, background = "rgb(135,174,237)", color = "black") %>%
        kableExtra::row_spec(row=1,extra_css=xtraCss3)
    
    return(list("famTable" = famTable, "grpTable" = grpTable))
}


SuppInfoTable <- function(dat, RGset, msetDat){
    suppinfo <- GetSuppInfo(dat, RGset, msetDat)
    suppinfo <- as.data.frame(t(suppinfo))
    rownames(suppinfo) <- NULL
    kgb <- c("striped", font_size = 9, bootstrap_options = c("bordered"), position = "float_left")
    totCol <- ncol(suppinfo)
    suppTab <- suppinfo %>%
        knitr::kable("html",c(booktabs = T, escape = F, linesep = ""), align='clc') %>%
        kableExtra::kable_styling(kgb, full_width = F, position="left") %>%
        column_spec(column = c(1:totCol), width = "200px") %>%
        column_spec(1, border_left = "3px solid white") %>%
        column_spec(5, border_right = "3px solid white") %>%
        kableExtra::row_spec(row = 0, font_size = 12, background="rgb(255, 235, 205)", color = "black") %>%
        kableExtra::row_spec(row = 0, background="rgb(255, 235, 205)") %>%
        kableExtra::row_spec(row = 1, font_size = 11)
    return(suppTab)
}


PrintGainLoss <- function(gnLss){
    knitr::asis_output('<p class="cnvdesc"></p> <hr class="dotted"><h4>CNV Gains and Loss</h4>')
    if (nrow(gnLss$gainDf) > 0) {knitr::asis_output(gnLss$gainTab)}
    if (nrow(gnLss$lossDf) > 0) {knitr::asis_output(gnLss$lossTab)}
    knitr::asis_output('<hr class="solid">')
}


GrabClassDetails <- function(refData, sclass) {
  sampleMatch <- refData$Molecular_subclass == sclass$predicted
  if (any(sampleMatch)) {
    theRefLi <- refData$description[sampleMatch]
  } else{
    theRefLi <- "No description availible."
  }
  return(theRefLi)
}
