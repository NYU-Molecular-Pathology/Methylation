#!/usr/bin/env Rscript
## ---------------------------
## Script name: ClassTables.R
## Purpose: source of global functions used for methylation report.Rmd knitting
## Date Created: Aug 5, 2021
## Date Last Modified: February 5, 2024
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

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
    if (!is_installed("UniD")) {
        unid_deps <- c("plyr", "mlr", "wateRmelon", "glmnet", "RPMM", "minfi", "impute",
                   "IlluminaHumanMethylationEPICmanifest", "IlluminaHumanMethylation450kmanifest")
        librarian::shelf(unid_deps, ask = F, verbose = T, warn.conflicts = F, quietly = T)
        install.packages("/Volumes/CBioinformatics/Methylation/classifiers/UniD", type = "source",
                         dependencies = c("Depends", "Imports", "LinkingTo"), repo = NULL, quiet = T,
                         ask = FALSE, lib = .libPaths()[1], update = F)
    }
    library("BiocManager")
    library("UniD")
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

GetSexMsetBa <- function(RGset, FFPE = NULL) {
    is450k <- RGset@annotation[[1]] == "IlluminaHumanMethylation450k"
    if (is450k) {
        library(verbose=F, warn.conflicts = F, quietly = T, package = "IlluminaHumanMethylation450kmanifest")
        library(verbose=F, warn.conflicts = F, quietly = T, package= "mnp.v11b4")
        Mset <- mnp.v11b4::MNPpreprocessIllumina(RGset)
        if (is.null(FFPE)) {FFPE <- mnp.v11b4::MNPgetFFPE(RGset)}
        Mset_ba <- mnp.v11b4::MNPbatchadjust(Mset, FFPE)
        if (FFPE == "Frozen") {
            Mset@preprocessMethod <- c(Mset_ba@preprocessMethod, FFPE_Frozen.mnp.adjustment = '0.11')
            }
        sex <- ifelse(mnp.v11b4::MNPgetSex(Mset)$predictedSex == "M", "Male", "Female")
    } else {
        library(verbose=F, warn.conflicts = F, quietly = T, package= "IlluminaHumanMethylationEPICmanifest")
        library(verbose=F, warn.conflicts = F, quietly = T, package= "mnp.v11b6")
        Mset <- mnp.v11b6::MNPpreprocessIllumina(RGset)
        Mset@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19")
        if (is.null(FFPE)) {FFPE <- mnp.v11b6::MNPgetFFPE(RGset)}
        Mset_ba <- mnp.v11b6::MNPbatchadjust(Mset, FFPE)
        detP <- minfi::detectionP(RGset)
        bs <- minfi::getBeta(Mset)
        sexEstimate <- as.data.frame(signif(sest::get.proportion_table(bs, detP), digits = 2))
        yest <- as.double(sexEstimate$`p.Y:(-18,-5]`) >= 0.75
        yest1 <- as.double(sexEstimate$`Y:(0,0.1]`) >= 0.12
        sex <- ifelse((yest == TRUE && yest1 == TRUE), "male", "female")
    }
    return(list(
        "sex" = sex,
        "Mset" = Mset,
        "Mset_ba" = Mset_ba,
        "FFPE" = FFPE
    ))
}


SpecialMNPpredict6 <- function (betas, calibrate = T, type = "response", MCF = F) {
    betas <- betas[match(rownames(mnp.v11b6::rf.pred$importance), rownames(betas)), drop = F]
    if (calibrate) {
        pred <- stats::predict(mnp.v11b6::rf.pred, t(betas), type = "prob")
        pred2 <- matrix(stats::predict(mnp.v11b6::calfit, newx = pred, type = "response", 
                                       s = 0.0006363791)[, , 1], nrow = nrow(pred))
        colnames(pred2) <- colnames(pred)
        colnames(pred2) <- mnp.v11b6::reflist[match(colnames(pred2), paste0(mnp.v11b6::reflist[[4]])), 2]
        rownames(pred2) <- rownames(pred)
        out <- pred2
        if (MCF) {
            pred3 <- lapply(cgroups, function(x) rowSums(pred2[,x, drop = F]))
            pred3 <- do.call(cbind, pred3)
            pred2 <- pred2[, -which(colnames(pred2) %in% unlist(cgroups)), drop = F]
            pred2 <- cbind(pred3, pred2)
            out <- pred2
        }
        if (type == "response") {
            out <- factor(colnames(pred2)[apply(pred2, 1, which.max)], levels = colnames(pred2))
            names(out) <- rownames(pred2)
        }
    } else {
        if (type == "response") {
            out <- stats::predict(mnp.v11b6::rf.pred, t(betas), type)
            out <- mnp.v11b6::reflist[match(out, paste0(mnp.v11b6::reflist[[4]])), 2]
            names(out) <- colnames(betas)
        }
        if (type == "prob") {
            out <- stats::predict(mnp.v11b6::rf.pred, t(betas), type)
            colnames(out) <- mnp.v11b6::reflist[match(colnames(out), paste0(mnp.v11b6::reflist[[4]])), 2]
        }
    }
    return(out)
}


GetFamilyProb <- function(is450k, Mset_ba, Mset){
    if(is450k == T){
        library(
            verbose = F,
            warn.conflicts = F,
            quietly = T,
            package = "IlluminaHumanMethylation450kmanifest"
        )
        library(verbose=F, warn.conflicts = F, quietly = T, package= "mnp.v11b4")
        loadNamespace("mnp.v11b4")
        gb$rf.pred <- mnp.v11b4::rf.pred
        gb$calfit <- mnp.v11b4::calfit
        gb$refset.center <- mnp.v11b4::refset.center
        gb$reflist <- mnp.v11b4::reflist
        betas <- minfi::getBeta(Mset_ba[, 1])  
        probs_mcf <- tryCatch(
            expr = {
                probs_mcf <- mnp.v11b4::MNPpredict_betas(betas, type = 'prob', MCF = TRUE)
                return(probs_mcf)
            },
            error = function(e) {
                message("Error at GetFamilyProb in Brain Classifier v11 prediction mnp.v11b4::MNPpredict: \n", e)
                message("Using MNPpredict(Mset[, 1]) instead of Mset_ba\n")
                probs_mcf <- mnp.v11b4::MNPpredict(Mset[, 1], type = 'prob', MCF = TRUE)
                return(probs_mcf)
            }
        )
        return(probs_mcf)
    }else{
        library(verbose=F, warn.conflicts = F, quietly = T, package= "IlluminaHumanMethylationEPICmanifest")
        library(verbose=F, warn.conflicts = F, quietly = T, package= "mnp.v11b6")
        loadNamespace("mnp.v11b6")
        gb$rf.pred <- mnp.v11b6::rf.pred
        gb$calfit <- mnp.v11b6::calfit
        gb$refset.center <- mnp.v11b6::refset.center
        gb$reflist <- mnp.v11b6::reflist
        probs_mcf <- tryCatch(
            expr = {
                probs_mcf <- mnp.v11b6::MNPpredict(Mset_ba[, 1], type = 'prob', MCF = TRUE)
                return(probs_mcf)
            },
            error = function(e) {
                message("Error at GetFamilyProb in Brain Classifier v11 function mnp.v11b6::MNPpredict:\n", e)
                message("Trying SpecialMNPpredict6...")
                betas <- minfi::getBeta(Mset_ba[, 1], type = "Illumina")
                probs_mcf <- SpecialMNPpredict6(betas, type = 'prob', MCF = TRUE)
                return(probs_mcf)
            }
        )
    }
}


GetMnpV11Prob <- function(Mset_ba) {
    probs <- tryCatch(
        expr = {
            probs <- mnp.v11b6::MNPpredict(Mset_ba[, 1], type = 'prob')
            return(probs)
        },
        error = function(e) {
            message("Error at GetMnpV11Prob in Brain Classifier v11 function mnp.v11b6::MNPpredict:\n", e)
            message("Trying SpecialMNPpredict6...")
            betas <- minfi::getBeta(Mset_ba[, 1], type = "Illumina")
            probs <- SpecialMNPpredict6(betas, type = 'prob')
            return(probs)
        }
    )
    return(probs)
}

LoadMnpManifest <- function(is450k){
    if(is450k==T){
        library(verbose = F, warn.conflicts = F, quietly = T, package = "IlluminaHumanMethylation450kmanifest")
        library(verbose = F, warn.conflicts = F, quietly = T, package = "mnp.v11b4")
        try(unloadNamespace("mnp.v11b6"), silent = T)
        try(detach_package("mnp.v11b6"), silent = T)
        invisible(loadNamespace("mnp.v11b4")) 
        invisible(requireNamespace("mnp.v11b4"))
        require(warn.conflicts = F, quietly = T, package = "mnp.v11b4")
    }else{
        library(verbose = F, warn.conflicts = F, quietly = T, package = "IlluminaHumanMethylationEPICmanifest")
        library(verbose = F, warn.conflicts = F, quietly = T, package = "mnp.v11b6")
        try(unloadNamespace("mnp.v11b4"), silent = T)
        try(detach_package("mnp.v11b4"), silent = T)
        invisible(loadNamespace("mnp.v11b6"))
        invisible(requireNamespace("mnp.v11b6"))
        require(warn.conflicts = F, quietly = T, package = "mnp.v11b6")
    }
}


GetFamScore <- function(out_class_family){
    fsco <- as.numeric(paste0(out_class_family$`Class Score`[1]))
    famVal <- NULL
    if (is.null(fsco)|is.na(fsco)) {
        warning("fsco value (family score) is NULL. Assigning value of 0.00")
        fsco <- 0.000
    }
    if (fsco >= 0.900) {famVal <- "Positive"}
    if (fsco < 0.900 & fsco > 0.300) {famVal <- "Indeterminate"}
    if (fsco <= 0.300) {famVal <- "Negative"}
    if (is.null(famVal)) {
        warning("family value (famVal) is NULL")
        famVal <- "Indeterminate"
    }
    return(famVal)
}


GetOutFamily <- function(is450k, Mset_ba, Mset) {
    probs_mcf <- GetFamilyProb(is450k, Mset_ba, Mset)
    oo_mcf <- order(probs_mcf, decreasing = T)
    eps <- 1e-3
    out_class_family <- probs_mcf[oo_mcf[1:5]]
    out_class_family <- cbind(
        round(pmax(pmin(out_class_family, 1 - eps), eps),3),
        colnames(probs_mcf)[oo_mcf][1:5]
    )
    colnames(out_class_family) <- c("Class Score","Methylation Family")
    out_class_family <- as.data.frame(out_class_family)
    famVal <- GetFamScore(out_class_family)
    out_class_family$Interpretation = c(famVal,"","","","")
    return(out_class_family)
}

detach_package <- function(pkg, character.only = FALSE) {
    if (!character.only) {
        pkg <- deparse(substitute(pkg))
    }
    search_item <- paste("package", pkg, sep = ":")
    while (search_item %in% search()){
        base::detach(search_item, unload = TRUE, character.only = TRUE)
    }
}


Get450kProb <- function(betas){
    betas <- betas[match(rownames(mnp.v11b4::rf.pred$importance), rownames(betas)), drop = F]
    pred <- predict(mnp.v11b4::rf.pred, t(betas), type = "prob")
    sidat <- rownames(pred)
    classes <- colnames(pred)
    currData <- predict(mnp.v11b4::calfit$glmnet.fit, newx = pred, type = "response", s = 0.000683)
    pred <- matrix(currData, nrow = nrow(pred))
    class(pred) <- "numeric"
    storage.mode(pred) <- "numeric"
    colnames(pred) <- classes
    rownames(pred) <- sidat
    return(pred)
}


GetProbData <- function(is450k, Mset_ba, Mset) {
    LoadMnpManifest(is450k)
    if(!exists('calfit')){
        load(file.path("/Volumes/CBioinformatics/Methylation/Methylation_classifier_v11b6/",
                       "mnp.v116/mnp.v11b6/data/rfpred.v11b6.RData"))
    }
    if(is450k==T) {
        betas <- minfi::getBeta(Mset_ba[, 1])
        probs <- Get450kProb(betas)
        return(probs)
    }else{
        probs <- GetMnpV11Prob(Mset_ba)
        return(probs)
    }
}


GetOutScore <- function(out) {
    out_score <- as.numeric(paste0(out$`Subgroup Score`[1]))
    subVal_int <- NULL
    if (is.null(out_score)|is.na(out_score)) {
        warning("out_score value (sub score) is NULL. Assigning value of 0.00")
        out_score <- 0.000
    }
    if (out_score >= 0.500) {subVal_int <- "Positive"}
    if (out_score < 0.500 & out_score >= 0.100) {subVal_int <- "Indeterminate"}
    if (out_score < 0.300) {subVal_int <- "Negative"}
    if (is.null(subVal_int)) {
        warning("subscore value (subVal_int) is NULL")
        subVal_int <- "Indeterminate"
    }
    return(subVal_int)
}


GetOutClass <- function(msetDat) {
    library(verbose=F, warn.conflicts = F, quietly = T, package= "knitr")
    library(verbose=F, warn.conflicts = F, quietly = T, package= "kableExtra")  
    Mset_ba <- msetDat$Mset_ba
    Mset <- msetDat$Mset
    is450k <- Mset_ba@annotation[["array"]] != "IlluminaHumanMethylationEPIC"
    if(is450k==T){
        try(unloadNamespace("mnp.v11b6"),T)
        library(verbose=F, warn.conflicts = F, quietly = T, package= "mnp.v11b4", mask.ok = T)
        requireNamespace("mnp.v11b4", quietly = T)
    }else{
        try(unloadNamespace("mnp.v11b4"),T)
        library(verbose=F, warn.conflicts = F, quietly = T, package= "mnp.v11b6", mask.ok = T)
        requireNamespace("mnp.v11b6", quietly = T)
    }
    probs <- gb$GetProbData(is450k, Mset_ba, Mset)
    oo <- base::order(probs, decreasing = T)
    eps <- 1e-3
    out <- probs[oo[1:5]]
    out <- cbind(round(pmax(pmin(out,1 - eps),eps),3),colnames(probs)[oo][1:5])
    colnames(out) <- c("Subgroup Score","Methylation Subgroup")
    
    if(is450k==T){
        colToMatch <- colnames(probs)[oo][1]
        idx <- match(colToMatch, mnp.v11b4::reflist[, 4])
        if(is.na(idx)){
            colToMatch <- stringr::str_replace_all(colToMatch, "_", "")
            newList <- stringr::str_replace_all( mnp.v11b4::reflist[, 4], "_", "")
            idx <- match(colToMatch, newList)
        }
    }else{
        idx <- match(colnames(probs)[oo][1], mnp.v11b6::reflist[,2])
    }
    stopifnot(!is.na(idx))
    out <- as.data.frame(out)
    subVal_int <- GetOutScore(out)
    out$Interpretation = c(subVal_int,"","","","")
    out_class_family <- GetOutFamily(is450k, Mset_ba, Mset)
    return(list("out"=out,"idx"=idx, "out_class_family"=out_class_family))
}

# GetV12score <- function(RGset, FFPE = NULL) {
#     library(verbose=F, warn.conflicts = F, quietly = T, package= "mnp.v12b6")
#     library(verbose=F, warn.conflicts = F, quietly = T, package= "IlluminaHumanMethylationEPICmanifest")
#     #    load("/Volumes/CBioinformatics/Methylation/mnp_v12_R-package/mnp.v12b6/R/sysdata.rda")
#     RGset <- RGset[,1]
#     Mset12 <- mnp.v12b6::MNPpreprocessIllumina(RGset)
#     if(is.null(FFPE)) FFPE <- mnp.v12b6::MNPgetFFPE(RGset)
#     Mset12_ba <- mnp.v12b6::MNPbatchadjust(Mset12,FFPE)
#     sex12 <- ifelse(mnp.v12b6::MNPgetSex(Mset12)$predictedSex=="M","Male","Female")
#     super <- mnp.v12b6::MNPpredict(Mset12_ba,MCF_level="superfamily")[,1:2]
#     fam  <- mnp.v12b6::MNPpredict(Mset12_ba,MCF_level="family")[,1:2]
#     class <- mnp.v12b6::MNPpredict(Mset12_ba,MCF_level="class")[,1:2]
#     sclass <- mnp.v12b6::MNPpredict(Mset12_ba,abbreviation=FALSE)[,1:2]
#     out <- rbind(super,fam,class,sclass)
#     rownames(out) <- c("Super Family","Family","Class","Subclass")
#     out$maxscore <- sapply(out$maxscore,function(x)round(pmax(pmin(x,1-1e-4),1e-4),4))
#     return(out)
# }

PrintScoreTable <- function(outV12, dat) {
    btso = c("bordered")
    kgb <- c("striped", font_size = 14, bootstrap_options = btso, position = "left")
    txtc = "text-align:center;"
    be = c(booktabs = T, escape = F, linesep = "")
    outTable12 <- outV12 %>% mutate_all(as.character) %>% knitr::kable("html", be, align = 'clc') %>%
        kableExtra::kable_styling(kgb, full_width = F, position = "left") %>%
        kableExtra::column_spec(column = 1, background = "palegreen", bold = T, extra_css = txtc)
    v12Data <- as.data.frame(t(outV12))
    #message(paste0(capture.output(v12Data), collapse="\n"))
    v12redcap <- as.data.frame(matrix(nrow=1))
    v12redcap$record_id <- paste0(dat$sampleID)
    v12redcap$Super <-  gsub(pattern=",", "", v12Data$`Super Family`[1])
    v12redcap$Super_score <- v12Data$`Super Family`[2]
    v12redcap$Class_predict <- gsub(pattern=",", "", v12Data$Class[1])
    v12redcap$Class_score <- v12Data$Class[2]
    v12redcap$Subgroup_predict <- gsub(pattern=",", "", v12Data$Subclass[1])
    v12redcap$Subgroup_score <- v12Data$Subclass[2]
    v12redcap <- v12redcap[1,-1]
    rownames(v12redcap) <- NULL
    v12FiOut <- paste0(dat$run_id, "_v12.csv")
    v12DirOut <- file.path(fs::path_home(), "Desktop", dat$run_id)
    redcsv <- file.path(v12DirOut, v12FiOut)
    if(!dir.exists(v12DirOut)){dir.create(v12DirOut)}
    if (file.exists(redcsv)) {
        dfRedcap = as.data.frame(read.csv(redcsv, header = T, row.names = NULL))
        redDF <- rbind(dfRedcap, v12redcap)
        row.names(redDF) = NULL
    } else{
        redDF <- v12redcap
    }
    write.csv(x = redDF,file = redcsv, row.names = F, quote=F)
    return(outTable12)
}


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


GetCNVTables <- function(dra){
    gainDf <- t(data.frame(Gains = c(rownames(dra[(dra$Gain == T),]))))
    lossDf <- t(data.frame(Loss = c(rownames(dra[(dra$Loss == T),]))))
    btso = c("bordered")
    kgb <- c("striped",font_size = 14, bootstrap_options = btso, position = "left")
    kgh <- c(booktabs = T, escape = F, linesep = "")
    extra_css1 = "border-radius:0px;border-width:1px;border-style:solid;border-color:rgb(192,192,192);"
    txtc = "text-align:center;"
    gainTab <-
        gainDf %>% knitr::kable("html", kgh, align = 'clc') %>%
        kableExtra::kable_styling(kgb, full_width = F, position="left") %>%
        kableExtra::row_spec(row = 1, extra_css = extra_css1) %>%
        kableExtra::column_spec(column = 1, background = "palegreen", extra_css = txtc)
    lossTab <-
        lossDf %>% knitr::kable("html", kgh, align = 'clc') %>%
        kableExtra::kable_styling(kgb, full_width = F, position="left") %>%
        kableExtra::row_spec(row = 1, extra_css = extra_css1) %>%
        kableExtra::column_spec(column = 1, background = "#CD5C5C", extra_css = txtc)
    return(list("gainDf"=gainDf, "lossDf"=lossDf, "gainTab"=gainTab, "lossTab"=lossTab))
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

# GetMgmtPlot <- function(Mset_raw){
#     be = c(booktabs = T, escape = F, linesep = "")
#     btso = c("bordered")
#     txtc = "text-align:center;"
#     kgb <- c("striped",font_size = 14, bootstrap_options = btso, position = "left")
#     plotmgmt <- as.data.frame(mnp.v12b6::MNPpredict_mgmt(Mset_raw), row.names = NULL)
#     newVals <- sapply(plotmgmt[1,2:5],function(x)plyr::round_any(x,.0001))
#     plotmgmt[1,2:5] <- newVals
#     mgmtPlot <- plotmgmt %>% mutate_all(as.character) %>%
#         knitr::kable("html", be, align = 'clc') %>%
#         kableExtra::kable_styling(kgb, full_width = F, position = "left")
#     return(list("mgmtVal" = plotmgmt, "mgmtPlot" = mgmtPlot))
# }

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


PrintClassTable <- function(outList,gitPath){
    classTables <- gb$GetClassProbTables(outList)
    if(length(classTables) == 2){
        knitr::asis_output(classTables$famTable)
        knitr::asis_output(classTables$grpTable)
    }else{
        knitr::asis_output(classTables)
    }
    htmltools::includeHTML(file.path(gitPath,"InterpretKey.html"))
    knitr::asis_output("<h4>Methylation Class Description</h4>")
}


PrintGainLoss <- function(gnLss){
    knitr::asis_output('<p class="cnvdesc"></p> <hr class="dotted"><h4>CNV Gains and Loss</h4>')
    if (nrow(gnLss$gainDf) > 0) {knitr::asis_output(gnLss$gainTab)}
    if (nrow(gnLss$lossDf) > 0) {knitr::asis_output(gnLss$lossTab)}
    knitr::asis_output('<hr class="solid">')
}


GetRefList <- function(Mset_ba, outList) {
    if (Mset_ba@annotation[["array"]] != "IlluminaHumanMethylationEPIC") {
        return(I(as.character(mnp.v11b4::reflist[outList$idx,7])))
    } else{
        return(I(as.character(mnp.v11b6::reflist[outList$idx,7])))
    }
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
