gb <- globalenv(); assign("gb", gb)
supM <- function(objTing){return(suppressMessages(suppressWarnings(objTing)))}
assign(x = "supM", value = supM, envir = .GlobalEnv)

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
    if(!require("librarian")){install.packages("librarian", dependencies=T, quiet=T)}
    if(!require("devtools")){install.packages("devtools", dependencies=T)}
    if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager", dependencies = T)}
    if(!requireNamespace("UniD", quietly = TRUE)) {
        install.packages(
            "/Volumes/CBioinformatics/Methylation/UniD", type = "source", dependencies = T, repo = NULL
            )
    }
    librarian::shelf(pkgs, ask=F, verbose=F, warn.conflicts = F, quietly = T)
    require("needs", quietly = T, warn.conflicts=F)
    library(verbose=F, warn.conflicts = F, quietly = T, package = "knitr")
    knitr::opts_chunk$set(optsLi)
    knitr::opts_knit$set(optsLi)
    knitr::opts_current$set(optsLi)
    knitr::opts_knit$set(root.dir = getwd())
    options(width = 300, scipen = 5)
    library(verbose=F, warn.conflicts = F, quietly = T, package = "pander")
    panderOptions('table.alignment.default', "left")
    require("compiler")
#    invisible(supM(compiler::enableJIT(3)))
    invisible(supM(compiler::compilePKGS(enable = TRUE)))
 #   supM(compiler::setCompilerOptions(suppressAll = TRUE, optimize = 3))
    
knitr::opts_chunk$set(chunkOpts)
supM(library(verbose=F, warn.conflicts = F, quietly = T, package="ggplot2"))
}


LoadReportPkgs(pkgs, optsLi, chunkOpts)

GetSexMsetBa <- function(RGset, FFPE=NULL){
        is450k <- RGset@annotation[[1]] == "IlluminaHumanMethylation450k"
    if (is450k) {
        library(verbose=F, warn.conflicts = F, quietly = T, package = "IlluminaHumanMethylation450kmanifest")
        library(verbose=F, warn.conflicts = F, quietly = T, package= "mnp.v11b4")
        Mset <- mnp.v11b4::MNPpreprocessIllumina(RGset)
        if (is.null(FFPE)) {FFPE <- mnp.v11b4::MNPgetFFPE(RGset)}
        Mset_ba <- mnp.v11b4::MNPbatchadjust(Mset, FFPE)
        if (FFPE == "Frozen") {Mset@preprocessMethod <-
            c(Mset_ba@preprocessMethod, FFPE_Frozen.mnp.adjustment = '0.11')}
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
    return(list("sex"=sex, "Mset"= Mset, "Mset_ba"=Mset_ba, "FFPE"=FFPE))
}

GetFamilyProb <- function(is450k, Mset_ba, Mset){
    library(verbose=F, warn.conflicts = F, quietly = T, package= "IlluminaHumanMethylationEPICmanifest")
    library(verbose=F, warn.conflicts = F, quietly = T, package= "mnp.v11b6")
    tryCatch(
        expr = {
            probs_mcf <- mnp.v11b6::MNPpredict(Mset_ba[, 1], type = 'prob', MCF = TRUE)
        },
        error = function(e) {
            message("Error occured at Brain Classifier v11 prediction:")
            message("Using MNPpredict(Mset[, 1]) instead of Mset_ba\n")
            probs_mcf <- mnp.v11b6::MNPpredict(Mset[, 1], type = 'prob', MCF = TRUE)
        }
    )
    
    return(probs_mcf)
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

GetOutFamily <- function(is450k, Mset_ba, Mset){
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

GetProbData <- function(is450k, Mset_ba, Mset) {
    library(verbose = F, warn.conflicts = F, quietly = T, package = "IlluminaHumanMethylationEPICmanifest")
    library(verbose = F, warn.conflicts = F, quietly = T, package = "mnp.v11b6")
    try(unloadNamespace("mnp.v11b4"), silent = T)
    try(detach_package("mnp.v11b4"), silent = T)
    invisible(loadNamespace("mnp.v11b6"))
    invisible(requireNamespace("mnp.v11b6"))
    if(!exists('calfit')){
        load(file.path(
            "/Volumes/CBioinformatics/Methylation/Methylation_classifier_v11b6/",
            "mnp.v116/mnp.v11b6/data/rfpred.v11b6.RData"))    
        }
    tryCatch(
        expr = {
            probs <- mnp.v11b6::MNPpredict(Mset_ba[, 1], type = 'prob')
        },
        error = function(e) {
            message("Error occured at Brain Classifier v11 prediction:")
            message("Using Mset instead of Mset_ba:\n", e)
            probs <- mnp.v11b6::MNPpredict(Mset[, 1], type = 'prob')
        }
    )
    return(probs)
}

GetOutScore <- function(out){
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

GetOutClass <- function(Mset_ba, Mset){
        is450k <- Mset_ba@annotation[["array"]] == "IlluminaHumanMethylationEPIC"
    library(verbose=F, warn.conflicts = F, quietly = T, package= "knitr")
    library(verbose=F, warn.conflicts = F, quietly = T, package= "kableExtra")
    library(verbose=F, warn.conflicts = F, quietly = T, package= "mnp.v11b6", mask.ok = T)
    requireNamespace("mnp.v11b6", quietly = T)
    probs <- gb$GetProbData(is450k, Mset_ba, Mset)
    oo <- base::order(probs, decreasing = T)
    eps <- 1e-3
    out <- probs[oo[1:5]]
    out <- cbind(round(pmax(pmin(out,1 - eps),eps),3),colnames(probs)[oo][1:5])
    colnames(out) <- c("Subgroup Score","Methylation Subgroup")
    idx <- match(colnames(probs)[oo][1], mnp.v11b6::reflist[,2])
    stopifnot(!is.na(idx))
    out <- as.data.frame(out)
    subVal_int <- GetOutScore(out)
    out$Interpretation = c(subVal_int,"","","","")
    out_class_family <- GetOutFamily(is450k, Mset_ba, Mset)
    return(list("out"=out,"idx"=idx, "out_class_family"=out_class_family))
}

GetV12score <- function(RGset, FFPE=NULL){
    library(verbose=F, warn.conflicts = F, quietly = T, package= "mnp.v12b6")
    library(verbose=F, warn.conflicts = F, quietly = T, package= "IlluminaHumanMethylationEPICmanifest")
#    load("/Volumes/CBioinformatics/Methylation/mnp_v12_R-package/mnp.v12b6/R/sysdata.rda")
    RGset <- RGset[,1]
    Mset12 <- mnp.v12b6::MNPpreprocessIllumina(RGset)
    if(is.null(FFPE)) FFPE <- mnp.v12b6::MNPgetFFPE(RGset)
    Mset12_ba <- mnp.v12b6::MNPbatchadjust(Mset12,FFPE)
    sex12 <- ifelse(mnp.v12b6::MNPgetSex(Mset12)$predictedSex=="M","Male","Female")
    super <- mnp.v12b6::MNPpredict(Mset12_ba,MCF_level="superfamily")[,1:2]
    fam  <- mnp.v12b6::MNPpredict(Mset12_ba,MCF_level="family")[,1:2]
    class <- mnp.v12b6::MNPpredict(Mset12_ba,MCF_level="class")[,1:2]
    sclass <- mnp.v12b6::MNPpredict(Mset12_ba,abbreviation=FALSE)[,1:2]
    out <- rbind(super,fam,class,sclass)
    rownames(out) <- c("Super Family","Family","Class","Subclass")
    out$maxscore <- sapply(out$maxscore,function(x)round(pmax(pmin(x,1-1e-4),1e-4),4))
    return(out)
}

PrintScoreTable <- function(outV12) {
    btso = c("bordered")
    kgb <- c("striped",font_size = 14, bootstrap_options = btso, position = "left")
    txtc = "text-align:center;"
    be = c(booktabs = T, escape = F, linesep = "")
    outTable12 <- outV12 %>% mutate_all(as.character) %>% knitr::kable("html", be, align = 'clc') %>%
        kableExtra::kable_styling(kgb, full_width = F, position = "left") %>%
        kableExtra::column_spec(column = 1, background = "palegreen", bold = T, extra_css = txtc)
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
    suppinfo <-
        c(
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
                           
GetMgmtPlot <- function(Mset_raw){
    be = c(booktabs = T, escape = F, linesep = "")
    btso = c("bordered")
    txtc = "text-align:center;"
    kgb <- c("striped",font_size = 14, bootstrap_options = btso, position = "left")
    plotmgmt <- as.data.frame(mnp.v12b6::MNPpredict_mgmt(Mset_raw), row.names = NULL)
    newVals <- sapply(plotmgmt[1,2:5],function(x)plyr::round_any(x,.0001))
    plotmgmt[1,2:5] <- newVals
    mgmtPlot <- plotmgmt %>% mutate_all(as.character) %>% 
        knitr::kable("html", be, align = 'clc') %>%
        kableExtra::kable_styling(kgb, full_width = F, position="left")
    return(list("mgmtVal" = plotmgmt, "mgmtPlot" = mgmtPlot))
}

GetClassProbTables <- function(out_class_family=NULL, outList=NULL){
    stopifnot(!is.null(outList) & !is.null(out_class_family))
    out <- outList$out
    xtraCss1="border-radius:0px;border-width:1px;border-style:solid;border-color:rgb(192,192,192);"
    xtraCss2="border-radius:0px;border-width:3px;border-style:solid;border-color:rgb(26,105,16);"
    xtraCss3="border-radius:0px;border-width:3px;border-style:solid;border-color:rgb(16,28,105);"
    txtc = "text-align:center;"
    btso = c("bordered")
    kgb <- c("striped",font_size = 14, bootstrap_options = btso, position = "left")
    kgh <- c(booktabs = T, escape = F, linesep = "")
    
    famTable <- out_class_family %>% 
    	knitr::kable("html",kgh,align='clc') %>%
    	kableExtra::kable_styling(kgb, full_width = F, position="float_left") %>%
    	kableExtra::column_spec(column=c(1,2),extra_css=xtraCss1) %>%
        kableExtra::column_spec(column = 2, background = "rgb(204,255,204)", extra_css = txtc) %>%
        kableExtra::row_spec(row = 0, font_size = 16, background = "rgb(127,217,126)", color = "black") %>%
    	kableExtra::row_spec(row=1,extra_css=xtraCss2)
    			
    grpTable <- out %>% 
    	knitr::kable("html",kgh,align='clc') %>%
        kableExtra::kable_styling(kgb, full_width = F) %>%
    	kableExtra::column_spec(column=c(1,2), extra_css=xtraCss1) %>%
    	kableExtra::column_spec(column=2,background="rgb(204,230,255)", extra_css=txtc) %>%
    	kableExtra::row_spec(row=0,font_size=16,background="rgb(135,174,237)", color="black") %>%
    	kableExtra::row_spec(row=1,extra_css=xtraCss3)
    
    return(list("famTable"=famTable,"grpTable"=grpTable))
}
                      
SuppInfoTable <- function(suppinfo){
    suppinfo <- as.data.frame(t(suppinfo))
    rownames(suppinfo) <- NULL
    kgb <- c("striped",font_size = 9, bootstrap_options = c("bordered"), position = "float_left")
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
