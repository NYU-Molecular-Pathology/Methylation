gb <- globalenv(); assign("gb", gb)
require("compiler")
invisible(supM(compiler::enableJIT(3)))
invisible(supM(compiler::compilePKGS(enable = TRUE)))
supM(compiler::setCompilerOptions(suppressAll = TRUE, optimize = 3))

GetSexMsetBa <- function(is450k, RGset, FFPE=NULL){
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
    if (is450k==T) {
        library(verbose=F, warn.conflicts = F, quietly = T, package= "IlluminaHumanMethylation450kmanifest")
        library(verbose=F, warn.conflicts = F, quietly = T, package= "mnp.v11b4")

        tryCatch(
            expr = {
                load("/Volumes/CBioinformatics/Methylation/Methylation_classifier_v11b4/mnp.v11b4/data/rfpred.v11b4.RData")
                probs_mcf <- mnp.v11b4::MNPpredict(Mset_ba[, 1], type='prob',MCF=TRUE)
            },
            error = function(e) {
                message("Error caught at mnp.v11b4::MNPpredict(), trying uncalibrated...")
                probs_mcf <- mnp.v11b4::MNPpredict(Mset_ba[, 1], type='prob',MCF=TRUE,calibrate = FALSE)
            }
        )

    } else {
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
    }
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

GetProbData <- function(is450k, Mset_ba, Mset){
    if (is450k==T) {
        library(verbose=F, warn.conflicts = F, quietly = T, package= "IlluminaHumanMethylation450kmanifest")
        library(verbose=F, warn.conflicts = F, quietly = T, package= "mnp.v11b4")
        load("/Volumes/CBioinformatics/Methylation/Methylation_classifier_v11b4/mnp.v11b4/data/rfpred.v11b4.RData")
        tryCatch(
            expr = {
                probs <- mnp.v11b4::MNPpredict(Mset_ba[, 1], type = 'prob')
            },
            error = function(e) {
                message("Error caught at mnp.v11b4::MNPpredict(Mset_ba[, 1], type = 'prob'):\n")
                message(e)
                message("\nNow trying mnp.v11b6::MNPpredict(Mset[, 1], type = 'prob')")
                probs <- mnp.v11b4::MNPpredict(Mset[, 1], type = 'prob')
            }
        )
    } else {
        library(verbose=F, warn.conflicts = F, quietly = T, package= "IlluminaHumanMethylationEPICmanifest")
        library(verbose=F, warn.conflicts = F, quietly = T, package= "mnp.v11b6")
        load("/Volumes/CBioinformatics/Methylation/Methylation_classifier_v11b6/mnp.v116/mnp.v11b6/data/rfpred.v11b6.RData")
        tryCatch(
            expr = {
                probs <- mnp.v11b6::MNPpredict(Mset_ba[, 1], type = 'prob')
            },
            error = function(e) {
                message("Error occured at Brain Classifier v11 prediction: using Mset instead of Mset_ba:\n")
                message(e)
                probs <- mnp.v11b6::MNPpredict(Mset[, 1], type = 'prob')
            }
        )
    }
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

GetOutClass <- function(is450k, Mset_ba, Mset){
    probs <- GetProbData(is450k, Mset_ba, Mset)
    oo <- order(probs, decreasing = T)
    eps <- 1e-3
    out <- probs[oo[1:5]]
    out <- cbind(round(pmax(pmin(out,1 - eps),eps),3),colnames(probs)[oo][1:5])
    colnames(out) <- c("Subgroup Score","Methylation Subgroup")
    if(is450k){
        idx <- match(colnames(probs)[oo][1], mnp.v11b4::reflist[,2])
    }else{
        idx <- match(colnames(probs)[oo][1], mnp.v11b6::reflist[,2])
    }
    stopifnot(!is.na(idx))
    out <- as.data.frame(out)
    message(out)
    subVal_int <- GetOutScore(out)
    out$Interpretation = c(subVal_int,"","","","")
    return(list("out"=out,"idx"=idx))
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

GetClassProbTables <- function(out_class_family, out){
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
