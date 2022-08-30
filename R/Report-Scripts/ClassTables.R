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

GetOutFamily <- function(is450k, Mset_ba, Mset){
    if (is450k==T) {
        library(verbose=F, warn.conflicts = F, quietly = T, package= "IlluminaHumanMethylation450kmanifest")
        library(verbose=F, warn.conflicts = F, quietly = T, package= "mnp.v11b4")
        tryCatch(
            expr = {
                probs_mcf <- mnp.v11b4::MNPpredict(Mset_ba[, 1], type='prob',MCF=TRUE)
            },
            error = function(e) {
                message("Error caught at mnp.v11b4::MNPpredict(), trying mnp.v11b6::MNPpredict()...")
                message(e)
                library(verbose=F,warn.conflicts=F, quietly=T, package="mnp.v11b6")
                probs_mcf <- mnp.v11b6::MNPpredict(Mset[, 1], type = 'prob', MCF = TRUE)
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
    oo_mcf <- order(probs_mcf, decreasing = T)
    eps <- 1e-3
    out_class_family <- probs_mcf[oo_mcf[1:5]]
    out_class_family <- cbind(
        round(pmax(pmin(out_class_family, 1 - eps), eps),3), colnames(probs_mcf)[oo_mcf][1:5])
    colnames(out_class_family) <- c("Class Score","Methylation Family")
    out_class_family <- as.data.frame(out_class_family)
    fsco <- as.numeric(paste0(out_class_family$`Class Score`[1]))
    famVal <- NULL
    if (is.null(fsco)|is.na(fsco)) {warning("fsco value (family score) is NULL. Assigning value of 0.00")
        fsco <- 0.000
    }
    if (fsco >= 0.900) {famVal <- "Positive"}
    if (fsco < 0.900 & fsco > 0.300) {famVal <- "Indeterminate"}
    if (fsco <= 0.300) {famVal <- "Negative"}
    
    if (is.null(famVal)) {warning("family value (famVal) is NULL")
        famVal <- "Indeterminate"
    }
    out_class_family$Interpretation = c(famVal,"","","","")
    return(out_class_family)
}

GetOutClass <- function(is450k, Mset_ba, Mset){
    if (is450k==T) {
        library(verbose=F, warn.conflicts = F, quietly = T, package= "IlluminaHumanMethylation450kmanifest")
        library(verbose=F, warn.conflicts = F, quietly = T, package= "mnp.v11b4")
        tryCatch(
            expr = {
                probs <- mnp.v11b4::MNPpredict(Mset_ba[, 1], type = 'prob')
            },
            error = function(e) {
                message("Error caught at mnp.v11b4::MNPpredict(Mset_ba[, 1], type = 'prob'):\n")
                message(e)
                message("\nNow trying mnp.v11b6::MNPpredict(Mset[, 1], type = 'prob')...")
                library(verbose=F,warn.conflicts = F, quietly = T, package="mnp.v11b6")
                probs <- mnp.v11b6::MNPpredict(Mset[, 1], type = 'prob')
            }
        )
    } else {
        library(verbose=F, warn.conflicts = F, quietly = T, package= "IlluminaHumanMethylationEPICmanifest")
        library(verbose=F, warn.conflicts = F, quietly = T, package= "mnp.v11b6")
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
    
    oo <- order(probs, decreasing = T)
    eps <- 1e-3
    out <- probs[oo[1:5]]
    out <- cbind(round(pmax(pmin(out,1 - eps),eps),3),colnames(probs)[oo][1:5])
    colnames(out) <- c("Subgroup Score","Methylation Subgroup")
    idx <- match(colnames(probs)[oo][1],mnp.v11b6::reflist[,2])
    out <- as.data.frame(out)
    out_score <- as.numeric(paste0(out$`Subgroup Score`[1]))
    subVal_int <- NULL
    if (is.null(out_score)|is.na(out_score)) {
        warning("out_score value (sub score) is NULL.  Assigning value of 0.00")
        out_score <- 0.000
    }
    if (out_score >= 0.500) {subVal_int <- "Positive"}
    if (out_score < 0.500 & out_score >= 0.100) {subVal_int <- "Indeterminate"}
    if (out_score < 0.300) {subVal_int <- "Negative"}
    if (is.null(subVal_int)) {
        warning("subscore value (subVal_int) is NULL")
        subVal_int <- "Indeterminate"
    }
    out$Interpretation = c(subVal_int,"","","","")
    return(list("out"=out,"idx"=idx))
}

GetV12score <- function(RGset, FFPE=NULL){
    library(verbose=F, warn.conflicts = F, quietly = T, package= "mnp.v12b6")
    library(verbose=F, warn.conflicts = F, quietly = T, package= "IlluminaHumanMethylationEPICmanifest")
    load("/Volumes/CBioinformatics/Methylation/mnp_v12_R-package/mnp.v12b6/R/sysdata.rda")
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
