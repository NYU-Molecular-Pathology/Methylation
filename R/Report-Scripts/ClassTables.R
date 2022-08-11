GetOutFamily <- function(is450k, Mset_ba){
    if (is450k==T) {
        library(verbose=F,warn.conflicts = F, quietly = T, package="mnp.v11b4")

        probs_mcf <- mnp.v11b4::MNPpredict(Mset_ba[, 1], type='prob',MCF=TRUE)
    } else {
        library(verbose=F,warn.conflicts=F, quietly=T, package="mnp.v11b6")
        tryCatch(
            expr = {

                probs_mcf <- mnp.v11b6::MNPpredict(Mset_ba[, 1], type = 'prob', MCF = TRUE)
            },
            error = function(e) {
                message("Error occured at Brain Classifier v11 prediction:")
                message("Value of Mset_ba is: ", paste0(Mset_ba[, 1]))
                message("Using MNPpredict(Mset[, 1]) instead of Mset_ba:\n")
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
    out_class_family$Interpretation = c(famVal,"","","","")
    return(out_class_family)
}

GetOutClass <- function(is450k, Mset_ba){
    if (is450k==T) {
        library(verbose=F,warn.conflicts = F, quietly = T, package="mnp.v11b4")
        probs <- mnp.v11b4::MNPpredict(Mset_ba[, 1], type = 'prob')
    } else {
        library(verbose=F,warn.conflicts=F, quietly=T, package="mnp.v11b6")
        tryCatch(
            expr = {
                probs <- mnp.v11b6::MNPpredict(Mset_ba[, 1], type = 'prob')
            },
            error = function(e) {
                message("Error occured at Brain Classifier v11 prediction:")
                message("Value of Mset_ba is: ", paste0(Mset_ba[, 1]))
                message("Using MNPpredict(Mset[, 1]) instead of Mset_ba:\n")
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
    return(out)
}