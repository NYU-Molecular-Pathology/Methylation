## title: "sesameEthnicity"
## author: "Jonathan Serrano"
## date: "9/15/2021"
# https://bioconductor.org/packages/release/bioc/vignettes/sesame/inst/doc/sesame.html
# https://doi.org/10.1093/nar/gky691

rqPk <- function(pk) {
    if (!require(pk, character.only = T)) {
        install.packages(
            pk,
            warn.conflicts = F,
            ask = F,
            dependencies = T,
            Ncpus = 5
        )
    } else{
        library(
            pk,
            warn.conflicts = F,
            mask.ok = T,
            character.only = T,
            quietly = T
        )
    }
}


SwapSesNames <- function(sesameSet, targets) {
    idatName <- if(is.null(names(sesameSet))) {colnames(sesameSet)} else{names(sesameSet)}
    matchIndex <- match(idatName, targets$SentrixID_Pos)
    newNamesId <- targets[matchIndex, "Sample_Name"]
    if (is.null(names(sesameSet))) {
        colnames(sesameSet) <- newNamesId
    } else{
        names(sesameSet) <- newNamesId
    }
    return(sesameSet)
}

# Install Packages ------------------------------------------------------------------------
rqPk("devtools"); rqPk("pacman")
if(!require("sesameData")){devtools::install_github("zwdzwd/sesameData"); sesameData::sesameDataCache()}
if(!require("sesame")){devtools::install_github("zwdzwd/sesame")}
pacman::p_load("parallel","BiocParallel","sesame","sesameData","fs")
#mft <- sesame::sesameAnno_get("EPIC/EPIC.hg19.manifest.tsv.gz")
#mft <- sesame::sesameAnno_get("EPIC/homo_sapiens.tsv.gz")
#mft <- sesame::sesameAnno_get("Anno/EPIC/EPIC.hg19.manifest.tsv.gz")
#annoS <- sesame::sesameAnno_get("Anno/EPIC/EPIC.hg19.typeI_overlap_b151.rds")
#EPIC.hg19.manifest "EH3670"

GenerateSesameTab <- function(sampleSheet, outFiName="sesameOutput.csv"){
    mcc = parallel::detectCores()-2
    targets <- read.csv(sampleSheet, row.names=NULL)
    ssets = mclapply(sesame::searchIDATprefixes(getwd()), sesame::readIDATpair, mc.cores = mcc)
    ssets <- SwapSesNames(sesameSet = ssets, targets)
    betas = sesame::openSesame(getwd(), BPPARAM=BiocParallel::MulticoreParam(mcc))
    betas <- SwapSesNames(sesameSet = betas, targets)
    sesame_sex <- mclapply(ssets, sesame::inferSex, platform="EPIC", mc.cores = mcc)
    karyotype_sex <- mclapply(ssets, sesame::inferSexKaryotypes, mc.cores = mcc)
    predicted_ethnicity <- mclapply(ssets, sesame::inferEthnicity, mc.cores = mcc)
    ages <- mclapply(betas, sesame::predictAgeHorvath353, mc.cores = mcc)
    # tissues <- mclapply(betas, inferTissue,platform="EPIC", mc.cores = mcc)
    sesameOutput <- data.frame(
        karyotype_sex = cbind(karyotype_sex),
        predicted_ethnicity = cbind(predicted_ethnicity),
        sesame_sex = cbind(sesame_sex)
        #sesame_ages = cbind(ages)
    )
    sesameOutput$sesameID = rownames(sesameOutput)
    rownames(targets) <- targets$Sample_Name
    reOrdered <- targets[rownames(sesameOutput)]
    sesameOut <- cbind(sesameOutput, reOrdered)
    sesameOut <- as.data.frame(sesameOut,fix.empty.names=T)
    sesameOut <- apply(sesameOut,2,as.character)
    write.table(sesameOut, file = file.path(fs::path_home(),"Desktop",outFiName), quote=F, sep=",", row.names=F)
    return(sesameOut)
    #knitr::kable(sesameOut)
}
