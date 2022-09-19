gb <- globalenv(); assign("gb", gb)
SilentLoadLib <- function(pkg){
    stopifnot(suppressPackageStartupMessages(library(
        verbose = F,
        warn.conflicts = F,
        quietly = T,
        character.only = T,
        package = pkg,
        logical.return = T
    )))
}
mealPkgs <- c("MEAL", "MultiDataSet",  "minfi", "ggplot2")
invisible(lapply(mealPkgs, SilentLoadLib))

GetMethSet <- function(RGSet){
    Mset <- minfi::preprocessIllumina(RGSet)
    meth <- mapToGenome(ratioConvert(Mset))
    rowData(meth) <- getAnnotation(meth)[, -c(1:3)]
    meth <- dropMethylationLoci(meth) ## Remove probes measuring SNPs
    meth <- dropLociWithSnps(meth) ## Remove probes with SNPs
    meth <- meth[!apply(getBeta(meth), 1, function(x) any(is.na(x))), ] ## Remove probes with NAs
    set.seed(0)  ## Select a subset of samples
    meth <- meth[sample(nrow(meth), 100000), ]
}

GetMealGene <- function(X, resAdj, listOfGenes){
    geneName <- paste0(listOfGenes[X])
    return(MEAL::getGeneVals(resAdj, geneName, genecol = "UCSC_RefGene_Name", fNames = c("chromosome", "start")))
}

GetNewGeneVals <- function(geneDataLi, resAdj) {
    listOfGenes <- geneDataLi[,"Gene"]
    X = 1:length(listOfGenes)
    
    geneValDf <- lapply(
        X,
        FUN = function(X, resAdj, listOfGenes) {
            geneName <- paste0(listOfGenes[X])
            theRow <-
                MEAL::getGeneVals(
                    resAdj,
                    geneName,
                    genecol = "UCSC_RefGene_Name",
                    fNames = c("chromosome", "start")
                )
            return(as.data.frame(theRow))
        },
        resAdj,
        listOfGenes
    )
    
    newValues <- data.frame()
    for (x in geneValDf) {
        newValues <- rbind(newValues, x)
    }
    
    for (gn in 1:length(newValues$UCSC_RefGene_Name)) {
        newName <-
            unlist(stringr::str_split(newValues$UCSC_RefGene_Name[gn], pattern = ";"))
        newName <- paste(unique(newName), collapse = ";")
        newValues[gn, "UCSC_RefGene_Name"] <- newName
    }
    newValues$Region <- ""
    for (rgn in 1:length(geneDataLi$Gene)) {
        theName <- geneDataLi$Gene[rgn]
        for (xName in 1:length(newValues$UCSC_RefGene_Name)) {
            theGene <-
                unlist(stringr::str_split(newValues$UCSC_RefGene_Name[xName], ";"))
            if (any(grepl(pattern = theName, theGene)) == T) {
                newValues[xName, "Region"] <- geneDataLi$Site[rgn]
            }
        }
    }
    return(newValues)
}
