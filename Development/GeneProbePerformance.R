#!/usr/bin/env Rscript
## ---------------------------
## Script name: GenePerformProbes.R
## Purpose: Select probes for island regions of an input gene and test probe performance
## Author: Jonathan Serrano
## Copyright (c) NYULH 2023
## ---------------------------

# Load or Install required packages
LibLoad <- function(pkg){
    libOpt <- list(pkg, quietly = T, logical.return = T, warn.conflicts = F)
    suppressPackageStartupMessages(suppressWarnings(do.call(library, libOpt)))
}
CheckReq <- function(pkg){
    return(suppressWarnings(!require(pkg, character.only = T, warn.conflicts = F)))
}
BioCinst <- function(pkg) {
    if (CheckReq(pkg)) {do.call(BiocManager::install, list(pkg, update = F, ask = F))}
    return(LibLoad(pkg))
}
CheckCran <- function(pkg) {
    if(CheckReq(pkg)){do.call(install.packages, list(pkg, dependencies=T))}
    return(LibLoad(pkg))
}

pkgLis <- c("utils", "grDevices", "stringr", "BiocManager", "ggplot2", "pals", "gridExtra")
bioPkg <- c("minfi", "minfiDataEPIC", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
pkgLoad <- unlist(lapply(pkgLis, CheckCran))
bioLoad <- unlist(lapply(bioPkg, BioCinst))
stopifnot(all(pkgLoad) & all(bioLoad))

GetIslandProbes <- function(anno, geneName){
    gene_probes <- anno[grepl(anno$UCSC_RefGene_Name, pattern=geneName),]
    island_probes <- as.data.frame(gene_probes[gene_probes$Relation_to_Island == "Island",])
    shortName <- lapply(lapply(stringr::str_split(island_probes$UCSC_RefGene_Name, ";", simplify = F), unique), sort)
    toKeep <- unlist(lapply(shortName, FUN=function(i){any(i==geneName)}))
    island_probes <- island_probes[toKeep,]
    shortName <- shortName[toKeep]
    island_probes$UCSC_RefGene_Name <- unlist(lapply(shortName, paste, collapse = ";"))
    message(paste(geneName, "probes found:"))
    island_probes <- island_probes[,c("chr", "pos", "UCSC_RefGene_Name", "Relation_to_Island" )]
    message(paste0(capture.output(island_probes), collapse="\n"))
    return(island_probes)
}


GenerateProbePlot <- function(plot_data, plotCols){
    mainTitle <- "CG Probe Performance on Illumina EPIC Methylation Arrays for"
    txtTheme <- ggplot2::theme(axis.text.x = element_text(angle = 90, face="bold", color="black"),
                               legend.position="none")
    plot_data <- stack(plot_data)
    means <- aggregate(values ~ ind, plot_data, mean)
    medians <- aggregate(values ~ ind, plot_data, median)
    thePlot <-
        ggplot2::ggplot(plot_data, aes(x = ind, y = values, fill = ind)) +
        geom_boxplot(alpha=0.8, outlier.colour="red",  outlier.size=4) +
        scale_fill_manual(values=plotCols) +
        ggplot2::geom_dotplot(binaxis='y', stackdir = 'center', dotsize = 0.5, binwidth=0.004) +
        geom_label(data= means, aes(label = paste("Mean:",round(values, 3)), y = -0.02),
                   size = 3, alpha=0.5, color="darkred", fill="white", fontface = "bold") +
        geom_label(data = medians, aes(label = paste("Median:", round(values, 3)), y = -0.05),
                   size = 3,alpha=0.8, color="darkblue", fill="white", fontface = "bold") +
        ggplot2::theme_bw() + txtTheme + labs(x="Probe Names", y= "Beta Values") +
        ggplot2::ggtitle(paste(mainTitle, geneName, "Islands"),)
    return(thePlot)
}

# Function to plot cg probe performance for gene islands
PlotIslandProbes <- function(beta_values, anno, geneName="BRCA1") {
    imgName <- paste(geneName, "probes_CNV.png", sep = "_")
    if(!file.exists(file.path(getwd(), imgName))){
        island_probes <- GetIslandProbes(anno, geneName)
        beta_subset <- beta_values[rownames(island_probes),]
        plot_data <- as.data.frame(t(beta_subset))
        allCols <- paste(pals::polychrome(36))
        allCols <- allCols[-2]
        plotCols <- if(ncol(plot_data) > 2) {sample(allCols, ncol(plot_data))} else{c("darkred", "royalblue")}
        ggTable <- as.data.frame(t(island_probes))
        mytheme <-
            gridExtra::ttheme_default(
                base_size = 10,
                colhead = list(bg_params = list(fill = plotCols, alpha = 0.8)),
                core = list(padding = unit(c(1.5, 1.5), "mm"))
            )
        if(ncol(ggTable) > 10){
            ggTable1 <- ggTable[, 1:10]
            ggTable2 <- ggTable[, 10:ncol(ggTable)]
        }
        if(ncol(ggTable)<=10){
        tbl <- tableGrob(ggTable, theme = mytheme)
        }else{
            tbl1 <- tableGrob(ggTable1, theme = mytheme)
            mytheme2 <- gridExtra::ttheme_default(
                base_size = 10, colhead = list(bg_params = list(
                    fill = plotCols[10:ncol(ggTable)], alpha = 0.8)),
                    core = list(padding = unit(c(1.5, 1.5), "mm"))
                )
            tbl2 <- tableGrob(ggTable2, theme = mytheme2)
        }
        thePlot <- GenerateProbePlot(plot_data, plotCols)
        png(file.path(getwd(), imgName), width = 17, height = 10, units="in", res=300)
        if(ncol(ggTable)>10){
            gridExtra::grid.arrange(thePlot, tbl1, tbl2, nrow = 3, as.table = T, heights=c(5, 1, 1))
        }else{
            gridExtra::grid.arrange(thePlot, tbl, nrow = 2, as.table = T, heights=c(5, 1))
        }
        invisible(dev.off())

        }
}

#rgset <- minfiDataEPIC::RGsetEPIC
targets <- as.data.frame(read.csv("samplesheet.csv"))
rgset <- minfi::read.metharray.exp(getwd(), targets, verbose = T, force = T)
mset <- minfi::preprocessRaw(rgset)
beta_values <- minfi::getBeta(mset)
anno <- minfi::getAnnotation(rgset)

myGenes <- c("BRCA1", "BRCA2", "ATM", "RAD51C")

for (gn in myGenes) {PlotIslandProbes(beta_values, anno, gn)}
