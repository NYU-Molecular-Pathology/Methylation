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
    return(suppressWarnings(!requireNamespace(pkg, character.only = T, warn.conflicts = F)))
}
BioCinst <- function(pkg) {
    if (CheckReq(pkg)) {do.call(BiocManager::install, list(pkg, update = F, ask = F))}
    return(LibLoad(pkg))
}
CheckCran <- function(pkg) {
    if(CheckReq(pkg)){do.call(install.packages, list(pkg, dependencies=T))}
    return(LibLoad(pkg))
}

pkgLis <- c("utils", "grDevices", "stringr", "BiocManager", "ggplot2", "ggpubr")
bioPkg <- c("minfi", "minfiDataEPIC", "IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
pkgLoad <- unlist(lapply(pkgLis, CheckCran))
bioLoad <- unlist(lapply(bioPkg, BioCinst))
stopifnot(all(pkgLoad) & all(bioLoad))

# Function to plot cg probe performance for gene islands
PlotIslandProbes <- function(rgset, geneName="BRCA1") {
    mset <- preprocessRaw(rgset)
    beta_values <- getBeta(mset)
    anno <- getAnnotation(rgset)
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
    beta_subset <- beta_values[rownames(island_probes),]
    plot_data <- as.data.frame(t(beta_subset))
    plotCols <- if(ncol(plot_data) > 2) {pals::brewer.dark2(ncol(plot_data))} else{c("darkred", "royalblue")}
    txtTheme <- ggplot2::theme_bw() + ggplot2::theme(axis.text.x = element_text(angle = 90, face="bold", color=plotCols))
    plot_data <- stack(plot_data)
    means <- aggregate(values ~ ind, plot_data, mean)
    medians <- aggregate(values ~ ind, plot_data, median)
    ggTable <- as.data.frame(t(island_probes))
    mytheme <- gridExtra::ttheme_default(base_size = 10, colhead=list(bg_params=list(fill=plotCols, alpha=0.8)),
                                         core = list(padding = unit(c(5, 1.5), "mm")))
    tbl <- tableGrob(ggTable, theme = mytheme)
    thePlot <-
        ggplot2::ggplot(plot_data, aes(x = ind, y = values, fill = ind)) +
        geom_boxplot(alpha=0.8) + scale_fill_manual(values=plotCols) +
        ggplot2::geom_dotplot(binaxis='y', stackdir = 'center', dotsize = 0.5, binwidth=0.003) +
        geom_label(data= means, aes(label = paste("Mean:",round(values, 3)), y = -0.01),
                   size = 3, alpha=0.5, color="darkred", fill="white", fontface = "bold") +
        geom_label(data = medians, aes(label = paste("Median:", round(values, 3)), y = 0),
                   size = 3,alpha=0.8, color="darkblue", fill="white", fontface = "bold") +
        txtTheme + theme(legend.position="none") + labs(x="Probe Names", y= "Beta Values") +
        ggplot2::ggtitle(paste("CG Probe Performance on Illumina EPIC Methylation Arrays for", geneName, "Islands"),)

    return(
        grid.arrange(thePlot, tbl, nrow = 2, as.table = T, heights=c(5,1))
        )
}

rgset <- minfiDataEPIC::RGsetEPIC
myGenes <- c("BRCA1", "BRCA2", "ATM", "RAD51C")

for (gn in myGenes) {
    print(PlotIslandProbes(rgset, gn))
}

