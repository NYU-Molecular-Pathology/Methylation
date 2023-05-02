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

pkgLis <- c("utils", "grDevices", "stringr", "BiocManager")
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
    shortName <- lapply(stringr::str_split(island_probes$UCSC_RefGene_Name, ";", simplify = F), unique)
    island_probes$UCSC_RefGene_Name <- unlist(lapply(shortName, paste, collapse = ";"))
    message(paste(geneName, "probes found:"))
    message(paste0(capture.output(island_probes[,c("chr", "pos", "UCSC_RefGene_Name", "Relation_to_Island" )]), collapse="\n"))
    beta_subset <- beta_values[rownames(island_probes),]
    plot_data <- as.data.frame(t(beta_subset))
    return(boxplot(
        plot_data, names = rownames(island_probes), col = palette.colors(n=ncol(nrow(plot_data))),
        main = paste("CG Probe Performance on Illumina EPIC Methylation Arrays for", geneName, "Islands"),
        xlab = "Probe Names", ylab = "Beta Values", las = 2, cex.axis = 0.8
    ))
}

rgset <- minfiDataEPIC::RGsetEPIC
myGenes <- c("BRCA1", "BRCA2", "ATM", "RAD51C")

for (gn in myGenes) {
    PlotIslandProbes(rgset, gn)
}

