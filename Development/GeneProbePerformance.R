#!/usr/bin/env Rscript
## ---------------------------
## Script name: GenePerformProbes.R
## Purpose: Select probes for island regions of an input gene and test probe performance
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
library(minfi)
library(minfiDataEPIC)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Function to plot cg probe performance for gene islands
PlotIslandProbes <- function(rgset, geneName="BRCA1") {
    mset <- preprocessRaw(rgset)
    beta_values <- getBeta(mset)
    anno <- getAnnotation(rgset)
    gene_probes <- anno[grepl(anno$UCSC_RefGene_Name, pattern=geneName),]
    island_probes <- gene_probes[gene_probes$Relation_to_Island == "Island",]
    beta_subset <- beta_values[rownames(island_probes),]
    plot_data <- as.data.frame(cbind(island_probes, t(beta_subset)))
    scopeRange <- c((ncol(island_probes)+1):ncol(plot_data))
    thePlot <- boxplot(
        plot_data[, scopeRange], names = rownames(island_probes),
        main = paste("CG Probe Performance on Illumina EPIC Methylation Arrays for", geneName, "Islands"),
        xlab = "Probe Names", ylab = "Beta Values", las = 2, cex.axis = 0.8
        )
    return(thePlot)
}

rgset <- minfiDataEPIC::RGsetEPIC
PlotIslandProbes(rgset, "BRCA1")
