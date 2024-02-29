#!/usr/bin/env Rscript
## ---------------------------
## Script name: Params_output.R
## Purpose: Source global parameters for research report output file naming
## Date Created: May 17, 2022
## Date Last Modified: February 29, 2024
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

gb <- globalenv(); assign("gb", gb)
# Default Output file names  -----------------------------------------------
td <- format(Sys.Date(),"%b%d") # today's date (i.e., "Jun17")
rgOut <- file.path(gb$runDir, "data", paste0(td, "_RGset.Rdata"))
rawBetaFi <- file.path(gb$runDir, "data", paste0(td, "_betaRaw.Rdata"))
pValsOutFi <- file.path(gb$runDir, "data", paste0(td, "_DetPvals.Rdata"))
mbfile <- file.path(gb$runDir, "data", paste0(td, "_mSetSq.beta.Rdata"))
unbetaVariance <- file.path(gb$runDir, "data", paste0(td, "_unbetaVariance.RData"))
combatOut <- file.path(gb$runDir, "data", paste0(td, "_combatBetas.Rdata"))

# Default Output directory names  -----------------------------------------------------
tsneOutU = "figures/tsne/unsupervised/"
tsneOutS = "figures/tsne/supervised/"
hmOutU = "figures/heatmaps/unsupervised/"
hmOutS = "figures/heatmaps/supervised/"
cnvSegDir = file.path(gb$runDir, "data", "CNV_segments")
supbetaOut <- file.path(gb$runDir, "data", paste0(td, "_supervisedBetas"))
annotFi <- file.path(gb$runDir, "data", paste0(td, "_annotations.rds"))

# FUNCTION: Checks if directory exists, if not then creates recursively
CheckDirCreate <- function(pathLocation){
    dataOutDir <- file.path(gb$runDir, pathLocation)
    if (!dir.exists(dataOutDir)) {dir.create(dataOutDir, recursive = T)}
}

# Create OutPut Directories  ----------------------------------------------------------
CheckDirCreate("data")
CheckDirCreate("idats")
CheckDirCreate("figures/mds")
CheckDirCreate("csv")
CheckDirCreate(tsneOutU)
if (gb$supervisedRun) {
    CheckDirCreate(tsneOutS)
}
if (gb$genHeatMaps == T) {
    CheckDirCreate(hmOutU)
    if (gb$supervisedRun) CheckDirCreate(hmOutS)
}
if (gb$genCNchunk == T) {
    CheckDirCreate("figures/cnv/")
    CheckDirCreate(cnvSegDir)
}
if (gb$genPairwise == T) {CheckDirCreate("/figures/diffmean/")}
if (gb$genPathChunk == T) {
    CheckDirCreate("figures/pathway/")
    CheckDirCreate("figures/cluster/")
    CheckDirCreate("figures/heatmaps/gene_cluster/")
}

# Segments Copy Number Output Files
segFile <- file.path(cnvSegDir, paste0(td, "_segfile.csv"))
seg_clust_file <- file.path(cnvSegDir, paste0(td, "_seg_clusters.txt"))
gsetFile <- file.path(gb$runDir, "data", paste0(td, "_gsetfunnorm.rds"))
gBetaFile <- file.path(gb$runDir, "data", paste0(td, "_gsetbeta.rds"))
ClusfiNam <- file.path(gb$runDir, "data", paste0(td, "_dmp.csv"))

# Pathway Analysis Output Files
if (gb$genPathChunk) {
    siteSpecific <- as.data.frame(readxl::read_excel(gb$genesInputFi))
    your_genes <- siteSpecific$Gene
    geneListIn <- file.path(gb$runDir, "data", paste0(td,"_Genes_Pathway.csv"))
    pathCsvOut <- file.path(gb$runDir, "data", paste0(td,"-signaling_pathway.csv"))
}

# Plot Title Names  -----------------------------------------------
clusType <- c("Unsupervised", "Supervised")
if (gb$supervisedRun == F) {
    clusType <- clusType[1]
}

varProbes <- c(100, 1000, 10000) # Which top Variance probes to pull i.e. c(100, 1000)
hmVarProbes <- c(100, 1000, 3000)

gb$tsne_titles <- gb$generateTitles(clusType, topTitle = as.character(varProbes), gb$titleMain)
topN = 1000
topVar = 10000

mdsTitle <- paste("Top", topN, "Common", "mSet Sq.beta", "MDS plot")

assignOpts <- function() {
    geneClusPath = "figures/heatmaps/gene_cluster/"

    commonOpts <- list(echo=FALSE, error=FALSE, message=FALSE, warning=FALSE, results='asis')
    figDefOpts <- c(commonOpts, list(fig.keep='all'))

    knitr::opts_template$set(htmlasis = c(figDefOpts, list(cache=FALSE)))

    knitr::opts_template$set(mdsopts = c(figDefOpts, list(
        fig.height=7, fig.width=12, out.width='100%', dpi=350,
        fig.path = "figures/mds/", include=TRUE)))

    knitr::opts_template$set(tsnechunk = c(figDefOpts, list(
        dpi=350, out.height="650px", out.width="3650px")))

    knitr::opts_template$set(hmopts = c(figDefOpts, list(fig.width=18, fig.height=22, dpi=350)))

    knitr::opts_template$set(cn_opts_1 = c(figDefOpts, list(
        fig.height=8, fig.width=15, fig.asp=0.75, fig.path = "figures/cnv/")))

    knitr::opts_template$set(gencnv = c(figDefOpts, list(fig.path="figures/cnv/")))

    knitr::opts_template$set(cn_opts_2 = c(figDefOpts, list(
        fig.height=10, fig.width=16, dpi=300, out.width='100%', fig.path='figures/cnv/')))

    knitr::opts_template$set(clustprof = c(commonOpts, list(
        fig.path="figures/cluster/", dpi=300, fig.height=5, fig.width=8, out.width='800px')))

    knitr::opts_template$set(pathhm = c(figDefOpts, list(
        fig.path=geneClusPath, fig.align='left', fig.width=18, fig.height=22, dpi=350)))

    knitr::opts_template$set(genepath = c(figDefOpts, list(
        fig.path = "figures/pathway/", fig.height=7, fig.width=12, out.width='100%', dpi=350)))
}


gb$cssLink <- "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/Research/Research.css"
