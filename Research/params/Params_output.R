#!/usr/bin/env Rscript
## Script name: Params_output.R
## Purpose: Source global parameters for research report output file naming
## Date Created: May 17, 2022
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024

gb <- globalenv(); assign("gb", gb)
gb$cssLink <-
    "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/Research/Research.css"

# Top Variance probe(s) to pull ----------------------------------------------
varProbes <- c(100, 1000, 10000) # Vector of top N probes to plot T-sne
hmVarProbes <- c(100, 1000, 3000) # Vector of top N probes to plot HeatMaps
topN = 1000    # Top N probes to use for MDS plot
topVar = 10000 # Top Variance probes to pull

# Plot Title Names ------------------------------------------------------------
clusType <- c("Unsupervised", "Supervised")
if (gb$supervisedRun == F) clusType <- clusType[1]

# Generates output plot titles using Params_input.R titleMain and probes
gb$tsne_titles <- gb$generateTitles(
    clusType, topTitle = as.character(varProbes), gb$titleMain)

# MDS plot title
mdsTitle <- paste("Top", topN, "Common", "mSet Sq.beta", "MDS plot")

# Default Output file names --------------------------------------------------
samsheet <- "samplesheet.csv" # default output primary sample sheet name/path
td <- format(Sys.Date(),"%b%d") # today's date (i.e., "Jun17")
gb$rgOut <- file.path(gb$runDir, "data", paste0(td, "_RGset.Rdata"))
rawBetaFi <- file.path(gb$runDir, "data", paste0(td, "_betaRaw.Rdata"))
pValsOutFi <- file.path(gb$runDir, "data", paste0(td, "_DetPvals.Rdata"))
mbfile <- file.path(gb$runDir, "data", paste0(td, "_mSetSq.beta.Rdata"))
unbetaVariance <- file.path(gb$runDir, "data", paste0(td, "_unbetaVar.RData"))
combatOut <- file.path(gb$runDir, "data", paste0(td, "_combatBetas.Rdata"))

# Default Output Paths -------------------------------------------------------
tsneOutU = "figures/tsne/unsupervised/"
tsneOutS = "figures/tsne/supervised/"
hmOutU = "figures/heatmaps/unsupervised/"
hmOutS = "figures/heatmaps/supervised/"
supbetaOut <- file.path("data", paste0(td, "_supervisedBetas"))
annotFi <- file.path("data", paste0(td, "_annotations.rds"))

# Create OutPut Directories ---------------------------------------------------
CheckDirCreate <- function(path_dirs) {
    for (dir_name in path_dirs) {
        out_dir <- file.path(gb$runDir, dir_name)
        if (!dir.exists(out_dir)) dir.create(out_dir, recursive = T)
    }
}

CheckDirCreate(c("data", "idats", "figures/mds", "csv", tsneOutU))
if (gb$supervisedRun == T) CheckDirCreate(tsneOutS)
if (gb$genHeatMaps == T) {
    CheckDirCreate(hmOutU)
    if (gb$supervisedRun) CheckDirCreate(hmOutS)
}

if (gb$genCNchunk == T) {
    cnvSegDir <- file.path("data", "CNV_segments")
    cnvPath <- file.path("figures", "cnv", "interactive_cnv")
    CheckDirCreate(c(cnvSegDir, cnvPath))
}

if (gb$genPairwise == T) CheckDirCreate("/figures/diffmean/")
if (gb$genPathChunk == T) {
    CheckDirCreate(c("figures/pathway/", "figures/cluster/"))
    CheckDirCreate("figures/heatmaps/gene_cluster/")
}

# Segments Copy Number Output Files -------------------------------------------
if (gb$genCNchunk == T) {
    segFile <- file.path(cnvSegDir, paste0(td, "_segfile.csv"))
    seg_clust_file <- file.path(cnvSegDir, paste0(td, "_seg_clusters.txt"))
}
gsetFile <- file.path(gb$runDir, "data", paste0(td, "_gsetfunnorm.rds"))
gBetaFile <- file.path(gb$runDir, "data", paste0(td, "_gsetbeta.rds"))
ClusfiNam <- file.path(gb$runDir, "data", paste0(td, "_dmp.csv"))

# Pathway Analysis Output Files -----------------------------------------------
if (gb$genPathChunk) {
    siteSpecific <- as.data.frame(readxl::read_excel(gb$genesInputFi))
    stopifnot("Gene" %in% colnames(siteSpecific))
    your_genes <- siteSpecific$Gene
    geneListIn <- file.path(gb$runDir, "data", paste0(td,"_Genes_Pathway.csv"))
    pathCsvOut <- file.path(gb$runDir, "data", paste0(td,"-signal_paths.csv"))
}

# Chunk Header Options --------------------------------------------------------
assignOpts <- function() {
    geneClusPath = "figures/heatmaps/gene_cluster/"

    commonOpts <- list(
        echo = FALSE, error = FALSE, message = FALSE, warning = FALSE,
        results = 'asis'
        )

    figDefOpts <- c(commonOpts, list(fig.keep = 'all'))

    knitr::opts_template$set(mdsopts = c(figDefOpts, list(
        fig.height = 7, fig.width = 12, out.width = '100%', dpi = 350,
        fig.path = "figures/mds/", include = TRUE)))

    knitr::opts_template$set(tsnechunk = c(figDefOpts, list(
        dpi = 350, out.height = "650px", out.width = "3650px")))

    knitr::opts_template$set(hmopts = c(figDefOpts, list(
        fig.width = 18, fig.height = 22, dpi = 350)))

    knitr::opts_template$set(cn_opts1 = c(figDefOpts, list(
        fig.height = 8, fig.width = 15, fig.asp = 0.75,
        fig.path = "figures/cnv/", include = FALSE)))

    knitr::opts_template$set(gencnv = c(figDefOpts, list(
        fig.path = "figures/cnv/")))

    knitr::opts_template$set(cn_opts2 = c(figDefOpts, list(
        fig.height = 10, fig.width = 16, dpi = 300, out.width = '100%',
        fig.path = 'figures/cnv/', include = TRUE)))

    knitr::opts_template$set(clustprof = c(commonOpts, list(
        fig.path = "figures/cluster/", dpi = 300,
        fig.height = 5, fig.width = 8, out.width = '800px')))

    knitr::opts_template$set(pathhm = c(figDefOpts, list(
        fig.path = geneClusPath, fig.align = 'left',
        fig.width = 18, fig.height = 22, dpi = 350)))

    knitr::opts_template$set(genepath = c(figDefOpts, list(
        fig.path = "figures/pathway/",
        fig.height = 7, fig.width = 12, out.width = '100%', dpi = 350)))
}
