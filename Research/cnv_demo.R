#!/usr/bin/env Rscript
## Script name: conumee2_demo.R
## Purpose: Demos how to use conumee2 to annotate genes and save genomeplots
## Date Created: July 30, 2025
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025

# Helper function to silently load packages
quite_load <- function(pkg) suppressWarnings(suppressPackageStartupMessages(library(
    pkg, character.only = TRUE, logical.return = TRUE)))

options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Helper to install 'pak'
install_pak <- function() {
    pak_arch <- sprintf("https://r-lib.github.io/p/pak/stable/%s/%s/%s",
                        .Platform$pkgType, R.Version()$os, R.Version()$arch)
    tryCatch(
        install.packages("pak", repos = pak_arch, dependencies = T, ask = F),
        error = function(e) install.packages("pak", dependencies = T, ask = F)
    )
}

# Load/install packages; prioritizes BiocManager
ensure_packages <- function(..., bioconductor_pkgs = NULL) {
    pkgs <- c(...)
    all_pkgs_to_check <- unique(c(pkgs, bioconductor_pkgs))
    installed <- rownames(installed.packages())

    # Install BiocManager if needed for Bioconductor packages
    if (!is.null(bioconductor_pkgs) && length(bioconductor_pkgs) > 0) {
        if (!"BiocManager" %in% installed) {
            message("Installing 'BiocManager'...")
            install.packages("BiocManager", dependencies = TRUE, ask = FALSE)
        }
        quite_load("BiocManager")
        missing_bioc_pkgs <- setdiff(bioconductor_pkgs, installed)
        if (length(missing_bioc_pkgs) > 0) {
            message(paste("Installing Bioc packages:", paste(missing_bioc_pkgs, collapse = ", ")))
            BiocManager::install(missing_bioc_pkgs, update = FALSE, ask = FALSE)
        }
    }

    # Load all required packages and report status
    load_results <- sapply(all_pkgs_to_check, quite_load)
    for (i in seq_along(load_results)) {
        pkg_name <- names(load_results)[i]
        status <- ifelse(load_results[i], "TRUE", "FALSE")
        message(paste0(pkg_name, " loaded... ", status))
    }

    stopifnot(all(load_results)) # Ensure all packages were successfully loaded
    message("All packages loading complete.")
}

# Package lists - Now includes org.Hs.eg.db for gene symbol mapping
CRAN_PACKAGES <- c(
    "data.table", "RColorBrewer", "lattice", "nor1mix",
    "beanplot", "reshape", "MASS", "quadprog", "GEOquery",
    "stats", "grDevices", "graphics", "utils"
)

BIOCONDUCTOR_PACKAGES <- c(
    "methods", "BiocGenerics", "GenomicRanges", "SummarizedExperiment",
    "Biostrings", "bumphunter", "S4Vectors", "GenomeInfoDb",
    "Biobase", "IRanges", "siggenes", "limma", "preprocessCore",
    "illuminaio", "DelayedMatrixStats", "mclust", "genefilter", "nlme",
    "DelayedArray", "HDF5Array", "BiocParallel", "ExperimentHub",
    "minfi", "minfiDataEPIC", "FlowSorted.Blood.EPIC",
    "GenomicFeatures", "AnnotationDbi", "TxDb.Hsapiens.UCSC.hg19.knownGene",
    "org.Hs.eg.db"
)

# Ensure packages are installed and loaded
ensure_packages(CRAN_PACKAGES, bioconductor_pkgs = BIOCONDUCTOR_PACKAGES)


# Replace "conumee2" with "conumee2" to match the new package name
if (!"conumee2" %in% rownames(installed.packages())) {
    devtools::install_github("hovestadtlab/conumee2", subdir = "conumee2")
}

# ------------------------ User Input ------------------------
GENES_TO_HIGHLIGHT <- c("TP53","EGFR", "BRCA1", "BRCA2")

# 1. Load public EPIC example data from minfi (case samples) and controls
rgSet_case <- minfiDataEPIC::RGsetEPIC

hub <- ExperimentHub::ExperimentHub()
rgSet_control <- hub[["EH1136"]]

# 2. Preprocess cases and control data
mSet_case <- minfi::preprocessIllumina(rgSet_case, bg.correct = TRUE,
                                       normalize = "controls")
mSet_control <- minfi::preprocessIllumina(rgSet_control, bg.correct = TRUE,
                                          normalize = "controls")

# 3. Create separate CNV.data objects for query and reference
cnv_data <- conumee2::CNV.load(input = mSet_case)
cnv_ref <- conumee2::CNV.load(input = mSet_control)


# 4. (Optional) Create annotation object for detail regions
check_gene_anno <- function(GENES_TO_HIGHLIGHT) {
    if (length(GENES_TO_HIGHLIGHT) == 0) {
        return(NULL)
    }
    # Map SYMBOL → ENTREZID and remove missing values
    gene_ids <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
        keys = GENES_TO_HIGHLIGHT, column = "ENTREZID",
        keytype = "SYMBOL", multiVals = "first"
    )
    gene_ids <- gene_ids[!is.na(gene_ids)]

    if (length(gene_ids) == 0) {
        warning("None of the specified gene symbols were found.")
        return(NULL)
    }

    # Load all known genes and subset by Entrez IDs
    txdb <-
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    all_genes <- suppressMessages(GenomicFeatures::genes(txdb))
    genes_of_interest <-
        all_genes[S4Vectors::mcols(all_genes)$gene_id %in% gene_ids]

    if (length(genes_of_interest) == 0) {
        warning("No genomic regions matched the requested genes.")
        return(NULL)
    }

    # Annotate each range with its SYMBOL
    gene_entrez <- S4Vectors::mcols(genes_of_interest)$gene_id
    gene_symbol <- names(gene_ids)[match(gene_entrez, gene_ids)]
    S4Vectors::mcols(genes_of_interest)$name <- gene_symbol
    S4Vectors::mcols(genes_of_interest)$gene_id <- gene_symbol

    # Build and return the annotation object
    anno <- conumee2::CNV.create_anno(array_type = "EPIC", genome = "hg19",
         detail_regions = genes_of_interest, chrXY = TRUE
    )

    return(anno)
}

anno <- check_gene_anno(GENES_TO_HIGHLIGHT)


names(cnv_data@intensity) <- colnames(mSet_case)
# 5. Fit and segment CNVs
cnv_fit <- conumee2::CNV.fit(query = cnv_data, ref = cnv_ref, anno = anno)
cnv_bin <- conumee2::CNV.bin(cnv_fit)
cnv_seg <- conumee2::CNV.segment(cnv_bin)

#  Call CNV.detail if you want to add in detail regions to the CNV
if (!is.null(anno)) {
    cnv_detail <- conumee2::CNV.detail(cnv_seg)
    print(cnv_detail)
}

# 6. Generate genome-wide CNV plot
cnv_detail@name <- names(cnv_detail@seg$p)

for (sam_name in 1:length(cnv_detail@name)) {
    sam_obj <- cnv_detail[sam_name]
    output_png_name <- paste(sam_name, "cnv_genome_plot.pdf", sep = "_")
    png(
        filename = output_png_name,
        width = 1200,
        height = 800,
        res = 150
        )
    conumee2::CNV.genomeplot(sam_obj,
                             main = paste(sam_name, "CNV Genome Plot"))
    dev.off()
}
