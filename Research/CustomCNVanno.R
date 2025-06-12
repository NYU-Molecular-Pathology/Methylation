#!/usr/bin/env Rscript
## Script name: CustomCNVanno.R
## Purpose: Functions for annotating custom gene details to CNV genomeplot
## Date Created: August 10, 2023
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025


# Global Paths to Reference Files ---------------------------------------------
manifest_path = "/Volumes/CBioinformatics/Methylation/Manifests"

# User Input: Project Directory and Sample Sheet -----------------------------
PROJ_DIR <- "/Volumes/CBioinformatics/Methylation/Clinical_Runs/pull_redcap_idats"
OUTPUT_DIR <- file.path(PROJ_DIR, "CNV_annotated_plots")
sam_sheet <- "samplesheet_rds.csv"
geneNames <- c("BRCA1", "FOXA1", "TP53")
min_probes = 4  # Minimum number of probes to bin for genes

if (!dir.exists(OUTPUT_DIR)) {dir.create(OUTPUT_DIR, recursive = TRUE)}
setwd(PROJ_DIR)

snapshot_date <- "2025-05-01"
bioc_version   <- "3.19"
options(
    repos = c(
        CRAN = sprintf("https://cran.microsoft.com/snapshot/%s/", snapshot_date),
        BioC = sprintf("https://packagemanager.rstudio.com/bioconductor/%s@%s",
                       bioc_version, snapshot_date)
    )
)

pkgs <- c(
    "devtools",
    "sesame",
    "minfi",
    "S4Vectors",
    "GenomicRanges",
    "IRanges",
    "dplyr",
    "stringr",
    "biomaRt",
    "AnnotationDbi",
    "plotly",
    "ggplot2",
    "ExperimentHub",
    "sesameData",
    "rmarkdown",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    "htmlwidgets",
    "plotly",
    "reticulate"
)

# Helper function: Check if a package is not installed
not_installed <- function(pkgName) return(!pkgName %in% rownames(installed.packages()))

# Helper function: Install pak package if not already installed
install_pak <- function() {
    tryCatch(
        install.packages(
            "pak", repos = sprintf("https://r-lib.github.io/p/pak/stable/%s/%s/%s",
                                   .Platform$pkgType, R.Version()$os, R.Version()$arch)),
        error = function(e) {
            install.packages("pak", ask = FALSE, dependencies = TRUE,
                             repos = "https://packagemanager.rstudio.com/all/latest"
            )
        })
}

# Function to ensure required packages are installed and loaded
ensure_packages <- function(pkgs) {
    installed_pkgs <- rownames(installed.packages())
    missing_pkgs <- setdiff(pkgs, installed_pkgs)
    if (length(missing_pkgs) > 0) {
        message(missing_pkgs)
        if (!"pak" %in% installed_pkgs) install_pak()
        library("pak")
        for (pkg in missing_pkgs) {
            tryCatch(
                pak::pkg_install(pkg, ask = FALSE),
                error = function(e) {
                    install.packages(pkg, ask = FALSE, dependencies = TRUE)
                }
            )
        }
    }
    stopifnot(all(sapply(pkgs, function(pkg) {
        suppressPackageStartupMessages(library(
            pkg, mask.ok = TRUE, character.only = TRUE, logical.return = TRUE
        ))
    })))
}

if (not_installed("devtools")) install.packages("devtools", ask = F, dependencies = T)
library("devtools")
ensure_packages(pkgs)

if (not_installed("conumee2")) {
    devtools::install_github("hovestadtlab/conumee2", subdir = "conumee2", upgrade = "always")
}
ensure_packages("conumee2")

# If no usable Python is found, bootstrap a Miniconda env named "r-reticulate"
if (!reticulate::py_available(initialize = FALSE)) {
    reticulate::install_miniconda()
    reticulate::conda_install("r-reticulate", "python-kaleido")
    reticulate::conda_install("r-reticulate", "plotly", channel = "plotly")
    reticulate::use_condaenv("r-reticulate", required = TRUE)
}


# Helper functions for plotly pip save html widget debug
pip_install <- function(pkg) {
    reticulate::py_run_string(
        paste0("import sys; import subprocess;", " ", "subprocess.check_call",
               "([sys.executable, '-m', 'pip', 'install', '", pkg, "'])")
    )
}

# Helper to load Python module using reticulate for kaleido needed for plotly image export
load_module <- function(pkg) {
    if (!reticulate::py_module_available(pkg)) pip_install(pkg)
    reticulate::py_run_string(paste("import", pkg))
    reticulate::import(pkg)
}

reticulate::py_run_string("import sys")
load_module("kaleido")
load_module("plotly")


# Helper function: detect array type and genome build from an RGset
detect_array_info <- function(RGset) {
    array_type <- "EPIC"
    if (RGset@annotation[["array"]] == "IlluminaHumanMethylation450k") {
        array_type <- "450k"
    }
    if (RGset@annotation[["array"]] == "IlluminaHumanMethylationEPICv2") {
        array_type <- "EPICv2"
    }
    genome <- ifelse(array_type == "EPICv2", "hg38", "hg19")
    array_info <- list(array = array_type, genome = genome)
    message(paste("Array type:", array_info$array, "Genome build:", array_info$genome))
    return(array_info)
}


# Helper function: retrieve gene ranges (e.g., TSS/promoter region) for a given gene symbol.
get_gene_ranges <- function(geneSymbols, array_info) {

    if (array_info$genome != "hg19") {
        host_url <- "https://grch37.ensembl.org"  # Use the GRCh37 archive for hg19
    } else {
        host_url <- "https://www.ensembl.org"  # Use the main Ensembl site for hg38
    }

    all_marts <- biomaRt::listMarts(host = host_url)
    mart_name <- all_marts$biomart[1]
    mart_set <- all_marts$version[1]
    
    mart <- biomaRt::useMart(
        biomart = mart_name,
        dataset = "hsapiens_gene_ensembl",
        host = host_url,
        version = mart_set
    )

    geneInfo <- biomaRt::getBM(
        attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand"),
        filters = "hgnc_symbol",
        values = geneSymbols,
        mart = mart
    )
    if (nrow(geneInfo) == 0) return(NULL)
    gr <- GenomicRanges::GRanges(
        seqnames = paste0("chr", geneInfo$chromosome_name),
        ranges = IRanges::IRanges(start = geneInfo$start_position, end = geneInfo$end_position),
        strand = geneInfo$strand,
        symbol = geneInfo$hgnc_symbol
    )
    gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = "coarse")
    S4Vectors::mcols(gr)$name <- S4Vectors::mcols(gr)$symbol
    return(gr)
}


# Helper function: for each gene, count overlapping probes
count_gene_probes <- function(gene_gr, probe_gr) {
    hits <- GenomicRanges::findOverlaps(gene_gr, probe_gr)
    hits <- as(hits, "Hits")
    length(unique(S4Vectors::subjectHits(hits)))
}

# Function: Aligns query and reference CNV objects, ensuring common probes
alignCNVObjects <- function(query_obj, ref_obj, anno_obj) {
    qProbs <- rownames(query_obj@intensity)
    query_int_df <- as.data.frame(
        do.call(cbind, lapply(query_obj@intensity, function(x) data.frame(x))))
    rownames(query_int_df) <- qProbs

    rProbs <- rownames(ref_obj@intensity)
    ref_int_df   <- as.data.frame(
        do.call(cbind, lapply(ref_obj@intensity, function(x) data.frame(x))))
    rownames(ref_int_df) <- rProbs

    query_obj@intensity <- query_int_df
    ref_obj@intensity   <- ref_int_df

    # Remove duplicate rows in intensity data frames
    query_int <- query_obj@intensity
    ref_int <- ref_obj@intensity
    query_int <- query_int[!duplicated(rownames(query_int)), , drop = FALSE]
    ref_int <- ref_int[!duplicated(rownames(ref_int)), , drop = FALSE]
    query_obj@intensity <- query_int
    ref_obj@intensity <- ref_int

    # Extract probe IDs
    query_probes <- rownames(query_int)
    ref_probes <- rownames(ref_int)

    anno_probes <- anno_obj@probes@elementMetadata@listData[["IlmnID"]]
    if (is.null(anno_probes)) {
        anno_probes <- names(anno_obj@probes)
    } else {
        names(anno_obj@probes) <-
            anno_obj@probes@elementMetadata@listData[["IlmnID"]]
    }
    common_probes <- base::Reduce(intersect,
                                  list(query_probes, ref_probes, anno_probes))
    if (length(common_probes) == 0)
        stop("No common probes found among query, reference, and annotation.")
    query_obj@intensity <- query_obj@intensity[common_probes, , drop = FALSE]
    ref_obj@intensity <- ref_obj@intensity[common_probes, , drop = FALSE]
    anno_obj@probes <- anno_obj@probes[common_probes]
    list(query = query_obj, ref = ref_obj, anno = anno_obj)
}

# Function: Alters apearence of plot points, colors, and arrows in the output plot
modify_markers <- function(p_interactive, gene_gr) {
    gene_order <- p_interactive[["x"]][["data"]][[5]][["text"]]
    is_found <- gene_order %in% gene_gr$symbol
    orig_m <- p_interactive[["x"]][["data"]][[5]][["marker"]][["size"]]
    marks <- rep(orig_m, length(gene_order))
    marks[is_found] <- 12
    p_interactive[["x"]][["data"]][[5]][["marker"]][["size"]] <- marks

    orig_c <- p_interactive[["x"]][["data"]][[5]][["marker"]][["color"]]
    new_cols <- rep(orig_c, length(gene_order))
    new_cols[is_found] <- "yellow"
    p_interactive[["x"]][["data"]][[5]][["marker"]][["color"]] <- new_cols

    x_points <- as.numeric(paste(p_interactive[["x"]][["data"]][[5]]$x))
    y_points <- as.numeric(paste(p_interactive[["x"]][["data"]][[5]]$y))

    p_interactive <- p_interactive %>% plotly::layout(
        annotations = lapply(seq_along(x_points), function(i) {
            list(
                x = x_points[i], y = y_points[i],
                text = p_interactive[["x"]][["data"]][[5]][["text"]][i],
                showarrow = TRUE,
                ax = 20, ay = -30, arrowhead = 4, arrowwidth = 0.8,
                font = list(size = 12, color = "black")
            )
        })
    )

    return(p_interactive)
}

# Function: Process gene ranges and check for probe counts
process_gene_ranges <- function(gene_gr, anno) {
    if (is.null(gene_gr)) {
        warning("No gene ranges could be retrieved for the supplied gene names.")
        return("No gene ranges found for input gene list.")
    }

    probe_gr <- anno@probes
    gene_probe_counts <- sapply(split(gene_gr, gene_gr$symbol),
                                count_gene_probes, probe_gr = probe_gr)

    gene_anno_msg <- ""
    for (gene in names(gene_probe_counts)) {
        nprobes <- gene_probe_counts[gene]
        if (nprobes == 0) {
            gene_anno_msg <- paste0(gene_anno_msg, "No matches found for gene ", gene,
                                    " on the methylation array. ")
        } else if (nprobes < min_probes) {
            gene_anno_msg <- paste0(gene_anno_msg, "Caution: only ", nprobes,
                                    " probes annotated for gene ", gene, ". ")
        }
    }

    return(gene_anno_msg)
}

# Helper function: Generate CNV object from aligned query and reference objects
get_cnv_obj <- function(align_obj, sample_name) {
    cnv_obj <- conumee2::CNV.fit(query = align_obj$query,
                                 ref = align_obj$ref,
                                 anno = align_obj$anno)
    cnv_obj <- conumee2::CNV.bin(cnv_obj)
    cnv_obj <- conumee2::CNV.detail(cnv_obj)
    cnv_obj <- conumee2::CNV.segment(cnv_obj)
    cnv_obj@name <- sample_name
    return(cnv_obj)
}

# Function: Generates the final plot with custom gene annotations
get_final_plot <- function(cnv_obj, anno, gene_gr) {

    geneNames <- unique(gene_gr$symbol)

    p_interactive <- conumee2::CNV.plotly(cnv_obj)
    if (!is.null(gene_gr)) {
        p_interactive <- modify_markers(p_interactive, gene_gr)
    }

    if (length(geneNames) <= 5) {
        final_title <- paste(cnv_obj@name, "CNV Genome Plot With Genes –",
                             paste(geneNames, collapse = ", "))
    } else {
        final_title <- "CNV Genome Plot With Custom Gene Annotations"
    }

    gene_anno_msg <- process_gene_ranges(gene_gr, anno)
    if (nchar(gene_anno_msg) > 0) {
        final_title <- paste0(final_title, "\nNote: ", gene_anno_msg)
    }
    p_interactive <- plotly::layout(p_interactive, title = list(text = final_title))
    return(p_interactive)
}

# Function: loads pre-computed reference files for CNV analysis
get_reference <- function(array_info) {
    if (array_info$genome == "hg19") {
        ref_path <- file.path(manifest_path, "ref_v1_hg19.rds")
    } else {
        ref_path <- file.path(manifest_path, "ref_v2_hg38.rds")
    }
    if (file.exists(ref_path)) {
        message("Loading reference file...")
        ref <- readRDS(ref_path)
        return(ref)
    }

    if (array_info$array == "EPICv2") {
        message("Attempting to load EPICv2 reference from sesameData.")
        available_refs <- sesameData::sesameDataList()
        if ("EPICv2.8.SigDF" %in% available_refs$Title) {
            ref_list <- sesameData::sesameDataGet("EPICv2.8.SigDF")
            ref_mat <- do.call(cbind, lapply(ref_list, sesame::totalIntensities))
            ref_df <- as.data.frame(ref_mat, check.names = FALSE)
            ref <- conumee2::CNV.load(ref_df)
            return(ref)
        } else {
            stop("EPICv2.8.SigDF not available in sesameData::sesameDataList()")
        }
    }

    message("Using default reference from minfiData")
    dt_pkgs <- c("minfiDataEPIC", "minfiData")
    missing_pkgs <- sapply(pkgs, not_installed)
    if (any(missing_pkgs)) pak::pkg_install(pkgs[missing_pkgs], ask = FALSE)

    if (array_info$array == "EPIC") {
        library("minfiDataEPIC")
        data("MsetEPIC", package = "minfiDataEPIC", envir = environment())
        ref <- conumee2::CNV.load(get("MsetEPIC", envir = environment()))
    } else {
        library("minfiData")
        data("MsetEx", package = "minfiData", envir = environment())
        ref <- conumee2::CNV.load(get("MsetEx", envir = environment()))
    }

    return(ref)
}

# Helper function: Use minfi to get RGset and set genome build
get_rgset <- function(sam_data) {
    RGset <- minfi::read.metharray(sam_data$Basename, force = TRUE)
    if (RGset@annotation[["array"]] == "IlluminaHumanMethylationEPICv2") {
        RGset@annotation[["annotation"]] <- "20a1.hg38"
        RGset@annotation[["Genome"]] <- "hg38"
    } else{
        RGset@annotation[["Genome"]] <- "hg19"
    }
    return(RGset)
}

# Helper function: Save output files - HTML & PNG plots to current directory
save_output <- function(cnv_obj, p_interactive) {
    html_file <- paste0(cnv_obj@name, "_cnv_plot_anno.html")
    out_html <- file.path(OUTPUT_DIR, html_file)
    if (!file.exists(html_file)) {
        message("Saving file:\n", out_html)
        htmlwidgets::saveWidget(p_interactive, out_html)
    }

    png_file <- paste0(cnv_obj@name, "_cnv_plot_anno.png")
    out_png <- file.path(OUTPUT_DIR, png_file)
    if (!file.exists(png_file)) {
        message("Saving file:\n", file.path(getwd(), png_file))
        plotly::save_image(p_interactive, out_png, width = "1600", height = "900", scale = 2.5)
    }
}


# Main function: process RGset, perform CNV analysis, annotate input genes, plots generate plots
annotate_cnv_with_genes <- function(sam_data, geneNames, ref = NULL, doXY = TRUE) {

    RGset <- get_rgset(sam_data)
    array_info <- detect_array_info(RGset)
    gene_gr <- get_gene_ranges(geneNames, array_info)

    Mset <- minfi::preprocessIllumina(RGset, normalize = "controls")
    query_cnv <- conumee2::CNV.load(Mset)

    con_path <- file.path(path.package("conumee2"), "data")
    exclude_regions <- get(load(file.path(con_path, "exclude_regions.rda")))
    load(file.path(con_path, "exclude_regions.rda"))

    anno <- conumee2::CNV.create_anno(
        array_type = array_info$array,
        chrXY = doXY,
        genome = array_info$genome,
        exclude_regions = exclude_regions,
        detail_regions = gene_gr
    )

    aligned <- alignCNVObjects(query_cnv, ref, anno)
    colnames(aligned$ref@intensity) <- colnames(ref@intensity)

    cnv_obj <- get_cnv_obj(aligned, sam_data$Sample_Name)

    p_interactive <- get_final_plot(cnv_obj, anno, gene_gr)

    save_output(cnv_obj, p_interactive)
    #list(plot = p_interactive, cnv_analysis = cnv_obj)
}

# Main function: loop through samples in the sample sheet and call annotate_cnv_with_genes
loop_samples <- function(PROJ_DIR, sam_sheet, geneNames, ref = NULL) {
    targets <- read.csv(sam_sheet)

    if (is.null(ref)) {
        ref <- get_reference(array_info)
    }

    for (samRow in 1:nrow(targets)) {
        sam_data <- targets[samRow, ]
        annotate_cnv_with_genes(sam_data, geneNames, ref = ref)
    }
}

# Example usage:
loop_samples(PROJ_DIR, sam_sheet, geneNames)
