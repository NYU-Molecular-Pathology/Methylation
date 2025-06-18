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
merged_name = NULL # Optional argument if two or more genes combined annotation
min_probes = 4  # Minimum number of probes to bin for genes
doXY = TRUE

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
    "fs",
    "readr",
    "S4Vectors",
    "IRanges",
    "GenomicRanges",
    "AnnotationDbi",
    "biomaRt",
    "ExperimentHub",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    "sesameData",
    "sesame",
    "minfi",
    "stringr",
    "dplyr",
    "ggplot2",
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
        suppressWarnings(suppressPackageStartupMessages(library(
            pkg, mask.ok = TRUE, character.only = TRUE, logical.return = TRUE
        )))
    })))
}

if (not_installed("devtools")) install.packages("devtools", ask = F, dependencies = T)
library("devtools")
if (not_installed("BiocManager")) install.packages("BiocManager", ask = F, dependencies = T)
ensure_packages(pkgs)

if (not_installed("conumee2")) {
    devtools::install_github("hovestadtlab/conumee2", subdir = "conumee2", upgrade = "always")
}
ensure_packages("conumee2")

# If no usable Python is found, bootstrap a Miniconda env named "r-reticulate"
if (!reticulate::py_available(initialize = FALSE)) {
    if (is.null(reticulate::miniconda_path())) try(reticulate::install_miniconda(), TRUE)
    reticulate::use_condaenv("r-reticulate", required = TRUE)
    if (!reticulate::py_module_available("kaleido")) {
        reticulate::conda_install("r-reticulate", "python-kaleido")
    }
    if (!reticulate::py_module_available("plotly")) {
        reticulate::conda_install("r-reticulate", "plotly", channel = "plotly")
    }
}


# Helper functions for python module pip installation
pip_install <- function(pkg) {
    reticulate::py_run_string(
        paste0("import sys; import subprocess;", " ", "subprocess.check_call",
               "([sys.executable, '-m', 'pip', 'install', '", pkg, "'])")
    )
}

# Helper loads Python module with reticulate: kaleido is needed for plotly image export
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
get_gene_ranges <- function(geneNames, array_info) {

    if (array_info$genome == "hg19") {
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
        values = geneNames,
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


# Helper function: for each gene, counts overlapping probes
count_gene_probes <- function(gene_gr, probe_gr) {
    hits <- GenomicRanges::findOverlaps(gene_gr, probe_gr)
    hits <- as(hits, "Hits")
    length(unique(S4Vectors::subjectHits(hits)))
}

# Aligns query and reference CNV objects, ensuring common probes
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

# Alters appearance of plot points, colors, and arrows in the final plot
modify_markers <- function(p_interactive, gene_gr) {
    gene_order <- p_interactive[["x"]][["data"]][[5]][["text"]]
    is_found <- gene_order %in% gene_gr$symbol
    orig_m <- p_interactive[["x"]][["data"]][[5]][["marker"]][["size"]]
    marks <- rep(orig_m, length(gene_order))
    marks[is_found] <- 12
    p_interactive[["x"]][["data"]][[5]][["marker"]][["size"]] <- marks

    orig_c <- p_interactive[["x"]][["data"]][[5]][["marker"]][["color"]]
    new_cols <- rep(orig_c, length(gene_order))
    new_cols[is_found] <- "#fdbc4b"
    p_interactive[["x"]][["data"]][[5]][["marker"]][["color"]] <- new_cols

    x_points <- as.numeric(paste(p_interactive[["x"]][["data"]][[5]]$x))
    y_points <- as.numeric(paste(p_interactive[["x"]][["data"]][[5]]$y))

    p_interactive <- p_interactive %>% plotly::layout(
        annotations = lapply(seq_along(x_points), function(i) {
            list(
                x = x_points[i], y = y_points[i],
                text = p_interactive[["x"]][["data"]][[5]][["text"]][i],
                showarrow = TRUE, arrowsize  = 0.7,
                ax = 30, ay = -60, arrowhead = 3, arrowwidth = 2.5,
                bgcolor     = "rgba(255,255,255,0.7)",
                bordercolor = "#0a4209",
                borderpad   = 4,
                font = list(family = "Arial Black", size = 14, color = "black")
            )
        })
    )

    return(p_interactive)
}

# Processes gene ranges and checks probe counts
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

# Helper function: Generates CNV object from aligned query and reference objects
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

# Generates the final plot with custom gene annotations
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

# Loads pre-computed normal sample reference files for CNV analysis
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
    message("Saving file:\n", out_html)
    filesDir <- paste0(tools::file_path_sans_ext(basename(out_html)), "_files")

    htmlwidgets::saveWidget(p_interactive, out_html, libdir = filesDir)
    outFiles <- file.path(dirname(out_html), filesDir)
    try(unlink(filesDir, recursive = TRUE, force = TRUE), silent = TRUE)

    png_file <- paste0(cnv_obj@name, "_cnv_plot_anno.png")
    out_png <- file.path(OUTPUT_DIR, png_file)
    message("Saving file:\n", file.path(getwd(), png_file))
    plotly::save_image(p_interactive, out_png, width = "1000", height = "600",
                       scale = 2.0)
}


merge_gene_regions <- function(gene_gr, merged_name) {
    chrs <- unique(as.character(GenomeInfoDb::seqnames(gene_gr)))
    if (length(chrs) != 1) {
        stop("Cannot merge genes on multiple chromosomes.")
    }
    merged_start <- min(GenomicRanges::start(gene_gr))
    merged_end   <- max(GenomicRanges::end(gene_gr))
    gene_gr <- GenomicRanges::GRanges(
        seqnames = chrs,
        ranges   = IRanges::IRanges(start = merged_start, end = merged_end),
        strand   = "*"
    )
    S4Vectors::mcols(gene_gr)$symbol <- merged_name
    S4Vectors::mcols(gene_gr)$name <- merged_name
    return(gene_gr)
}


con_path        <- file.path(path.package("conumee2"), "data")
exclude_regions <- get(load(file.path(con_path, "exclude_regions.rda")))
load(file.path(con_path, "exclude_regions.rda"))

annotate_cnv_with_genes <- function(sam_data, geneNames, merged_name = NULL, refs = NULL) {

    RGset <- get_rgset(sam_data)
    array_info <- detect_array_info(RGset)
    gene_gr <- get_gene_ranges(geneNames, array_info)

    if (!is.null(merged_name)) {
        gene_gr <- merge_gene_regions(gene_gr, merged_name)
    }

    Mset <- minfi::preprocessIllumina(RGset, normalize = "controls")
    query_cnv <- conumee2::CNV.load(Mset)
    anno <- conumee2::CNV.create_anno(
        array_type = array_info$array,
        chrXY = doXY,
        genome = array_info$genome,
        exclude_regions = exclude_regions,
        detail_regions = gene_gr
    )

    ref <- NULL
    ref <- if (is.null(refs)) {get_reference(array_info)} else {
        switch(array_info$genome, hg19 = refs$hg19, hg38 = refs$hg38)
    }

    aligned <- alignCNVObjects(query_cnv, ref, anno)
    colnames(aligned$ref@intensity) <- colnames(ref@intensity)
    cnv_obj <- get_cnv_obj(aligned, sam_data$Sample_Name)
    p_interactive <- get_final_plot(cnv_obj, anno, gene_gr)

    save_output(cnv_obj, p_interactive)
    #list(plot = p_interactive, cnv_analysis = cnv_obj)
}


read_zipfile <- function(zip_url) {
    tmp <- tempfile(fileext = ".zip")
    download.file(zip_url, tmp, mode = "wb")
    csv_name <- unzip(tmp, list = TRUE)$Name[1]
    man_epic <- suppressWarnings(readr::read_csv(unz(tmp, csv_name), skip = 7,
                                                 show_col_types = FALSE))
    unlink(tmp)
    man_epic <- as.data.frame(man_epic)
    return(man_epic)
}


load_manifests <- function() {
    illumina_url <- "https://webdata.illumina.com/downloads/productfiles"
    man_450k <- suppressWarnings(readr::read_csv(
        skip = 7, file.path(illumina_url, "humanmethylation450",
                            "humanmethylation450_15017482_v1-2.csv"),
        show_col_types = FALSE))
    man_450k <- as.data.frame(man_450k)

    zip_url <- file.path(illumina_url, "methylationEPIC",
                         "infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip")
    man_epic <- read_zipfile(zip_url)

    zip_url <- file.path(illumina_url, "methylationEPIC",
                         "infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip")
    man_epic2 <- read_zipfile(zip_url)

    array_manifests <- list("450k" = man_450k, "EPIC" = man_epic, "EPICv2" = man_epic2)

    return(array_manifests)
}

# Returns list of genes with enough probes to plot and saves coverage summary
evaluate_gene_coverage <- function(geneNames, array_manifests, exclude_probes = NULL) {
    tile_size = 200
    summary_df <- data.frame(Gene = geneNames, stringsAsFactors = FALSE)
    detail_list <- list()

    for (array in names(array_manifests)) {
        df <- array_manifests[[array]]
        if (!is.null(exclude_probes)) {
            df <- df[!df$IlmnID %in% exclude_probes, ]
        }
        df <- df[!is.na(df$CHR) & !is.na(df$MAPINFO), ]

        probe_gr <- GenomicRanges::GRanges(
            seqnames = paste0("chr", df$CHR),
            ranges   = IRanges::IRanges(start = df$MAPINFO, width = 1),
            probe_id = df$IlmnID,
            Island   = df$Relation_to_UCSC_CpG_Island,
            Region   = df$UCSC_RefGene_Group
        )

        per_array <- lapply(geneNames, function(g) {
            gene_gr <- get_gene_ranges(
                g, list(array = array, genome = ifelse(array == "EPICv2", "hg38", "hg19"))
            )
            if (is.null(gene_gr)) {
                return(data.frame(
                    Gene           = g,
                    Array          = array,
                    nProbes        = 0L,
                    pctTiles       = 0,
                    tilesCovered   = 0L,
                    nTiles         = 0L,
                    GeneRegion     = NA_character_,
                    ProbeIDs       = I(list(character())),
                    CpG_Islands    = I(list(character())),
                    RefGene_Groups = I(list(character())),
                    stringsAsFactors = FALSE
                ))
            }

            # tile the gene into fixed-width bins
            gene_width <- sum(GenomicRanges::width(gene_gr))
            nTiles     <- ceiling(gene_width / tile_size)
            tiles      <- GenomicRanges::tile(gene_gr, n = nTiles)[[1]]

            # count overlaps per tile
            hits_tiles    <- GenomicRanges::countOverlaps(tiles, probe_gr)
            tilesCovered  <- sum(hits_tiles > 0)
            pctTiles      <- round(tilesCovered / nTiles * 100, 1)
            hits          <- GenomicRanges::findOverlaps(gene_gr, probe_gr)
            idx           <- unique(S4Vectors::subjectHits(hits))
            nProbes       <- length(idx)
            probe_ids     <- df$IlmnID[idx]
            islands       <- unique(df$Relation_to_UCSC_CpG_Island[idx])
            regions       <- unique(df$UCSC_RefGene_Group[idx])
            region_string <- paste0(as.character(GenomeInfoDb::seqnames(gene_gr)),
                                    ":", GenomicRanges::start(gene_gr),
                                    "-", GenomicRanges::end(gene_gr))

            data.frame(
                Gene           = g,
                Array          = array,
                nProbes        = nProbes,
                pctTiles       = pctTiles,
                tilesCovered   = tilesCovered,
                nTiles         = nTiles,
                GeneRegion     = region_string,
                ProbeIDs       = if (length(probe_ids) > 0) paste(probe_ids, collapse = ";") else "",
                CpG_Islands    = if (length(islands)  > 0) paste(islands,    collapse = ";") else "",
                RefGene_Groups = if (length(regions)  > 0) paste(regions,    collapse = ";") else "",
                stringsAsFactors = FALSE
            )

        })

        per_array_df <- do.call(rbind, per_array)
        summary_df[[paste0(array, "_nProbes")]] <- per_array_df$nProbes
        summary_df[[paste0(array, "_pctTiles")]] <- per_array_df$pctTiles
        detail_list[[array]] <- per_array_df
    }

    write.table(
        summary_df,
        file      = file.path(OUTPUT_DIR, "gene_coverage_summary.tsv"),
        sep       = "\t",
        row.names = FALSE,
        quote     = FALSE
    )
    detail_df <- do.call(rbind, detail_list)

    write.table(
        detail_df,
        file      = file.path(OUTPUT_DIR, "gene_coverage_details.tsv"),
        sep       = "\t",
        row.names = FALSE,
        quote     = FALSE
    )

    # filter genes
    maxP <- apply(summary_df[, grep("_nProbes$", names(summary_df))],
                  1, max, na.rm = TRUE)
    good <- summary_df$Gene[maxP >= min_probes]
    bad  <- summary_df$Gene[maxP <  min_probes]
    if (length(bad)) warning("Excluded (low coverage): ", paste(bad, collapse = ", "))

    final_out <- invisible(list(
        goodGenes = good,
        summary   = summary_df,
        details   = detail_df
    ))
    return(final_out$good)
}


# Helper function: Loads pre-computed reference files if available
load_references <- function(manifest_path) {
    if (!dir.exists(manifest_path)) return(NULL)
    message("Loading ref_v1_hg19.rds")
    hg19_ref <- readRDS(file.path(manifest_path, "ref_v1_hg19.rds"))
    message("Loading ref_v2_hg38.rds")
    hg38_ref <- readRDS(file.path(manifest_path, "ref_v2_hg38.rds"))
    return(list(hg19 = hg19_ref, hg38 = hg38_ref))
}


# MAIN: Load references if available and iterate through sample sheet to render plots ------------
loop_samples <- function(manifest_path, sam_sheet, geneNames, merged_name) {
    targets <- read.csv(sam_sheet)
    refs <- load_references(manifest_path)
    array_manifests <- load_manifests()
    geneNames <- evaluate_gene_coverage(geneNames, array_manifests)
    for (samRow in 1:nrow(targets)) {
        sam_data <- targets[samRow, ]
        message(paste0(capture.output(sam_data), collapse = "\n"))
        annotate_cnv_with_genes(sam_data, geneNames, merged_name, refs = refs)
    }
}


loop_samples(manifest_path, sam_sheet, geneNames, merged_name)
