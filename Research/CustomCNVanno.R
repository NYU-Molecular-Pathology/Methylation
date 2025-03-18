#!/usr/bin/env Rscript
## Script name: CustomCNVanno.R
## Purpose: Functions for taking a list of input genes and output custom CNV annotation for Conumee
## Date Created: August 10, 2023
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025

# Example usage:
my_idat <- "/path/to/file.idat"
geneNames <- c("BRCA1", "FOXA1", "TP53")

manifest_path = "/Volumes/CBioinformatics/Methylation/Manifests"
ref_path <- "/Volumes/CBioinformatics/Methylation/Manifests/ref_v2_hg38.rds"
pkgs <- c(
    "devtools",
    "conumee2",
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
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
)

not_installed <- function(pkgName) {
    return(!pkgName %in% rownames(installed.packages()))
}

if ( !not_installed("conumee2.0")) {
    remove.packages("conumee2.0")
    devtools::install_github("hovestadtlab/conumee2", subdir = "conumee2")
}

missing_pkgs <- sapply(pkgs, not_installed)

if (not_installed("pak")) install.packages("pak", ask = F, dependencies = T)
if (any(missing_pkgs)) pak::pkg_install(pkgs[missing_pkgs], ask = FALSE)
stopifnot(all(sapply(pkgs, requireNamespace, quietly = FALSE)))


# Helper function: detect array type and genome build from an RGset.
detect_array_info <- function(RGset) {
    ann <- minfi::getAnnotation(RGset)
    array_type <- if (grepl("450k", RGset@annotation[["array"]], ignore.case = TRUE)) {
        "450k"
    } else if (grepl("EPICv2", RGset@annotation[["array"]], ignore.case = TRUE)) {
        "EPICv2"
    } else {
        "EPIC"
    }
    genome <- if (grepl("EPICv2", RGset@annotation[["array"]], ignore.case = TRUE)) "hg38" else "hg19"
    list(array = array_type, genome = genome)
}


liftOver_Mset_manifest <- function(mset, v2Man) {
    mapping_vec <- setNames(v2Man$EPICv1_Loci, v2Man$IlmnID)

    keep <- mset@NAMES %in% v2Man$IlmnID
    mset <- mset[keep, ]
    orig <- mset@NAMES
    new_names <- ifelse(orig %in% names(mapping_vec), mapping_vec[orig], orig)
    mset@NAMES <- new_names
    rownames(mset) <- new_names
    mset@annotation["annotation"] <- "ilm10b4.hg19"
    mset@annotation["IlluminaHumanMethylationEPIC"] <- "ilm10b4.hg19"

    return(mset)
}


# Helper function: retrieve gene ranges (e.g., TSS/promoter region) for a given gene symbol.
get_gene_ranges <- function(geneSymbols) {
    mart <- biomaRt::useMart(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = "hsapiens_gene_ensembl",
        host = "https://www.ensembl.org"
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
    ref_probes   <- rownames(ref_int)
    anno_probes <- anno_obj@probes@elementMetadata@listData[["IlmnID"]]
    names(anno_obj@probes) <-
        anno_obj@probes@elementMetadata@listData[["IlmnID"]]
    common_probes <- base::Reduce(intersect,
                                  list(query_probes, ref_probes, anno_probes))
    if (length(common_probes) == 0)
        stop("No common probes found among query, reference, and annotation.")
    query_obj@intensity <- query_obj@intensity[common_probes, , drop = FALSE]
    ref_obj@intensity <- ref_obj@intensity[common_probes, , drop = FALSE]
    anno_obj@probes <- anno_obj@probes[common_probes]
    list(query = query_obj, ref = ref_obj, anno = anno_obj)
}


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
    return(p_interactive)
}


process_gene_ranges <- function(gene_gr, anno, min_probes = 4) {
    if (is.null(gene_gr)) {
        warning("No gene ranges could be retrieved for the supplied gene names.")
        return("No gene ranges found for input gene list.")
    }

    probe_gr <- anno@probes
    gene_probe_counts <- sapply(split(gene_gr, gene_gr$symbol), count_gene_probes, probe_gr = probe_gr)

    gene_anno_msg <- ""
    for (gene in names(gene_probe_counts)) {
        nprobes <- gene_probe_counts[gene]
        if (nprobes == 0) {
            gene_anno_msg <- paste0(gene_anno_msg,
                                    "No matches found for gene ",
                                    gene,
                                    " on the methylation array. ")
        } else if (nprobes < min_probes) {
            gene_anno_msg <- paste0(gene_anno_msg,
                                    "Caution: only ",
                                    nprobes,
                                    " probes annotated for gene ",
                                    gene,
                                    ". ")
        }
    }

    return(gene_anno_msg)
}


get_cnv_obj <- function(align_obj) {
    cnv_obj <- conumee2::CNV.fit(query = align_obj$query,
                                 ref = align_obj$ref,
                                 anno = align_obj$anno)
    cnv_obj <- conumee2::CNV.bin(cnv_obj)
    cnv_obj <- conumee2::CNV.detail(cnv_obj)
    cnv_obj <- conumee2::CNV.segment(cnv_obj)
    return(cnv_obj)
}


get_final_plot <- function(cnv_obj, geneNames, gene_anno_msg, gene_gr) {
    p_interactive <- conumee2::CNV.plotly(cnv_obj)
    if (!is.null(gene_gr)) {
        p_interactive <- modify_markers(p_interactive, gene_gr)
    }

    if (length(geneNames) <= 5) {
        final_title <- paste("CNV Genome Plot With Genes –",
                             paste(geneNames, collapse = ", "))
    } else {
        final_title <- "CNV Genome Plot With Custom Gene Annotations"
    }

    if (nchar(gene_anno_msg) > 0) {
        final_title <- paste0(final_title, "\nNote: ", gene_anno_msg)
    }
    p_interactive <- plotly::layout(p_interactive, title = list(text = final_title))
    return(p_interactive)
}


get_reference <- function(array_info) {
    if (file.exists(ref_path)) {
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
    } else{
        library("minfiData")
        data("MsetEx", package = "minfiData", envir = environment())
        ref <- conumee2::CNV.load(get("MsetEx", envir = environment()))
    }

    return(ref)
}



# Main function: process RGset, perform CNV analysis, annotate requested genes,
# and generate an interactive genome plot with gene annotation.
annotate_cnv_with_genes <- function(my_idat, geneNames,
                                                 ref = NULL, min_probes = 4) {

    RGset <- minfi::read.metharray(my_idat, force = TRUE)

    if (RGset@annotation[["array"]] == "IlluminaHumanMethylationEPICv2") {
        RGset@annotation[["annotation"]] <- "20a1.hg38"
        RGset@annotation[["Genome"]] <- "hg38"
    } else{
        RGset@annotation[["Genome"]] <- "hg19"
    }

    array_info <- detect_array_info(RGset)
    message("Detected array type: ", array_info$array,
            " and genome build: ",array_info$genome)

    gene_gr <- get_gene_ranges(geneNames)

    library("conumee2")
    Mset <- minfi::preprocessIllumina(RGset, normalize = "controls")
    query_cnv <- conumee2::CNV.load(Mset)

    con_path <- file.path(path.package("conumee2"), "data")
    load(file.path(con_path, "exclude_regions.rda"))

    anno <- conumee2::CNV.create_anno(
        array_type = array_info$array,
        chrXY = TRUE,
        genome = array_info$genome,
        exclude_regions = exclude_regions,
        detail_regions = gene_gr
    )

    if (is.null(ref)) {
        ref <- get_reference(array_info)
    }

    aligned <- alignCNVObjects(query_cnv, ref, anno)
    colnames(aligned$ref@intensity) <- colnames(ref@intensity)

    cnv_obj <- get_cnv_obj(aligned)

    gene_anno_msg <- process_gene_ranges(gene_gr, anno, min_probes)
    p_interactive <- get_final_plot(cnv_obj, geneNames, gene_anno_msg, gene_gr)

    list(plot = p_interactive, cnv_analysis = cnv_obj)
}


# Example usage:
result <- annotate_cnv_with_genes(my_idat, geneNames)
result$plot
