#!/usr/bin/env Rscript
## Script name: CustomCNVanno.R
## Purpose: Functions for taking a list of input genes and output custom CNV annotation for Conumee
## Date Created: August 10, 2023
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025

gb <- globalenv()

my_idat <- "/path/to/myfiles"

pkgs <- c(
    "devtools",
    "conumee2",
    "sesame",
    "minfi",
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
    "IlluminaHumanMethylationEPICv2.anno.20a1.hg38",
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

if (not_installed("pak")) {
    install.packages("pak", ask = F, dependencies = T)
}


if (any(missing_pkgs)) {
    needs_install <- pkgs[missing_pkgs]
    pak::pkg_install(needs_install)
}

sesameData::sesameDataCacheAll()

stopifnot(all(sapply(pkgs, requireNamespace, quietly = TRUE)))

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
    return(gr)
}


# Helper function: for each gene, count overlapping probes in the CNV annotation's probe GRanges.
count_gene_probes <- function(gene_gr, probe_gr) {
    hits <- GenomicRanges::findOverlaps(gene_gr, probe_gr)
    hits <- as(hits, "Hits")  # convert to a plain Hits object
    length(unique(S4Vectors::subjectHits(hits)))
}



alignCNVObjects <- function(query_obj, ref_obj, anno_obj) {
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
    names(anno_obj@probes) <- anno_obj@probes@elementMetadata@listData[["IlmnID"]]
    common_probes <- base::Reduce(intersect, list(query_probes, ref_probes, anno_probes))
    if (length(common_probes) == 0)
        stop("No common probes found among query, reference, and annotation.")
    query_obj@intensity <- query_obj@intensity[common_probes, , drop = FALSE]
    ref_obj@intensity <- ref_obj@intensity[common_probes, , drop = FALSE]
    anno_obj@probes <- anno_obj@probes[common_probes]
    list(query = query_obj, ref = ref_obj, anno = anno_obj)
}


# Main function: process RGset, perform CNV analysis, annotate requested genes,
# and generate an interactive genome plot with gene annotation.
custom_cnv_plot_with_gene_annotation <- function(RGset, geneNames, ref = NULL, minProbeThreshold = 4) {
    # Detect array type and genome build.
    array_info <- detect_array_info(RGset)
    message("Detected array type: ", array_info$array, " and genome build: ", array_info$genome)
    library("conumee2")
    # Preprocess the queried RGset using minfi.
    Mset <- minfi::preprocessIllumina(RGset, normalize = "controls")
    query_cnv <- conumee2::CNV.load(Mset)

    # Load built-in exclude_regions and detail_regions.
    con_path <- file.path(path.package("conumee2"), "data")
    load(file.path(con_path, "exclude_regions.rda"))  # loads object 'exclude_regions'
    load(file.path(con_path, if (array_info$genome == "hg38") "detail_regions.hg38.rda" else "detail_regions.rda"))

    ref <- NULL
    anno <- conumee2::CNV.create_anno(
        bin_minprobes = 15,
        bin_minsize = 50000,
        bin_maxsize = 5e6,
        array_type = array_info$array,
        chrXY = TRUE,
        genome = array_info$genome,
        exclude_regions = exclude_regions,
        detail_regions = detail_regions
    )

    if (is.null(ref)) {
        if (array_info$array == "EPICv2") {
            message("No reference provided; attempting to load EPICv2 reference from sesameData.")
            available_refs <- sesameData::sesameDataList()
            if ("EPICv2.8.SigDF" %in% available_refs$Title) {
                ref_list <- sesameData::sesameDataGet("EPICv2.8.SigDF")
                ref_mat <- do.call(cbind, lapply(ref_list, totalIntensities))
                ref_df <- as.data.frame(ref_mat, check.names = FALSE)
                ref <- conumee2::CNV.load(ref_df)
            } else {
                stop("EPICv2.8.SigDF is not available. Please update sesameData and cache all data, or choose an alternative reference.")
            }
        } else {
            message("No reference provided; using default reference from minfiData ('MsetEx').")
            data("MsetEx", package = "minfiData", envir = environment())
            ref <- conumee2::CNV.load(get("MsetEx", envir = environment()))
        }
    }

    qProbs <- rownames(query_cnv@intensity)
    query_int_df <- as.data.frame(do.call(cbind, lapply(query_cnv@intensity, function(x) data.frame(x))))
    rownames(query_int_df) <- qProbs

    rProbs <- rownames(ref@intensity)
    ref_int_df   <- as.data.frame(do.call(cbind, lapply(ref@intensity, function(x) data.frame(x))))
    rownames(ref_int_df) <- rProbs

    query_cnv@intensity <- query_int_df
    ref@intensity   <- ref_int_df

    aligned <- alignCNVObjects(query_cnv, ref, anno)

    colnames(aligned$ref@intensity) <- names(ref_list)
    cnv_obj <- conumee2::CNV.fit(query = aligned$query, ref = aligned$ref, anno = aligned$anno)
    cnv_obj <- conumee2::CNV.bin(cnv_obj)
    cnv_obj <- conumee2::CNV.detail(cnv_obj)
    cnv_obj <- conumee2::CNV.segment(cnv_obj)

    # Grab the gene ranges if possible
    gene_gr <- get_gene_ranges(geneNames)

    if (is.null(gene_gr)) {
        warning("No gene ranges could be retrieved for the supplied gene names.")
        gene_annotation_message <- "No gene ranges found for input gene list."
    } else {
        probe_gr <- anno@probes
        gene_probe_counts <- sapply(split(gene_gr, gene_gr$symbol), count_gene_probes, probe_gr = probe_gr)

        gene_annotation_message <- ""
        for (gene in names(gene_probe_counts)) {
            nprobes <- gene_probe_counts[gene]
            if (nprobes == 0) {
                gene_annotation_message <- paste0(gene_annotation_message,
                                                  "No matches found for gene ", gene, " on the methylation array. ")
            } else if (nprobes < minProbeThreshold) {
                gene_annotation_message <- paste0(gene_annotation_message,
                                                  "Caution: only ", nprobes, " probes annotated for gene ", gene, ". ")
            }
        }
    }

    p_interactive <- conumee2::CNV.plotly(cnv_obj)

    if (!is.null(gene_gr)) {

        gene_df <- as.data.frame(gene_gr)
        gene_df$midpoint <- (gene_df$start + gene_df$end) / 2
        p <- p_interactive +
            ggplot2::geom_vline(data = gene_df,
                                ggplot2::aes(xintercept = midpoint),
                                color = "purple", linetype = "dashed", size = 0.5) +
            ggplot2::geom_text(data = gene_df,
                               ggplot2::aes(x = midpoint, y = Inf, label = symbol),
                               angle = 90, vjust = -0.5, hjust = 1, size = 3, color = "purple")
    }
    final_title <- paste("CNV Genome Plot With Genes –", paste(geneNames, collapse = ", "))
    if (nchar(gene_annotation_message) > 0) {
        final_title <- paste0(final_title, "\nNote: ", gene_annotation_message)
    }
    p_interactive <- plotly::layout(p_interactive, title = list(text = final_title))
    list(plot = p_interactive, cnv_analysis = cnv_obj)
}

# Example usage:
geneNames <- c("BRCA1", "FOXA1", "TP53")
RGset <- minfi::read.metharray(my_idat, force = T)
if (RGset@annotation[["array"]] == "IlluminaHumanMethylationEPICv2") {
    RGset@annotation[["annotation"]] <- "20a1.hg38"
    RGset@annotation[["Genome"]] <- "hg38"
} else{
    RGset@annotation[["Genome"]] <- "hg19"
}

result <- custom_cnv_plot_with_gene_annotation(RGset, myGenes)
result$plot  # Interactive plotly object
