#!/usr/bin/env Rscript
## ---------------------------
## Script name: CustomCNVanno.R
## Purpose: Functions for taking a list of input genes and output custom CNV annotation for Conumee
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------
gb <- globalenv(); assign("gb", gb)
options("install.packages.compile.from.source" = "No")
options("install.packages.check.source" = "no")

pkgs <- c(
    "dplyr",
    "conumee",
    "minfi",
    "org.Hs.eg.db",
    "GenomicRanges",
    "rtracklayer",
    "AnnotationHub",
    "clusterProfiler",
    "pathview",
    "GenomeInfoDb",
    "IRanges",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
    "Homo.sapiens",
    "GenomicRanges",
    "AnnotationDbi",
    "biomaRt",
    "ggplot2",
    "plotly",
    "tictoc",
    "tidyr"
)

unlist(lapply(pkgs, librarian::shelf))

library("mnp.v11b6")
stopifnot(all(sapply(pkgs, requireNamespace, quietly = TRUE)))

GetOverlaps <- function(mycoords.gr){
    return(suppressMessages(suppressWarnings(
        IRanges::subsetByOverlaps(GenomicRanges::genes(TxDb.Hsapiens.UCSC.hg19.knownGene, single.strand.genes.only = T), mycoords.gr)
    )))
}


GetGeneRanges <- function(myGenes) {
  entrezIDs <- AnnotationDbi::select(org.Hs.eg.db, keys = myGenes, keytype = "SYMBOL", columns = "ENTREZID")
  if (any(is.na(entrezIDs$ENTREZID))) {
      missing <- myGenes[is.na(match(myGenes, entrezIDs$SYMBOL))]
      stop("Gene(s) not found in the database: ",
           paste(missing, collapse = ", "))
  }
  mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  genesInfo <-
      biomaRt::getBM(
          attributes = c(
              "entrezgene_id",
              "chromosome_name",
              "start_position",
              "end_position"
          ),
          filters = "entrezgene_id",
          values = entrezIDs$ENTREZID,
          mart = mart
      )
  granges <- GenomicRanges::GRanges(
      seqnames = genesInfo$chromosome_name,
      ranges = IRanges::IRanges(
          start = genesInfo$start_position,
          end = genesInfo$end_position
      ),
      name = entrezIDs$SYMBOL
  )
  granges@elementMetadata@listData[["thick"]] = granges@ranges
  granges@ranges@NAMES <- paste(genesInfo$entrezgene_id)
  return(as.data.frame(granges))
}


GetGenesListRange <- function(grangesDF){
    mycoords.gr = GenomicRanges::makeGRangesFromDataFrame(grangesDF)
    grNames <- names(mycoords.gr)
    gene_list <- suppressMessages(GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene,single.strand.genes.only=T))[grNames,]
    entrezIDs <- gene_list@elementMetadata@listData[["gene_id"]]
    geneNames <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = entrezIDs, keytype = "ENTREZID", columns = "SYMBOL")
    gene_list@elementMetadata@listData[["name"]] = geneNames$SYMBOL
    gene_list@elementMetadata@listData[["gene_id"]] = NULL
    gene_list@elementMetadata@listData[["thick"]] = gene_list@ranges
    return(gene_list)
}


GetMnpAnno <- function (chiptype, gb) {
    if (chiptype == "EPIC") {
        require("mnp.v11b6")
        path <- file.path(path.package("mnp.v11b6"), "ext")
        return(
            load(file.path(path, "IlluminaArrayDBconumee_annotation_EPIC_B4.2017-06-07.RData"), envir = gb)
            )
    } else {
        require("mnp.v11b4")
        path <- file.path(path.package("mnp.v11b4"), "ext")
        return(
            load(file.path(path, "IlluminaArrayDBconumee_annotation_EPIC_B4.2017-06-07.RData", envir = gb))
            )
    }
}


MakeCustomAnno <- function (bin_minprobes = 15,
                                    bin_minsize = 50000,
                                    bin_maxsize = 5e+06,
                                    array_type = "EPIC",
                                    chrXY = T,
                                    exclude_regions = NULL,
                                    detail_regions = NULL) {
    object <- new("CNV.anno")
    object@date <- base::date()
    a1 <- formals()
    a2 <- as.list(match.call())[-1]
    getElt <- function(an) {
        if (is.element(an, names(a2))) {
            a2[[an]]
        }
        else {
            a1[[an]]
        }
    }
    object@args <- as.list(sapply(unique(names(c(
        a1, a2
    ))), getElt,
    simplify = F))
    if (is.null(array_type)) {
        array_type <- "450k"
    }
    if (!is.element(array_type, c("450k", "EPIC", "overlap"))) {
        stop("array_type must be on of 450k, EPIC, or overlap")
    } else {
        object@args[["bin_minsize"]] <- bin_minsize
        object@args[["array_type"]] <- array_type
        object@args[["chrXY"]] <- chrXY
    }
    if (chrXY) {
        object@genome <- data.frame(chr = paste("chr",
                                                c(1:22, "X", "Y"), sep = ""),
                                    stringsAsFactors = FALSE)
    } else {
        object@genome <- data.frame(chr = paste("chr", 1:22,
                                                sep = ""),
                                    stringsAsFactors = FALSE)
    }
    rownames(object@genome) <- object@genome$chr
    message("using genome annotations from UCSC")
    tbl.chromInfo <-
        conumee:::tbl_ucsc$chromInfo[match(object@genome$chr, conumee:::tbl_ucsc$chromInfo$chrom), "size"]
    object@genome$size <- tbl.chromInfo
    tbl.gap <-
        conumee:::tbl_ucsc$gap[is.element(conumee:::tbl_ucsc$gap$chrom, object@genome$chr),]
    object@gap <-
        sort(GRanges(
            as.vector(tbl.gap$chrom),
            IRanges(tbl.gap$chromStart +
                        1, tbl.gap$chromEnd),
            seqinfo = Seqinfo(object@genome$chr,
                              object@genome$size)
        ))
    tbl.cytoBand <-
        conumee:::tbl_ucsc$cytoBand[is.element(conumee:::tbl_ucsc$cytoBand$chrom,
                                               object@genome$chr),]
    pq <-
        sapply(split(tbl.cytoBand$chromEnd[grepl("p", tbl.cytoBand$name)],
                     as.vector(tbl.cytoBand$chrom[grepl("p", tbl.cytoBand$name)])),
               max)
    irangeObj <- IRanges(pq, pq)
    grangeObj <- GRanges(names(pq), irangeObj)
    object@genome$pq <-
        BiocGenerics::start(IRanges::resize(IRanges::subsetByOverlaps(object@gap, grangeObj), 1))

    probes450k <- probesEPIC <- GRanges()
    if (is.element(array_type, c("450k", "overlap"))) {
        message("getting 450k annotations")
        probes450k <-
            sort(minfi::getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19))
    }
    if (is.element(array_type, c("EPIC", "overlap"))) {
        message("getting EPIC annotations")
        probesEPIC <-sort(minfi::getLocations(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
    }
    if (array_type == "overlap") {
        probes <- sort(subsetByOverlaps(probes450k, probesEPIC))
    } else {
        probes <- unique(sort(c(probes450k, probesEPIC)))
    }
    probes <-
        probes[substr(names(probes), 1, 2) == "cg" &
                   is.element(as.vector(seqnames(probes)),
                              object@genome$chr)]
    object@probes <- sort(GRanges(as.vector(seqnames(probes)), ranges(probes), seqinfo = Seqinfo(object@genome$chr, object@genome$size)))
    message(" - ", length(object@probes), " probes used")
    if (!is.null(exclude_regions)) {
        message("importing regions to exclude from analysis")
        if (class(exclude_regions) == "GRanges") {
            object@exclude <- GRanges(
                as.vector(seqnames(exclude_regions)),
                ranges(exclude_regions),
                seqinfo = Seqinfo(object@genome$chr,
                                  object@genome$size)
            )
            values(object@exclude) <- values(exclude_regions)
            object@exclude <- sort(object@exclude)
        } else {
            object@exclude <- sort(rtracklayer::import(
                exclude_regions,
                seqinfo = Seqinfo(object@genome$chr, object@genome$size)
            ))
        }
    } else {
        object@exclude <- GRanges(seqinfo = Seqinfo(object@genome$chr, object@genome$size))
    }
    if (!is.null(detail_regions)) {
        message("importing regions for detailed analysis")
        if (class(detail_regions) == "GRanges") {
            object@detail <- GRanges(as.vector(seqnames(detail_regions)), ranges(detail_regions), seqinfo = Seqinfo(object@genome$chr, object@genome$size))
            if (any(grepl("name", names(values(detail_regions))))) {
                values(object@detail)$name <- values(detail_regions)[[grep("name", names(values(detail_regions)))[1]]]
            }
            if (any(grepl("IRanges", sapply(values(detail_regions), class)))) {
                values(object@detail)$thick <- values(detail_regions)[[grep("IRanges", sapply(values(detail_regions), class))[1]]]
            }
            object@detail <- sort(object@detail)
        }else {
            object@detail <- sort(rtracklayer::import(detail_regions, seqinfo = Seqinfo(object@genome$chr, object@genome$size)))
        }
        if (!is.element("name", names(values(object@detail)))) {
            stop("Detailed region bed file must contain name column")
        }
        if (!all(table(values(object@detail)$name) == 1)) {
            stop("Detailed region names must be unique")
        }
    } else {
        object@detail <- GRanges(seqinfo = Seqinfo(object@genome$chr, object@genome$size))
    }
    if (!is.element("thick", names(values(object@detail)))) {
        values(object@detail)$thick <- IRanges::resize(IRanges::ranges(object@detail), fix = "center", 1e+06)
    }
    message("creating bins")
    anno.tile <- conumee:::CNV.create_bins(hg19.anno = object@genome, bin_minsize = bin_minsize, hg19.gap = object@gap, hg19.exclude = object@exclude)
    message(" - ", length(anno.tile), " bins created")
    message("merging bins")
    object@bins <-  conumee:::CNV.merge_bins(
        hg19.anno = object@genome,
        hg19.tile = anno.tile,
        bin_minprobes = bin_minprobes,
        hg19.probes = object@probes,
        bin_maxsize = bin_maxsize
    )
    message(" - ", length(object@bins), " bins remaining")
    return(object)
}


merge_and_filter_genes <- function(input_directory, gene1 = "MTAP", gene2 = "CDKN2A/B") {
    csv_files <- grep("^genes_", list.files(input_directory, pattern = "\\.csv$"), value = T)

    if(!file.exists("all_gene_segments.csv")){
        all_gene_data <- do.call(rbind, lapply(
            file.path(input_directory, csv_files), read.csv, stringsAsFactors = F))
        write.csv(all_gene_data, file = "all_gene_segments.csv", row.names = FALSE)
    }else{
        all_gene_data <- as.data.frame(read.csv("all_gene_segments.csv"))
    }


    filtered_gene_data <- all_gene_data[all_gene_data$Gene %in% c(gene1, gene2), ]
    outFi = paste(gene1, gene2, "filtered.csv", sep = "_")
    outFi <- stringr::str_replace_all(outFi, .Platform$file.sep, "")
    write.csv(filtered_gene_data, file = outFi, row.names = FALSE)
    return(outFi)
}
