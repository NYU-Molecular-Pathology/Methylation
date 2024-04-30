#!/usr/bin/env Rscript
## ---------------------------
## Script name: CustomCNVanno.R
## Purpose: Functions for taking a list of input genes and output custom CNV annotation for Conumee
## Date Created: August 10, 2023
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


drop_missing_genes <- function(myGenes, entrezIDs) {
    is_missing <- is.na(entrezIDs$ENTREZID)
    missing <- myGenes[is_missing]
    message(crayon::bgRed("Gene(s) not found in the database:"))
    message(paste(missing, collapse = "\n"))
    missing_df <- data.frame(missing = missing)
    write.csv(missing_df, file = "missing_genes.csv", quote = F, row.names = F)
    myGenes <- myGenes[!is_missing]
    entrezIDs <- AnnotationDbi::select(org.Hs.eg.db, keys = myGenes, keytype = "SYMBOL", columns = "ENTREZID")
    rownames(entrezIDs) <- 1:nrow(entrezIDs)
    return(entrezIDs)
}


grab_gene_info <- function(entrezIDs) {
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
    genesInfo <- subset(genesInfo, nchar(genesInfo$chromosome_name) <= 2)
    dupes <- which(base::duplicated(genesInfo$entrezgene_id))
    if (any(dupes)) {
        message("DUPLICATES:")
        message(paste(genesInfo[dupes,"entrezgene_id"],collapse = "\n"))
    }
    return(genesInfo)
}


grab_granges <- function(genesInfo, entrezIDs){
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


GetGeneRanges <- function(myGenes) {
    entrezIDs <- AnnotationDbi::select(org.Hs.eg.db, keys = myGenes, keytype = "SYMBOL", columns = "ENTREZID")
    
    if (any(is.na(entrezIDs$ENTREZID))) {
        entrezIDs <- drop_missing_genes(myGenes, entrezIDs)
    }

    genesInfo <- grab_gene_info(entrezIDs)
    granges_df <- grab_granges(genesInfo, entrezIDs)
    
   return(granges_df)
}


GetGenesListRange <- function(grangesDF, array_type = "EPIC") {
    mycoords.gr = GenomicRanges::makeGRangesFromDataFrame(grangesDF)
    grNames <- names(mycoords.gr)
    if (array_type == "EPICv2") {
        #genomeType <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
        genomeType <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    } else{
        genomeType <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    }
    gene_list <- suppressMessages(GenomicFeatures::genes(genomeType, single.strand.genes.only = T))[grNames,]
    entrezIDs <- gene_list@elementMetadata@listData[["gene_id"]]
    geneNames <- AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db, keys = entrezIDs, keytype = "ENTREZID", columns = "SYMBOL")
    gene_list@elementMetadata@listData[["name"]] = geneNames$SYMBOL
    gene_list@elementMetadata@listData[["gene_id"]] = NULL
    gene_list@elementMetadata@listData[["thick"]] = gene_list@ranges
    return(gene_list)
}


GetMnpAnno <- function(array_type, gb) {
    if (array_type == "EPIC") {
        require("mnp.v11b6")
        pkgPath <- file.path(path.package("mnp.v11b6"), "ext")
        annoFi <- "IlluminaArrayDBconumee_annotation_EPIC_B4.2017-06-07.RData"
    }
    if (array_type == "450k") {
        require("mnp.v11b4")
        pkgPath <- file.path(path.package("mnp.v11b4"), "ext")
        annoFi <- "IlluminaArrayDBconumee_annotation_EPIC_B4.2017-06-07.RData"
    }
    if (array_type == "EPICv2") {
        require("mnp.v12epicv2")
        pkgPath <- file.path(path.package("mnp.v12epicv2"), "ext")
        annoFi <- "conumee_annotation_EPIC_v2.2023-02-08.RData"
    }
    return(load(file.path(pkgPath, annoFi), envir = gb))
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
        library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
        probes450k <-
            sort(minfi::getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19))
    }
    if (is.element(array_type, c("EPIC", "overlap"))) {
        message("getting EPIC annotations")
        library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
        probesEPIC <-sort(minfi::getLocations(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
    }
    if (is.element(array_type, c("EPICv2", "overlap"))) {
        message("getting EPIC annotations")
        library("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
        probesEPIC <-sort(minfi::getLocations(IlluminaHumanMethylationEPICv2anno.20a1.hg38::IlluminaHumanMethylationEPICv2anno.20a1.hg38))
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


MakeCustomAnno2 <- function(gene_range_li, array_type, gb = .GlobalEnv){
    if (array_type == "EPIC") {
        custom_anno <- conumee2.0::CNV.create_anno(
            array_type = array_type, detail_regions = c(gene_range_li, gb$annoEPICxy@detail))
    }
    if (array_type == "EPICv2") {
        require("mnp.v12epicv2")
        pkgPath <- file.path(path.package("mnp.v12epicv2"), "ext")
        annoFi <- "conumee_annotation_EPIC_v2.2023-02-08.RData"
        load(file.path(pkgPath, annoFi), envir = gb)
        gb$annoEPICv2_xy@detail@seqinfo@genome <- rep("hg19", 24)
        custom_anno <- conumee2.0::CNV.create_anno(
            array_type = array_type, detail_regions = c(gene_range_li, gb$annoEPICv2_xy@detail))
    }
    if (array_type == "450k") {
        custom_anno <- conumee2.0::CNV.create_anno(
            array_type = array_type, detail_regions = c(gene_range_li, gb$anno450k@detail))
    }
    return(custom_anno)
}


LoadSaveAnno <- function(gene_range_li, array_type = "EPIC", customAnnoFi= NULL, gb = .GlobalEnv){
    if (is.null(customAnnoFi)) {
        customAnnoFi <- file.path(getwd(), paste0(array_type, "_custom_anno.Rdata"))
    }
    if (!file.exists(customAnnoFi)) {
        custom_anno <- MakeCustomAnno2(gene_range_li, array_type, gb)
        gb$SaveObj(custom_anno, customAnnoFi)
    } else{
        custom_anno <- gb$LoadRdatObj(customAnnoFi)
    }
    return(custom_anno)
}


generateSeg2 <- function(a, b, c) {
    x <- conumee2.0::CNV.fit(a, b, c)
    x <- conumee2.0::CNV.bin(x)
    x <- conumee2.0::CNV.detail(x)
    x <- conumee2.0::CNV.segment(x)
    return(x)
}


customCNV2 <- function(Mset, samName = NULL, sex = NULL, customAnno = NULL) {
    gb <- globalenv()
    if (is.null(samName)) {samName <- colnames(Mset)[1]}

    chiptype <- minfi::annotation(Mset)[[1]]
    Rset <- minfi::ratioConvert(Mset, what = "both", keepCN = TRUE)
    cndata <- conumee2.0::CNV.load(Mset)

    if (chiptype == "IlluminaHumanMethylationEPIC") {
        require("mnp.v11b6")
        if (is.null(sex)) {
            sex <- ifelse(mnp.v11b6::MNPgetSex(Rset)$predictedSex == "M", "Male", "Female")
        }
        path <- file.path(path.package("mnp.v11b6"), "ext")
        load(file.path(path,"conumee_annotation_EPIC_B6.2019-11-29.RData"), envir = gb)
        load(file.path(path,"CNanalysis6_conumee_REF_M.2018-09-19.RData"), envir = gb)
        load(file.path(path,"CNanalysis6_conumee_REF_F.2018-09-19.RData"), envir = gb)
        sexRefData <- if (sex == "Male") {gb$refM_epic} else{gb$refF_epic}
        mainAnno <- gb$annoEPICxy
    }
    if (chiptype == "IlluminaHumanMethylation450k") {
        require("mnp.v11b4")
        if (is.null(sex)) {
            sex <- ifelse(mnp.v11b4::MNPgetSex(Rset)$predictedSex == "M", "Male", "Female")
        }
        path <- file.path(path.package("mnp.v11b4"), "ext")
        load(file.path(path,"CNanalysis4_conumee_ANNO.vh20150715.RData"), envir = gb)
        load(file.path(path,"CNanalysis4_conumee_REF-M.vh20150715.RData"), envir = gb)
        load(file.path(path,"CNanalysis4_conumee_REF-F.vh20150715.RData"), envir = gb)
        sexRefData <- if (sex == "Male") {gb$refM.data} else{gb$refF.data}
        mainAnno <- gb$annoXY
    }

    if (chiptype == "IlluminaHumanMethylationEPICv2") {
        require("mnp.v12epicv2")
        if (is.null(sex)) {
            Rset@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
            sex <- ifelse(mnp.v12epicv2::MNPgetSex(Rset)$predictedSex[[1]] == "M", "Male", "Female")
        }
        path <- file.path(path.package("mnp.v12epicv2"), "ext")
        load(file.path(path,"conumee_annotation_EPIC_v2.2023-02-08.RData"), envir = gb)
        load(file.path(path,"CNanalysis6_conumee_REF_M.2018-09-19.RData"), envir = gb)
        load(file.path(path,"CNanalysis6_conumee_REF_F.2018-09-19.RData"), envir = gb)
        sexRefData <- if (sex == "Male") {gb$refM_epic} else{gb$gb$refF_epic}
        mainAnno <- gb$annoEPICv2_xy
        # Extracting all probe names
        sexRows <- rownames(sexRefData@intensity)
        mainProbes <- names(mainAnno@probes)
        customProbes <- names(customAnno@probes)
        cndataProbes <- names(cndata@intensity)

        # Finding common probe names across all three
        commonProbes <- intersect(intersect(sexRows, mainProbes), customProbes)
        commonProbes <- intersect(commonProbes, cndataProbes)

        # Filtering each object to keep only rows/columns with common probe names
        sexRefData@intensity <- sexRefData@intensity[commonProbes, ]
        mainAnno@probes <- mainAnno@probes[commonProbes, ]
        customAnno@probes <- customAnno@probes[commonProbes, ]
        cndata@intensity <- cndata@intensity[commonProbes, ]
    }
    
    if (!is.null(customAnno)) {
        return(generateSeg2(cndata, sexRefData, customAnno))
    } else{
        return(generateSeg2(cndata, sexRefData, mainAnno))
    }
}


get_cnv_data <- function(Mset, myGenes){
    library("conumee2.0")
    grangesDF <- gb$GetGeneRanges(myGenes)
    gene_range_li <- gb$GetGenesListRange(grangesDF)
    array_type <- switch(
        Mset@annotation[["array"]],
        IlluminaHumanMethylationEPICv2 = "EPICv2",
        IlluminaHumanMethylationEPIC = "EPIC",
        IlluminaHumanMethylation450k = "450k"
    )
    gb$GetMnpAnno(array_type, gb)
    custom_anno <- gb$LoadSaveAnno(gene_range_li, array_type)
    cnv_seg_data <- gb$customCNV2(Mset, samName = NULL, sex = NULL, customAnno = custom_anno)
    return(cnv_seg_data)
}

