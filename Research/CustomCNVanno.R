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

#library("mnp.v11b6")
stopifnot(all(sapply(pkgs, requireNamespace, quietly = TRUE)))

GetOverlaps <- function(mycoords.gr){
    return(suppressMessages(suppressWarnings(
        IRanges::subsetByOverlaps(GenomicRanges::genes(TxDb.Hsapiens.UCSC.hg19.knownGene, single.strand.genes.only = T), mycoords.gr)
    )))
}


v2Manifest <- "IlluminaHumanMethylationEPICv2manifest"
if (!requireNamespace(v2Manifest, quietly = T)) {
    devtools::install_github(file.path("mwsill", v2Manifest), upgrade = "always")
    devtools::install_github("mwsill/minfi", upgrade = "always")
}


grab_mset <- function(target1){
    RGSet <- minfi::read.metharray(target1$Basename, force = T)
    colnames(RGSet) <- target1$Sample_Name
    #RGSet@annotation[["array"]] <- "IlluminaHumanMethylationEPIC"
    Mset <- mnp.v12epicv2::MNPpreprocessIllumina(RGSet, normalize = "controls")
}


get_epic_manifest <- function(epic_vers){

    if (epic_vers == 1) {
        man_db <-
            IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19
    } else{
        man_db <-
            IlluminaHumanMethylationEPICv2anno.20a1.hg38::IlluminaHumanMethylationEPICv2anno.20a1.hg38
    }

    anno <- minfi::getAnnotation(man_db)
    man_df <- as.data.frame(anno@listData)
    row.names(man_df) <- anno@rownames
    return(man_df)
}

# Define a function to process gene names
process_genes <- function(gene_names) {
    if (length(gene_names) > 0) {
        split_names <- stringr::str_split(gene_names, ";")[[1]]
        gene_names <- paste(unique(split_names), collapse = ";")
    }
    return(gene_names)
}


get_gene_manifest <- function(epic_vers) {
    man_df <- get_epic_manifest(epic_vers)

    new_man_df <- man_df %>%
        dplyr::mutate(UCSC_RefGene_Name = purrr::map(UCSC_RefGene_Name, process_genes))

    if ("GencodeV41_Name" %in% colnames(man_df)) {
        new_man_df <- new_man_df %>%
            dplyr::mutate(GencodeV41_Name = purrr::map(GencodeV41_Name, process_genes))
    } else{
        new_man_df <- new_man_df %>%
            dplyr::mutate(GencodeCompV12_NAME = purrr::map(GencodeCompV12_NAME, process_genes))
    }

    new_man_df <- new_man_df %>% mutate_all(as.character)
    return(new_man_df)
}


filter_by_gene <- function(gene_names, new_man_df) {
    new_man_df %>%
        mutate(match = sapply(str_split(UCSC_RefGene_Name, ";"), function(genes_row)
            any(genes_row %in% gene_names))) %>%
        filter(match)
}


get_filtered_man <- function(manifest_out, myGenes, epic_vers = 1){
    if (!file.exists(manifest_out)) {
        new_man_df <- get_gene_manifest(epic_vers)
        filtered_df <- filter_by_gene(gene_names = myGenes, new_man_df)
        gb$SaveObj(filtered_df, manifest_out)
    } else{
        filtered_df <- gb$LoadRdatObj(manifest_out)
    }
    return(filtered_df)
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


drop_dupe_genes <- function(genesInfo, dupes){
    message("DUPLICATES:")
    the_extras <- genesInfo[dupes,"entrezgene_id"]
    extra <- which(genesInfo$entrezgene_id %in% the_extras)
    extra_df <- data.frame(genesInfo[extra, ])
    extra_msg <- paste(capture.output(extra_df), collapse = '\n')
    message(extra_msg)
    genesInfo <- genesInfo %>% group_by(entrezgene_id) %>% slice(n())
    rownames(genesInfo) <- 1:nrow(genesInfo)
    extra_msg <- paste(capture.output(data.frame(genesInfo[genesInfo$entrezgene_id %in% the_extras,])), collapse = '\n')
    message("Kept Ranges:\n", extra_msg)
    return(genesInfo)
}


grab_gene_info <- function(entrezIDs) {
    mart <- biomaRt::useMart(
        "ensembl", dataset = "hsapiens_gene_ensembl", version = "Ensembl Genes 111"
        )

    genesInfo <- biomaRt::getBM(
        attributes = c("entrezgene_id", "chromosome_name", "start_position", "end_position"),
        filters = "entrezgene_id", values = entrezIDs$ENTREZID, mart = mart
    )

    genesInfo <- subset(genesInfo, nchar(genesInfo$chromosome_name) <= 2)
    dupes <- base::duplicated(genesInfo$entrezgene_id)

    if (any(dupes)) {
        genesInfo <- drop_dupe_genes(genesInfo, dupes)
    }

    return(genesInfo)
}


get_anno_overlap <- function(granges){
    anno <- minfi::getAnnotation(
        IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19
    )
    anno_granges <- GRanges(
        seqnames = Rle(anno$chr),
        ranges = IRanges(start = anno$pos, end = anno$pos),
        strand = if ("strand" %in% colnames(anno)) Rle(anno$strand) else Rle("*")
    )
    
    seqlevelsStyle(granges) <- "UCSC"
    seqlevelsStyle(anno_granges) <- "UCSC"
    
    contained <- anno_granges %within% granges
    contained_positions <- anno_granges[contained]
    overlaps <- findOverlaps(granges, contained_positions)
    overlapping_indices <- queryHits(overlaps)
    granges <- granges[unique(overlapping_indices)]
    return(granges)
}


grab_granges <- function(genesInfo, entrezIDs){
    genesInfo$entrezgene_id <- as.character(genesInfo$entrezgene_id)
    entrezIDs$ENTREZID <- as.character(entrezIDs$ENTREZID)
    mergedData <- dplyr::inner_join(genesInfo, entrezIDs, by = c("entrezgene_id" = "ENTREZID"))
    
    granges <- GenomicRanges::GRanges(
        seqnames = Rle(as.character(mergedData$chromosome_name)),
        ranges = IRanges::IRanges(start = mergedData$start_position, end = mergedData$end_position),
        gene_id = mergedData$entrezgene_id,
        symbol = mergedData$SYMBOL
    )
    
    granges <- get_anno_overlap(granges)
    names(granges) <- granges$gene_id
    grangesDF <- as.data.frame(granges)

    return(grangesDF)
}


GetGeneRanges <- function(myGenes) {
    entrezIDs <- AnnotationDbi::select(org.Hs.eg.db, keys = myGenes, keytype = "SYMBOL", columns = "ENTREZID")
    
    if (any(is.na(entrezIDs$ENTREZID))) {
        entrezIDs <- drop_missing_genes(myGenes, entrezIDs)
    }

    genesInfo <- grab_gene_info(entrezIDs)
    grangesDF <- grab_granges(genesInfo, entrezIDs)
    return(grangesDF)
}


GetGenesListRange <- function(grangesDF, array_type = "EPIC") {
    mycoords.gr <- GenomicRanges::makeGRangesFromDataFrame(grangesDF, keep.extra.columns = T)
    grNames <- names(mycoords.gr)
    genomeType <- GenomicFeatures::makeTxDbFromUCSC(genome = "hg19")
    gene_list <- suppressMessages(GenomicFeatures::genes(genomeType, single.strand.genes.only = T))

    grNames <- grNames[grNames %in% names(gene_list)]
    gene_list <- gene_list[grNames, ]

    entrezIDs <- gene_list@elementMetadata@listData[["gene_id"]]
    geneNames <- AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db, keys = entrezIDs, keytype = "ENTREZID", columns = "SYMBOL"
    )

    gene_list@elementMetadata@listData[["name"]] = geneNames$SYMBOL
    gene_list@elementMetadata@listData[["thick"]] = gene_list@ranges
    seqlevelsStyle(gene_list) <- "UCSC"
    gene_range_li <- keepStandardChromosomes(gene_list, pruning.mode = "coarse")
    return(gene_range_li)
}


GetMnpAnno <- function(array_type, gb) {
    if (array_type == "EPIC") {
        #require("mnp.v11b6")
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

    chromNames <- paste(gene_range_li@seqnames@values)
    add_chrXY <-any(stringr::str_detect(chromNames, "chrX|chrY"))

    if (array_type == "EPIC") {
        custom_anno <- conumee2.0::CNV.create_anno(array_type = array_type,
                                                   detail_regions = gene_range_li,
                                                   chrXY = add_chrXY)
    }

    if (array_type == "EPICv2") {
        require("mnp.v12epicv2")
        pkgPath <- file.path(path.package("mnp.v12epicv2"), "ext")
        annoFi <- "conumee_annotation_EPIC_v2.2023-02-08.RData"
        load(file.path(pkgPath, annoFi), envir = gb)
        gb$annoEPICv2_xy@detail@seqinfo@genome <- rep("hg19", 24)

        custom_anno <- conumee2.0::CNV.create_anno(
            array_type = array_type,
            detail_regions = c(gene_range_li, gb$annoEPICv2_xy@detail),
            chrXY = add_chrXY
        )
    }

    if (array_type == "450k") {
        custom_anno <- conumee2.0::CNV.create_anno(
            array_type = array_type, detail_regions = c(gene_range_li, gb$anno450k@detail),
            chrXY = add_chrXY)
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
        #require("mnp.v11b6")
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

