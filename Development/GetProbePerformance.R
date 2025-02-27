#!/usr/bin/env Rscript
## ---------------------------
## Script name: GetProbePerformance.R
## Purpose: Input a list of genes and output list of optimally performing probe names
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------


# Load or Install required packages
LibLoad <- function(pkg){
    libOpt <- list(pkg, quietly = T, logical.return = T, warn.conflicts = F)
    suppressPackageStartupMessages(suppressWarnings(do.call(library, libOpt)))
}
CheckReq <- function(pkg){
    return(suppressWarnings(!require(pkg, character.only = T, warn.conflicts = F)))
}
BioCinst <- function(pkg) {
    if (CheckReq(pkg)) {do.call(BiocManager::install, list(pkg, update = F, ask = F))}
    return(LibLoad(pkg))
}
CheckCran <- function(pkg) {
    if(CheckReq(pkg)){do.call(install.packages, list(pkg, dependencies=T))}
    return(LibLoad(pkg))
}

pkgLis <- c("utils", "grDevices", "stringr", "BiocManager", "ggplot2", "BiocParallel", "reshape",
            "pals", "gridExtra", "fitdistrplus", "ggh4x", "dplyr", "purrr")
bioPkg <- c("minfi", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
pkgLoad <- unlist(lapply(pkgLis, CheckCran))
bioLoad <- unlist(lapply(bioPkg, BioCinst))
stopifnot(all(pkgLoad) & all(bioLoad))
# FUN: Returns unique BaseNames from idat file names
GetBaseNames <- function(filenames){unique(gsub("_(Red|Grn).idat(|\\.gz)$", "", filenames))}

# FUN: Generates the basenames in the input path
CreateBasenames <- function(path,recursive=FALSE) {
    grn.files <- list.files(path, pattern = "_Grn.idat(|\\.gz)$", recursive = recursive,
                            ignore.case = TRUE, full.names = TRUE)
    red.files <- list.files(path, pattern = "_Red.idat(|\\.gz)$", recursive = recursive,
                            ignore.case = TRUE, full.names = TRUE)
    GetBaseNames(c(grn.files, red.files))
}

# FUN: Gets the unique basenames from samples
NameFromBaseName <- function(basenames) {
    Sample_Name <- basename(basenames)
    if(any(duplicated(Sample_Name)))
    {
        warning("Some duplicated Sample_Name entries")
        Sample_Name <- make.unique(Sample_Name)
    }
    return(Sample_Name)
}

ReadMethParallel <- function(targets, verbose = TRUE, ...) {
    nworkers <- BiocParallel::bpworkers(BiocParallel::bpparam())
    if (nworkers <= 1) {
        stop("Did you registered a biocparallel back-end?")
    }
    y <- rep(1, ceiling(nrow(targets)/nworkers))
    for (i in 2:nworkers) {
        y <- c(y, rep(i, ceiling(nrow(targets) / nworkers)))
    }
    y <- y[1:nrow(targets)]
    jobs <- split(targets, y)

    fun <- function(x, ...) {
        ## Must be loaded on the worker nodes explicitly for BatchJobs!
        requireNamespace("minfi")
        requireNamespace("Biobase")
        read.metharray.exp(targets = x, ...)
    }

    message("Reading multiple idat-files in parallel...")
    res <- bplapply(jobs, FUN = fun, ...)
    if(verbose) {
        message(str(res))
    }
    message("Combining the RGsets to one big RGset...")
    rgSet <- res[[1]]

    for (i in 2:length(res)) {
        rgSet <- BiocGenerics::combine(rgSet, res[[i]])
    }

    message("DONE!")
    return(rgSet)
}


# FUN: Create sample sheet if an Illumina one isn't available
CreateSheetFromIdats <- function(idatPath, recursive = F, delim = "_") {
    pathEnd = substr(idatPath, nchar(idatPath), nchar(idatPath))
    basenames=CreateBasenames(path, recursive)
    if (length(basenames) == 0) {
        warning("No idat files found.")
        return(NULL)
    }
    dat <- data.frame(do.call(rbind, strsplit(basename(basenames), split=delim)), stringsAsFactors=F)
    if (ncol(dat) < 2) {
        warning(path," basenames do not appear to correspond to idat files")
    }
    sentrixpos <- grep("^R[0-9][0-9]C[0-9][0-9]$", as.character(unlist(dat[1,])))
    if(length(sentrixpos)==1) {
        temp <- do.call(rbind, strsplit(as.character(dat[,sentrixpos]), split="C"))
        dat$sentrix_row <- gsub("R", "", temp[,1])
        dat$sentrix_col <- temp[,2]
        dat <- dat[,-sentrixpos]
    }

    slidecol <- grep("^[0-9]{9}[0-9]*$", as.character(unlist(dat[1,])))
    if (length(slidecol) == 1) {
        colnames(dat)[slidecol] <- "Slide"
    }

    idcol <- which(apply(dat, 2, function(x) all(!duplicated(x))))
    if(length(idcol) >= 1) {
        Sample_Name <- dat[,idcol[1]]
        dat <- dat[,-idcol[1],drop=F]
    }
    else {
        Sample_Name <- NameFromBaseName(basenames)
    }

    samplesheet <- data.frame(Sample_Name = Sample_Name, Sex = NA, dat, Basename = basenames, stringsAsFactors=F)
    write.csv(samplesheet, file="samplesheet.csv", quote=F, row.names=F)
    return(samplesheet)
}


# FUN: reads .idat files and preprocess data
ReadProcessData <- function(idatPath) {
    targets <- CreateSheetFromIdats(idatPath)
    BiocParallel::register(BatchtoolsParam(workers = 8), default = T)
    RGset <- ReadMethParallel(targets, force = T)
    Mset <- minfi::preprocessQuantile(RGset)
    return(Mset)
}


GrabSnpProbes <- function(anno){
    snpCols <- c("Probe_rs" , "Probe_maf", "CpG_rs", "CpG_maf" , "SBE_rs" , "SBE_maf")
    snpAnno <- as.data.frame(anno@listData[snpCols])
    rownames(snpAnno) <- rownames(anno)
    snpAnno <- snpAnno %>% dplyr::filter_all(any_vars(!is.na(.)))
    snpProbes <- rownames(snpAnno)
    return(snpProbes)
}


GetIslandProbes <- function(anno, geneName, minProbes = 4) {
    gene_probes <- anno[grepl(anno$UCSC_RefGene_Name, pattern = geneName),]
    snpProbes <- GrabSnpProbes(anno)
    toKeep <- which(!rownames(gene_probes) %in% snpProbes)
    gene_probes <- gene_probes[toKeep,]
    #islandRegions <- gene_probes$Relation_to_Island == "Island"
    island_probes <- gene_probes %>% as.data.frame() %>%
        mutate(UCSC_RefGene_Name = str_split(UCSC_RefGene_Name, ";", simplify = F) %>%
                   map( ~ sort(unique(.x)))) %>%
        filter(map_lgl(UCSC_RefGene_Name, ~ any(.x == geneName))) %>%
        mutate(UCSC_RefGene_Name = map_chr(UCSC_RefGene_Name, ~ paste(.x, collapse = ";")))
    message("")
    message(paste(geneName, "probes found:"))
    island_probes <-
        island_probes[, c("chr", "pos", "UCSC_RefGene_Name", "Relation_to_Island", "UCSC_RefGene_Group")]
    csvFiName <- paste0(geneName, "_probes_used.csv")
    csvPath <- file.path(getwd(), "ProbesCSV")
    if (!dir.exists(csvPath)) {
        dir.create(csvPath)
    }
    message(paste0(capture.output(island_probes), collapse = "\n"))
    write.csv(
        island_probes,
        file = file.path(csvPath, csvFiName),
        quote = F,
        row.names = T
    )
    return(island_probes)
}


# FUN: identifies the best probes for genes input
FindBestProbes <- function(Mset, gene_list, threshold = 0.05) {
    anno <- getAnnotation(Mset)
    gene_anno <- NULL
    for (geneName in gene_list) {
        probeList <- GetIslandProbes(anno, geneName, minProbes = 4)
        if(is.null(gene_anno)){
            gene_anno <- probeList
        }else{
            gene_anno <- rbind(gene_anno, probeList)
        }
    }
    gene_probes <- rownames(gene_anno)
    beta_values <- getBeta(Mset)
    beta_subset <- beta_values[gene_probes,]
    signif_probes <- which(apply(beta_subset, 2, function(x) any(x > threshold)))
    best_probes <- gene_probes[signif_probes]
    write.csv(gene_anno[best_probes,], "BestCandidateProbes.csv", quote=F)
    return(best_probes)
}

# Main function
main <- function(gene_list, idatPath=NULL) {
    # Provide path to the folder containing the .idat files
    if(is.null(idatPath)){idatPath <- getwd()}
    Mset <- ReadProcessData(idatPath)
    # Define genes of interest
    best_probes <- FindBestProbes(Mset, gene_list)
    print(paste("The best candidate probes for the genes of interest are:", paste(best_probes, collapse = ", ")))
}

gene_list <- c("BRCA1", "BRCA2", "ATM", "RAD51C")
idatPath <- "/Volumes/CBioinformatics/Methylation/Clinical_Runs/Probe_performance/"
setwd(idatPath)
main(gene_list, idatPath)
