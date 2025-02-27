#!/usr/bin/env Rscript
## ---------------------------
## Script name: GetGeneRanges.R
## Purpose: source of global scripts imported for research methylation analysis
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------

LibLoad <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        suppressMessages(suppressWarnings(install.packages(pkg, dependencies = TRUE)))
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE, warn.conflicts = FALSE))
}

pkgs <- c("org.Hs.eg.db", "GenomicRanges", "AnnotationDbi", "biomaRt")
lapply(pkgs, LibLoad)
stopifnot(all(sapply(pkgs, requireNamespace, quietly = TRUE)))

GetGeneRanges <- function(myGenes) {
  entrezIDs <- AnnotationDbi::select(org.Hs.eg.db, keys = myGenes, keytype = "SYMBOL", columns = "ENTREZID")
  if (any(is.na(entrezIDs$ENTREZID))) {
      missing <- myGenes[is.na(match(myGenes, entrezIDs$SYMBOL))]
      stop("Gene(s) not found in the database: ",
           paste(missing, collapse = ", "))
  }
  mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  genes <-
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
      seqnames = genes$chromosome_name,
      ranges = IRanges::IRanges(
          start = genes$start_position,
          end = genes$end_position
      ),
      gene = entrezIDs$SYMBOL
  )
  return(granges)
}

#myGenes <- c("BRCA1", "BRCA2", "ATM", "RAD51C")
#granges <- GetGeneRanges(myGenes)
