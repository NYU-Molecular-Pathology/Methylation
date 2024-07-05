#!/usr/bin/env Rscript
## Script name: MakeIndelList.R
## Purpose: Filter out somatic variants for PACT consensus into a csv file
## Date Created: August 19, 2022
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024

library("base")
args <- commandArgs(TRUE)

# Parameters Input trailing command line --------------------------------------
args[1] -> pactRunName
args[2] -> concensusDir

DEFAULT_DIR <- "/Volumes/CBioinformatics/jonathan/pact/consensus/"
DEFAULT_OUT <- "/Volumes/molecular/MOLECULAR LAB ONLY/NYU PACT Patient Data/Results/Bioinformatics"

if (is.na(concensusDir)) {
    concensusDir <- paste0(DEFAULT_DIR, pactRunName, "_consensus")
}

deps <- c("Depends", "Imports", "LinkingTo")
if (!requireNamespace("stringr", quietly = TRUE)) {
    install.packages("stringr", dependencies = deps, ask = F, type = "binary")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr", dependencies = deps, ask = F, type = "binary")
}
library("stringr", warn.conflicts = F, quietly = T)
library("dplyr", warn.conflicts = F, quietly = T)

runYear <- stringr::str_split_fixed(pactRunName, "-", 3)[1, 2]
csvPath <- file.path(DEFAULT_OUT, paste0("20", runYear), pactRunName)

stopifnot(dir.exists(csvPath))

thisFile <- "somatic_variants.csv"
inputFile <- dir(path = csvPath, pattern = thisFile, full.names = T)

if (length(inputFile) == 0) {
    message("File not found:\n", thisFile)
    message("Ensure the file was copied here:\n", inputFile)
}

variantsData <- as.data.frame(read.csv(inputFile, strip.white = T))
varColumns <- c("Test_Number", "Gene.refGene", "Variant", "ExonicFunc.refGene")

exonicFilter <- stringr::str_detect(variantsData$ExonicFunc.refGene, "frame|delet|insert")
indelsList <- variantsData[exonicFilter, varColumns]
fix_genes <- stringr::str_replace_all(indelsList$Gene.refGene, ",", " ")
indelsList$Gene.refGene <- fix_genes

positions <- stringr::str_split_fixed(indelsList$Variant, ":", 3)[, 1:2]
indelsList$Position <- paste(positions[, 1], positions[, 2], sep = "_")
indelsList$nyu <- "Yes"
indelsList$philips <- ""
indelsList$Variant <- "SNV"
indelsList <- indelsList %>% dplyr::distinct(Position, .keep_all = T)

blank_row <- data.frame(
    Test_Number = "",
    Gene.refGene = "",
    Variant = "",
    ExonicFunc.refGene = "",
    Position = "",
    nyu = "",
    philips = ""
)

for (ngs in unique(variantsData$Test_Number)) {
    if (ngs %in% unique(indelsList$Test_Number) == F) {
        newSNVRow <- blank_row
        newSNVRow$Test_Number <- ngs
        newSNVRow$Gene.refGene <- "No Indels or FrameShifts called in-house"
        newSNVRow$Variant <- "SNV"
        indelsList <- rbind(indelsList, newSNVRow)
    }
    newCnvRow <- blank_row
    newCnvRow$Test_Number <- ngs
    newCnvRow$Position <- "NA"
    newCnvRow$Variant <- "CNV"
    indelsList <- rbind(indelsList, newCnvRow)
}

varsToCheck <- data.frame(
    "Test_Case" = indelsList$Test_Number,
    "Gene" = indelsList$Gene.refGene,
    "Mutation Type" = indelsList$ExonicFunc.refGene,
    "Other" = indelsList$Position,
    "In NYU" = indelsList$nyu,
    "In Philips" = indelsList$philips,
    "IGV" = '',
    "Comments" = '',
    "Variant" = indelsList$Variant
)

outPutFile <- file.path(concensusDir, paste0(pactRunName, "_desc.csv"))
write.csv(varsToCheck, file = outPutFile, quote = F, row.names = F)
message("Output file saved: '", outPutFile, "'")
