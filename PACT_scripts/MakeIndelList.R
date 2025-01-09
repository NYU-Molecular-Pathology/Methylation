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

somatic_file <- "somatic_variants.csv"
inputFile <- dir(path = csvPath, pattern = somatic_file, full.names = T)

if (length(inputFile) == 0) {
    message("File not found:\n", somatic_file)
    message("Ensure the file was copied here:\n", inputFile)
}

# Read the somatic variants input file
variantsData <- as.data.frame(read.csv(inputFile, strip.white = T))

# Filter out Strelka-only calls
filteredData <- variantsData %>%
    filter(!(Strelka == "Yes" & MuTect2 == "No" & LoFreqSomatic == "No"))

# Filter select columns and remove commas from rows
varColumns <- c(
    "Test_Number",
    "Tumor",
    "Normal",
    "Gene.refGene",
    "Variant",
    "DP",
    "AF",
    "MuTect2",
    "LoFreqSomatic",
    "ExonicFunc.refGene",
    "AAChange.refGene"
)
callsList <- filteredData[, varColumns]
fix_genes <- stringr::str_replace_all(callsList$Gene.refGene, ",", " ")
callsList$Gene.refGene <- fix_genes
fix_aa <- stringr::str_replace_all(callsList$AAChange.refGene, ",", "|")
callsList$AAChange.refGene <- fix_aa

# Paste gene positions separated by colons and create consensus column defaults
positions <- stringr::str_split_fixed(callsList$Variant, ":", 3)[, 1:2]
callsList$Position <- paste(positions[, 1], positions[, 2], sep = "_")
callsList$nyu <- "Yes"
callsList$philips <- ""
callsList$Variant <- "SNV"

# Arrange so that Mutect2 DP is used first when both callers are Yes
callsList <- callsList %>% arrange(desc(MuTect2 == "Yes")) %>% 
    distinct(Test_Number, Position, .keep_all = TRUE)

blank_row <- data.frame(
    Test_Number = "",
    Tumor = "",
    Normal = "",
    Gene.refGene = "",
    Variant = "",
    DP = "",
    AF = "",
    MuTect2 = "",
    LoFreqSomatic = "",
    ExonicFunc.refGene = "",
    AAChange.refGene = "",
    Position = "",
    nyu = "",
    philips = ""
)

# Iterate through NGS cases to ensure each is added if not in caller file
for (ngs in unique(variantsData$Test_Number)) {
    if (ngs %in% unique(callsList$Test_Number) == F) {
        newSNVRow <- blank_row
        newSNVRow$Test_Number <- ngs
        newSNVRow$Gene.refGene <- "No Indels or FrameShifts called in-house"
        newSNVRow$Variant <- "SNV"
        message("\nBinding new row:\n")
        message(paste0(capture.output(newSNVRow), collapse = "\n"))
        callsList <- rbind(callsList, newSNVRow)
    }
    newCnvRow <- blank_row
    newCnvRow$Test_Number <- ngs
    newCnvRow$Position <- "NA"
    newCnvRow$Variant <- "CNV"
    callsList <- rbind(callsList, newCnvRow)
}

# Generate final data frame with specific column names and column ordering
varsToCheck <- data.frame(
    "Test_Case" = callsList$Test_Number,
    "Tumor" = callsList$Tumor,
    "Normal" = callsList$Normal,
    "Gene" = callsList$Gene.refGene,
    "Mutation Type" = callsList$ExonicFunc.refGene,
    "Other" = callsList$Position,
    "In NYU" = callsList$nyu,
    "In Philips" = callsList$philips,
    "Depth" = callsList$DP,
    "AF" = callsList$AF,
    "MuTect2" = callsList$MuTect2,
    "LoFreqSomatic" = callsList$LoFreqSomatic,
    "IGV" = '',
    "Comments" = '',
    "AAChange" = callsList$AAChange.refGene,
    "Variant" = callsList$Variant
)

outPutFile <- file.path(concensusDir, paste0(pactRunName, "_desc.csv"))
write.csv(varsToCheck, file = outPutFile, quote = F, row.names = F)
message("Output file saved:\n", outPutFile)
