#!/usr/bin/env Rscript
## ---------------------------
## Script name: MakeIndelList.R
## Purpose: Filter out somatic variants for the PACT consensus report into a csv file
## Date Last Modified: January 16, 2024
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

library("base")
args <- commandArgs(TRUE)

# Parameters Input trailing commandline -------------------------------------------------------
args[1] -> pactRunName
args[2] -> concensusDir

if(is.na(concensusDir)){
    concensusDir <- paste0("/Volumes/CBioinformatics/jonathan/pact/consensus/", pactRunName, "_consensus")
}

deps <- c("Depends", "Imports", "LinkingTo")
if(!requireNamespace("stringr", quietly = TRUE)) {
    install.packages("stringr", dependencies = deps, ask = F)
}
if(!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr", dependencies = deps, ask = F)
}
library("dplyr", mask.ok = T, warn.conflicts = F, quietly = T)

runYear <- stringr::str_split_fixed(pactRunName, "-", 3)[1,2]
molecOut <- "/Volumes/molecular/MOLECULAR LAB ONLY/NYU PACT Patient Data/Results/Bioinformatics"
csvPath <- file.path(molecOut, paste0("20", runYear), pactRunName)

stopifnot(dir.exists(csvPath))

thisFile <- "somatic_variants.csv"
inputFile <- inputFile <- dir(path = csvPath, pattern = thisFile, full.names = T)

if(length(inputFile)==0){
    message("File not found: ", thisFile, "\n" , "Ensure the file was copied here:")
    message(inputFile)
}

variantsData <-as.data.frame(read.csv(inputFile, strip.white = T))

varColumns <- c("Test_Number", "Gene.refGene", "Variant", "ExonicFunc.refGene")
exonicFilter <- stringr::str_detect(variantsData$ExonicFunc.refGene, "frame|delet|insert")
indelsList <- variantsData[exonicFilter, varColumns]
indelsList$Gene.refGene <- stringr::str_replace_all(indelsList$Gene.refGene, ",", " ")

positions <- stringr::str_split_fixed(indelsList$Variant, ":", 3)[,1:2]
indelsList$Position <- paste(positions[,1], positions[,2], sep="_")
indelsList$nyu <- "Yes"
indelsList$philips <- ""
indelsList$Variant <- "SNV"
indelsList <- indelsList %>% dplyr::distinct(Position, .keep_all =T)

for(ngs in unique(variantsData$Test_Number)){
    if(!(ngs %in% unique(indelsList$Test_Number))){
        newSNVRow <- c(ngs, "No Indels or FrameShifts called in-house", "", "", "", "", "SNV")
        indelsList[nrow(indelsList) + 1,] = newSNVRow
    }
    newCnvRow <- c(ngs, "", "", "", "NA", "", "", "", "CNV")
    indelsList[nrow(indelsList) + 1,] = newCnvRow
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

#outPutFile <- file.path(fs::path_home(), "Desktop", paste0(pactRunName, "_desc.csv"))
outPutFile <- file.path(concensusDir, paste0(pactRunName, "_desc.csv"))
write.csv(varsToCheck, file=outPutFile, quote=F, row.names=F)
message("Output file saved: '", outPutFile, "'")
