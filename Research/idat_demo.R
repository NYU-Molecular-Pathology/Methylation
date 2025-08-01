#!/usr/bin/env Rscript
## Script name: idats_demo.R
## Purpose: Downloads public idat files to demonstrate minfi samplesheet
## Date Created: July 30, 2025
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025

## 1. Load/install packages
if (!"BiocManager" %in% rownames(installed.packages())) {
    install.packages("BiocManager", ask = FALSE, dependencies = TRUE)
}

if (!"minfi" %in% rownames(installed.packages())) {
    BiocManager::install("minfi", update = FALSE, ask = FALSE)
}
library("minfi")

PROJ_DIR <- "/Users/serraj10/Documents/conumee_test"

setwd(PROJ_DIR)
## 2. Download and unpack public EPIC IDATs for GSE42861
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE86831&format=file

ftp_url <-
    "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE86nnn/GSE86831/suppl/GSE86831_RAW.tar"
tar_file <- "GSE86831_RAW.tar"
options(timeout = 1200)
download.file(ftp_url, destfile = tar_file, mode = "wb", method = "libcurl")

# Decompress each file in place (removes .gz)
untar(tar_file, exdir = ".")
gz_files <- dir(getwd(), pattern = ".gz", full.names = TRUE)
lapply(gz_files, function(f) {R.utils::gunzip(f, overwrite = TRUE)})

## 3. Build sample sheet from file names
idat_files <- dir(getwd(), pattern = ".idat", full.names = TRUE)
basenames  <- unique(sub("_(Grn|Red)\\.idat$", "", idat_files, ignore.case = TRUE))
parts <- strsplit(basename(basenames), "_")
print(paste(parts, collapse = " "))

samples_df <- data.frame(
    Sample_Name = basename(basenames),
    Slide = sapply(parts, `[`, 1),
    Array = sapply(parts, `[`, 2),
    Basename    = basenames,
    stringsAsFactors = FALSE
)

print(samples_df)

write.csv(samples_df, file.path(getwd(), "samplesheet.csv"),
          row.names = FALSE, quote = FALSE)

## 4. Read into RGChannelSet
targets <- as.data.frame(read.csv("samplesheet.csv"))

rgSet <- minfi::read.metharray.exp(base = getwd(), targets = targets,
    verbose = TRUE, force = TRUE
)

rgSet
