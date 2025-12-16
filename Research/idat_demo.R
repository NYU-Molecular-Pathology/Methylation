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
    BiocManager::install("minfi", update = FALSE, ask = FALSE, dependencies = TRUE)
}
if (!"ComplexHeatmap" %in% rownames(installed.packages())) {
    BiocManager::install("ComplexHeatmap", update = FALSE, ask = FALSE, dependencies = TRUE)
}
if (!"matrixStats" %in% rownames(installed.packages())) {
    BiocManager::install("matrixStats", update = FALSE, ask = FALSE, dependencies = TRUE)
}
if (!"grDevices" %in% rownames(installed.packages())) {
    install.packages("grDevices", ask = FALSE, dependencies = TRUE)
}
library("minfi")
library("ComplexHeatmap")
library("matrixStats")
library("grDevices")

PROJ_DIR <- file.path(Sys.getenv("HOME"), "Documents", "conumee_test")
SAM_SHEET <- file.path(PROJ_DIR, "samplesheet.csv")
if (!dir.exists(PROJ_DIR)) dir.create(PROJ_DIR, recursive = TRUE)
setwd(PROJ_DIR)
## 2. Download and unpack public EPIC IDATs for GSE42861
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE86831&format=file

ftp_url <-
    "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE86nnn/GSE86831/suppl/GSE86831_RAW.tar"
geo_file <-
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE86831&format=file"
tar_file <- "GSE86831_RAW.tar"
options(timeout = 1200)
utils::download.file(geo_file, destfile = tar_file, mode = "wb")

# Decompress each file in place (removes .gz)
untar(tar_file, exdir = ".")
gz_files <- dir(getwd(), pattern = ".gz", full.names = TRUE)
lapply(gz_files, function(f) {R.utils::gunzip(f, overwrite = TRUE)})

## 3. Build sample sheet from file names
make_minfi_samplesheet <- function(idats_folder = NULL) {
    idats_folder <- ifelse(is.null(idats_folder), getwd(), idats_folder)
    idat_files <- dir(idats_folder, pattern = ".idat", full.names = TRUE)
    basenames  <- unique(sub("_(Grn|Red)\\.idat$", "",
                             idat_files, ignore.case = TRUE))
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

    write.csv(samples_df, file = SAM_SHEET, row.names = FALSE, quote = FALSE)
}

## 4. Read into RGChannelSet
# NOTE: You can add columns and anotations here to targets
make_minfi_samplesheet()
targets <- as.data.frame(read.csv(SAM_SHEET))

rgSet <- minfi::read.metharray.exp(base = getwd(), targets = targets,
                                   verbose = TRUE, force = TRUE
)

## ------------------------------------------------------------------
## 1) Compute detection P-values (probe x sample)
## ------------------------------------------------------------------
detP <- minfi::detectionP(rgSet)

## Align detP sample names to your targets (assume first column = sample IDs)
sampleNameCol <- colnames(targets)[1]
message("Sample Name Column is: ", sampleNameCol)
sample_ids <- rgSet@colData[[sampleNameCol]]
colnames(detP) <- sample_ids

## ------------------------------------------------------------------
## 2) Filter bad samples (mean detection P too high)
##    Typical thresholds: 0.05 (lenient) to 0.01 (stricter)
## ------------------------------------------------------------------
samplePvalThreshold <- 0.05
goodSamples <- colMeans(detP, na.rm = TRUE) < samplePvalThreshold

rgSet   <- rgSet[, goodSamples]
targets <- targets[goodSamples, , drop = FALSE]
detP    <- detP[, goodSamples, drop = FALSE]

## Refresh sample_ids after filtering
sample_ids <- rgSet@colData[[sampleNameCol]]
colnames(detP) <- sample_ids

## ------------------------------------------------------------------
## 3) Normalize / preprocess RGChannelSet -> (Genomic)MethylSet
##    Choose ONE method; common choices: preprocessNoob, preprocessFunnorm, preprocessQuantile
## ------------------------------------------------------------------
normMethod <- "noob"  # change to "funnorm" or "quantile" if desired

mSet <- if (identical(normMethod, "noob")) {
    suppressWarnings(minfi::preprocessNoob(rgSet))
} else if (identical(normMethod, "funnorm")) {
    suppressWarnings(minfi::preprocessFunnorm(rgSet))
} else if (identical(normMethod, "quantile")) {
    suppressWarnings(minfi::preprocessQuantile(rgSet))
} else {
    stop("normMethod must be one of: 'noob', 'funnorm', 'quantile'")
}

## ------------------------------------------------------------------
## 4) Filter poor-quality probes (require probe detection P < probePval in ALL samples)
## ------------------------------------------------------------------
probePval <- 0.01

## Ensure detP rows match the featureNames of the methylation object
detP_aligned <- detP[match(minfi::featureNames(mSet), rownames(detP)), , drop = FALSE]
keepProbes <- rowSums(detP_aligned < probePval, na.rm = FALSE) == ncol(mSet)
mSetClean <- mSet[keepProbes, ]

## ------------------------------------------------------------------
## 5) Add SNP annotation + optionally drop SNP-affected probes and sex-chr probes
## ------------------------------------------------------------------
removeSexSnps <- TRUE
if (isTRUE(removeSexSnps)) {
    ## Add SNP info first
    mSetClean <- minfi::addSnpInfo(mSetClean)

    ## Drop loci with SNPs at SBE and CpG (maf=0 removes any annotated SNP overlap)
    mSetClean <- minfi::dropLociWithSnps(mSetClean, snps = c("SBE", "CpG"), maf = 0)

    ## Drop chrX/chrY probes
    annot <- minfi::getAnnotation(mSetClean)
    sexProbes <- annot$Name[annot$chr %in% c("chrX", "chrY")]
    mSetClean <- mSetClean[!(rownames(mSetClean) %in% sexProbes), ]
}

## ------------------------------------------------------------------
## 6) Extract Beta values
## ------------------------------------------------------------------
betas <- minfi::getBeta(mSetClean)
colnames(betas) <- sample_ids

## betas is now a matrix: probes x samples
## Optional: quick sanity checks
stopifnot(is.matrix(betas))
stopifnot(all(is.finite(betas)))
stopifnot(all(betas >= 0 & betas <= 1))
dim(betas)
head(betas[, 1:min(3, ncol(betas)), drop = FALSE])

common_ids <- intersect(colnames(betas), targets[[sampleNameCol]])
if (length(common_ids) < 2L) {
    stop("Insufficient overlap between colnames(betas) and targets[[sampleNameCol]].")
}

betas <- betas[, common_ids, drop = FALSE]
targets <- targets[match(common_ids, targets[[sampleNameCol]]), , drop = FALSE]
stopifnot(identical(colnames(betas), targets[[sampleNameCol]]))

## ------------------------------------------------------------------
## 1) Adding in DEMO columns in `targets`
## (These are synthetic example annotations; replace with REAL metadata when available.)
## ------------------------------------------------------------------
n <- nrow(targets)

if (!("Group" %in% colnames(targets))) {
    ## Example: alternate Tumor/Normal
    targets$Group <- factor(rep(c("Tumor", "Normal"), length.out = n), levels = c("Tumor", "Normal"))
}

if (!("Batch" %in% colnames(targets))) {
    ## Example: 3 batches
    targets$Batch <- factor(paste0("Batch", ((seq_len(n) - 1L) %% 3L) + 1L))
}

if (!("Sex" %in% colnames(targets))) {
    ## Example: Made up F/M sample labels
    targets$Sex <- factor(sample(rep(c("F", "M"), length.out = n)),levels = c("F", "M"))
}

if (!("Age" %in% colnames(targets))) {
    ## Example: deterministic ages (useful for when you have a continuous annotation)
    targets$Age <- 30 + ((seq_len(n) - 1L) %% 25L)
}

## ------------------------------------------------------------------
## 2) Subset probes for visualization (beta matrices can be ~850k rows)
##    Use top variable probes across samples (fast and interpretable).
## ------------------------------------------------------------------
max_probes_to_plot <- 2000L # You can change this value to be the Top N variance probes
probe_var <- matrixStats::rowVars(betas, na.rm = TRUE)
ord <- order(probe_var, decreasing = TRUE)
keep_idx <- ord[seq_len(min(max_probes_to_plot, length(ord)))]
betas_plot <- betas[keep_idx, , drop = FALSE]

## Optional: center each probe for heatmap contrast (keeps values in roughly [-0.5, 0.5])
betas_plot_centered <- betas_plot - rowMeans(betas_plot, na.rm = TRUE)

## ------------------------------------------------------------------
## 3) Build ComplexHeatmap annotations from `targets`
## ------------------------------------------------------------------

## Create discrete age groups from the continuous 'Age' variable
targets$Age_Group <- cut(targets$Age,
                         breaks = c(29, 40, 50, 60),
                         labels = c("30-39", "40-49", "50-55"),
                         right = TRUE)


## IMPORTANT: keys in `anno_colors` must match the *annotation names* used below.
## In your code, if you used `Age = targets$Age_Group`, so the key must be "Age" (not "Age_Group"),
## unless you rename the annotation to Age_Group. https://r-charts.com/colors/
# Note: you can google hex color picker or install Rstudio hex color viewer
anno_colors <- list(
    Group = c("Tumor" = "black", "Normal" = "yellow"),
    Sex   = c("F" = "pink", "M" = "royalblue"),
    Age_Group   = c("30-39" = "#76c3e6", "40-49" = "#fc7808", "50-55" = "#a1111a")
)

## Optional: explicitly map Batch colors if you want "Batch" column show up
## (otherwise ComplexHeatmap will auto-assign; this is often what appears random colors)
if ("Batch" %in% colnames(targets)) {
    batch_levels <- levels(targets$Batch)
    batch_cols <- structure(grDevices::hcl.colors(length(batch_levels), "Dark 3"), names = batch_levels)
    anno_colors$Batch <- batch_cols
}

## ------------------------------------------------------------------
## 2A) Keep Batch annotation (with coloring)
## ------------------------------------------------------------------
col_ha_with_batch <- ComplexHeatmap::HeatmapAnnotation(
    Group = targets$Group,
    Batch = targets$Batch,
    Sex   = targets$Sex,
    Age   = targets$Age_Group,  # annotation name is "Age" notice colors do not match
    col   = anno_colors,
    annotation_name_side = "left"
)

ht1 <- ComplexHeatmap::Heatmap(
    betas_plot_centered,
    name = "Beta",
    top_annotation = col_ha_with_batch,
    column_split = targets$Group,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_rot = 45,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    use_raster = TRUE
)

ComplexHeatmap::draw(ht1, heatmap_legend_side = "right", annotation_legend_side = "right")

## ------------------------------------------------------------------
## 2B) Omit Batch annotation entirely (remove it from annotation + color list)
## ------------------------------------------------------------------
anno_colors_no_batch <- anno_colors
anno_colors_no_batch$Batch <- NULL

col_ha_no_batch <- ComplexHeatmap::HeatmapAnnotation(
    Group = targets$Group,
    Sex   = targets$Sex,
    Age_Group   = targets$Age_Group, # annotation name is "Age_Group" notice colors now match
    col   = anno_colors_no_batch,
    annotation_name_side = "left"
)

ht2 <- ComplexHeatmap::Heatmap(
    betas_plot_centered,
    name = "Beta",
    top_annotation = col_ha_no_batch,
    column_split = targets$Group,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_rot = 45,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    use_raster = TRUE
)

ComplexHeatmap::draw(ht2, heatmap_legend_side = "right", annotation_legend_side = "right")
