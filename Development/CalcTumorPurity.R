# Calculate tumor purity from Illumina EPIC microarray .idat files in R, use the 'estimate' package.

# Install R packages, 'estimate' package for tumor purity estimation and the 'minfi' package for reading .idat files:
if (!requireNamespace("BiocManager", quietly = T)){install.packages("BiocManager")}
if (!requireNamespace("minfi", quietly = T)){BiocManager::install("minfi")}
if (!requireNamespace("minfi", quietly = T)){BiocManager::install("estimate")}
#2. Load the 'minfi' and 'estimate' packages:
library("minfi")
library("estimate")
# Read the .idat files using the 'minfi' package:
base_path <- "path/to/your/idat/files"
targets <- read.metharray.sheet(base_path)
RGset <- read.metharray.exp(targets = targets)
# Preprocess the data and create a methylation M-value matrix:
Mset <- preprocessQuantile(RGset, normalize = "controls")
Mval_matrix <- getM(Mset)
# Convert the M-value matrix to a beta-value matrix:
Beta_matrix <- minfi::m2beta(Mval_matrix)
# Estimate tumor purity using the 'estimate' package:
purity_results <- estimate_purity(Beta_matrix, platform = "IlluminaHumanMethylationEPIC")
#The 'purity_results' object contains the tumor purity estimates for each sample. Extract the tumor purity values using:
tumor_purity_values <- purity_results$Purity
