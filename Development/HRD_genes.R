
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("conumee")
BiocManager::install("minfi")
BiocManager::install("minfiDataEPIC")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("ChAMP")
BiocManager::install("minfiDataEPIC")

library(ChAMP)
library(minfiDataEPIC)
library(conumee)
library(minfi)
library(minfiDataEPIC)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Load the example data "RGsetEPIC"
data("RGsetEPIC")

# Define gene list and cytobands
genes <- c("BRCA1", "BRCA2", "ATM", "RAD51C")
gene_symbols <- org.Hs.eg.db::select(org.Hs.eg.db, keys=genes, columns=c("SYMBOL", "TXCHROM", "TXSTART", "TXEND"), keytype="SYMBOL")
cyto_bands <- TxDb.Hsapiens.UCSC.hg19.knownGene::get_cytoband(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Preprocess the data and calculate CNV values
MsetEPIC <- minfi::preprocessIllumina(RGsetEPIC, bg.correct = TRUE, normalize = "controls")
beta_values <- minfi::getBeta(MsetEPIC)


# Load the example data "RGsetEPIC"
data("RGsetEPIC")

# Preprocess the data using ChAMP package
myLoad <- champ.load(RGsetEPIC)
myNorm <- champ.norm(myLoad)

# Calculate Log R Ratio (LRR) and B Allele Frequency (BAF)
# Preprocess the data using preprocessNoob from minfi package
noob_data <- minfi::preprocessNoob(RGset)

# Custom function to calculate Log R Ratio (LRR)
calculateLRR <- function(Mset) {
  green <- minfi::getGreen(Mset)
  red <- minfi::getRed(Mset)
  lrr <- log2(green + red)
  return(lrr)
}

# Custom function to calculate B Allele Frequency (BAF)
calculateBAF <- function(Mset) {
  green <- minfi::getGreen(Mset)
  red <- minfi::getRed(Mset)
  baf <- green / (green + red)
  return(baf)
}

# Calculate LRR and BAF values
lrr <- calculateLRR(noob_data)
baf <- calculateBAF(noob_data)

# Extract the CNV values for the specific genes
cnv_values_list <- list()

for (gene in gene_symbols$SYMBOL) {
  gene_info <- gene_symbols[gene_symbols$SYMBOL == gene, ]
  chromosome <- as.numeric(gsub("chr", "", gene_info$TXCHROM))
  start <- gene_info$TXSTART
  end <- gene_info$TXEND
  
  gene_cnv_values <- cnv_data[chromosome, start:end, , drop = FALSE]
  cnv_values_list[[gene]] <- gene_cnv_values
}

output_file <- "cnv_values_genes.csv"
cnv_values_df <- do.call(cbind, cnv_values_list)
write.csv(cnv_values_df, output_file, row.names = FALSE)

#### Simplified

library(conumee)
library(minfi)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

genes <- c("BRCA1", "BRCA2", "ATM", "RAD51C")
gene_symbols <- org.Hs.eg.db::select(org.Hs.eg.db, keys = genes, columns = c("SYMBOL", "TXCHROM", "TXSTART", "TXEND"), keytype = "SYMBOL")
cyto_bands <- TxDb.Hsapiens.UCSC.hg19.knownGene::get_cytoband(TxDb.Hsapiens.UCSC.hg19.knownGene)

cnv_values_list <- lapply(gene_symbols$SYMBOL, function(gene) {
  gene_info <- gene_symbols[gene_symbols$SYMBOL == gene, ]
  chromosome <- as.numeric(gsub("chr", "", gene_info$TXCHROM))
  start <- gene_info$TXSTART
  end <- gene_info$TXEND
  cnv_data[chromosome, start:end, , drop = FALSE]
})

output_file <- "cnv_values_genes.csv"
cnv_values_df <- do.call(cbind, cnv_values_list)
write.csv(cnv_values_df, output_file, row.names = FALSE)


