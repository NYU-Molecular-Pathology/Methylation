
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("minfi")
BiocManager::install("FlowSorted.Blood.EPIC")

library(minfi)
library(FlowSorted.Blood.EPIC)

# Set the path to the folder containing the .idat files
idat_path <- "path/to/your/idat/files"

# Read .idat files using minfi package
RGset <- minfi::read.metharray.exp(targets = minfi::getAnnotation(idat_path))

# Preprocess the data using preprocessNoob from minfi package
noob_data <- minfi::preprocessNoob(RGset)

# Obtain beta values
beta_values <- minfi::getBeta(noob_data)

# Load the normal reference set
data("FlowSorted.Blood.EPIC")
normal_reference <- FlowSorted.Blood.EPIC

# Estimate cell counts using the Houseman algorithm
cell_counts <- minfi::estimateCellCounts(beta_values, reference = normal_reference)

# Calculate tumor purity
tumor_purity <- 1 - cell_counts[,"Bcell"]
