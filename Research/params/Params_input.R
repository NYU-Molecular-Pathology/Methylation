#!/usr/bin/env Rscript
## ---------------------------
## Script name: Params_input.R
## Purpose: Source global parameters and flags for research report inputs
## Date Created: May 17, 2022
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

gb <- globalenv(); assign("gb", gb)

# Input Variable Parameters ---------------------------------------------------
xlsxFile <- "myExcelFile.xlsx" # Your main input file
gb$yourName <- "Your_Name_Here" # author name for report
gb$yourDate <- format(Sys.Date(), '%B %d, %Y')  # date for report
token = "REDCAP_API_METH_TOKEN" # REDCap API token string
sampleType <- "Sample Type Tumors" # Name of sample types
htmlTitle <- paste("Methylation Clusters of", sampleType) # Title of html file
runDir <- "/Volumes/CBioinformatics/path/to/Project_folder"
inputFi <- file.path(runDir, xlsxFile)

# T-SNE plot titles (note: add newline (\n) if title is long) -----------------
titleMain <- paste("Clustering of", sampleType, ":\nBy Top")
samSheetDir <- file.path(runDir, "csv", "samplesheet.csv")
idatPath <- file.path(runDir, "idats") # Default path to idat files

# Hard coded Parameters -------------------------------------------------------
gb$runDir <- runDir <- sub("/$", "", runDir)
gb$setKnitDir(runDir)
stopifnot(getwd() == runDir)
stopifnot(file.exists(xlsxFile))

# Functional Parameter Flags --------------------------------------------------
F -> mergeProbes  # Do idat files need to be merged 450K + EPIC?
F -> batchEffect  # Are the idats processed in different facilities?
F -> getNoob      # Do want to perform NOOB minfi preprocessing?
T -> needFi       # Create a minfi samplesheet and copy idat files?
T -> addGenesHm   # Annotate the gene/probe names on heatmap y-axis?
T -> supervisedRun # Will this run have supervised sample clustering?
T -> genHeatMaps   # Should beta value HeatMap plots be output?
F -> genCNchunk    # Will CNV analysis be performed?
F -> genPathChunk  # Will Pathway analysis be performed?
F -> genPairwise   # Will Pathway analysis be performed?
gb$superHeat <- supervisedRun & genHeatMaps

# Input worksheet column names to subset samples ------------------------------
NULL -> col_batchEffect -> batch_col # Institution names or batch column name
NULL -> col_arrayType # "ArrayType" if the samplesheet is mixed EPIC and 450K
NULL -> variable_to_filter # is there any variable to cluster samples by?

# Sample Sheet Variable Columns -----------------------------------------------
"UNKNOWN" -> blank_keywd # Replacement for blank/empty values in samplesheet
"RD-Number" -> col_samNames # Name of the column containing unique sample IDs
col_samNames <- gsub("-", ".", col_samNames) # readxl drops dashes to dots
"Dx" -> col_samGrp -> col_samTypes # Column name for grouping sample COLORS
"Dx" -> col_shapes # Column name for plot point SHAPES\SYMBOLS
NULL -> col_Grouping # Column to Group samples by sub-type (optional)
NULL ->  col_sentrix # header name containing the sentrix ID if availible
NULL -> col_other # Any additional column(s) to colorize/label
gb$names2Label <- NULL  # Vector of specific samples to label within plots
selectedVars <- unique(c(col_shapes, col_samTypes)) # All variables to label

# HeatMaps Columns splitting --------------------------------------------------
colSplitHm <- 3 # Default number of columns to split or separate in heatmaps
hideTopAnno <- FALSE # Should the top annotation of the HeatMap be hidden?

# Site-Specific Analysis ------------------------------------------------------
genesInputFi <- NULL #file.path(runDir,"site-specific_genes.xlsx")
pathwayGrp_Filter <- NULL # Filter to select specific sample groups in pathways
totPaths <- 1:5 # Total Number of top pathways to list and display

# Pair-wise Analysis ----------------------------------------------------------
GroupSubset <- "Sample_Grouping"
nameGrp1 <- "Sample_A"
nameGrp2 <- "Sample_B"
nameGrp3 <- "Sample_C"

if (xlsxFile == "samplesheet.csv" | xlsxFile == "samplesheet_og.csv") {
    stop("Your input xlsxFile cannot be named: ", xlsxFile)
}

if (token == "REDCAP_API_METH_TOKEN") {
    stop("You did not provide a REDCap API token in Params_input.R") 
}
if (batchEffect == T & is.null(col_batchEffect)){
    stop("If there is a batch effect, col_batchEffect cannot be NULL")
}
