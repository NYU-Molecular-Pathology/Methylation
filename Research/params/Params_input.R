#!/usr/bin/env Rscript
## ---------------------------
## Script name: Params_input.R
## Purpose: Source global parameters and flags for research report inputs
## Date Created: May 17, 2022
## Date Last Modified: February 26, 2024
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

gb <- globalenv(); assign("gb", gb)

# Input File Names  -----------------------------------------------
xlsxFile <- "my_samples.xlsx" # Your main input file
samsheet <- "samplesheet.csv" # primary sample sheet name/path

# Input Variable Parameters -----------------------------------------------
gb$yourName <- "Your Name" # author name for report
token = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX" # REDCap API token string
sampleType <- "Tumors" # Name of sample types
htmlTitle <- paste("Methylation Clusters of", sampleType) # Title of html file

# Input Directory Names -----------------------------------------------
runDir <- "/Users/serraj10/Documents/Local_Rprojects/Auto_test"
inputFi <- file.path(runDir, xlsxFile)
titleMain <- paste("Clustering of", sampleType, ":\nBy Top") # newline (\n) if title is long
samSheetDir <- file.path(runDir, "csv", samsheet)
idatPath <- file.path(runDir, "idats") # Default path to idat files

# Hard coded Parameters -------------------------------------------------------------------------
gb$setKnitDir(runDir)
stopifnot(getwd() == runDir)


# Functional Parameter Flags --------------------------------------------------------------------
F -> mergeProbes  # Do idat files need to be merged 450K + EPIC?
F -> batchEffect  # Are the idats processed in different facilities?
F -> getNoob      # Do want to perform NOOB minfi preprocessing?
T -> needFi       # Do you need to create a minfi samplesheet and copy idat files?
F -> addGenesHm   # Do you want to annotate the gene/probe names on heatmap y-axis?
F -> supervisedRun # Will this run have supervised sample clustering?
F -> genHeatMaps   # Should beta value HeatMap plots be output?
F -> genCNchunk    # Will CNV analysis be performed?
F -> genPathChunk  # Will Pathway analysis be performed?
F -> genPairwise   # Will Pathway analysis be performed?

# Input worksheet column names to subset samples ------------------------------------------------
NULL ->   col_batchEffect -> batch_col # Institution names column for batched idats
NULL -> col_arrayType # "ArrayType"
NULL -> variable_to_filter

# Sample Sheet Variable Columns -----------------------------------------------------------------

# Keyword to replace any blank or empty values in the samplesheet
"UNKNOWN" -> blank_keywd

# Column containing unique sample IDS
"record_id" -> col_samNames

# Column name in samplesheet grouping sample COLORS
"ethnicity" -> col_samGrp -> col_samTypes

# Column name for plot point SHAPES\SYMBOLS
"sex" -> col_shapes

# Column to Group samples by type i.e. t-sne cluster each subgroup/sample by type
NULL -> col_Grouping

# header name containing the sentrix ID
"Sentrix_ID" ->  col_sentrix

# Any additional column(s) in the samplesheet you would like to colorize/label
"diagnosis" -> col_other

# Vector of names of columns to target coloring or shapes
selectedVars <- unique(c(col_shapes, col_samTypes))

# Vector of specific sample IDs to label within the plots
gb$names2Label <- NULL

# Site-Specific Analysis ------------------------------------------------------------------------
genesInputFi <-  NULL #file.path(runDir,"site-specific_genes.xlsx")
pathwayGrp_Filter <- NULL # Filter to select specific sample groups in pathways
totPaths <- 1:5 # Total Number of top pathways to list and display

# Pair-wise Analysis ----------------------------------------------------------------------------
GroupSubset <- "Sample_Grouping"
nameGrp1 <- "Sample_A"
nameGrp2 <- "Sample_B"
nameGrp3 <- "Sample_C"
