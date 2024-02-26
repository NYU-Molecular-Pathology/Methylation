gb <- globalenv(); assign("gb", gb)
# Input File Names  -----------------------------------------------
xlsxFile <- "my_samples.xlsx" # Your main input file

# Input Variable Parameters
token = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" # REDCap API token string
sampleType <- "Tumors"
gb$yourName <- "Your Name"
htmlTitle <- paste("Methylation Clusters of", sampleType)

# Input Directory Names  -----------------------------------------------
runDir <- "/Users/Jonathan/Documents/Rprojects_local/Auto_test"
samsheet <- "samplesheet.csv" # primary sample sheet name/path
inputFi <- file.path(runDir, xlsxFile)
titleMain <- paste("Clustering of", sampleType, ":\nBy Top") # newline if title is too long (\n)
samSheetDir <- file.path(runDir,"csv","samplesheet.csv")
idatPath <- file.path(runDir,"idats")

# Hard coded Parameters -------------------------------------------------------------------------
gb$setKnitDir(runDir)
stopifnot(getwd() == runDir)
gb$CheckDirCreate("idats")

# Functional Parameter Flags --------------------------------------------------------------------
F -> mergeProbes  # Do idat files need to be merged 450K + EPIC?
F -> batchEffect  # Are the idats processed in different facilities?
F -> getNoob      # Do want to perform NOOB minfi preprocessing?
T -> needFi           # Do you need to create a minfi samplesheet and copy idat files?
F -> addGenesHm       # Do you want to annotate the gene/probe names on heatmap y-axis?
T -> supervisedRun    # Will this run have supervised sample clustering?
T -> generateCNVchunk # Will CNV analysis be performed?
F -> genPathwaychunk  # Will Pathway analysis be performed?
F -> genPairwise  # Will Pathway analysis be performed?

# Input worksheet column names to subset samples ------------------------------------------------
NULL ->   col_batchEffect -> batch_col # Institution names column for batched idats
NULL -> col_arrayType # "ArrayType"
NULL -> variable_to_filter

# Input Sheet Variable Columns
"UNKNOWN" -> blank_keywd    # Keyword to replace blank or empty values in the samplesheet
"record_id" -> col_samNames # Column containing unique sample IDS
"ethnicity" -> col_samGrp -> col_samTypes # Column name in samplesheet grouping sample type COLORS
"sex" -> col_shapes   # Column name for plot point SHAPES\SYMBOLS
NULL -> col_Grouping # Group samples by type i.e. t-sne cluster each subgroup/sample by type
"Sentrix_ID" ->  col_sentrix  # header name containing the sentrix ID
"diagnosis" -> col_other # any additional columns you would like to colorize/label in plots

selectedVars <- unique(c(col_shapes, col_samTypes)) # vector of names of columns you want to target
gb$names2Label <- NULL

# Site-Specific Analysis
genesInputFi <-  NULL #file.path(runDir,"site-specific_genes.xlsx")
pathwayGrp_Filter <- NULL
