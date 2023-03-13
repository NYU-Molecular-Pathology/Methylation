gb <- globalenv(); assign("gb", gb)
# Input File Names  -----------------------------------------------
proj.name <- "my_project"
inputFi <- "my_samples.xlsx"
samsheet <- "samplesheet_output.csv" # primary sample sheet name/path
genesInputFi <-  "site-specific_genes.xlsx"

# Input Directory Names  -----------------------------------------------
userNam <- paste0(Sys.info()[["user"]])
runDir <- file.path("","Users", userNam, "Documents", proj.name)

token = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" # REDCap API token string
inputFi <- file.path(runDir, xlsxFile)
samsheet <- "samplesheet.csv" # primary sample sheet name/path
titleMain <- paste("Clustering of", sampleType, ":\nBy Top") # newline if title is too long (\n)
samSheetDir <- file.path(runDir,"csv","samplesheet.csv")
inputFi <- file.path(runDir, xlsxFile)
idatPath <- file.path(runDir,"idats")

# Hard coded Parameters -------------------------------------------------------------------------
gb$setKnitDir(runDir)
stopifnot(getwd() == runDir)
gb$CheckDirCreate("idats")

# Functional Parameter Flags --------------------------------------------------------------------
F -> mergeProbes  # Do idat files need to be merged 450K + EPIC?
T -> batchEffect  # Are the idats processed in different facilities?
F -> getNoob      # Do want to perform NOOB minfi preprocessing?
T -> needFi           # Do you need to create a minfi samplesheet and copy idat files?
F -> addGenesHm       # Do you want to annotate the gene/probe names on heatmap y-axis?
F -> supervisedRun    # Will this run have supervised sample clustering?
F -> generateCNVchunk # Will CNV analysis be performed?
F -> genPathwaychunk  # Will Pathway analysis be performed?
F -> pairWiseSamples  # will perform pairwise analysis
F -> genPairwise  # Will Pathway analysis be performed?

# Input worksheet column names to subset samples ------------------------------------------------
"Institution" ->   col_batchEffect -> batch_col # Institution names column for batched idats
NULL -> col_arrayType # "ArrayType"

# Variable Columns
"RD_Number" ->     col_samNames # Column containing unique sample IDS
"Methylation_Class" -> col_samTypes # Column name in samplesheet grouping sample type COLORS
"Driver" ->    col_shapes   # Column name for plot point shapes
"Driver" -> samGroup -> col_samGrp #header for sample groups/subtypes SHAPES
"All Tumors" ->   sam.grp.type # Grouping of tumors by type if relevant i.e. list specific names or new
"Sentrix_ID" ->  col_sentrix  # header name containing the sentrix ID
NULL -> col_other # any additional columns you would like to colorize/label in plots

selectedVars <- c(col_shapes, col_samTypes) # vector of names of columns you want to target
gb$names2Label <- NULL

# Site-Specific Analysis
genesInputFi <-  NULL #file.path(runDir,"site-specific_genes.xlsx")

gb$tsneHead1 <- gb$GetColorShape(gb$col_samGrp, gb$col_samTypes)
gb$tsneHead2 <- gb$GetColorShape(gb$col_samGrp, gb$sam.grp.type)
gb$tsneHead3 <- gb$GetColorShape(gb$col_samTypes, gb$sam.grp.type)
