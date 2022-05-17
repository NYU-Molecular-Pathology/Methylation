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
gb$setKnitDir(runDir)

# Functional Plotting Parameter Flags   -----------------------------------------------
batchCorrect <- F # Do your samples need to be batch corrected?
needFi = T # if you need to create a minfi samplesheet and copy idats
geneNamesHeatMap <- F # Do you want to annotate the gene/probe names on heatmap y-axis?
supervisedRun <- F

# CVS sheet column names to subset samples  -----------------------------------------------
col_samNames <- "record_id" # column containing unique sample IDS
col_samTypes <- "subgroup" # groupings of sample kinds
col_samGrp <- col_samGroup <- "tumor_subgroup" # header for sample groups/subtypes
sam.grp.type <- "Organ"

# Additional CVS sheet column names -----------------------------------------------
col_sentrix <- "SentrixID_Pos" # header name containing the sentrix ID
col_other <- "Type" # any additional columns you would like to colorize/label in plots
selectedVars <- c("Type", "Tumor") # vector of names of columns you want to target
