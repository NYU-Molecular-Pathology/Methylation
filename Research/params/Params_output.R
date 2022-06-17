gb <- globalenv(); assign("gb", gb)
# Default Output file names  -----------------------------------------------
td <- format(Sys.Date(),"%b%d") # today's date Jun17
rgOut <- paste0(td,"_RGset.rds") # output name of the RG data set
rawBetaFi <- paste0(td,"_betaRaw.rds") # output name of the raw beta data set
mbfile <- paste0(td,"mSetSq.beta.rds")
unbetaVariance <- paste0(td,"unbetaVariance.rds")
supbetaOut <- paste0(td,"supervisedBetas.rds")
annotFi <- paste0(td,"annotations.rds")

# Segments Copy Number
segFile <- paste0(td, "_segfile.csv")
seg_clust_file <- paste0("data_seg_clusters_",td,".txt")

# Pathway Analysis Output
siteSpecific <- as.data.frame(readxl::read_excel(gb$genesInputFi))
your_genes <- siteSpecific$Gene
geneListIn <- "Genes_Pathway.csv"
pathCsvOut <- "signaling_pathway.csv"

# Plot Title Names  -----------------------------------------------
clusType <- c("Unsupervised", "Supervised")
if(gb$supervisedRun==F){clusType <- clusType[1]}
titleMain <- "Clustering of SAMPLENAME tumors:\nBy Type Top" # NOTE: newline if title is too long "\n"
varProbes <- c(100, 1000, 10000) # Which top Variance probes to pull i.e. c(100, 1000)
tsne_titles <- gb$generateTitles(clusType, topTitle=as.character(varProbes), titleMain)
topN = 1000
mdsTitle <- paste("Top", topN, "Common", "mSet Sq.beta", "MDS plot")