gb <- globalenv(); assign("gb", gb)
# Default Output file names  -----------------------------------------------
td <- format(Sys.Date(),"%b%d") # today's date Jun17
rgOut <- paste0(td, "_RGset.Rdata")
rawBetaFi <- paste0(td, "_betaRaw.Rdata")
pValsOutFi <- paste0(td, "_DetPvals.Rdata")
mbfile <- paste0(td, "_mSetSq.beta.Rdata")
unbetaVariance <- paste0(td, "_unbetaVariance.RData")
combatOut <- paste0(td, "_combatBetas.Rdata")

supbetaOut <- paste0(td, "_supervisedBetas")
annotFi <- paste0(td, "_annotations.rds")

# OutPut Directories
gb$CheckDirCreate("data")
gb$CheckDirCreate("figures")
gb$CheckDirCreate("figures/mds")
gb$CheckDirCreate("csv")
gb$CheckDirCreate("figures/tsne/")
gb$CheckDirCreate("figures/heatmaps/")
if(gb$generateCNVchunk==T){
    gb$CheckDirCreate("figures/cnv/")
}
if(gb$genPathwaychunk==T){
    gb$CheckDirCreate("figures/pathway/")
    gb$CheckDirCreate("figures/cluster/")
}
if(gb$genPairwise==T){
    gb$CheckDirCreate("/figures/diffmean/")
}

# Segments Copy Number Output Files
segFile <- paste0(td, "_segfile.csv")
seg_clust_file <- paste0(td, "_data_seg_clusters.txt")
gsetFile <- paste0(td, "_gsetfunnorm.rds")
gBetaFile <- paste0(td, "_gsetbeta.rds")
ClusfiNam <- paste0(td, "_dmp.csv")

# Pathway Analysis Output Files
#genesInputFi <- "yourFile.xlsx"
#siteSpecific <- as.data.frame(readxl::read_excel(gb$genesInputFi))
#your_genes <- siteSpecific$Gene
geneListIn <- paste0(td,"_Genes_Pathway.csv")
pathCsvOut <- paste0(td,"-signaling_pathway.csv")

# Plot Title Names  -----------------------------------------------
clusType <- c("Unsupervised", "Supervised")
if(gb$supervisedRun==F){clusType <- clusType[1]}

varProbes <- c(100, 1000, 10000) # Which top Variance probes to pull i.e. c(100, 1000)
tsne_titles <- gb$generateTitles(clusType, topTitle=as.character(varProbes), gb$titleMain)
topN = 1000
topVar = 10000
mdsTitle <- paste("Top", topN, "Common", "mSet Sq.beta", "MDS plot")
