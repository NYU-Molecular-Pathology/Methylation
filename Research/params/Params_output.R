gb <- globalenv(); assign("gb", gb)
# Default Output file names  -----------------------------------------------
td <- format(Sys.Date(),"%b%d") # today's date Jun17
rgOut <- file.path(gb$runDir, "data", paste0(td, "_RGset.Rdata"))
rawBetaFi <- file.path(gb$runDir, "data", paste0(td, "_betaRaw.Rdata"))
pValsOutFi <- file.path(gb$runDir, "data", paste0(td, "_DetPvals.Rdata"))
mbfile <- file.path(gb$runDir, "data", paste0(td, "_mSetSq.beta.Rdata"))
unbetaVariance <- file.path(gb$runDir, "data", paste0(td, "_unbetaVariance.RData"))
combatOut <- file.path(gb$runDir, "data", paste0(td, "_combatBetas.Rdata"))

supbetaOut <- file.path(gb$runDir, "data", paste0(td, "_supervisedBetas"))
annotFi <- file.path(gb$runDir, "data", paste0(td, "_annotations.rds"))

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
segFile <- file.path(gb$runDir, "data", paste0(td, "_segfile.csv"))
seg_clust_file <- file.path(gb$runDir, "data", paste0(td, "_data_seg_clusters.txt"))
gsetFile <- file.path(gb$runDir, "data", paste0(td, "_gsetfunnorm.rds"))
gBetaFile <- file.path(gb$runDir, "data", paste0(td, "_gsetbeta.rds"))
ClusfiNam <- file.path(gb$runDir, "data", paste0(td, "_dmp.csv"))

# Pathway Analysis Output Files
#genesInputFi <- "yourFile.xlsx"
#siteSpecific <- as.data.frame(readxl::read_excel(gb$genesInputFi))
#your_genes <- siteSpecific$Gene
geneListIn <- file.path(gb$runDir, "data", paste0(td,"_Genes_Pathway.csv"))
pathCsvOut <- file.path(gb$runDir, "data", paste0(td,"-signaling_pathway.csv"))

# Plot Title Names  -----------------------------------------------
clusType <- c("Unsupervised", "Supervised")
if(gb$supervisedRun==F){clusType <- clusType[1]}

varProbes <- c(100, 1000, 10000) # Which top Variance probes to pull i.e. c(100, 1000)
tsne_titles <- gb$generateTitles(clusType, topTitle=as.character(varProbes), gb$titleMain)
topN = 1000
topVar = 10000
mdsTitle <- paste("Top", topN, "Common", "mSet Sq.beta", "MDS plot")
