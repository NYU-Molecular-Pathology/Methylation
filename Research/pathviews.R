gb <- globalenv(); assign("gb", gb)
if(!require("librarian")){install.packages("librarian", dependencies=T, verbose=T, Ncpus = 6, quiet=T)}

# Setting US CRAN REPO
rlis = getOption("repos")
rlis["CRAN"] = "http://cran.us.r-project.org"
options(repos = rlis)

# List of all Packages
pkgs <- c(
  'clusterProfiler',
  'pathview',
  'org.Hs.eg.db',
  'DOSE',
  'enrichplot')

librarian::shelf(pkgs, ask=F)

library(clusterProfiler); require("clusterProfiler")
library(pathview)
library(org.Hs.eg.db)
library("DOSE"); require("DOSE")
library("enrichplot"); require("enrichplot")
data(geneList)

####### Pathview analysis using clusterprofiler ########
pathview_promoter <- function(file, pathwayid) {
  pathway = read.csv(file, header = T, sep = ",", row.names = 1)
  
  ids <- clusterProfiler::bitr(rownames(pathway),
              fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = "org.Hs.eg.db")
  
  geneList_promoter <- pathway[, 1]#promoter
  names(geneList_promoter) <- as.character(ids[, 2])
  geneList_promoter <- sort(geneList_promoter, decreasing = TRUE)
  
  pathview::pathview(
    gene.data  = geneList_promoter,
    pathway.id = pathwayid,
    species    = "hsa",
    limit      = list(gene = c(0, 1), cpd = 0.1),
    out.suffix = "promotor"
  )
}

pathview_body <- function(file, pathwayid) {
  pathway = read.csv(file, header = T, sep = ",", row.names = 1)
  ids <- clusterProfiler::bitr(rownames(pathway),
              fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = "org.Hs.eg.db")
  
  geneList_body <- pathway[, 2] # body
  names(geneList_body) <- as.character(ids[, 2])
  geneList_body <- sort(geneList_body, decreasing = TRUE)
  
  pathwayid_body <- pathview::pathview(
    gene.data  = geneList_body,
    pathway.id = pathwayid,
    species    = "hsa",
    limit      = list(gene = c(0, 1), cpd = 0.1),
    out.suffix = "body"
  )
}
