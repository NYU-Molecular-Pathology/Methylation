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

#' FUNCTION: writeMeansDmp outputs csv files of pathway gene promotor and body regions
#'
#' @param topPaths a lisst of dataframes of the top results from geneVals from `geneVals<-kk_final@result`
#' @param betas beta values dataframe
#' @param annot minfi annotations dataframe with UCSC_RefGene_Group or TSS regions
#' @param pathCsvOut default file suffix is "signaling_pathway.csv" written by hsa gene name
#'
#' @return CSV file with promotor and body gene mean values
#'
writeMeansDmp<-
function(topPaths, betas, annot, geneListIn= "Genes_Pathway.csv",
         pathCsvOut="signaling_pathway.csv"){
    genelist <- read.delim(geneListIn,header = T,sep = ",",row.names = NULL)
    for (geneRow in 1:nrow(topPaths)) {
        rap1 <- genelist[geneRow, ]
        beta <- as.matrix(betas)
        specific_genes <- as.data.frame(str_split(rap1$geneID, "/"))
        endFile <- paste0(rap1$ID[1], "_", pathCsvOut)
        data_final <- data.frame(Promoter = numeric(0), Body = numeric(0))
        message(endFile)
        for (rw in 1:nrow(specific_genes)) {
            i = specific_genes[rw, 1]
            dmp.sub <- subset(annot, grepl(i, annot$UCSC_RefGene_Name))
            sub.ucsc <- dmp.sub$UCSC_RefGene_Group
            tss200.sub <- subset(dmp.sub, grepl("TSS200|TSS1500", sub.ucsc))
            bod.probes <- subset(dmp.sub,grepl("Body", sub.ucsc))
            beta.promo <- subset(beta, rownames(beta) %in% tss200.sub$Name)
            beta.body <- subset(beta, rownames(beta) %in% bod.probes$Name)
            beta.promo$Promoter_Mean <- rowMeans(beta.promo)
            beta.body$Body_Mean <- rowMeans(beta.body)
            promo.mean <- mean(beta.promo$Promoter_Mean)
            body.mean <- mean(beta.body$Body_Mean)
            means.joined <- c(promo.mean, body.mean)
            means.joined_df <- data.frame(means.joined)
            data_final[i,] <- c(means.joined[1], means.joined[2])
        }
        write.csv(data_final, file = endFile, quote = F)
    }
}

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

splitByPathway <- function(pathDesc, geneVals) {
  selectPath <- geneVals$Description == pathDesc
  pathWayGenes <- data.frame(str_split(geneVals[selectPath, "geneID"], "/"))
  colnames(pathWayGenes) <- pathDesc
  return(pathWayGenes)
}
