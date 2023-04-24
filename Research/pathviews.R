#!/usr/bin/env Rscript
## ---------------------------
## Script name: pathviews.R
## Purpose: source of global scripts imported for research methylation pathview analysis
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------

gb <- globalenv(); assign("gb", gb)
formals(library)$quietly <- T
formals(library)$warn.conflicts <- F
formals(require)$warn.conflicts <- F
formals(install.packages)$dependencies <- T
formals(install.packages)$verbose <- T
formals(install.packages)$ask <- F 
if(!require("librarian")){install.packages("librarian", dependencies=T, verbose=T, Ncpus = 4, quiet=T)}

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

librarian::shelf(pkgs, ask=F, update_all = F, quiet = F)
library(clusterProfiler); require("clusterProfiler")
library(pathview)
library(org.Hs.eg.db)
library("DOSE"); require("DOSE")
library("enrichplot"); require("enrichplot")
#data(geneList)

if(!require("yulab.utils")){
  remotes::install_github("https://github.com/YuLab-SMU/yulab.utils")
}

subsetDmp <- function(annot, beta,i) {
  reg1 ="TSS200|TSS1500"
  reg2="Body"
  dmp.sub <- subset(annot, grepl(i, annot$UCSC_RefGene_Name))
  sub.ucsc <- dmp.sub$UCSC_RefGene_Group
  tss.sub <- subset(dmp.sub, grepl(reg1, sub.ucsc))
  bod.probes <- subset(dmp.sub, grepl(reg2, sub.ucsc))
  beta.promo <-subset(beta, rownames(beta) %in% tss.sub$Name)
  beta.body <-subset(beta, rownames(beta) %in% bod.probes$Name)
  beta.promo$Promoter_Mean <- rowMeans(beta.promo)
  beta.body$Body_Mean <- rowMeans(beta.body)
  beta.pb <- list('pro'=as.data.frame(beta.promo),'bod'=as.data.frame(beta.body))
  return(beta.pb)
}

#' FUNCTION: writeMeansDmp outputs csv files of pathway gene promotor and body regions
#'
#' @param topPaths a lisst of dataframes of the top results from geneVals from `geneVals<-kk_final@result`
#' @param betas beta values dataframe
#' @param annot minfi annotations dataframe with UCSC_RefGene_Group or TSS regions
#' @param pathCsvOut default file suffix is "signaling_pathway.csv" written by hsa gene name
#'
#' @return CSV file with promotor and body gene mean values
#'
writeMeansDmp <- function(topPaths,
                             betas,
                             annot,
                             geneListIn = "Genes_Pathway.csv",
                             pathCsvOut = "signaling_pathway.csv") {
  genelist <- read.delim(geneListIn, header = T, sep = ",", row.names = NULL)
  for (geneRow in 1:nrow(topPaths)) {
    rap1 <- genelist[geneRow,]
    beta <- as.data.frame(betas)
    specific_genes <- as.data.frame(str_split(rap1$geneID, "/"))
    endFile <- paste0(rap1$ID[1], "_", basename(pathCsvOut))
    endPath <- "./data"
    if (!dir.exists(endPath)) {
      dir.create(endPath, recursive = T)
    }
    endFile <- file.path(endPath, endFile)
    if (file.exists(endFile)) {
      message(endFile, " exists")
    } else{
      data_final <- data.frame(Promoter = numeric(0), Body = numeric(0))
      for (rw in 1:nrow(specific_genes)) {
        i = specific_genes[rw, 1]
        beta.pb <- subsetDmp(annot, beta, i)
        promo.mean <- mean(beta.pb$pro$Promoter_Mean)
        body.mean <- mean(beta.pb$bod$Body_Mean)
        avg.join <- c(promo.mean, body.mean)
        data_final[i,] <- c(avg.join[1], avg.join[2])
      }
      write.csv(data_final, file = endFile, quote = F)
    }
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
    kegg.dir = "./figures/pathway",
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
    kegg.dir = "./figures/pathway",
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

renameDots_og <- function(png1, txtIn){
  newName <-stringr::str_remove(png1,"\\.")
  file.rename(from=png1, to=newName)
  newName <- file.path(getwd(),newName)
  imageCap<- paste0("![",txtIn,"](",newName,"){width=90%}")
  return(imageCap)
}

renameDots <- function(png1, txtIn){
  newName <-stringr::str_remove(png1,"\\.")
  newName <- file.path(getwd(),"figures","pathway",newName)
  imageCap<- paste0("![",txtIn,"](",newName,"){width=90%}")
  return(imageCap)
}

rmDupeAnno2 <- function(geneBetas,out.fi) {
  geneRows <- rownames(geneBetas)
  geneRows <-unlist(lapply(1:length(geneRows), function(rw){
    spltRw <- stringr::str_split_fixed(geneRows[rw], pattern = "_", 2)
    rowGene <- spltRw[1, 1]
    if (rowGene == "") {rowGene = "NO.GENE.NAME"}
    rowGene <- paste(unique(stringr::str_split(rowGene, ";")[[1]]), collapse = ";")
    return(paste(rowGene, spltRw[1, 2], sep = "_"))
  }))
  toWrite <- as.data.frame(stringr::str_split_fixed(geneRows,"_",2))
  colnames(toWrite)<- c("Gene_Names", "Probe_Name")
  write.csv(toWrite, out.fi, quote=F,row.names=F)
  rownames(geneBetas)<- geneRows
  return(geneBetas)
}

GetMappedIds <- function(gene_char_unique){
   mappedIDs1 <- gb$supM(AnnotationDbi::mapIds(
       org.Hs.eg.db, keys = gene_char_unique$GeneNames, 
       column = "ENTREZID", keytype = "SYMBOL",
       multiVals = "filter", fuzzy = T))
    mappedIDs2 <- gb$supM(AnnotationDbi::mapIds(
        org.Hs.eg.db, keys = gene_char_unique$UCSCrefseq, 
        column = "ENTREZID", keytype = "REFSEQ",
        multiVals = "filter", fuzzy = T))
    mappedIDs <- c(mappedIDs1, mappedIDs2)
    mappedIDs <- mappedIDs[!is.na(mappedIDs)]
    mappedIDs <- mappedIDs[!duplicated(mappedIDs)]
    mappedIDs <- unique(mappedIDs)
    return(mappedIDs)
}

GetkkDotPlot <- function(termNames, enrichType="NCG", qVa=0.2, pVa=0.05){
  data(geneList, package="DOSE")  
  kkTtl <- paste("Top", enrichType, 
                   "Enriched Pathways qValue <", qVa ,"& pValue <", pVa )
    if(enrichType == "NCG"){
      kk <- DOSE::enrichNCG(names(termNames), qvalueCutoff = qVa, pvalueCutoff = pVa, universe = NULL)
    }else{
      kk <- clusterProfiler::enrichKEGG(gene = names(termNames), organism = 'hsa', universe = NULL,
                                        qvalueCutoff = qVa, pvalueCutoff = pVa)
    }
    totalQ <- length(which(kk@result$qvalue < qVa & kk@result$pvalue < pVa))
    if(totalQ < 10){
      warning(paste("Too few results with kk@result$qvalue <", qVa, "& pvalue", pVa))
      newQval = qVa + 0.1
      newPval = pVa + 0.1
      if(newPval > 1 | newQval > 1){
        cat("\n\n")
        print(enrichplot::dotplot(kk, title = kkTtl))
        return(kk)
      }
      message(paste("Trying higher p and q values:", newQval, "&", newPval))
      GetkkDotPlot(termNames, enrichType, qVa = newQval, pVa = newPval)
      }else{
        cat("\n\n")
        cat(paste("###", "Dotplot", enrichType, "Enriched Paths", "\n\n"))
        try(print(enrichplot::dotplot(kk, title = kkTtl, showCategory=10) +
                    theme(text = element_text(size = 12))), silent=T)
        cat("\n\n")
        cat(paste("###", "BarPlot", enrichType, "Enriched Paths", "\n\n"))
        try(print(graphics::barplot(kk) + ggplot2::ggtitle(kkTtl) +
                    theme(text=element_text(size=12))), silent = T)
        cat("\n\n")
        return(kk)
        }
}


PrintKkDotPlots <- function(termNames){
  edo <- GetkkDotPlot(termNames, enrichType = "NCG", qVa= 0.2, pVa = 0.05)
  kk <- GetkkDotPlot(termNames, enrichType = "KEGG", qVa= 0.2, pVa = 0.05)
  return(list("edo" = edo, "kk" = kk))
}


entrz2kegg <- function(gene_char_unique) {
    if (!require("ggstar")) {install.packages("ggstar", ask = F)}
    library("org.Hs.eg.db")
    library("clusterProfiler")
    library("enrichplot")
    if (!require("DOSE")) {remotes::install_github("GuangchuangYu/DOSE")}
    library("DOSE")
    requireNamespace("clusterProfiler")
    requireNamespace("enrichplot")
    mappedIDs <- GetMappedIds(gene_char_unique)
    mappedIDs <- unique(mappedIDs)
    termNames <- gb$supM(
      AnnotationDbi::mapIds(org.Hs.eg.db, keys = mappedIDs, column = "SYMBOL", 
                            keytype = "ENTREZID", multiVals = "filter", fuzzy = T)
      )
    kk <- PrintKkDotPlots(termNames)
    return(kk)
}


PrintPathways <- function(pathWayGenes){
    dtOpts <- list(scrollX = T, scrollY=T, info = F, autoWidth = F,
                   pageLength = 10, rownames=F, lengthChange = T, searchable = T)
    cat("\n\n")
    cat("## Top Pathways Genes List Tables {.tabset}")
    cat("\n\n")
    for(pthwy in pathWayGenes$Description) {
        pgenes <- gb$splitByPathway(pthwy, pathWayGenes)
        cat(paste0("### **", colnames(pgenes),"** \n\n"))
        dtTab <- htmltools::tagList(DT::datatable(
            pgenes, rownames = F, options=dtOpts, height = "120%", width="100%")
            )
        print(htmltools::tagList(dtTab))
        cat("\n\n")
    }
    cat("#")
}

PrintDotBarUpset <- function(kk_final){
    cat("### DotPlot \n\n")
    enrichplot::dotplot(kk_final, showCategory=10)
    cat("\n\n")
    cat("### BarPlot \n\n")
    barplot(kk_final, showCategory = 10, beside = T) 
    cat("\n\n")
    cat("### Upset Plot \n\n")
    enrichplot::upsetplot(kk_final)
    cat("\n\n")
}


GetHsaPath <- function(topPaths, pathCsvOut){
    endFile <- paste0(topPaths$ID,"_", basename(pathCsvOut))
    hsaOutDir <- file.path(".", "figures", "pathway")
    if (!dir.exists(hsaOutDir)) {dir.create(hsaOutDir)}
    hsaOutFi <- file.path(hsaOutDir, endFile)
}


LoopHSAfiles <- function(hsaOutFi){
    for (pathFi in hsaOutFi) {
    hsaPath <- stringr::str_split_fixed(basename(pathFi), "_", 2)[[1]]
    hsaPath <- file.path(".", "figures", "pathway", hsaPath)
    if (file.exists(pathFi)) {
        suppressMessages(gb$pathview_promoter(pathFi, hsaPath))
        suppressMessages(gb$pathview_body(pathFi, hsaPath))
    } else{
        message(pathFi, " does not exist.")
    }
}
}


LoopSaveHsaPng <- function(pathWayGenes, pathCsvOut){
  endFile <- paste0(pathWayGenes$ID,"_",pathCsvOut)
  for(i in 1:length(endFile)){
    hsaPath <- stringr::str_split_fixed(endFile[i],"_",2)[[1]]  
    message(hsaPath)
    png1 <- paste0(hsaPath,".body.png")
    png2 <- paste0(hsaPath,".promotor.png")
    if(file.exists(png1)){
    imageCap <- gb$renameDots(png1, "Gene Body Mean Values")
    cat(paste0("## ", pathWayGenes$Description[i], "\n\n"))
    cat(imageCap)
    cat("\n\n")
    imageCap <- gb$renameDots(png2,"Gene Promotor Mean Values")
    cat(imageCap)
    cat("\n\n")
    }
  }
}


WritePathVals <- function(geneVals, geneListIn, pathNum=c(1:35)){
    # Sort lowest Pvalues and lowest qvalue
    message("Min p-value: ", min(geneVals$pvalue))
    cat("\n\n")
    print(htmltools::tagList(DT::datatable(geneVals[order(geneVals$pvalue), 1:ncol(geneVals)], options = list(rownames=F))))
    cat("\n\n")
    topPaths <- geneVals[order(geneVals$qvalue),]
    if(nrow(topPaths) >=10){
      topPaths <- topPaths[pathNum,] # take top 10 pathways
    }
    pathWayGenes <- as.data.frame(topPaths)
    pathOutFi <- file.path(getwd(),"figures","pathway", geneListIn)
    if(!dir.exists(dirname(pathOutFi))){dir.create(dirname(pathOutFi))}
    write.csv(pathWayGenes, file = pathOutFi, row.names = F, quote = F)
    return(pathWayGenes)
}

GetKeggGeneVals <- function(RGSet, targets, nameGrp=NULL, gb) {
    if(is.null(nameGrp)){
      nameGrp <- unique(targets[,gb$col_samGrp])[1]
    }
    allUniqueGenes <- file.path(getwd(), "data", "gene_char_unique.rds")
    if(!file.exists(allUniqueGenes)){
      gset.funnorm <- gb$grabGsetFun(RGSet, targets, gb)
      condition <- minfi::pData(gset.funnorm)[,gb$col_samGrp] == nameGrp
      gene_char_unique <- gb$getDmpData(gb$ClusfiNam, gset.funnorm, condition, gb)  
      saveRDS(gene_char_unique, allUniqueGenes)
    }else{
      gene_char_unique <- readRDS(allUniqueGenes)
    }
    kkData <- gb$entrz2kegg(gene_char_unique)
    kk_final <- DOSE::setReadable(kkData$kk, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID")
    geneVals <- kk_final@result
    cat("\n\n")
    if (any(geneVals$pvalue < 0.05)) {
        geneVals <- geneVals[geneVals$pvalue < 0.05, ]
    } else{
        cat("# No P-values below 0.05! \n\n")
    }
    return(list("geneVals" = geneVals, "kkData" = kkData))
}


# Get the current figure size in pixels:
# get_w <- function() {
#   with(knitr::opts_current$get(c("fig.width", "dpi", "fig.retina")),
#        fig.width*dpi/fig.retina)
# }

# get_h <- function() {
#   with(knitr::opts_current$get(c("fig.height", "dpi", "fig.retina")),
#        fig.height*dpi/fig.retina)
# }

