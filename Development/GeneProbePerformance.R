#!/usr/bin/env Rscript
## ---------------------------
## Script name: GenePerformProbes.R
## Purpose: Select probes for island regions of an input gene and test probe performance
## Author: Jonathan Serrano
## Copyright (c) NYULH 2023
## ---------------------------

# Global input variables
myGenes <- c("BRCA1", "BRCA2", "ATM", "RAD51C")
setwd("/Volumes/CBioinformatics/Methylation/Clinical_Runs/Probe_performance")

# Load or Install required packages
LibLoad <- function(pkg){
    libOpt <- list(pkg, quietly = T, logical.return = T, warn.conflicts = F)
    suppressPackageStartupMessages(suppressWarnings(do.call(library, libOpt)))
}
CheckReq <- function(pkg){
    return(suppressWarnings(!require(pkg, character.only = T, warn.conflicts = F)))
}
BioCinst <- function(pkg) {
    if (CheckReq(pkg)) {do.call(BiocManager::install, list(pkg, update = F, ask = F))}
    return(LibLoad(pkg))
}
CheckCran <- function(pkg) {
    if(CheckReq(pkg)){do.call(install.packages, list(pkg, dependencies=T))}
    return(LibLoad(pkg))
}

pkgLis <- c("utils", "grDevices", "stringr", "BiocManager", "ggplot2",
            "pals", "gridExtra", "fitdistrplus", "ggh4x", "dplyr", "purrr")
bioPkg <- c("minfi", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
pkgLoad <- unlist(lapply(pkgLis, CheckCran))
bioLoad <- unlist(lapply(bioPkg, BioCinst))
stopifnot(all(pkgLoad) & all(bioLoad))

GrabSnpProbes <- function(anno){
    snpCols <- c("Probe_rs" , "Probe_maf", "CpG_rs", "CpG_maf" , "SBE_rs" , "SBE_maf")
    snpAnno <- as.data.frame(anno@listData[snpCols])
    rownames(snpAnno) <- rownames(anno)
    snpAnno <- snpAnno %>% filter_all(any_vars(!is.na(.)))
    snpProbes <- rownames(snpAnno)
    return(snpProbes)
}


GetIslandProbes <- function(anno, geneName, minProbes = 4) {
    gene_probes <- anno[grepl(anno$UCSC_RefGene_Name, pattern = geneName),]
    snpProbes <- GrabSnpProbes(anno)
    toKeep <- which(!rownames(gene_probes) %in% snpProbes)
    gene_probes <- gene_probes[toKeep,]
    #islandRegions <- gene_probes$Relation_to_Island == "Island"
    island_probes <- gene_probes %>% as.data.frame() %>%
        mutate(UCSC_RefGene_Name = str_split(UCSC_RefGene_Name, ";", simplify = F) %>%
                   map( ~ sort(unique(.x)))) %>%
        filter(map_lgl(UCSC_RefGene_Name, ~ any(.x == geneName))) %>%
        mutate(UCSC_RefGene_Name = map_chr(UCSC_RefGene_Name, ~ paste(.x, collapse = ";")))
    message("")
    message(paste(geneName, "probes found:"))
    probeCols <- c("chr", "pos", "UCSC_RefGene_Name", "Relation_to_Island", "UCSC_RefGene_Group")
    island_probes <- island_probes[, probeCols]
    csvFiName <- paste0(geneName, "_probes_used.csv")
    csvPath <- file.path(getwd(), "ProbesCSV")
    if (!dir.exists(csvPath)) {
        dir.create(csvPath)
    }
    message(paste0(capture.output(island_probes), collapse = "\n"))
    write.csv(
        island_probes, file = file.path(csvPath, csvFiName), quote = F, row.names = T
        )
    return(island_probes)
}


GenerateProbePlot <- function(plot_data, plotCols, geneName) {
    mainTitle <-
        "CG Probe Performance on Illumina EPIC Methylation Arrays for"
    txtTheme <-
        ggplot2::theme(
            axis.text.x = element_text(
                angle = 90,
                face = "bold",
                color = "black"
            ),
            legend.position = "none"
        )
    plot_data <- stack(plot_data)
    means <- aggregate(values ~ ind, plot_data, mean)
    medians <- aggregate(values ~ ind, plot_data, median)
    thePlot <-
        ggplot2::ggplot(plot_data, aes(x = ind, y = values, fill = ind)) +
        geom_boxplot(
            alpha = 0.8,
            outlier.colour = "red",
            outlier.size = 4
        ) +
        scale_fill_manual(values = plotCols) +
        ggplot2::geom_dotplot(
            binaxis = 'y',
            stackdir = 'center',
            dotsize = 0.5,
            binwidth = 0.004
        ) +
        geom_label(
            data = means,
            aes(label = paste("Mean:", round(values, 3)), y = -0.02),
            size = 3,
            alpha = 0.5,
            color = "darkred",
            fill = "white",
            fontface = "bold"
        ) +
        geom_label(
            data = medians,
            aes(label = paste("Median:", round(values, 3)), y = -0.05),
            size = 3,
            alpha = 0.8,
            color = "darkblue",
            fill = "white",
            fontface = "bold"
        ) +
        ggplot2::theme_bw() + txtTheme + labs(x = "Probe Names", y = "Beta Values") +
        ggplot2::ggtitle(paste(mainTitle, geneName, "Islands"),)
    return(thePlot)
}


SaveDistPNG <- function(geneName, plot2) {
    imgName <- paste(geneName, "probes_dist.png", sep = "_")
    imgPath <- file.path(getwd(), "BetaDist")
    if (!dir.exists(imgPath)) {
        dir.create(imgPath)
    }
    message(paste0("\nSaving file:\n"), file.path(imgPath, imgName))
    png(
        file.path(imgPath, imgName),
        width = 17,
        height = 10,
        units = "in",
        res = 300
    )
    print(plot2)
    invisible(dev.off())
}


PlotBetaDist <- function(betasMat, island_probes, geneName) {
    pltTtl <- paste("EPIC Methylation Beta Values Distributions for",
              geneName, "Probes")
    betaValSub <- data.frame(Beta_Value = as.vector(betasMat), Probe_Name = rownames(betasMat)) %>%
        filter(Probe_Name %in% rownames(island_probes)) %>%
        mutate(Probe_Name = paste(Probe_Name, island_probes[Probe_Name, "Relation_to_Island"],
                                  island_probes[Probe_Name, "UCSC_RefGene_Group"], sep = "_"))

    plot2 <- ggplot(betaValSub, aes(x = Beta_Value)) +
        geom_histogram(bins = 50, fill = "blue", alpha = 0.6, aes(y = after_stat(density))) +
        stat_theodensity() +
        labs(title = pltTtl, x = "Beta Values", y = "Frequency") +
        theme_minimal() + facet_wrap( ~ Probe_Name, scales = "free_y") +
        theme(
            strip.text.x = element_text(size = 10),
            panel.spacing = unit(.05, "lines"),
            panel.border = element_rect(color = "black", fill = NA,linewidth = 1),
            strip.background = element_rect(color = "black", linewidth = 1)
        )
    SaveDistPNG(geneName, plot2)
}

SaveBoxPlot <- function(thePlot, ggTable, tbl1, tbl2 = NULL, geneName) {
    imgName <- paste(geneName, "probes_CNV.png", sep = "_")
    imgPath <- file.path(getwd(), "BetaBoxplots")
    if (!dir.exists(imgPath)) {
        dir.create(imgPath)
    }
    message(paste0("\nSaving file:\n"), file.path(imgPath, imgName))
    png(
        file.path(imgPath, imgName),
        width = 17,
        height = 10,
        units = "in",
        res = 300
    )
    if (ncol(ggTable) > 10) {
        gridExtra::grid.arrange(
            thePlot,
            tbl1,
            tbl2,
            nrow = 3,
            as.table = T,
            heights = c(5, 1, 1)
        )
    } else{
        gridExtra::grid.arrange(
            thePlot,
            tbl1,
            nrow = 2,
            as.table = T,
            heights = c(5, 1)
        )
    }
    invisible(dev.off())
}

GrabGGtable <- function(island_probes, plotCols) {
    ggTable <- as.data.frame(t(island_probes))
    mytheme <- gridExtra::ttheme_default(
        base_size = 10,
        colhead = list(bg_params = list(fill = plotCols, alpha = 0.8)),
        core = list(padding = unit(c(1.5, 1.5), "mm"))
    )
    if (ncol(ggTable) > 10) {
        ggTable1 <- ggTable[, 1:10]
        ggTable2 <- ggTable[, 10:ncol(ggTable)]
    }
    if (ncol(ggTable) <= 10) {
        tbl1 <- tableGrob(ggTable, theme = mytheme)
        tbl2 = NULL
    } else{
        tbl1 <- tableGrob(ggTable1, theme = mytheme)
        mytheme2 <- gridExtra::ttheme_default(
            base_size = 10,
            colhead = list(bg_params = list(
                fill = plotCols[10:ncol(ggTable)], alpha = 0.8
            )),
            core = list(padding = unit(c(1.5, 1.5), "mm"))
        )
        tbl2 <- tableGrob(ggTable2, theme = mytheme2)
    }
    return(list(ggTable=ggTable, tbl1=tbl1, tbl2=tbl2))

}

# Function to plot cg probe performance for gene islands
PlotIslandProbes <- function(beta_values, anno, geneName) {

    island_probes <- GetIslandProbes(anno, geneName)
    beta_subset <- beta_values[rownames(island_probes),]
    plot_data <- as.data.frame(t(beta_subset))
    allCols <- paste(pals::polychrome(36))[-2]
    plotCols <-  c("darkred", "royalblue")
    if (ncol(plot_data) > 2) {
        plotCols <- sample(allCols, ncol(plot_data), replace = T)
    }
    ggTab <- GrabGGtable(island_probes, plotCols)
    #thePlot <- GenerateProbePlot(plot_data, plotCols, geneName)
    #SaveBoxPlot(thePlot, ggTable=ggTab[["ggTable"]], tbl1=ggTab[["tbl1"]], tbl2=ggTab[["tbl2"]], geneName)
    PlotBetaDist(beta_values, island_probes, geneName)
}


#rgset <- minfiDataEPIC::RGsetEPIC
targets <- as.data.frame(read.csv(file.path(getwd(), "samplesheet.csv")))
targets$Basename <- file.path(getwd(), targets$SentrixID_Pos)
rgset <- minfi::read.metharray.exp(getwd(), targets, verbose = T, force = T)
snpData <- minfi::getSnpInfo(rgset)
mset <- minfi::preprocessRaw(rgset)
beta_values <- minfi::getBeta(mset)
anno <- minfi::getAnnotation(rgset)

for (gn in myGenes) {
    PlotIslandProbes(beta_values, anno, gn)
}

# function to plot beta values 
plot_beta_values <- function(Mset, best_probes, sample_groups) {
  beta_values <- getBeta(Mset)
  data <- data.frame(Probe = rep(best_probes, each = length(sample_groups)),
                     Beta_Value = beta_values[best_probes, ],
                     Sample_Group = rep(sample_groups, times = length(best_probes)))

  ggplot(data, aes(x = Probe, y = Beta_Value, color = Sample_Group)) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Probe", y = "Beta Value", color = "Sample Group")
}

# With a line added:

plot_beta_values <- function(Mset, best_probes, sample_groups) {
  beta_values <- getBeta(Mset)
  data <- data.frame(Sample = rep(1:nrow(Mset), each = length(best_probes)),
                     Probe = rep(best_probes, each = nrow(Mset)),
                     Beta_Value = c(beta_values[best_probes, ]),
                     Sample_Group = rep(sample_groups, times = length(best_probes)))

  ggplot(data, aes(x = Probe, y = Beta_Value, group = Sample, color = Sample_Group)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("A" = "red", "B" = "blue", "C" = "green")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Probe", y = "Beta Value", color = "Sample Group")
}

# Main function
main <- function() {
  # Provide path to the folder containing the .idat files
  idat_files_path <- "path/to/your/idat/files/folder"

  Mset <- read_and_preprocess_data(idat_files_path)
  best_probes <- identify_best_probes(Mset, gene_list)

  # Define sample groups
  sample_groups <- c(rep("A", 12), rep("B", 8), rep("C", 9))

  plot_beta_values(Mset, best_probes, sample_groups)
}
