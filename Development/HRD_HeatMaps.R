#!/usr/bin/env Rscript
## ---------------------------
## Script name: HRD_HeatMaps.R
## Purpose: source of global scripts imported for research methylation probe analysis
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------

gb <- globalenv(); assign("gb", gb)

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

pkgLis <- c("utils", "grDevices", "stringr", "BiocManager", "ggplot2", "pals", "gridExtra", "fitdistrplus", "ggh4x", "dplyr", "purrr")
bioPkg <- c("minfi", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "BiocParallel", "Biobase")

pkgLoad <- unlist(lapply(pkgLis, CheckCran))
bioLoad <- unlist(lapply(bioPkg, BioCinst))
stopifnot(all(pkgLoad) & all(bioLoad))

mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/Research"
script.list <- c("tableMaker.R", "TsnePlotter.R", "getRGsetBetas.R", "hmFunctions.R")
scripts <- file.path(mainHub, script.list)

invisible(lapply(scripts, function(i){suppressPackageStartupMessages(devtools::source_url(i))}))

gb$LoadHeatMapLibs()
library("tidyverse")
require("tidyverse")


GrabSnpProbes <- function(anno){
    snpCols <- c("Probe_rs" , "Probe_maf", "CpG_rs", "CpG_maf" , "SBE_rs" , "SBE_maf")
    snpAnno <- as.data.frame(anno@listData[snpCols])
    rownames(snpAnno) <- rownames(anno)
    snpAnno <- snpAnno %>% dplyr::filter_all(dplyr::any_vars(!is.na(.)))
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
    island_probes <- island_probes[, c("chr", "pos", "UCSC_RefGene_Name", "Relation_to_Island", "UCSC_RefGene_Group")]
    island_probes$UCSC_RefGene_Group <- stringr::str_split(island_probes$UCSC_RefGene_Group, ";", simplify=T)
    island_probes$UCSC_RefGene_Group <- unlist(lapply(X=1:nrow(island_probes), function(X){
        return(paste(collapse=";", sep=";", unique(island_probes$UCSC_RefGene_Group[X,])))}))
    csvPath <- file.path(getwd(), "ProbesCSV")
    csvFiName <- file.path(csvPath, paste0(geneName, "_probes_used.csv"))
    if (!dir.exists(csvPath)) {dir.create(csvPath, recursive=T)}
    island_probes_sort <- island_probes[order(island_probes$pos),]
    island_probes_sort$ProbeNames <- rownames(island_probes_sort)
    message(paste0(capture.output(island_probes_sort), collapse = "\n"))
    write.csv(island_probes_sort, file = csvFiName, quote = F, row.names = F)
    message("File Saved:", csvFiName)
    return(island_probes_sort)
}


GenerateProbePlot <- function(plot_data, plotCols, geneName) {
    mainTitle <-
        "CG Probe Performance on Illumina EPIC Methylation Arrays for"
    txtTheme <- ggplot2::theme(axis.text.x = element_text(angle = 90, face = "bold", color = "black"), legend.position = "none")
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


SaveDistPNG <- function(geneName, plot2, folderNam = "All_Tissues_Single Gene") {
    imgName <- paste(geneName, "probes_dist.png", sep = "_")
    imgPath <- file.path(getwd(), "BetaDist", folderNam)
    if (!dir.exists(imgPath)) {dir.create(imgPath, recursive = T)}
    message(paste0("\nSaving file:\n"), file.path(imgPath, imgName))
    png(file.path(imgPath, imgName), width = 17, height = 10, units = "in", res = 300)
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
        ggh4x::stat_theodensity() +
        labs(title = pltTtl, x = "Beta Values", y = "Frequency") +
        theme_minimal() + facet_wrap( ~ Probe_Name, scales = "free_y") +
        theme(
            strip.text.x = element_text(size = 10),
            panel.spacing = unit(.05, "lines"),
            panel.border = element_rect(color = "black", fill = NA,linewidth = 1),
            strip.background = element_rect(color = "black", linewidth = 1),
            strip.text =  element_text(size = 8)
        )
    SaveDistPNG(geneName, plot2, "Probe_distributions")
}


SaveBoxPlot <- function(thePlot, ggTable, tbl1, tbl2 = NULL, geneName) {
    imgName <- paste(geneName, "probes_CNV.png", sep = "_")
    imgPath <- file.path(getwd(), "BetaBoxplots")
    
    if (!dir.exists(imgPath)) {
        dir.create(imgPath)
    }
    
    message(paste0("\nSaving file:\n"), file.path(imgPath, imgName))
    png(file.path(imgPath, imgName), width = 17, height = 10, units = "in", res = 300)
    
    if (ncol(ggTable) > 10) {
        gridExtra::grid.arrange(thePlot, tbl1, tbl2, nrow = 3, as.table = T, heights = c(5, 1, 1))
    } else{
        gridExtra::grid.arrange(thePlot, tbl1, nrow = 2, as.table = T, heights = c(5, 1))
    }
    invisible(dev.off())
}

GrabGGtable <- function(island_probes, plotCols) {
    ggTable <- as.data.frame(t(island_probes))
    mytheme <- gridExtra::ttheme_default(
        base_size = 10,
        colhead = list(bg_params = list(fill = plotCols, alpha = 0.8)),
        core = list(padding = grid::unit(c(1.5, 1.5), "mm"))
    )
    if (ncol(ggTable) > 10) {
        ggTable1 <- ggTable[, 1:10]
        ggTable2 <- ggTable[, 10:ncol(ggTable)]
    }
    if (ncol(ggTable) <= 10) {
        tbl1 <-  gridExtra::tableGrob(ggTable, theme = mytheme)
        tbl2 = NULL
    } else{
        tbl1 <- gridExtra::tableGrob(ggTable1, theme = mytheme)
        mytheme2 <- gridExtra::ttheme_default(
            base_size = 10,
            colhead = list(bg_params = list(
                fill = plotCols[10:ncol(ggTable)], alpha = 0.8
            )),
            core = list(padding = grid::unit(c(1.5, 1.5), "mm"))
        )
        tbl2 <- gridExtra::tableGrob(ggTable2, theme = mytheme2)
    }
    return(list(ggTable=ggTable, tbl1=tbl1, tbl2=tbl2))

}

# Function to plot cg probe performance for gene islands
PlotIslandProbes <- function(beta_values, anno, geneName) {
    island_probes <- GetIslandProbes(anno, geneName)
    beta_subset <- beta_values[rownames(island_probes), ]
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

GetBaseNames <- function(filenames){unique(gsub("_(Red|Grn).idat(|\\.gz)$", "", filenames))}

# FUN: Generates the basenames in the input path
CreateBasenames <- function(path,recursive=FALSE) {
    grn.files <- list.files(path, pattern = "_Grn.idat(|\\.gz)$", recursive = recursive,
                            ignore.case = TRUE, full.names = TRUE)
    red.files <- list.files(path, pattern = "_Red.idat(|\\.gz)$", recursive = recursive,
                            ignore.case = TRUE, full.names = TRUE)
    GetBaseNames(c(grn.files, red.files))
}

# FUN: Gets the unique basenames from samples
NameFromBaseName <- function(basenames) {
    Sample_Name <- basename(basenames)
    if(any(duplicated(Sample_Name))) {
        warning("Some duplicated Sample_Name entries")
        Sample_Name <- make.unique(Sample_Name)
    }
    return(Sample_Name)
}

ReadMethParallel <- function(targets, verbose = TRUE) {
    nworkers <- BiocParallel::bpworkers(BiocParallel::bpparam())
    if (nworkers <= 1) {
        stop("Did you registered a biocparallel back-end?")
    }
    y <- rep(1, ceiling(nrow(targets) / nworkers))
    for (i in 2:nworkers) {
        y <- c(y, rep(i, ceiling(nrow(targets) / nworkers)))
    }
    y <- y[1:nrow(targets)]
    jobs <- split(targets, y)
    fun <- function(x) {
        requireNamespace("minfi")
        requireNamespace("Biobase")
        return(try(minfi::read.metharray.exp(targets = x, verbose = TRUE, force = TRUE), silent = T))
    }
    message("Reading multiple idat-files in parallel...")
    res <- BiocParallel::bplapply(jobs, FUN = fun)
    message("Combining the RGsets to one big RGset...")
    rgSet <- res[[1]]
    for (i in 2:length(res)) {
        rgSet <- BiocGenerics::combine(rgSet, res[[i]])
    }
    message("DONE!")
    return(rgSet)
}


# FUN: Create sample sheet if an Illumina one isn't available
CreateSheetFromIdats <- function(idatPath, recursive = F, delim = "_") {
    basenames=CreateBasenames(path, recursive)
    if (length(basenames) == 0) {
        warning("No idat files found.")
        return(NULL)
    }
    dat <- data.frame(do.call(rbind, strsplit(basename(basenames), split=delim)), stringsAsFactors=F)
    if (ncol(dat) < 2) {
        warning(path, " basenames do not appear to correspond to idat files")
    }
    sentrixpos <- grep("^R[0-9][0-9]C[0-9][0-9]$", as.character(unlist(dat[1,])))
    if(length(sentrixpos)==1) {
        temp <- do.call(rbind, strsplit(as.character(dat[,sentrixpos]), split="C"))
        dat$sentrix_row <- gsub("R", "", temp[,1])
        dat$sentrix_col <- temp[,2]
        dat <- dat[,-sentrixpos]
    }
    slidecol <- grep("^[0-9]{9}[0-9]*$", as.character(unlist(dat[1,])))
    
    if (length(slidecol) == 1) {
        colnames(dat)[slidecol] <- "Slide"
    }
    
    idcol <- which(apply(dat, 2, function(x) all(!duplicated(x))))
    
    if(length(idcol) >= 1) {
        Sample_Name <- dat[,idcol[1]]
        dat <- dat[,-idcol[1],drop=F]
    }
    else {
        Sample_Name <- NameFromBaseName(basenames)
    }

    samplesheet <- data.frame(Sample_Name = Sample_Name, Sex = NA, dat, Basename = basenames, stringsAsFactors=F)
    write.csv(samplesheet, file="samplesheet.csv", quote=F, row.names=F)
    return(samplesheet)
}


# FUN: reads .idat files and preprocess data
ReadProcessData <- function(idatPath = NULL, targets=NULL) {
    invisible(gc(verbose=F))
    if(is.null(idatPath)){
        idatPath <- getwd()
    }
    if(is.null(targets)){
        targets <- CreateSheetFromIdats(idatPath)
    }
    if(file.exists("RGset.rds")){
        RGset <- readRDS("RGset.rds")
    }else{
        BiocParallel::register(BiocParallel::BatchtoolsParam(workers = 8), default = T)
        RGset <- tryCatch(ReadMethParallel(targets),
                          error = function(e) {
                              RGset <- minfi::read.metharray.exp(targets = targets, verbose = T, force = T)
                              return(RGset)
                              }
                          )
        saveRDS(RGset,file="RGset.rds")
    }
    if(file.exists("Mset.rds")){
        Mset <- readRDS("Mset.rds")
    }else{
        Mset <- minfi::preprocessQuantile(RGset)
        saveRDS(Mset, file="Mset.rds")
    }

    anno <- minfi::getAnnotation(RGset)
    return(list(Mset = Mset, anno = anno, RGset = RGset))
}


PlotTissueSample <- function(betaValSub, gn) {
    tissueTypes <- unique(betaValSub$Tissue)
    txtTheme <- ggplot2::theme(axis.text.x = element_text(angle = 90, color = "black"))
    for (tissue in tissueTypes) {
        toKeep <- betaValSub$Tissue == tissue
        tissueSubset <- betaValSub[toKeep, ]

        row.names(tissueSubset) <- 1:nrow(tissueSubset)
        pltTtl3 <- paste("Beta Values for", gn, "Probes in", tissue)
        plot2 <-
            ggplot(tissueSubset, aes(x = Probe_Name, y = Beta_Value)) +
            geom_point(aes(Probe_Name, Beta_Value, color = SampleName), show.legend = F) +
            geom_line(aes(color = SampleName, group = SampleName), show.legend = F) +
            theme_minimal(base_size = 6) + txtTheme +
            labs(title = pltTtl3, x = "Probes", y = "Beta Values") +
            theme(panel.border = element_rect(fill = NA, colour = "grey40")) +
            ylim(0, 1.0) #+ facet_wrap( ~ SampleName, ncol = 4) + theme(strip.text = element_text(size = 12))

        SaveDistPNG(paste0(gn, "_", tissue), plot2, folderNam = "AllSamplesInGene")

        plot3 <-
            ggplot(tissueSubset, aes(x = Probe_Name, y = Beta_Value)) +
            geom_point(aes(Probe_Name, Beta_Value, color = SampleName), show.legend = F) +
            geom_line(aes(color = SampleName, group = SampleName), show.legend = F) +
            theme_minimal(base_size = 6) + txtTheme +
            labs(title = pltTtl3, x = "Probes", y = "Beta Values") +
            theme(panel.border = element_rect(fill = NA, colour = "grey40")) +
            ylim(0, 1.0) + facet_wrap( ~ SampleName, ncol = 4) + theme(strip.text = element_text(size = 12))

        SaveDistPNG(paste0(gn, "_", tissue), plot3, folderNam = file.path("perSample",paste0(tissue,"_AllSamples")))

        plot4 <-
            ggplot(tissueSubset, aes(x = Probe_Name, y = Beta_Value)) +
            geom_point(aes(Probe_Name, Beta_Value, color = SampleName), show.legend = F) +
            geom_line(aes(color = SampleName, group = SampleName), show.legend = F) +
            theme_minimal(base_size = 6) + txtTheme +
            labs(title = pltTtl3, x = "Probes", y = "Beta Values") +
            theme(panel.border = element_rect(fill = NA, colour = "grey40")) +
            ylim(0, 1.0)

        SaveDistPNG(paste0(gn, "_", tissue), plot4, folderNam = file.path("All_Tissue_Per_Gene", gn))
    }
}


SubsetBetaMat <- function(beta_values, cgProbes, targets) {
    levelOrder <- paste(cgProbes)

    betasMat <- as.matrix(t(beta_values[cgProbes,]))

    var_probes <- apply(beta_values[cgProbes,], 1, var)
    if(length(var_probes)>=10){
        topSamples <- var_probes[1:10]
    }else{
        topSamples <- var_probes
    }
    betasMat <- betasMat[, names(topSamples)]
    betaGroups <- targets[rownames(betasMat), "organ"]
    betaValSub <- data.frame(
        Beta_Value = as.vector(betasMat),
        SampleName = rownames(betasMat),
        Probe_Name = colnames(betasMat),
        Tissue = as.vector(betaGroups)
    )
    betaValSub$Probe_Name <- factor(betaValSub$Probe_Name, levels = levelOrder, ordered = T)
    return(betaValSub)
}



SaveAllProbesPlot <- function(gn, betaValSub) {
    txtTheme <- ggplot2::theme(axis.text.x = element_text(angle = 90, color = "black"))
    pltTtl <- paste("EPIC Methylation Beta Values for", gn, "Probes")
    plot2 <- ggplot(betaValSub, aes(x = Probe_Name, y = Beta_Value)) +
        geom_point(aes(Probe_Name, Beta_Value, color = Tissue), show.legend = T) + txtTheme +
        labs(title = pltTtl, x = "Probes", y = "Beta Values")
    SaveDistPNG(paste0(gn, "_organs"), plot2)
}


GetProbesPlotData <- function(myGenes, msetAnno, targets) {
    beta_values <- minfi::getBeta(msetAnno$Mset)
    colnames(beta_values) <- targets$Sample_Name
    rownames(targets) <- targets$Sample_Name

    for (gn in myGenes) {
        csvPath <- file.path(getwd(), "ProbesCSV")
        csvFiName <- file.path(csvPath, paste0(gn, "_probes_used.csv"))
        if(!file.exists(csvFiName)){
            PlotIslandProbes(beta_values, msetAnno$anno, gn)
        }
        probeData <- as.data.frame(read.csv(csvFiName))

        positNames <- probeData[order(probeData$pos),]
        cgProbes <- positNames[,"ProbeNames"]
        betaValSub <- SubsetBetaMat(beta_values, cgProbes, targets)
        SaveAllProbesPlot(gn, betaValSub)

        plot5 <-
            ggplot(betaValSub, aes(x = Probe_Name, y = Beta_Value)) +
            geom_point(aes(Probe_Name, Beta_Value, color = SampleName), show.legend = F) +
            geom_line(aes(color = SampleName, group = SampleName), show.legend = F) +
            theme_minimal(base_size = 6) + txtTheme +
            labs(title = paste("All Samples For", gn), x = "Probes", y = "Beta Values") +
            theme(panel.border = element_rect(fill = NA, colour = "grey40")) +
            ylim(0, 1.0)

        SaveDistPNG(paste0(gn), plot5, folderNam = "All_Samples_Per_Gene")

        PlotTissueSample(betaValSub, gn)
    }
}


# organInfo <- file.path(getwd(), "raw_organs_data.csv")
# targx <- as.data.frame(read.csv(organInfo))
# targy <- as.data.frame(read.csv(file.path(getwd(), "samplesheet_og.csv")))
# targets <- base::merge(targx, targy, by.x="Sample_Name", all=T)
# targets <- as.data.frame(read.csv(file.path(getwd(), "samplesheet.csv")))
# targets$Basename <- file.path(getwd(), targets$SentrixID_Pos)
# targets$Basename <- file.path(idatPath, targets$SentrixID_Pos)
# write.csv(targets, "samplesheet.csv", quote=F, row.names=F)
# targets <- as.data.frame(read.csv(file.path(getwd(), "samplesheet.csv")))

CombineAllProbes <- function(myGenes, betas, msetAnno, csvPath=NULL, csvSuffix = "_probes_used.csv", islandFilter=NULL) {
    allProbes <- NULL
    probesDF <- matrix(ncol = length(myGenes), nrow = 1000, dimnames = list(1:1000, myGenes))
    if(is.null(csvPath)){
        csvPath <- file.path(getwd(), "ProbesCSV")
    }
    for (gn in myGenes) {
        csvFiName <- file.path(csvPath, paste0(gn, csvSuffix))
        if(!file.exists(csvFiName)){
            PlotIslandProbes(betas, msetAnno$anno, gn)
        }
        probeData <- as.data.frame(read.csv(csvFiName))
        if(!is.null(islandFilter)){
            probeData <- probeData[probeData$Relation_to_Island==islandFilter,]
        }
        positNames <- probeData[order(probeData$pos),]
        cgProbes <- positNames[,"ProbeNames"]
        filler <- 1000-length(cgProbes)
        probesDF[,gn] <- c(cgProbes, rep(NA, filler))
        if(is.null(allProbes)){
            allProbes <- cgProbes
        }else{
            allProbes <- c(allProbes, cgProbes)
        }
    }
    return(list(allProbes= allProbes, probesDF=probesDF))
}

gb$drawHeatMap <- function(yourHeatMap) {
    return(
        ComplexHeatmap::draw(
            yourHeatMap,
            merge_legend = F,
            ht_gap = unit(5.5, "cm"),
            heatmap_legend_side = "left",
            annotation_legend_side = "left",
            column_title_gp=grid::gpar(fontsize=10),
            padding = unit(c(20, 20, 20, 20), "mm")
        )
    )
}


gb$GetProbeAverage <- function(csvColumns, betas, pathwayName){
    avgBetas <- GetEmptyDf(csvColumns,betas)
    for (geneCol in 1:ncol(csvColumns)) {
        probeList <- csvColumns[,geneCol]
        for (sam in 1:dim(betas)[2]) {
            probeBetas <- betas[rownames(betas) %in% probeList, sam]
            sampleNam <- colnames(betas)[sam]
            geneNam <- colnames(csvColumns)[geneCol]
            probeAvg <-mean(probeBetas)
            avgBetas[geneNam,sampleNam] <- probeAvg
        }
    }
    pathwayName <- paste(stringr::str_split(pathwayName, " ", simplify = T), collapse="_")
    outFile <- paste(pathwayName, "avgBetas_per_gene.csv" , sep = "_")
    outPath <- file.path(gb$outDir, outFile)
    avgBetas$XXX <- rownames(avgBetas)
    write.csv(avgBetas, file = outPath, row.names=F, na="", quote=F)
    rwNm = NULL
    rwNm <- which(colnames(avgBetas) == "XXX")
    row.names(avgBetas) <- avgBetas[, rwNm]
    avgBetas <- avgBetas[, -rwNm]
    avgBetas <- as.matrix(avgBetas)
    avgBetas <- na.omit(avgBetas)
    colnames(avgBetas) <- gsub(".", "-",  colnames(avgBetas), fixed = TRUE)
    return(avgBetas)
}


GetFontSize <- function(avgBetas, geneNamesHeatMap){
    toLabRows <- ifelse(nrow(avgBetas)<=50, T, F)
    rowsFont = 12
    if(toLabRows==F){
        toLabRows <- T
        warning(paste("Total Number of row labels is >= 50, decreasing row_names_gp font size!"))
        rowsFont = 10
    }
    if(toLabRows==F & nrow(avgBetas)<=100 & geneNamesHeatMap == T){
        toLabRows <- T
        warning(paste("Total Number of row labels is under 100, decreasing row_names_gp font size!"))
        rowsFont = 9
    }
    if(toLabRows==F & nrow(avgBetas)<=100 & geneNamesHeatMap == T){
        warning(paste("Total Number of row labels is over 100 turning off row labels"))
        toLabRows <- F
    }
    if(geneNamesHeatMap == T){
        geneNamesHeatMap <- toLabRows
    }
    return(list(geneNamesHeatMap=geneNamesHeatMap, rowsFont=rowsFont))
}


tierBetas <- function(betas, col_sentrix, RGSet, batchCorrect = F, getSuper = F, topVar = 1:10000, getAll=F) {
    rgLi <- RGSet@colData@listData
    selectSams <- rgLi[[col_sentrix]][rgLi[["Sample_ID"]] %in% colnames(betas)]
    newRg <- RGSet[, RGSet@colData@rownames %in% selectSams]

    if (batchCorrect == T & getSuper == T) {
        superbetas <- gb$batchCorrectBs(betas, newRg, gb$batchCol)
        return(superbetas)
    }
    if (batchCorrect == T & getSuper == F){
        unBetas <- gb$batchCorrectBs(betas, newRg , gb$batchCol)
        return(unBetas)
    }
    if (getSuper == T) {
        superbetas <- gb$getSupervise(betas, newRg, topVar)
        return(superbetas)
    }
    if(getAll==T){
        return(betas)
    } else{
        unBetas <- gb$takeTopVariance(betas, topVar)
        return(unBetas)
    }
}


gb$GetAvgGeneHeatMap <- function(avgBetas, titleValue, ha, geneNamesHeatMap = T, colSplt = NULL, rwsplt = NULL, ...){
    col_fun2 <- circlize::colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("#00008BFF", "#ADD8E6FF", "#FFFFFFFF", "#FF0000FF", "#FF0000FF"))
    titleOfPlot <- paste("Heatmap of", titleValue)

    fontDat <- GetFontSize(avgBetas, geneNamesHeatMap)
    geneNamesHeatMap <- fontDat[[1]]
    rowsFont <- fontDat[[2]]

    colTall <- ncol(avgBetas)*unit(5, "mm")
    if(ncol(avgBetas)*5 <= 20){
        colTall <- colTall + unit(45, "mm")
    }
    rowTall <- ifelse(fontDat[[1]] == T, 5, 2)
    if(nrow(avgBetas)<=10){
        rowTall <- rowTall*2
    }
    hmTopNumbers <- ComplexHeatmap::Heatmap(
        avgBetas,
        col = gb$col_fun2,  ## Define the color scale
        cluster_columns = F,  ## Cluster the columns
        cluster_rows = T,
        show_column_names = T,  ## Show the Column Names (which is sample #)
        column_names_gp = gpar(fontsize = 12),  ## Column Name Size
        show_row_names = fontDat[[1]],  ## Show Row names (which is probes)
        row_names_side = "left",
        row_title_side = "left",
        row_names_gp = gpar(fontsize = rowsFont, labels=rownames(avgBetas)),
        row_title_gp = gpar(fontsize = 12, fontface = "bold"),
        show_row_dend = F,
        show_column_dend = F,
        column_order = colnames(avgBetas),
        width = colTall,
        height = nrow(avgBetas)*unit(rowTall, "mm"),
        use_raster=T,
        show_heatmap_legend = T,
        top_annotation = ha,
        column_title = titleOfPlot,
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        raster_device = "CairoPNG",
        raster_quality = 1,
        raster_resize_mat = TRUE,
        heatmap_legend_param = list(
            title = "Beta Value",
            labels_gp = gpar( fontsize = 10),
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            legend_direction = "vertical",
            heatmap_legend_side = "left", annotation_legend_side = "left",
            legend_height =  unit(2, "in"),
            at = c(0, 0.25, 0.5, 0.75, 1),
            labels = c("0", "0.25", "0.5", "0.75", "1")
        ),
        column_split = colSplt,
        row_split = rwsplt
    )
    #hm@ht_list[[1]]@matrix_color_mapping@levels <- c(0, 0.25, 0.5, 0.75, 1)
    gb$GetHmDimensions(hmTopNumbers)
    return(
        gb$drawHeatMap(hmTopNumbers)
        )
}


GetCsvGeneColumns <- function(pathwayName, z, outDir=NULL){
    if(is.null(outDir)){
        outDir <- file.path(getwd(),"csv")
    }
    outFile1 <- paste(pathwayName,"genes_columns_withMissing.csv", sep = "_")

    if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}
    write.csv(z, file = file.path(outDir, outFile1), row.names=F, na="")
    z <- as.data.frame(z)
    toDrop <- c(is.na(z[1,]))
    z <- z[,!toDrop]
    outFile2 <- paste(pathwayName,"genes_columns.csv", sep = "_")
    write.csv(z, file = file.path(outDir, outFile2), row.names=F, na="")
    csvColumns <- as.data.frame(read.csv(file.path(outDir, outFile2)))
    return(csvColumns)
}


SaveAvgBetaPlot <- function(gb, outDir, targets, avgBetas, titleValue, preFix = "HRD_Genes", topVar = "TopVarianceProbes") {
    if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}
    ha <- NULL
    ha <- gb$GetAvgBetaAnno(targets, avgBetas, gb$selectedVars)
    ha@anno_list[[2]] <- NULL
    hm <- gb$GetAvgGeneHeatMap(avgBetas, titleValue, ha, geneNamesHeatMap = T)
    gb$SaveHmPng(preFix, fi_suffix = "samples", hm, topvar = topVar, outDir = NULL)
}


SaveTopVarProbes <- function(myGenes, hmData, probesDF, targets, gb){
    varProbeColumns <- matrix(nrow=10, ncol=length(myGenes), dimnames=list(1:10, myGenes))
    unBetas <- hmData$unBetas[,targets[order(targets$organ),"Sample_Name"]]
    titleValue <- paste("Gene Methylation Status for:\n", "Top Variance Probes")
    for(gn in myGenes){
        selecGene <- probesDF[, gn]
        selecGene <- paste(selecGene[!is.na(selecGene)])
        nameOrder <- targets[order(targets$organ),"Sample_Name"]
        unBetas2 <- unBetas[selecGene, nameOrder]
        unBetaSub <- gb$takeTopVariance(unBetas2, 1:10)
        varProbeColumns[,gn] <- rownames(unBetaSub)
    }
    csvColumns <- gb$GetCsvGeneColumns(pathwayName = "TopProbeGenes", varProbeColumns)
    csvColumns <- as.data.frame(csvColumns)
    #avgBetas <- gb$GetProbeAverage(csvColumns, betas=unBetas, "AllGenesProbes")
    allprob <- NULL

    for(x in colnames(csvColumns)){
        if(is.null(allprob)){
            allprob <- csvColumns[, x]
        }else{
            allprob <- c(csvColumns[, x],allprob)
        }

    }
    unBetas3 <-unBetas[allprob, ]
    SaveAvgBetaPlot(gb, gb$outDir, targets, avgBetas=unBetas3, titleValue, preFix="HRD_Genes", topVar="TopVarProbes3")

}


PlotGeneAvgBeta <- function(myGenes, hmData, probesDF, targets, gb) {
    unBetas <- hmData$unBetas
    for(gn in myGenes){
        titleValue <- paste("Gene Methylation Status for:\n", gn)
        selecGene <- probesDF[, gn]
        selecGene <- paste(selecGene[!is.na(selecGene)])
        nameOrder <- targets[order(targets$organ),"Sample_Name"]
        unBetaSub <- unBetas[selecGene,nameOrder]
        unBetas2 <- gb$takeTopVariance(unBetaSub, 1:10)
        ha <- NULL
        ha <- gb$GetAvgBetaAnno(targets, unBetas2, gb$selectedVars)
        ha@anno_list[[2]] <- NULL
        hm <- gb$GetAvgGeneHeatMap(unBetas2, titleValue, ha, geneNamesHeatMap = T)
        gb$SaveHmPng(gn, fi_suffix = "All_Samples", hm, topvar = "Top10_probes_rowclust", outDir = NULL)
    }
}

