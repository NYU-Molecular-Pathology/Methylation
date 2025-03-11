#!/usr/bin/env Rscript
## Script name: QC-Scripts.R
## Purpose: source of global scripts imported for QC methylation analysis
## Date Created: June 13, 2022
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025

gb <- globalenv(); assign("gb", gb)
options("install.packages.compile.from.source" = "No")
options("install.packages.check.source" = "no")
options(warn = -1)
# Load/install missing pacakges without asking
supM <- function(pk){return(suppressPackageStartupMessages(suppressWarnings(pk)))}

# Setting US CRAN REPO
options(repos = c(getOption("repos"), CRAN = "http://cran.us.r-project.org"))

pkgs <- c(
    "knitr",
    "kableExtra",
    "magick",
    "webshot",
    "plyr",
    "ggplot2",
    "reshape2",
    "data.table",
    "DT",
    "plotly",
    "MethylAid",
    "minfi",
    "scales",
    "htmltools",
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylationEPICmanifest",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    "Biobase",
    "RColorBrewer",
    "limma",
    "ggfortify",
    "Rtsne",
    "qdapTools",
    "gplots",
    "readxl",
    "stringr",
    "ggrepel",
    "Polychrome",
    "tinytex",
    "gridExtra",
    "rmarkdown",
    "BiocParallel",
    "pals",
    "grid",
    "grDevices",
    "magrittr",
    "dplyr"
)


not_installed <- function(pkgName) {
    return(!pkgName %in% rownames(installed.packages()))
}

if (not_installed("pak")) {
    install.packages("pak", dependencies = T, ask = F, type = "binary")
}
library("pak")

missing_pkgs <- sapply(pkgs, not_installed)

if (any(missing_pkgs)) {
    to_install <- names(missing_pkgs[missing_pkgs == T])
    try(pak::pkg_install(to_install, ask = F), silent = T)
}

if (not_installed("librarian")) {
    try(pak::pkg_install("librarian", ask = F), silent = T)
}

library("librarian")
not_loaded <- setdiff(pkgs, loadedNamespaces())

supM(librarian::shelf(not_loaded, ask = F))

# FUN: Increases vertical spacing between legend keys
draw_key_polygon3 <- function(data, params, size) {
    lwd <- min(data$size, min(size) / 4)
    grid::rectGrob(
        width = grid::unit(1.5, "npc"),
        height = grid::unit(1.5, "npc"),
        gp = grid::gpar(
            col = data$colour,
            fill = scales::alpha(data$fill, data$alpha),
            lty = data$linetype,
            lwd = lwd * .pt,
            linejoin = "mitre"
        )
    )
}
# register new key drawing fun, global & persistent during R session
GeomBar$draw_key = draw_key_polygon3
set.seed(123)

## Generate Plot labels  ------------------
makeLabels <- function(totNum, xName, yName, plotName, thePlot) {
    legendLabel <- thePlot + xlab(xName) + ylab(yName) + ggtitle(plotName) +
        theme(
            plot.background = element_blank(),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10),
            legend.key = element_rect(size = 2),
            axis.title = element_text(size = 14),
            legend.direction = "horizontal",
            legend.position = "bottom",
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            plot.margin = grid::unit(c(1, 1, 1, 1), 'cm'),
            legend.spacing.x = unit(0.20, 'cm'),
            legend.spacing.y = unit(0.65, 'cm'),
            legend.spacing = unit(0.20, 'cm'),
            plot.title = element_text(color = "navy", size = 26, face = "bold"),
            plot.subtitle = element_text(color = "black", size = 12),
            plot.caption = element_text(color = "darkgreen",face = "italic",size = 12)
        ) +
        labs(label = "", color = "") +
        guides(colour = guide_legend(
            title = "Samples",
            override.aes = list(shape = 19),
            byrow = TRUE, ncol = 9), fill = guide_legend(show.legend = F)
        ) + coord_cartesian(clip = 'off')
    legendLabel <- legendLabel + theme(legend.position="none")
    return(legendLabel)
}

## Replace any NA Values in the plot with fail values
ReplaceNAorNull <- function(dParam, xincept, yincept){
    if(any(is.na(dParam$x))){
        dParam$x[is.na(dParam$x)] <- xincept - 1
    }
    if(any(is.na(dParam$y))){
        dParam$y[is.na(dParam$y)] <- yincept - 1
    }
    return(dParam)
}


## Generate Plots for Probes --------------------------------------
plotParams <- function(totNum, dParam, xincept, yincept) {
    dParam = dParam$final_data
    dParam <- ReplaceNAorNull(dParam, xincept, yincept)
    dParam$Sample_Name = paste(dParam$Sample_Name, dParam$MP_num, sep = "\n")
    plot.colours <- glasbey()[1:(length(dParam$x))]


    thePlot <-
        ggplot(dParam, aes(x = dParam[, 2], y = dParam[, 3],
                           color = dParam$Sample_Name, label = dParam$Sample_Name
        ), show.legend = F) +
        scale_color_manual(values = plot.colours) +
        geom_point(shape = 19, size = 5, alpha = 0.8) +
        theme_bw() +
        scale_fill_manual(values = plot.colours) +
        guides(fill = guide_legend(show.legend = F))
    if (yincept == 0){
        thePlot <- thePlot + ggrepel::geom_label_repel(
            aes(
                label = dParam$Sample_Name,
                size = 4,
                label.padding = 0.15,
                label.size = 0.15,
                x = dParam[, 2],
                y = dParam[, 3],
                show.legend = F,
                fill = dParam$Sample_Name,
                colour = scales::alpha(c("black"), 1.0), inherit.aes = F, legend=F
            ),
            fontface = 'bold',
            colour = scales::alpha(c("black"), 1.0),
            alpha = 0.50,
            show.legend = F,
            segment.alpha = 0.30,
            segment.size = 0.75,
            direction = "both",
            max.overlaps = Inf,
            min.segment.length = 0.015,
            color = "black",
            label.size = 0.4,
            size = 4,
            nudge_x = ifelse(thePlot$data$x > mean(thePlot$data$x), -0.15, 0),
            nudge_y = ifelse(thePlot$data$y > mean(thePlot$data$y), -0.05, 0),
            point.padding = unit(0.25, "lines"),
            label.r = unit(0.5, "lines"),
            force = 12,
            max.iter = 10000) + theme(legend.position="none")
        thePlot <- thePlot +
            geom_vline(xintercept = xincept, linetype = 'dashed', colour = "red", inherit.aes = F) +
            coord_cartesian(clip="off") +
            expand_limits(x = min(dParam[,2]), y = max(dParam[,3])*0.25) +
            guides(fill = guide_legend(show.legend = F)) +
            theme(legend.position = "none") +
            scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
            scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
            guides(fill = guide_legend(show.legend = F))
    }
    if (yincept != 0) {
        samlab <- c(str_replace(dParam$Sample_Name, c("\n"), "_"))
        thePlot <- thePlot +
            geom_hline(yintercept = yincept, linetype = 'dashed', colour = "red", inherit.aes = F) +
            annotate("text", x = dParam[, 2], y = 0.5, label = samlab, angle = 90, size = 4) +
            scale_x_continuous(breaks = scales::pretty_breaks()) +
            guides(fill = guide_legend(show.legend = F)) +
            theme(legend.position = "none") + ylim(0.00, 1.00)
    }
    thePlot <- thePlot + theme(legend.position = "none")
    return(thePlot)
}


## Get Plot Specific Probes --------------------------------------
getProbes <- function(probeName) {
    gb <- globalenv(); assign("gb", gb)
    qcProbes <- list(
        BSI = "^BISULFITE CONVERSION I$",
        BSII = "^BISULFITE CONVERSION II$",
        EC = "^EXTENSION$",
        SPI = "^SPECIFICITY I$",
        HYB = "^HYBRIDIZATION$",
        NP = "^NON-POLYMORPHIC$",
        SPII = "^SPECIFICITY II$",
        TR = "^TARGET REMOVAL$",
        SC = "^STAINING$",
        NC = "^NEGATIVE$"
    )
    stopifnot(!is.null(gb$sdata))
    d <- gb$sdata@plotdata
    d <- d[grepl(qcProbes[probeName], d$Type),]
    return(d)
}

## Merges Dataframe Values --------------------------------------
mergeDF <- function(df, tg) {
    mdf <- merge(df, tg, by = "row.names", suffixes = c("", ".y"))
    return(mdf)
}

## MU data --------------------------------------
rotateData <- function(data, columns) {
    data[, columns] <-
        c(0.5 * (data[, columns[1]] + data[, columns[2]]), data[, columns[1]] - data[, columns[2]])
    return(data)
}

## Subset Custom Data --------------------------------------
GetNewDt <- function(data_final){
    dataList <- list(final_data = data_final)
    names(dataList) <- c("final_data")
    return(dataList)
}

## Get Data Values for Plot --------------------------------------
getData <- function(theD, exGrn, exRed, cutoff, targets){
    dGrn <- theD[(exGrn), c(1:5, 7)]
    x <- tapply((dGrn$IntGrn), dGrn$Samples, mean)
    is.na(x) <- !is.finite(x)
    dRed <- theD[as.array(exRed), c(1:6)]
    df <- data.frame(x, y = tapply(dRed$IntRed, dRed$Samples, mean))
    mdf <- mergeDF(df,tg=targets)
    plot_data <- gb$rotateData(mdf,columns = c("x", "y"))
    fdt <- plot_data; ott <- subset(plot_data, plot_data$x <= cutoff)
    return(list(final_data = fdt,outlier = ott))
}

get.hc.dat <- function(targets) {
    d <- getProbes("HYB")
    d <- d[order(d$Samples),]
    hiD <- grepl("High", d$ExtendedType)
    loD <- grepl("Low", d$ExtendedType)
    x <- 0.5 * (d$IntGrn[hiD] + d$IntGrn[loD])
    y <- d$IntGrn[hiD] - d$IntGrn[loD]
    df <- data.frame(x = x, y = y, row.names = d$Samples[hiD])
    mdf <- mergeDF(df,tg = targets)
    return(mdf)
}

get.bs.dat <- function(targets) {
    bsD = getProbes(probeName = "BSI")
    BSvals <- getData(theD = bsD, exGrn = grepl("C1|C2|C3", bsD$ExtendedType),
                      exRed = grepl("C4|C5|C6", bsD$ExtendedType), cutoff = 10, targets)
    return(BSvals)
}

get.op.dat <- function(targets) {
    newD <- getProbes(probeName = "NP")
    OPvals <- getData(theD = newD, exGrn = newD$ExtendedType %in% c("NP (C)", "NP (G)"),
                      exRed = newD$ExtendedType %in% c("NP (A)", "NP (T)"), cutoff = 11, targets)
    return(OPvals)
}

get.dp.dat <- function(gb){
    dpDt <- gb$sdata@DPfreq
    df <- data.frame(x = 1:length(dpDt), y = dpDt, row.names = names(dpDt))
    dat.dp <- gb$mergeDF(df, tg = gb$targets)
    return(dat.dp)
}

swm <- function(funObj){return(suppressMessages(suppressWarnings(funObj)))}

SetKnitPath <- function(runPath, baseDir){
    system(paste("cd", runPath))
    knitr::opts_knit$set(runPath) # rprojroot::find_rstudio_root_file()
    knitr::opts_knit$set(root.dir=runPath)
    message("\nUsing the following output Directory:\n", baseDir)
    message("\nUsing the following Knit Directory:\n", runPath)
}

GetSheetNamePath <- function(params, baseDir) {
    if (is.null(params$sheetNamePath)) {
        sheetNamePath <- list.files(baseDir, "_samplesheet.csv", full.names = T)
        if (length(sheetNamePath) > 1) {
            warning(">1 samplesheet in the folder")
        }
        sheetNamePath <- sheetNamePath[1]
    } else{
        sheetNamePath = params$sheetNamePath
    }
    if (length(sheetNamePath) == 0) {
        warning("No samplesheet found:\n", baseDir)
        samSheetDir <-
            dir(getwd(), "samplesheet.csv", full.names = T)
        sheetName <- paste0(params$runID, "_samplesheet.csv")
        newOut <-
            file.path(fs::path_home(), "Desktop", params$runID, sheetName)
        fs::file_copy(path = samSheetDir, new_path = newOut)
        sheetNamePath <- file.path(baseDir, sheetName)
    }
    if (is.na(sheetNamePath)) {
        sheetNamePath <- dir(getwd(), "samplesheet.csv", full.names = T)[1]
    }
    message(paste0("Sample sheet name is: ", sheetNamePath))
    return(sheetNamePath)
}

GetTotalPairs <- function(fixerrors) {
    # Checks here if you have odd number of arrays used in makeLabels function
    if (nrow(fixerrors) < 8) {
        warning("Less than 8 samples are run, disabling BiocParallel")
        totNum <- 1
    } else{
        # total arrays for biocParrallel workers
        totNum <- nrow(fixerrors) / 8
    }
    return(totNum)
}

genSumFail <- function(failPlot, failTex, sf = gb$summaryFail) {
    ftx = paste("<h4>The following samples **Failed", failTex, "QC:**</h4>")
    if (length(failPlot) > 0) {
        sf <- paste(sf, paste(ftx, paste0(failPlot, "</li>\n"), sep = "\n<li>"))
    }
    return(sf)
}


GetSummaryTab <- function(mnpOutTb){
    tableHeader <-
        c("RD-number", "B-Number", "TM-number", "Methylation Class", "Classifier Score",
          "Subgroup", "Subscore", "MGMT Status")
    gCol <- ifelse(mnpOutTb$mgmt_status == "methylated", "red", "blue")
    gCel <- ifelse(mnpOutTb$classifier_value < 0.90, "red", "green")
    bNum =  mnpOutTb$b_number
    tableSum <- mnpOutTb %>%
        dplyr::mutate(
            classifier_value = cell_spec(classifier_value, "html", color = gCel, bold = T),
            mgmt_status = cell_spec(mgmt_status, "html", color = gCol)) %>%
        dplyr::select(
            record_id, b_number, tm_number, classifier_score, classifier_value,
            subgroup, subgroup_score, mgmt_status
        ) %>%
        kable(format = "html", booktabs = T, escape = F,
              linesep = "", align = "c", col.names = tableHeader) %>%
        kable_styling("striped", position="left") %>%
        row_spec(which(grepl('control',bNum)),bold=T, color="white",background="orange") %>%
        row_spec(which(grepl('low',bNum)), bold=T, color="white",background="salmon") %>%
        row_spec(which(grepl('_',bNum)),bold=T, color="white",background="salmon") %>%
        kable_styling(latex_options="scale_down")
    tableSum <- tableSum %>% kable_styling(position="left")
    return(tableSum)
}


GetFailedSams <- function(mnpOutTb){
    rNum = mnpOutTb$record_id
    low_vals <- mnpOutTb$classifier_value < 0.90
    if(any(low_vals)){
        lowClassVals <- paste0("**", rNum[low_vals], "**")
        lowValScores <- mnpOutTb$classifier_value[low_vals]
        lowScoring <-  paste(lowClassVals, lowValScores, sep = ", ")
        failedSams <- unlist(lapply(lowScoring, function(x){
            paste("<li>", x,"</li>\n")
        }))
        return(failedSams)
    }else{
        return(NULL)
    }
}


GetControlSam <- function(mnpOutTb){
    tm <- mnpOutTb$tm_number
    controlSample <- paste0(tm[(which(grepl('control', mnpOutTb$b_number)))])
    if (length(controlSample) == 0) {
        controlSample <- paste0(tm[(which(grepl('control', mnpOutTb$record_id)))])
    }
    controlSample <- paste("PC Control:", paste(controlSample))
    if(!is.null(controlSample)) {
        conSam <- paste0('<p><h3 style="color:#3CB371">', controlSample, "</h3></p>")
        return(cat(conSam, sep = "\n\n"))
    }
}

GetNotesData <- function(xlsmSheet){
    library("knitr");library("kableExtra");library("dplyr")
    columnCss <- "border-width:2px;border-style:solid;background-color:rgb(255,250,205);border-color:rgb(105,105,105);width:900px;"
    rwcss <- "color:white;text-align:center;font-weight:bold;border-color:rgb(105,105,105);background-color:rgb(139,69,19);"
    noteData <- as.data.frame(readxl::read_excel(xlsmSheet, 2, "M1:M7", col_types = c("text")), row.names=NULL)
    noteData <- noteData[c(!is.na(noteData[, 1]) & noteData[, 1] != 0), 1]  %>% as.data.frame
    colnames(noteData) <- "Worksheet Notes"
    noteData <- noteData %>%
        kable("html", booktabs = T, escape = F, linesep = "") %>%
        kable_styling("striped", full_width = F, font_size = 16, position = "left") %>%
        column_spec(column = c(1), extra_css = columnCss) %>%
        row_spec(row = 0, font_size = 24, extra_css = rwcss)
    noteData <- noteData %>% kable_styling(position="left")
    return(noteData)
}


GetFixedDf <- function(sheetNamePath, runPath){
    fixerrors <- read.csv(sheetNamePath, strip.white = T)
    basePaths <- file.path(runPath, fixerrors$SentrixID_Pos)
    fixerrors$Basename <- basePaths
    return(fixerrors)
}


GrabTargetsDf <- function(runPath, fixerrors, samSheet = "samplesheet.csv"){
    targets <- minfi::read.metharray.sheet(runPath, samSheet, verbose=F)
    basePaths <- file.path(runPath, fixerrors$SentrixID_Pos)
    targets$basenames <- basePaths # path/to/idat/files
    return(targets)
}

CheckRedCsv <- function(runID) {
    redCsv <- file.path(fs::path_home(), "Desktop", runID, paste0(runID, "_Redcap.csv"))
    if (!file.exists(redCsv)) {
        warning(paste0("File not found - ", redCsv))
        warning("Classifier Summary Table will not generate without _Redcap.csv file")
    }else{
        return(redCsv)
    }
}

CheckParamRunID <- function(params, gb){
    defaultRunID <- paste0(basename(getwd()))
    if (is.null(params$runID)) {runID <- defaultRunID} else{runID = params$runID}
    stopifnot(!is.null(runID)); assign("runID", runID, envir = gb)
    gb$runID <- runID
    return(gb$runID)
}

CheckParamBaseDir <- function(params, runID){
    defaultBaseDir = file.path(fs::path_home(), "Desktop", runID)
    if (is.null(params$baseDir)) {baseDir = defaultBaseDir} else{baseDir <- paste0(params$baseDir)}
    if(!dir.exists(baseDir)){baseDir <- getwd()}
    return(baseDir)
}

CheckParamKnitDir <- function(params, runID){
    defaultRunPath = file.path("/Volumes/CBioinformatics/Methylation/Clinical_Runs", runID)
    if(is.null(params$knitDir)){runPath <- defaultRunPath}else{runPath <- params$knitDir}
    return(runPath)
}


MsgFailedSams <- function(failedSams) {
    if (length(failedSams) > 0) {
        failText = "\n\n<h4>The following samples did not classify with a score >0.90:</h4>\n\n"
        cat(failText, sep = "\n\n")
        cat(failedSams, sep = "")
        cat("", sep = "\n\n")
    } else{
        cat("All samples classified with positive matches", sep = "\n\n")
    }
}


MsgFailedQCs <- function(summaryFail) {
    cat('<p class="comment" style="width:80%;"><ul>')
    if (length(summaryFail)>0) {
        cat(summaryFail, sep = "\n\n")
    } else {
        cat("All samples **Passed each QC Critera**", sep = "\n\n")
    }
    cat('</ul>')
}


# FUNCTION: Filters out columns in Plot data and re-names
FilterColNames <- function(outData, xName, yName) {
    selectedColumns <- c("RunID", "Sample_Name", "DNA_Number", "MP_num", "x", "y")
    newColumnNames <- c("RunID", "RD-number", "B-number", "TM-number", as.character(xName), as.character(yName))

    dataSubset <- if ("final_data" %in% names(outData)) outData[["final_data"]] else outData

    filteredData <- dataSubset[, selectedColumns]
    colnames(filteredData) <- newColumnNames

    if (xName == "Samples") {filteredData <- filteredData[, -5]}
    return(filteredData)
}

# FUNCTION: Merges common columns between two plot data frames
Combine_QC_data <- function(data_1, data2){
    return(merge(data_1, data2, by = c("RunID", "RD-number", "B-number", "TM-number")))
}

# FUNCTION: Saves each qc metric x and y values to a csv file
SaveQCmetrics <- function(gb, dat.mu, dat.op, dat.bs, dat.hc, dat.dp) {
    fileOut <- paste(gb$runID, "qc_data.csv", sep = "_")
    fileOutDir <- file.path(gb$runPath, fileOut)

    outData_MU <- FilterColNames(dat.mu, "Log2sqrt(M*U)", "Log2(M/U)")
    outData_OP <- FilterColNames(dat.op, "log2sqrt(R*G)", "log2(R/G)")

    final_data <- Combine_QC_data(outData_MU, outData_OP)

    outData_BS <- FilterColNames(dat.bs, "BS_log2sqrt(R*G)", "BS_log2(R/G)")
    final_data <- Combine_QC_data(final_data, outData_BS)

    outData_HC <- FilterColNames(dat.hc, "log2sqrt(H*L)", "log2(H/L)")
    final_data <- Combine_QC_data(final_data, outData_HC)

    outData_DP <- FilterColNames(dat.dp, "Samples", "Pvalue")
    final_data <- Combine_QC_data(final_data, outData_DP)

    write.csv(final_data, file = fileOutDir, row.names = F, quote = F)
}
