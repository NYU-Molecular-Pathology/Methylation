gb <- globalenv(); assign("gb", gb)
options("install.packages.compile.from.source" = "No")
options("install.packages.check.source" = "no")
options(warn = -1)
# Load/install missing pacakges without asking
supM <- function(pk){return(suppressPackageStartupMessages(suppressWarnings(pk)))}

# Setting US CRAN REPO
rlis = getOption("repos")
rlis["CRAN"] = "http://cran.us.r-project.org"
options(repos = rlis)

pkgs <-
    c(
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

if(!require("librarian")){
    install.packages("librarian", dependencies=T, verbose=T, Ncpus = 4, quiet=T)}
librarian::shelf(pkgs, ask=F)

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
            override.aes = list(fill = NA, shape = 19),
            byrow = TRUE, ncol = 9)
            ) + coord_cartesian(clip = 'off') 
    return(legendLabel)
}
## Generate Plots for Probes --------------------------------------
plotParams <- function(totNum, dParam, xincept, yincept) {
    dParam = dParam$final_data
    dParam$Sample_Name = paste(dParam$Sample_Name, dParam$MP_num, sep = "\n")
    plot.colours <- glasbey()[1:(length(dParam$x))]
    thePlot <-
        ggplot(
            dParam,
            aes(
                x = dParam[, 2],
                y = dParam[, 3],
                color = dParam$Sample_Name,
                label = dParam$Sample_Name
            ),
            show.legend = F
        ) +
        scale_color_manual(values = plot.colours) +
        geom_point(shape = 19,
                   size = 5,
                   alpha = 0.8) +
        guides(color = guide_legend(nrow = totNum)) + theme_bw() +
        guides(color = guide_legend(ncol = 8)) +
        scale_fill_manual(values = plot.colours) 
  
  if (yincept == 0){
      thePlot <- thePlot + ggrepel::geom_label_repel(
          aes(
              label = dParam$Sample_Name,
              size = 3,
              fill = dParam$Sample_Name,
              colour = scales::alpha(c("black"), 1.0)
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
          size = 3,
          nudge_x = ifelse(thePlot$data$x > mean(thePlot$data$x), -0.15, 0),
          nudge_y = ifelse(thePlot$data$y > mean(thePlot$data$y), -0.05, 0),
          point.padding = unit(0.25, "lines"),
          label.r = unit(0.5, "lines"),
          force = 12,
          max.iter = 10000
      )
      
      thePlot <- thePlot +
          geom_vline(xintercept = xincept, linetype = 'dashed', colour = "red") +
          coord_cartesian(clip="off") +
          expand_limits(x=min(dParam[,2]), y=max(dParam[,3])*0.25) +
          guides(fill = guide_legend(show.legend = F)) + 
          theme(legend.position = "none") + 
          scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
          scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
    }
  if (yincept != 0) {
      samlab <- c(str_replace(dParam$Sample_Name, c("\n"), "_"))
      thePlot <- thePlot +
          geom_hline(yintercept = yincept,
                     linetype = 'dashed',
                     colour = "red") +
          annotate(
              "text",
              x = dParam[, 2],
              y = (dParam[, 3]) - 0.0075,
              label = samlab,
              angle = 90,
              size = 3
          ) +
          scale_x_continuous(breaks = scales::pretty_breaks()) +
          guides(fill = guide_legend(show.legend = F)) +
          theme(legend.position = "none")
      }
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
	mdf <- merge(df,tg,by = "row.names",suffixes = c("", ".y"))
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

genSumFail <- function(failPlot, failTex, sf = summaryFail) {
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
    lowClassVals <- paste0("**",rNum[mnpOutTb$classifier_value < 0.90],"**")
    lowValScores <- mnpOutTb$classifier_value[mnpOutTb$classifier_value < 0.90]
    lowScoring <-  paste(lowClassVals, lowValScores, sep = ", ")
    failedSams <- unlist(lapply(lowScoring, function(x){paste("<li>", x,"</li>\n")}))
    return(failedSams)
}

GetControlSam <- function(mnpOutTb){
    tm <- mnpOutTb$tm_number
    controlSample <- paste0(tm[(which(grepl('control', mnpOutTb$b_number)))])
    if (length(controlSample) == 0) {
        controlSample <- paste0(tm[(which(grepl('control', mnpOutTb$record_id)))])
    }
    controlSample <- paste("PC Control:", paste(controlSample))
    return(controlSample)
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


GrabTargetsDf <- function(runPath, fixerrors){
    targets <- minfi::read.metharray.sheet(runPath, "samplesheet.csv")
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


