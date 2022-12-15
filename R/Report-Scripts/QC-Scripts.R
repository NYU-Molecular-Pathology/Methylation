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
              y = (dParam[, 3]) - 0.007,
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
