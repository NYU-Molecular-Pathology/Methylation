#!/usr/bin/env Rscript
## Script name: KnitScripts.R
## Purpose: Source of global scripts used by the PACT_consensus.Rmd file
## Date Created: August 10, 2022
## Version: 1.0.1
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025

gb <- globalenv()
assign("gb", gb)


# Ensures the brew command can be seen in the PATH environment
fix_brew_path <- function() {
    brew_paths <- c("/opt/homebrew/bin/brew", "/usr/local/bin/brew")
    brew_dir <- dirname(brew_paths[file.exists(brew_paths)])
    current_path <- strsplit(Sys.getenv("PATH"), ":")[[1]]
    updated_path <- unique(c(brew_dir, current_path))
    Sys.setenv(PATH = paste(updated_path, collapse = ":"))
    
    renviron_file <- file.path(Sys.getenv("HOME"), ".Renviron")
    path_entry <- paste0('PATH="', paste(updated_path, collapse = ":"), '"')
    renviron_content <- character(0)
    if (file.exists(renviron_file)) {
        renviron_content <- readLines(renviron_file, warn = FALSE)
    }
    
    if (!any(grepl("^PATH=", renviron_content))) {
        new_content <- c(renviron_content, path_entry)
    } else {
        new_content <- sub("^PATH=.*", path_entry, renviron_content)
    }
    
    if (!identical(new_content, renviron_content)) {
        writeLines(new_content, renviron_file)
    }
}


# Check if Brew and Latex is installed
CheckBrewLatex <- function() {
    tryCatch(
        system("which brew", intern = T),
        warning = function(ww) {
            system('/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"')
        }
    )
    
    fix_brew_path()
    
    if (!"fontconfig" %in% system("brew list", intern = TRUE)) {
        system("brew install fontconfig")
    }
    
    tryCatch(
        system("which latex", intern = T),
        warning = function(ww) {
            system("brew install latex")
        }
    )
    
    tryCatch(
        system("which pdflatex", intern = T),
        warning = function(ww) {
            system("brew install --cask basictex")
        }
    )
}

if (Sys.info()[['sysname']] == "Darwin") {
    invisible(CheckBrewLatex())
}

# Librarian shelf loads or installs package from CRAN, BioConductor, & GitHub
if (!require("librarian")) {
    install.packages(
        "librarian",
        dependencies = T,
        verbose = T,
        Ncpus = 4,
        quiet = T
    )
}

pkgs <-
    c(
        "animation",
        "tidyverse",
        "knitr",
        "kableExtra",
        "stringr",
        "readr",
        "ggplot2",
        "data.table",
        "plotly",
        "DT",
        "datasets",
        "pdftools",
        "htmltools",
        "forcats",
        "stringr",
        "dplyr",
        "purrr",
        "tidyr",
        "tibble",
        "stats",
        "graphics",
        "grDevices",
        "utils",
        "datasets",
        "methods",
        "base",
        "magrittr",
        "crayon",
        "devtools"
    )

pdfToolsPkgs <-
    c(
        'Rcpp',
        'qpdf',
        'askpass',
        'curl',
        'sys',
        'rappdirs',
        'digest',
        'brio',
        'callr',
        'cli',
        'desc',
        'ellipsis',
        'evaluate',
        'jsonlite',
        'lifecycle',
        'pkgload',
        'praise',
        'processx',
        'ps',
        'R6',
        'rlang',
        'waldo',
        'withr',
        'rprojroot',
        'glue',
        'fs',
        'diffobj',
        'fansi',
        'rematch2',
        'pillar',
        'pkgconfig',
        'vctrs',
        'utf8',
        'png',
        'testthat'
    )

librarian::shelf(
    pkgs,
    ask = F,
    update_all = F,
    quiet = F,
    dependencies = T
)
librarian::shelf(
    pdfToolsPkgs,
    ask = F,
    update_all = F,
    quiet = F,
    dependencies = T
)

library("kableExtra")

bold_line <- function() {
  return(cat('<hr style="border:none; ' , 'height:3px; ' ,
        'background-color:#000; ' , 'width:100%;" />\n'))
}

CheckMntDirs <- function(critialMnts, params) {
    outDir <- file.path(params$workDir, paste0(params$pactName, "_consensus"))
    failMount <-
        lapply(critialMnts, function(mnt) {
            ifelse(!dir.exists(mnt), return(T), return(F))
        })
    if (any(failMount == T)) {
        toFix <- paste(critialMnts[which(failMount == T)])
        cat("PATH does not exist, ensure drive is mounted:",
            toFix,
            sep = "\n")
        message(paste(
            "\nTry executing:\nsystem('diskutil mountdisk",
            toFix,
            "')"
        ))
        stopifnot(!any(failMount == T))
    } else{
        message("All network Drive paths are accessible")
    }
    if (!dir.exists(outDir)) {
        dir.create(outDir)
    }
    knitr::opts_knit$set(root.dir = outDir)
    system(paste("cd", outDir))
    setwd(outDir)
    message("Working directory changed to:\n", outDir)
}

StopMissingFile <- function(fiNam, isDir = F) {
    if (isDir == F) {
        if (!file.exists(fiNam)) {
            stop(fiNam, " does not exist in the working directory!")
        }
    } else{
        if (!dir.exists(fiNam)) {
            stop(fiNam,
                 " does not exist as a directory! Check directory path name.")
        }
    }
}


get_ngs_path <- function(outPath, sam) {
    outFileDir <- list.dirs(path = outPath, full.names = T)
    found_folder <- stringr::str_detect(outFileDir, sam)
    if (any(found_folder)) {
        sam_folder <- outFileDir[found_folder]
        return(sam_folder)
    } else {
        return(NULL)
    }
}


CheckFiExist <- function(pactName, philipsFtp) {
    methSheet <- paste0(pactName, "_MethylMatch.xlsx")
    qcTsv <- paste0(pactName, "-QC.tsv")
    descrip <- paste0(pactName, "_desc.csv")
    samsheet <- list.files('.', "demux-samplesheet.csv", T)[1]
    if (!file.exists(methSheet)) {
        message(crayon::bgRed("No methylation sheet for this run!"))
    }
    StopMissingFile(samsheet)
    StopMissingFile(qcTsv)
    StopMissingFile(descrip)
    StopMissingFile(philipsFtp, T)
}


MoveFacets <- function(){
  facet_pngs <- list.files(getwd(), pattern = "*.FACETS.cnv.plot.png", full.names = T)
  facet_dirs <- file.path(getwd(), "facet_pngs")
  if (!dir.exists(facet_dirs)) try(dir.create(facet_dirs, showWarnings = F), T)
  if (length(facet_pngs) > 0) {
    try(fs::file_move(facet_pngs, facet_dirs), TRUE)
  }
}


CheckMethPaths <- function(methData, methSheet) {
    for (i in 1:length(methData$`Report Path`)) {
        currPath <- methData$`Report Path`[i]
        currSplit <- stringr::str_split_fixed(currPath, "/", 11)[1, ]
        if (stringr::str_detect(currSplit[10], "MGDM") == F) {
            next
        }
        runYear <- stringr::str_split_fixed(currSplit[10], "-", 2)[1, 1]
        runYear <- paste0("20", runYear)
        currSplit[9] <- runYear
        newPath <- paste(currSplit, collapse = "/")
        methData[i, "Report Path"] <- newPath
    }
    
    checkPaths <- stringr::str_replace_all(
        methData$`Report Path`,
        "smb://shares-cifs.nyumc.org/apps/acc_pathology",
        "/Volumes"
    )
    checkPaths <- checkPaths[checkPaths != ""]
    anyPathsFalse <- file.exists(checkPaths) == F
    
    if (any(anyPathsFalse)) {
        message("Some paths need to be fixed 'Report Path' column of in ",
                methSheet,
                ":")
        cat("Paths don't exist:\n")
        cat(paste(checkPaths[anyPathsFalse], collapse = "\n"))
    }
    stopifnot(all(file.exists(checkPaths)))
    return(methData)
}


FixMethDataNA <- function(methData, pactName) {
    samCsv <- dir(getwd(), "demux-samplesheet.csv", full.names = T)
    wsData <- as.data.frame(read.csv(samCsv, skip = 19))
    toFix <- which(is.na(methData$Test_Number))
    qcData <- gb$ReadQcFile(pactName)
    for (xRow in toFix) {
        nSam <- methData$accession_number[xRow]
        ngsRow <- which(stringr::str_detect(wsData$Specimen_ID, pattern = nSam))[[1]]
        methData$Test_Number[xRow] <- wsData$Test_Number[ngsRow]
    }
    return(methData)
}


GetMethDf <- function(pactName) {
    methSheet <- paste0(pactName, "_MethylMatch.xlsx")
    if (!file.exists(methSheet)) {
        return(NULL)
    }
    methData <- as.data.frame(readxl::read_excel(methSheet))
    methData <- methData[methData$report_complete == "YES", ]
    if (any(is.na(methData$Test_Number))) {
        methData <- FixMethDataNA(methData, pactName)
    }
    if (nrow(methData) > 0) {
        methData <- CheckMethPaths(methData, methSheet)
    }
    return(methData)
}


GetSamList <- function(pactName, listType = 1) {
    colFltr <- c("Test_Number", "Specimen_ID", "Tumor_Content")
    if (listType == 2) {
        colFltr <- c(colFltr, "Paired_Normal")
    }
    samsheet <- list.files('.', "demux-samplesheet.csv", T)[1]
    message("Reading file: ", samsheet)
    samList <- read.csv(samsheet, skip = 19)[, colFltr]
    isNGS <- !is.na(stringr::str_extract(samList$Test_Number, 'NGS'))
    toKeep <- samList$Test_Number != "0" & isNGS
    if (any(toKeep == F)) {
        message(
            crayon::bgRed(
                "The following cases have no NGS number and are being excluded:"
            ),
            "\n"
        )
        message(paste(utils::capture.output(as.data.frame(samList[!toKeep, ])), collapse = "\n"))
        samList <- samList[toKeep, ]
    }
    stopifnot(nrow(samList) > 2)
    return(samList)
}


GrabSamples <- function(samList) {
    sam <- unique(samList$Test_Number)
    samples <- sam[sam != 0 & !is.na(sam) & sam != "0"]
    return(samples)
}


ReadQcFile <- function(pactName) {
    qcTsv <- paste0(pactName, "-QC.tsv")
    qcData <- read.delim(file = qcTsv)
    row.names(qcData) <- qcData$Row
    return(qcData)
}


makePdfTab <- function(pdfFi, cnvTab, outDir) {
    ngsPng <- file.path(outDir, "ngs_cnv_png", paste0(pdfFi, ".png"))
    cat(
        "<span style='color: red;'>",
        "**Note**: CNV Amplifications are listed first,",
        "followed other Philips CNV abberations in table below this image\n\n"
    )
    cat(paste0("![In-house Facets Plot (Above)](", ngsPng, "){width=100%}"))
    cat("\n\n")
    if (nrow(cnvTab) > 0) {
        newTa <- knitr::kable(cnvTab,
                              row.names = F,
                              "html",
                              align = "c")
        newTa <-
            kableExtra::kable_styling(
                newTa,
                bootstrap_options = c("bordered", "condensed"),
                full_width = F,
                position = "left"
            )
        newTa <- kableExtra::column_spec(newTa, 1:5, width = "3cm")
        cat("#### Philips Data Dump CNV table:\n\n")
        print(newTa)
        cat("\n\n")
    } else{
        cat("\n\nNo CNV Amplifications >=5 in Philips data dump for this case\n\n")
        cat("\n\n")
    }
    
}


makeRdTab <- function(rdNumb, sam) {
    cnvFi <- file.path(".", "methCNV", paste0(rdNumb, "_cnv.png"))
    altTxt <- paste("![Methylation", rdNumb, sam, "CNV Plot](")
    cat(paste0(altTxt, cnvFi, "){height=400px}"))
    cat("\n\n")
}

makeSpecimenTab <- function(objDat) {
    newTa <- knitr::kable(objDat, row.names = F, "html")
    newTa <- kableExtra::kable_styling(
        newTa,
        bootstrap_options = c("condensed"),
        full_width = T,
        position = "left"
    )
    if (ncol(objDat) > 3) {
        newTa <- kableExtra::column_spec(newTa, 1:5, width = "2cm")
    } else{
        newTa <- kableExtra::column_spec(newTa, 1:ncol(objDat), width = "1cm")
    }
    print(newTa)
    cat("\n\n")
}


PrintHotspotTable <- function(objDat) {
    cat(
        "<span style='color: red;'>",
        "**Note**: The Hotspot calls below are from the unfiltered VCF files.  Any called (YES) will be listed first</span>\n\n"
    )
    
    # Reorder the data frame so that rows with "YES" are at the top
    rows_with_yes <- rowSums(objDat[, c("Strelka", "LoFreqSomatic", "Mutect")] == "YES") > 0
    sorted_objDat <- objDat[order(-rows_with_yes), ]
    
    dtOpts <- list(
        scrollX = T,
        scrollY = T,
        info = F,
        autoWidth = F,
        pageLength = 25,
        lengthChange = T,
        searchable = T
    )
    
    datTab <- DT::datatable(
        sorted_objDat,
        style = "bootstrap",
        rownames = FALSE,
        options = dtOpts,
        height = "120%",
        width = "120%"
    )
    
    # Apply the color to cells with "YES" in the specified columns
    datTab <- datTab %>%
        formatStyle(
            columns = c("Strelka", "LoFreqSomatic", "Mutect"),
            backgroundColor = styleEqual("YES", "#64a463"),
            target = 'cell'
        )
    
    # format "Common COSMIC Variant" column in italics
    variantCol <- which(stringr::str_detect(names(sorted_objDat), "Variant"))
    datTab <- datTab %>% formatStyle(variantCol, fontStyle = "italic")
    
    dtTab <- htmltools::tagList(datTab)
    
    print(dtTab)
    cat("\n\n")
}


makeQCTab <- function(objDat) {
    qcHeader <- c(
        "Sample",
        "Mapped Reads QC",
        "Deduplicated Reads QC",
        "Deduplication Rate QC",
        "Targets <50 Coverage",
        "Avg & Med Coverage",
        "Overlapped HOMO SNP"
    )
    colnames(objDat) <- qcHeader
    newTa <- knitr::kable(objDat, row.names = F, "html")
    newTa <- kableExtra::kable_styling(
        newTa,
        bootstrap_options = c("condensed"),
        full_width = F,
        position = "left"
    )
    newTa <- kableExtra::column_spec(newTa, 1:5, width = "5cm")
    print(newTa)
    cat("\n\n")
}

makeDefaultDt <- function(objDat) {
    dtOpts <- list(
        scrollX = T,
        scrollY = T,
        info = F,
        autoWidth = F,
        pageLength = 100,
        rownames = T,
        lengthChange = T,
        searchable = T
    )
    dtTab <- htmltools::tagList(
        DT::datatable(
            objDat,
            style = "bootstrap",
            rownames = F,
            options = dtOpts,
            height = "120%",
            width = "120%"
        )
    )
    print(dtTab)
    cat("\n\n")
}

CheckTabName <- function(tabNam) {
    if (grepl("Indels|INDEL", tabNam)) {
        cat(
            "<span style='color: red;'>",
            "All Variant calls have the following filters:",
            "</span>\n\n",
            "<span style='color: black;'>",
            "**Tumor freq >= 5%**, **Normal freq <2%**,",
            " and **Tumor Depth >50**",
            "</span>\n\n"
        )
    }
}

MakeTabColor <- function(tabNam) {
    color1 = "## <span style='color: blue;'> **"
    color2 = "## <span style='color: purple;'> **"
    if (stringr::str_detect(tabNam, "Philips") == T) {
        cat(paste0(color1, tabNam, "**</span> \n\n"))
    } else{
        cat(paste0(color2, tabNam, "**</span> \n\n"))
    }
}


ColorTable <- function(df) {
    df[] <- lapply(names(df), function(col) {
        if (col %in% c("Same")) {
            ifelse(
                df$Same == "No",
                kableExtra::cell_spec(df[[col]], "html", background = "#F06C78"),
                df[[col]]
            )
        } else {
            df[[col]]
        }
    })
    as.data.frame(df, stringsAsFactors = FALSE)
}


makeColorfulTab <- function(objDat) {
    cat(
        "<span style='color: red;'>",
        "All Variant calls have the following filters:",
        "</span>\n\n",
        "<span style='color: black;'>",
        "**Tumor freq >= 5%**, **Normal freq <2%**,",
        " and **Tumor Depth >50**",
        "</span>\n\n"
    )
    
    if (nrow(objDat) == 0) {
        cat("No calls after filtering")
    } else {
        objDat <- ColorTable(objDat)
        knitr::kable(objDat, "html", escape = FALSE, row.names = FALSE) %>%
            kableExtra::kable_styling(bootstrap_options = c("striped", "hover")) %>%
            kableExtra::row_spec(0, extra_css = "font-weight: bold;") %>%
            kableExtra::scroll_box(height = "100%", width = "100%") %>%
            print()
    }
    cat("\n\n")
}


MakeRegularTab <- function(tabNam, objDat) {
    if (!is.null(objDat) &
        stringr::str_detect(tabNam, "Methylation") == F) {
        if (tabNam == "QC") {
            makeQCTab(objDat)
        } else{
            CheckTabName(tabNam)
            if (tabNam == "In-House Somatic Variant Calls") {
                makeColorfulTab(objDat)
            } else{
                makeDefaultDt(objDat)
            }
            
        }
    } else{
        cat("\n\nNo additonal results for this case yet\n\n")
        cat("\n\n")
    }
}


# Prints out the table tab headings and tab information based on object params
makeDT <- function(tabNam,
                   objDat,
                   pdfFi = NULL,
                   rdNumb = NULL,
                   sam = NULL,
                   outDir = NULL) {
    MakeTabColor(tabNam)
    if (!is.null(pdfFi)) {
        return(makePdfTab(pdfFi, cnvTab = objDat, outDir))
    }
    if (!is.null(rdNumb)) {
        return(makeRdTab(rdNumb, sam))
    }
    if (stringr::str_detect(tabNam, "Info") == T) {
        return(makeSpecimenTab(objDat))
    } else{
        return(MakeRegularTab(tabNam, objDat))
    }
}


makeDT <- function(tabNam,
                   objDat,
                   pdfFi = NULL,
                   rdNumb = NULL,
                   sam = NULL,
                   outDir = NULL) {
    MakeTabColor(tabNam)
    if (!is.null(pdfFi)) {
        return(makePdfTab(pdfFi, cnvTab = objDat, outDir))
    }
    if (!is.null(rdNumb)) {
        return(makeRdTab(rdNumb, sam))
    }
    if (stringr::str_detect(tabNam, "Hotspots") == T) {
        return(PrintHotspotTable(objDat))
    }
    if (stringr::str_detect(tabNam, "Info") == T) {
        return(makeSpecimenTab(objDat))
    } else{
        return(MakeRegularTab(tabNam, objDat))
    }
}


MakeHStab <- function(sam, hsDat, samList, outDir) {
    ngsRows <- samList$Test_Number == sam
    ts_number <- samList$Specimen_ID[samList$Tumor_Content != 0 &
                                         ngsRows]
    hotspot_tsv <- file.path(outDir, 'hotspots', paste0(ts_number, "_Hotspots.tsv"))
    hsMain <- hsDat[hsDat$Test_Number == sam, ]
    names(hsMain)[names(hsMain) == 'Variant'] <- 'Common COSMIC Variant'
    makeDT(tabNam = "Hotspots", objDat = hsMain, NULL, NULL, sam)
}

# Adds MSI and TMB info under QC tab
make_msi_tmb_tab <- function(tmb_row, msi_row) {
    tmb_cols <- c("SampleID", "VariantCaller")
    msi_cols <- c("SampleID", "Status", "TumorType")
    
    idx2drop_tmb <- which(colnames(tmb_row) %in% tmb_cols) * -1
    
    newTa1 <- knitr::kable(tmb_row[, idx2drop_tmb], row.names = F, "html")
    tmbTable <- kableExtra::kable_styling(
        newTa1,
        bootstrap_options = c("condensed"),
        full_width = F,
        position = "left"
    )
    cat("\n*TMB DATA*\n")
    print(tmbTable)
    cat("\n\n")
    
    idx2drop_msi <- which(colnames(msi_row) %in% msi_cols) * -1
    
    newTa2 <- knitr::kable(msi_row[, idx2drop_msi], row.names = F, "html")
    msiTable <- kableExtra::kable_styling(
        newTa2,
        bootstrap_options = c("condensed"),
        full_width = F,
        position = "left"
    )
    cat("\n*MSI DATA*\n")
    print(msiTable)
    cat("\n\n")
}

# Parses the philips CNV abberations export from Philips dump csv file
parsePhilipsCn <- function(cnvInfo, sam, cnvTab) {
    cnvInfo$Test_Case <- sam
    cnvInfo$Variant <- "CNV"
    cnvInfo$Mutation.Type <- cnvInfo$AberrationType
    cnvInfo$Other <- cnvInfo$Chrom
    cnvInfo$Comments <- substr(cnvInfo$CopyNumber, 1, 3)
    cnvInfo$In.NYU <- cnvInfo$In.Philips <- NA
    cnvInfo$IGV <- "Copy Number"
    cnvTab <- cnvInfo
    dropMut <- cnvTab$Mutation.Type != "Hemizygous Loss"
    cnvTab <- cnvTab[dropMut, ]
    cnvTab <- cnvTab[cnvTab$Mutation.Type != "", ]
    cnvTab <- cnvTab[order(cnvTab$Mutation.Type), ]
    return(cnvTab)
}


logMissingData <- function(sam) {
    log_file <- paste(gb$pactName, "philips", "missing.txt", sep = "_")
    log_path <- file.path(Sys.getenv("HOME"), "Desktop", log_file)
    base::writeLines(sam, log_path)
}


parseCNV <- function(outPath, sam, cnvTab) {
    sam_folder <- get_ngs_path(outPath, sam)
    cnvFi <- file.path(outPath, sam, "aberration_cnv.csv")
    if (!is.null(sam_folder)) {
        cnvFi <- file.path(sam_folder, "aberration_cnv.csv")
    }
    
    if (length(cnvFi) > 1) {
    info <- file.info(cnvFi)
    newest <- cnvFi[which.max(info$mtime)]
    cnvFi <- newest
    }
    
    if (!file.exists(cnvFi)) {
        warning("File does not exist:\n", cnvFi, "\nSkipping Philips...")
        cnvTab$Gene <- "No Philips Data Dump CSV file availible, check ISPM for CNV"
        try(logMissingData(sam), silent = T)
        return(cnvTab)
    }
    cnvInfo <- as.data.frame(readr::read_csv(cnvFi, col_types = readr::cols(.default = readr::col_character())))
    if (nrow(cnvInfo) > 0) {
        cnvTab <- parsePhilipsCn(cnvInfo, sam, cnvTab)
    }
    return(cnvTab)
}


getDumpFiles <- function(outPath, sam, cnvTab, philipsFtp = "/Volumes/molecular/Molecular/Philips_SFTP") {
    sam_pat <- paste0(sam, "_")
    samFiles <- dir(
        path = philipsFtp,
        pattern = sam_pat,
        all.files = T,
        full.names = T
    )
    if (length(samFiles) > 1) {
        newest_file <- which.max(file.info(samFiles)$mtime)
        samFiles <- samFiles[newest_file]
    }
    if (length(samFiles) == 0) {
        message("No data dump found for pattern: ", sam_pat)
        return(cnvTab)
    }
    samZip <- basename(samFiles)
    dumpDir <- file.path(philipsFtp, samZip)
    destDir <- file.path(outPath, paste0(basename(samZip)))
    if (file.exists(dumpDir) & !file.exists(destDir)) {
        try(fs::file_copy(path = dumpDir, new_path = destDir),
            silent = T)
        unzip(zipfile = destDir, exdir = outPath)
    }
    cnvTab <- parseCNV(outPath, sam, cnvTab)
    return(cnvTab)
}

# Appends any additional cnv abberations from data dump
checkDataDump <- function(sam, cnvTab) {
    if (nrow(cnvTab) == 0) {
        cnvTab <- as.data.frame(matrix(
            nrow = 1,
            ncol = length(cnvTab),
            dimnames = list(1, names(cnvTab))
        ))
        cnvTab$Test_Case <- sam
    }
    outPath <- file.path(getwd(), "zipfiles")
    if (!dir.exists(outPath)) {
        dir.create(outPath)
    }
    cnvTab <- getDumpFiles(outPath, sam , cnvTab)
    drpCol <- which(stringr::str_detect(colnames(cnvTab), "In.|IGV"))
    cnvTab <- cnvTab[, -drpCol]
    names(cnvTab)[names(cnvTab) == 'Comments'] <- 'Copy Number'
    names(cnvTab)[names(cnvTab) == 'Other'] <- 'Chromosome'
    cnvTab$`Copy Number` <- as.numeric(cnvTab$`Copy Number`)
    cnvTab <- cnvTab[!is.na(cnvTab$Test_Case), ]
    return(cnvTab)
}


checkTumorPdf <-  function(samList, outDir) {
    pngOutDir <- file.path(outDir, "cnvpng") # output copy of cnvPNG files
    if (!dir.exists(pngOutDir)) {
        dir.create(pngOutDir)
    }
    pdfList <- list.files(
        outDir,
        pattern = "*.pdf",
        recursive = T,
        full.names = T
    )
    if (length(pdfList) == 0) {
        message("No PDFs found in current directory:\n", getwd())
        message(
            "Make sure this Runs Facet PDFs are availible in the directory:\n",
            "/Volumes/molecular/Molecular/REDCap/cnv_facets\n",
            "or copy facetpdf files to working directory"
        )
        stopifnot(length(pdfList) != 0)
    }
    toDrop <- samList$Paired_Normal != ""
    tumors <- samList[toDrop, ] # drop controls/normals
    row.names(tumors) <- 1:nrow(tumors) # re-number rows
    if (length(pdfList) != nrow(tumors)) {
        message(
            "Number of sample rows ",
            nrow(tumors),
            " does not equal length of pdf files: ",
            length(pdfList),
            "\nCheck if any are missing:\n"
        )
        message("Tumors:")
        print(tumors$Specimen_ID)
        message("Files:")
        print(pdfList)
        stopifnot(length(pdfList) == nrow(tumors))
    }
    return(tumors)
}

# Converts FACETS PDF cnv facets to PNG format renamed as NGS name-------------
convert.plots <- function(tumors, pdfList) {
    stopifnot(class(tumors) == "data.frame")
    ngsOrder <- unlist(lapply(
        X = tumors$Paired_Normal,
        FUN = function(X) {
            base::which(stringr::str_detect(pdfList, X))
        }
    ))
    ngsOrder = base::which(sapply(tumors$Specimen_ID, grepl, pdfList), arr.ind = T)[, "row"]
    for (X in 1:nrow(tumors)) {
        currCase <- tumors$Paired_Normal[X]
        fileFind <- stringr::str_detect(pdfList, pattern = currCase)
        pdfFile <- pdfList[fileFind]
        pngName <- paste0(tumors$Test_Number[X], ".png")
        message("PDF file: ", basename(pdfFile))
        message("Matching NGS: ",
                tumors$Test_Number[X],
                " ",
                tumors$Paired_Normal[X])
        outputFile <- file.path("cnvpng", pngName)
        if (!file.exists(outputFile)) {
            tryCatch(
                expr = {
                    suppressWarnings(pdftools::pdf_convert(
                        pdfFile,
                        filenames = outputFile,
                        dpi = 300
                    ))
                },
                error = function(e) {
                    message(e,
                            "\nTry running:\nbrew install fontconfig --universal")
                    system("brew update --auto-update")
                }
            )
        }
    }
}

# Copies Methylation CNV PNG files from Molecular to CBioinformatics
GetMethCnv <- function(params, methDir) {
    methData <- gb$GetMethDf(params$pactName)
    outDir <- file.path(params$workDir, paste0(params$pactName, "_consensus"))
    pngDir <- file.path(outDir, "methCNV") # output copy of meth png files
    methSamples <- paste(methData$record_id, "cnv.png", sep = "_")
    methSamples <- methSamples[methData$report_complete == "YES"]
    if (!dir.exists(pngDir)) {
        dir.create(pngDir)
    }
    if (length(methSamples) == 0) {
        return(message("No Methylation PNGs on this run!"))
    }
    to_copy <- methSamples[!file.exists(file.path(pngDir, methSamples))]
    
    if (length(to_copy) == 0)
        return(message("All files have already been copied."))
    
    for (mSam in to_copy) {
        methPath <- file.path(methDir, mSam)
        methOut <- file.path(pngDir, mSam)
        
        if (file.exists(methPath)) {
            tryCatch(
                fs::file_copy(methPath, methOut),
                error = function(e) {
                    message("Failed to copy: ", basename(methPath))
                    message("To: ", methOut, "\nError: ", e$message)
                }
            )
        } else {
            new_path <- file.path(fs::path_home(), "Desktop", basename(methPath))
            if (!file.exists(new_path)) {
                stop("PNG file not found:\n", methPath)
            } else {
                fs::file_copy(new_path, methOut)
                try(fs::file_copy(new_path, methPath), silent = T)
            }
        }
    }
}

# Checks if the facets pdfs have been converted to png ------------------------------------
CopyPdfsPngs <- function(params) {
    outDir <- file.path(params$workDir, paste0(params$pactName, "_consensus"))
    samList <- gb$GetSamList(params$pactName, 2)
    pdfDir <- file.path(outDir, "FACETpdfs") # input facet pdf directory
    tumors <- checkTumorPdf(samList, outDir)
    pdfList <- list.files(pattern = "*.pdf", recursive = T)
    pdfList <- pdfList[!grepl("^0_", pdfList)]
    if (length(pdfList) > 0) {
        convert.plots(tumors, pdfList)
        if (!dir.exists(pdfDir)) {
            dir.create(pdfDir)
        }
        try(fs::file_move(pdfList, pdfDir), silent = T)
    }
}


CheckHotspots <- function(params) {
    outDir <- file.path(params$workDir, paste0(params$pactName, "_consensus"))
    hotspot_tsv <- file.path(outDir, paste0(params$pactName, "_Hotspots.tsv"))
    if (file.exists(hotspot_tsv)) {
        message("Hotspots exist")
    } else{
        message("Hotspots Missing:\n", hotspot_tsv)
    }
}


# Generates the link to the BAM file based on PACT run ------------------------------------
makeBamLink <- function(sam, pactID) {
    samsheet <- list.files(".", "demux-samplesheet.csv", T)[1]
    samCSVfi <- as.data.frame(read.csv(samsheet, skip = 19))
    allNames <- samCSVfi[, 1]
    rowNGS <- which(samCSVfi$Test_Number == sam)
    bamFiNam <- paste0(allNames[rowNGS], ".dd.ra.rc.bam")
    bamFiUrl <- file.path("https://genome.med.nyu.edu/external/clinpathlab",
                          pactID,
                          bamFiNam)
    cat(paste0("```\nIGV Load from URL: ", bamFiUrl, "\n```\n\n"))
}


# Generates a new Sample Tabbed row in html ------------------------------------------------
makeNewTab <- function(sam, samList, qcData, pactID) {
    currSam <- samList[samList$Test_Number == sam, "Specimen_ID"]
    cat(' <div class="boxed"> ')
    cat(paste0("\n\n# **", sam, "** {.tabset}", "\n\n"))
    cat(paste0("(", currSam[1], ")\n\n"))
    makeBamLink(sam, pactID)
    currQC <- grepl(paste(currSam, collapse = "|"), qcData$Sample)
    qcTab <- qcData[currQC, 2:ncol(qcData)]
    if (nrow(qcTab) == 0) {
        warning(sam, " is missing from qcTsv file")
    } else{
        makeDT("QC", qcTab)
    }
}


GetMethMatchInfo <- function(methData, theRd) {
    rdDtRows <- methData[methData$record_id == theRd, ]
    cat("\n\n")
    newTa <- knitr::kable(rdDtRows, row.names = F, "html")
    newTa <- kableExtra::kable_styling(
        newTa,
        bootstrap_options = c("bordered", "condensed"),
        full_width = F,
        position = "left"
    )
    newTa <- kableExtra::column_spec(newTa, 1:5, width = "5cm")
    cat("#### Methylation Full Report Link:")
    cat("\n\n")
    txtLink <- paste0("[", theRd, "](", rdDtRows$`Report Path`, ")")
    cat(txtLink)
    cat("\n\n")
    pathLink <- file.path(
        stringr::str_split_fixed(dirname(rdDtRows$`Report Path`), "acc_pathology", 2)[1, 2],
        basename(rdDtRows$`Report Path`)
    )
    cat("The link should lead to the file here:", pathLink)
    cat("\n\n")
    print(newTa)
    cat("\n\n")
}


# Generates a new Sample Tabbed row in html  -------------------------------
makeMethTab <- function(sam, methCn, methData) {
    if (nrow(methCn) == 0) {
        methCn <- NULL
    }
    if (any(sam %in% methData$Test_Number) == F) {
        cat("## **No Methylation** \n\n")
        cat("\n\nNo additional data table for this sample tab\n\n")
    }
    if (any(sam %in% methData$Test_Number) == T) {
        currRD <- methData$record_id[methData$Test_Number == sam]
        passFail <- methData$qc_passed[methData$Test_Number == sam]
        if (length(currRD) > 1) {
            message(sam,
                    " has more than one RD-number:\n",
                    paste(currRD, sep = " ", collapse = " "))
            for (rd in 1:length(currRD)) {
                rdTab <- paste0("Methylation_", rd)
                if (!is.na(passFail[rd])) {
                    if (passFail[rd] == "No") {
                        rdTab <- paste0(rdTab, "_QC_FAILED")
                        curr_html <- methData$`Report Path`[methData$Test_Number == sam][rd]
                        new_html <- stringr::str_replace_all(curr_html,
                                                             ".html",
                                                             "_QC_FAILED.html")
                        methData$`Report Path`[methData$Test_Number == sam][rd] <- new_html
                        methData$`Report Link`[methData$Test_Number == sam][rd] <- basename(new_html)
                    }
                }
                makeDT(rdTab, methCn, pdfFi = NULL, currRD[rd], sam)
                GetMethMatchInfo(methData, currRD[rd])
            }
        } else{
            rdTab <- "Methylation"
            if (!is.na(passFail)) {
                if (passFail == "No") {
                    rdTab <- paste0(rdTab, "_QC_FAILED")
                    curr_html <- methData$`Report Path`[methData$Test_Number == sam]
                    new_html <- stringr::str_replace_all(curr_html, ".html", "_QC_FAILED.html")
                    methData$`Report Path`[methData$Test_Number == sam] <- new_html
                    methData$`Report Link`[methData$Test_Number == sam] <- basename(new_html)
                }
            }
            makeDT(rdTab, methCn, NULL, currRD, sam)
            GetMethMatchInfo(methData, currRD)
        }
    }
}


# Filters out Philips Abberations Tab -----------------------------
filterAbberations <- function(snvOut) {
    snvOut$tumor_freq <- readr::parse_number(snvOut$tumor_freq)
    snvOut <- snvOut[!is.na(snvOut$tumor_freq), ]
    snvOut <- snvOut[snvOut$tumor_freq >= 5, ]
    fltrStr <- stringr::str_detect(snvOut$AberrationType, pattern = "frame|delet|insert|indel")
    snvOut <- snvOut[fltrStr, ]
    drpFlt <- as.numeric(sub("%", "", snvOut$normal_freq)) < 2.0
    snvOut <- snvOut[drpFlt, ]
    drpFlt <- as.numeric(sub("%", "", snvOut$tumor_freq)) >= 5.0
    snvOut <- snvOut[drpFlt, ]
    drpFlt <- snvOut$tumor_dp >= 50
    snvOut <- snvOut[drpFlt, ]
    return(snvOut)
}


printSnvs <- function(snvCsv) {
    snvCsvCol <- c(
        "HGNC_gene",
        "AberrationType",
        "Coordinate",
        "tumor_freq",
        "normal_freq",
        "tumor_dp",
        "normal_dp",
        "THERAPY_AVAILABILITY",
        "HGVSp_Short",
        "SomaticStatus"
    )
    snvOutDf <- as.data.frame(read.delim(snvCsv, sep = ","))[, snvCsvCol]
    snvOut <- filterAbberations(snvOutDf)
    makeDT("Philips Filtered Indels", snvOut)
    cat("\n\n")
    return(snvOutDf)
}

printSpecInfo <- function(samCsv) {
    samOut <- as.data.frame(read.delim(samCsv, sep = ","))[, c(1, 3:13)]
    makeDT("Philips Specimen Info", samOut)
    cat("\n\n")
}

printDiagInfo <- function(diagCsv) {
    diagCsvCols <- c(
        "Diagnostic Order",
        "Orderer Name",
        "Orderer NPI",
        "Patient Name",
        "Gender",
        "Diagnosis Code",
        "Diagnosis",
        "Diagnosis for interpretation code",
        "Diagnosis for interpretation",
        "Tumor Stage"
    )
    samOut <- as.data.frame(read.csv(diagCsv, header = F))
    filterCol <- paste0(samOut[1, ]) %in% diagCsvCols
    samOut <- samOut[2, filterCol]
    colnames(samOut) <- diagCsvCols
    rownames(samOut) <- NULL
    makeDT("Philips Diagnostic Info", samOut)
    cat("\n\n")
}

# Renders error when csv file fails to be parsed ------------------------------------------------------------
PrintParseErr <- function(csvPath, tabTxt) {
    eTxt = "\n\n## **Failed to read Philips CSV File"
    snvErr = paste(eTxt, tabTxt, "**\n`", csvPath, "`\n\n")
    return(cat(snvErr))
}

# Makes Abberations Tab ------------------------------------------------------------
makeAbTab <- function(sam, philipsFtp = "/Volumes/molecular/Molecular/Philips_SFTP") {
    dumpDir <- get_ngs_path(outPath = file.path(getwd(), "zipfiles"), sam)
    
    if (is.null(dumpDir)) {
        td <- format(Sys.Date(), "%Y-%m-%d")
        zipFiN <- file.path(philipsFtp, paste0(sam, "_", td, ".zip"))
        cat("\n\n## **No Philips Data Dump**\n\n")
        cat("Data dump not found:\n")
        return(cat(paste0(zipFiN, "\n\n")))
    }
    
    snvCsv <- file.path(dumpDir, "aberration_snv.csv")
    
    if (length(snvCsv) > 1) {
    info <- file.info(snvCsv)
    newest <- snvCsv[which.max(info$mtime)]
    snvCsv <- newest
    }
    
    if (file.exists(snvCsv)) {
        philipsIndels <- tryCatch(
            printSnvs(snvCsv),
            error = function(e) {
                PrintParseErr(snvCsv, "SNV")
            }
        )
        return(philipsIndels)
    } else{
        return(NULL)
    }
}


makeBlankRow <- function(sam, snvDt) {
    nonMutant <- data.frame(matrix(ncol = ncol(snvDt), nrow = 1))
    colnames(nonMutant) <- colnames(snvDt)
    nonMutant[1, ] <- "None"
    nonMutant$Test_Case <- sam
    snvDt <- rbind(snvDt, nonMutant)
    return(snvDt)
}


GrabHotspots <- function(params) {
    outDir <- file.path(params$workDir, paste0(params$pactName, "_consensus"))
    hsTsv <- file.path(outDir, paste0(params$pactName, "_Hotspots.tsv"))
    if (file.exists(hsTsv)) {
        hsDat <- as.data.frame(read.csv(hsTsv, sep = '\t'))
    } else{
        message("No file found here:", hsTsv)
        hsDat <- NULL
    }
    return(hsDat)
}


CheckMissedSam <- function(sam, snvDt, pactName) {
    if (!(sam %in% snvDt$Test_Case)) {
        message(sam, " is missing from your Description File:")
        message(pactName, "_desc.csv")
        message("Adding sample as a blank row")
        snvDt <- makeBlankRow(sam, snvDt)
    }
    return(snvDt)
}


cleanTabCols <- function(snvTab) {
    snvTab <- subset(snvTab, select = -Variant)
    snvTab <- subset(snvTab, select = -IGV)
    snvTab <- subset(snvTab, select = -Comments)
    snvTab <- subset(snvTab, select = -MuTect2)
    snvTab <- subset(snvTab, select = -LoFreqSomatic)
    snvTab <- snvTab %>% rename(Position = Other)
    snvTab$Position <- gsub("_", ":", snvTab$Position)
    return(snvTab)
}


cleanPhilTab <- function(philipsIndels) {
    philipsIndels <- philipsIndels %>% rename(Position = Coordinate,
                                              Gene = HGNC_gene ,
                                              Mutation.Type = AberrationType)
    return(philipsIndels)
}


compare_philips <- function(snvTab, philipsIndels) {
    snvTab$Same <- "Yes"
    snvTab$In.Philips <- "No"
    
    genCols <- c("Gene", "Position")
    matched_rows <- logical(nrow(snvTab))
    
    for (inhouse in seq_len(nrow(snvTab))) {
        row <- snvTab[inhouse, ]
        house_gene <- row[["Gene"]]
        genes <- unlist(stringr::str_split(house_gene, " "))
        
        for (phil_row in seq_len(nrow(philipsIndels))) {
            indel_row <- philipsIndels[phil_row, ]
            if (any(sapply(genes, function(gene) all(c(gene, row["Position"]) == indel_row[genCols])))) {
                matched_rows[inhouse] <- TRUE
                break
            }
        }
    }
    
    snvTab$Same[!matched_rows] <- "No"
    snvTab$In.Philips[matched_rows] <- "Yes"
    combTab <- snvTab[order(snvTab$Same == "No", decreasing = TRUE), ]
    
    return(combTab)
}


MakeVAFtab <- function(sam) {
    if (dir.exists("VAF_plots")) {
        vaf_pngs <- dir(
            path = file.path(getwd(), "VAF_plots"),
            pattern = ".png",
            full.names = T
        )
        sam_find <- grepl(pattern = sam, x = vaf_pngs)
        
        if (any(sam_find)) {
            MakeTabColor("VAF Plots")
            cat(
                "Below are VAF distribution plots for MuTect2 and Strelka.\n\nBoth callers may not be plotted if one does not have enough data after filtering.\n\n"
            )
            foundPngs <- vaf_pngs[sam_find]
            
            for (vafPlot in foundPngs) {
                vaf_png_file <- file.path(".", "VAF_plots", paste0(basename(vafPlot)))
                caller <- stringr::str_split_fixed(basename(vafPlot), "_", 4)[1, 4]
                callerName <- stringr::str_replace_all(caller, ".png", "")
                altTxt <- paste("![VAF Distribution for", sam, callerName, "](")
                cat(paste0(altTxt, vaf_png_file, "){height=400px}"))
            }
        } else{
            cat("\n\n")
            #cat(paste(sam, "did not have enough VAF data to plot after filtering"))
        }
    } else{
        cat("\n\n")
        #cat("No VAF Plots were generated for this run")
    }
    cat("\n\n")
}


CreateVariantsTabs <- function(philipsIndels, snvTab) {
    snvTab <- cleanTabCols(snvTab)
    if (!is.null(philipsIndels)) {
        philipsIndels <- cleanPhilTab(philipsIndels)
        if (nrow(snvTab) == 0) {
            makeDT("In-House Somatic Variant Calls", objDat = snvTab)
        } else {
            combTab <- compare_philips(snvTab, philipsIndels)
            columns_to_front <- c("Same")
            combTab <- combTab[, c(columns_to_front, setdiff(names(combTab), columns_to_front))]
            makeDT("In-House Somatic Variant Calls", objDat = combTab)
        }
    } else{
        if (nrow(snvTab) > 0) {
            snvTab$In.Philips <- "NA"
            snvTab$Same <- "Pending"
            columns_to_front <- c("Same")
            snvTab <- snvTab[, c(columns_to_front, setdiff(names(snvTab), columns_to_front))]
        }
        makeDT("In-House Somatic Variant Calls", objDat = snvTab)
    }
}


makeOncoKbTab <- function(tabNam, objDat) {
    gb$MakeTabColor(tabNam)
  dtOpts <- list(
    info = F,
    autoWidth = F,
    lengthChange = F,
    searching = F,
    paging = FALSE
  )
    
    if (nrow(objDat) == 0) {
      cat("No matching OncoKB Level 1 genes called in-house")
    }else {
      dtTab <- htmltools::tagList(
          DT::datatable(
              objDat,
              style = "default",
              rownames = FALSE,
              options = dtOpts
          )
      )
      
      cat("OncoKB Level 1 genes called from In-House Somatic Variant tab:")
      cat("\n\n")
      print(dtTab)
    }
    cat("\n\n")
}


read_merge_metrics <- function(metrics_dir, fileType ) {
    filePatt <- paste0("_", fileType, ".*\\.csv$")
    normal_paths <- list.files(
        path      = metrics_dir,
        pattern   = filePatt,
        full.names = TRUE
    )
    normal_list <- lapply(normal_paths, function(filepath) {
        df <- read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)
        vec <- df[2]
        rownames(vec) <- df[[1]]
        return(vec)
    })
    normals_mat <- do.call(cbind, normal_list)
    normals_df <- as.data.frame(normals_mat, stringsAsFactors = FALSE)
    return(normals_df)
}


get_metrics_dir <- function(run_id) {
  runYear <- paste0("20", substring(run_id, 1, 2))
  MOL_DIR <- "/Volumes/molecular/Molecular/NGS607"
  OUTPUT_DIR <- file.path(MOL_DIR, runYear, run_id, "output")
  METRICS_DIR <- file.path(OUTPUT_DIR, "CollectHsMetrics")
  OUT_FOLDER <- file.path(OUTPUT_DIR, "HS_Metrics_CSV")
  return(OUT_FOLDER)
}


bold_center <- function(txt_in) {
  msg <- paste0(
    '<div style="text-align:left;">',
    '<span style="font-weight:bold; font-size:18px; ; color:black;">', txt_in, '</span>',
    '</div>\n\n'
  )
  return(cat(msg))
}


make_hs_metrics_tab <- function(sam, norm_metrics, tumor_metrics, tumorSams) {
  gb$MakeTabColor("HS Metrics")
  
  normals <- rownames(norm_metrics)
  tumors <- rownames(tumor_metrics)
  sam_match <- which(sam == tumorSams$Test_Number)  
  if (length(sam_match) > 0) {
    sam_row <- tumorSams[sam_match, ]

    tum_match <- which(stringr::str_detect(sam_row$Sample_ID, tumors))
    norm_match <- which(stringr::str_detect(sam_row$Paired_Normal, normals))

    norm_row <- norm_metrics[norm_match, ]
    norm_mean_sd <- norm_metrics[c("Mean:", "SD:"), ]
    norms <- rbind(norm_row, norm_mean_sd)
    
    old_cols <- colnames(norms)
    stringr::str_replace_all(old_cols, c("_" = " ", "normals" = "")) -> new_cols
    colnames(norms) <- new_cols
    
    tum_row <- tumor_metrics[tum_match, ]
    tumor_mean_sd <- tumor_metrics[c("Mean:", "SD:"), ]
    tums <- rbind(tum_row, tumor_mean_sd)
    
    old_cols <- colnames(tums)
    stringr::str_replace_all(old_cols, c("_" = " ", "tumors" = "")) -> new_cols
    colnames(tums) <- new_cols
    
    kableTab1 <- norms %>%
      knitr::kable(format = "html", row.names = TRUE) %>%
      kableExtra::kable_styling(full_width = FALSE) %>%
      kableExtra::row_spec(1, background = "lightblue")
    cat("\n\n")
    bold_center("Normal HS Metric Values")
    cat("\n\n"); print(kableTab1); cat("\n\n")

    kableTab2 <- tums %>%
      knitr::kable(format = "html", row.names = TRUE) %>%
      kableExtra::kable_styling(full_width = FALSE) %>%
      kableExtra::row_spec(1, background = "lightblue")

    cat('<hr style="border:none; ' , 'height:3px; ' ,
        'background-color:#000; ' , 'width:100%;" />\n')
      
    cat("\n\n")
    bold_center("Tumor HS Metric Values")
    cat("\n\n"); print(kableTab2); cat("\n\n")

  } else {
    cat("\n\n"); cat("No HS Metrics for this Sample"); cat("\n\n")
  }
}


LoopSampleTabs <- function(params) {
    pactName <- params$pactName
    methData <- gb$GetMethDf(params$pactName)
    qcData <- gb$ReadQcFile(pactName)
    samList <- gb$GetSamList(pactName, 2)
    toDrop <- grepl("^0_", samList$Paired_Normal)

    run_id <- read.csv("demux-samplesheet.csv", skip = 19)$Run_Number[1]
    metrics_dir <- get_metrics_dir(run_id)
    norm_metrics <- read_merge_metrics(metrics_dir, "normals")
    tumor_metrics <- read_merge_metrics(metrics_dir, "tumors")
    
    if (any(toDrop)) {
        ngs_drop <- samList$Test_Number[toDrop]
        toKeep <- !samList$Test_Number %in% ngs_drop
        samList <- samList[toKeep, ]
    }
    
    samsPaired <- names(table(samList$Test_Number))[table(samList$Test_Number) > 1]
    samList <- samList[samList$Test_Number %in% samsPaired, ]
    
    samples <- gb$GrabSamples(samList)
    hsDat <- gb$GrabHotspots(params)
    snvDt <- read.csv(paste0(pactName, "_desc.csv"))
    outDir <- file.path(params$workDir, paste0(pactName, "_consensus"))
    demux_sh <- as.data.frame(read.csv("demux-samplesheet.csv", skip = 19))
    toKeep <- demux_sh$Paired_Normal != ""
    tumorSams <- demux_sh[toKeep, ]
    rownames(tumorSams) <- NULL
    
    tmb_tsv <- "./TMB_MSI/annotations.paired.tmb.validation.2callers.tsv"
    msi_tsv <- "./TMB_MSI/msi_validation.tsv"
    onco_wsh <- "/Volumes/CBioinformatics/jonathan/pact/OncoKB-Level1-Genes-Feb2025.xlsx"
    
    tmb_dat <- as.data.frame(read.table(tmb_tsv, sep = "\t", header = T))
    msi_dat <- as.data.frame(read.table(msi_tsv, sep = "\t", header = T))
    onco_dat <- as.data.frame(readxl::read_xlsx(onco_wsh))
    
    onco_genes <- unique(onco_dat$Gene)
    
    for (sam in samples) {
        message("SAMPLE: ", sam)
        makeNewTab(sam, samList, qcData, pactName)
        sam_match <- which(sam == tumorSams$Test_Number)
        
        if (length(sam_match > 0)) {
            sam_id <- tumorSams$Sample_ID[sam_match]
            tmb_match <- tmb_dat$SampleID == sam_id
            msi_match <- msi_dat$SampleID == sam_id
            tmb_row <- tmb_dat[tmb_match, ]
            msi_row <- msi_dat[msi_match, ]
            make_msi_tmb_tab(tmb_row, msi_row)
        }
        
        snvDt <- CheckMissedSam(sam, snvDt, pactName)
        samRows <- snvDt$Test_Case == sam
        cnvTab <- snvDt[samRows & snvDt$Variant == "CNV", ]
        cnvTab <- checkDataDump(sam, cnvTab)
        
        snvTab <- snvDt[samRows & snvDt$Variant == "SNV", ]
        philipsIndels <- makeAbTab(sam)
        CreateVariantsTabs(philipsIndels, snvTab)
        
        cnvTab <- subset(cnvTab, select = -Variant)
        makeDT("CNV", cnvTab, pdfFi = sam, outDir = outDir)
        
        if (!is.null(methData)) {
            methCn <- methData[methData$Test_Number == sam, ]
            makeMethTab(sam, methCn, methData)
        }
        
        if (!is.null(hsDat)) {
            MakeHStab(sam, hsDat, samList, outDir)
        }
        cat("\n\n")
        MakeVAFtab(sam)
        cat("\n\n")
        
        onco_df <- as.data.frame(
          matrix(nrow = 0, ncol = ncol(onco_dat),
                 dimnames = list(NULL, colnames(onco_dat)))
          )
        
        if (any(snvTab$Gene %in% onco_genes)) {
          genes_found <- unique(snvTab$Gene[snvTab$Gene %in% onco_genes])
          onco_sub <- onco_dat[onco_dat$Gene %in% genes_found,]
          rownames(onco_sub) <- NULL
          onco_df <- onco_sub
        }
        
        makeOncoKbTab("OncoKB Level 1 Genes", onco_df)
        
        make_hs_metrics_tab(sam, norm_metrics, tumor_metrics, tumorSams)
        
        cat(' </div> ')
        
    }
}


loadHtmlTag <- function() {
    require("tidyverse")
    library("plotly")
    require("plotly")
    require("ggplot2")
    library("ggplot2")
    library("DT")
    htmltools::tagList()
    par(ask = F)
    devAskNewPage(ask = F)
    return(htmltools::tagList(
        DT::datatable(cars, options = list(
            pageLength = 100, width = "100%"
        )),
        plotly::ggplotly(ggplot2::ggplot())
    ))
}


RenamePngs <- function(tumors, pngList, ngsPngDir) {
    stopifnot(class(tumors) == "data.frame")
    ngsOrder <- base::which(sapply(tumors$Specimen_ID, grepl, pngList), arr.ind = T)[, "row"]
    for (X in 1:nrow(tumors)) {
        currCase <- tumors$Paired_Normal[X]
        fileFind <- stringr::str_detect(pngList, pattern = currCase)
        current_png <- pngList[fileFind]
        isInValid <- length(current_png) == 1
        if (!isInValid) {
            message("File png not found for sample: ", currCase)
            message(
                "Check Samplesheet CSV file and compare with the demux-samplesheet.csv used"
            )
        }
        stopifnot(length(current_png) == 1)
        new_pngName <- paste0(tumors$Test_Number[X], ".png")
        message("Current PNG file:\n", basename(current_png))
        message("Matching NGS name: ",
                tumors$Test_Number[X],
                "\n",
                tumors$Paired_Normal[X],
                "\n")
        renamedPng <- file.path(ngsPngDir, new_pngName)
        if (!file.exists(renamedPng)) {
            file.copy(current_png, renamedPng)
        }
    }
}


FixPngList <- function(tumors, pngList) {
    message("Mismatch in tumor samples and FACET PNG files!")
    message("Sample rows: ",
            nrow(tumors),
            ", PNG files: ",
            length(pngList))
    message("Check if any are missing:\n")
    message(paste0(capture.output(as.data.frame(
        tumors$Specimen_ID
    )), collapse = "\n"))
    message(paste0(capture.output(data.frame(
        PNG_Files = basename(pngList)
    )), collapse = "\n"))
    normLi <- paste(tumors$Paired_Normal, collapse = "|")
    toKeep <- unlist(lapply(
        pngList,
        FUN = function(X) {
            stringr::str_detect(X, pattern = normLi)
        }
    ))
    message("Dropping the following:")
    print(pngList[!toKeep])
    stopifnot(any(toKeep))
    pngList <- pngList[toKeep]
    stopifnot(length(pngList) == nrow(tumors))
    return(pngList)
}

ListNonMatching <- function(pngList, tumors, pngOutDir) {
    if (length(pngList) < nrow(tumors)) {
        nonContainedIDs <-
            sapply(tumors$Specimen_ID, function(id) {
                !any(grepl(id, basename(pngList)))
            })
        if (any(nonContainedIDs)) {
            nonMatchedIDs <- tumors$Specimen_ID[nonContainedIDs]
            message("The following samples are missing from: ", pngOutDir)
            message(paste0(capture.output(nonMatchedIDs), collapse = "\n"))
        }
        pngList <- FixPngList(tumors, pngList)
        return(pngList)
    } else{
        nonMatchingFiles <- sapply(basename(pngList), function(baseName) {
            !any(sapply(tumors$Specimen_ID, function(id) {
                grepl(id, baseName)
            }))
        })
        if (any(nonMatchingFiles)) {
            nonMatchedFiles <- pngList[nonMatchingFiles]
            message("The following PNGs are extra in: ", pngOutDir)
            message(paste0(capture.output(nonMatchedFiles), collapse = "\n"))
            return(pngList)
        } else{
            baseNames <- basename(pngList)
            matches <- expand.grid(
                Specimen_ID = tumors$Specimen_ID,
                FileName = baseNames,
                stringsAsFactors = FALSE
            )
            matched <- sapply(1:nrow(matches), function(i) {
                grepl(matches$Specimen_ID[i], matches$FileName[i])
            })
            duplicates <- matches[matched, ]
            duplicates$count <- ave(duplicates$Specimen_ID,
                                    duplicates$Specimen_ID,
                                    FUN = length)
            duplicateMatches <- duplicates[duplicates$count > 1, ]
            message(
                paste(
                    "The sample(s)",
                    paste(unique(
                        duplicateMatches$Specimen_ID
                    )),
                    "has duplicate matching CNV Facets PNG files:"
                )
            )
            message(paste(unique(duplicateMatches$FileName), collapse = "\n"))
            return(pngList)
        }
    }
}


CheckTumorPngs <- function(samList, outDir) {
    pngOutDir <- file.path(outDir, "cnvpng") # output copy of cnvPNG files
    if (!dir.exists(pngOutDir)) {
        dir.create(pngOutDir)
    }
    
    pngList <- list.files(
        pngOutDir,
        pattern = "*.png",
        recursive = T,
        full.names = T
    )
    pngList <- pngList[!grepl("^0_", basename(pngList))]
    
    if (length(pngList) == 0) {
        message("No CNV PNG Facets found in current directory:\n",
                pngOutDir)
        message(
            "Check if the Facet PNG files are available in the output directory:\n",
            "/Volumes/molecular/Molecular/REDCap/cnv_facets/ or the Desktop"
        )
        stopifnot(length(pngList) != 0)
    }
    
    tumors <- samList[samList$Paired_Normal != "", ] # drop controls/normals
    dupedNgs <- duplicated(tumors$Test_Number)
    
    if (any(dupedNgs)) {
        stop(c(
            "There are duplicated NGS numbers:\n",
            paste(tumors$Test_Number[dupedNgs], collapse = "\n")
        ))
    }
    row.names(tumors) <- 1:nrow(tumors) # correct row numbering
    
    if (length(pngList) != nrow(tumors)) {
        pngList <- ListNonMatching(pngList, tumors, pngOutDir)
    }
    
    ngsPngDir <- file.path(outDir, "ngs_cnv_png")
    if (!dir.exists(ngsPngDir)) {
        dir.create(ngsPngDir)
    }
    RenamePngs(tumors, pngList, ngsPngDir)
}


CopyCnvPngs <- function(params) {
    outDir <- file.path(params$workDir, paste0(params$pactName, "_consensus"))
    samList <- gb$GetSamList(params$pactName, 2)
    toDrop <- grepl("^0_", samList$Paired_Normal)
    if (any(toDrop)) {
        ngs_drop <- samList$Test_Number[toDrop]
        toKeep <- !samList$Test_Number %in% ngs_drop
        samList <- samList[toKeep, ]
    }
    CheckTumorPngs(samList, outDir)
}


animation::ani.options(autobrowse = FALSE)
options(width = 1600)
knitr::opts_chunk$set(echo = FALSE)

if (!require("xfun")) {
    install.packages("xfun", ask = FALSE)
}

knitr::opts_template$set(input_opts = list(
    echo = FALSE,
    message = TRUE,
    warning = FALSE,
    error = TRUE,
    include = FALSE
))

knitr::opts_template$set(
    tab_opts = list(
        echo = FALSE,
        warning = FALSE,
        message = FALSE,
        error = TRUE,
        fig.keep = 'all',
        include = TRUE,
        fig.align = "left",
        fig.show = "asis",
        out.width = '100%',
        results = 'asis'
    )
)

GetBamUrl <- function(params) {
    return(file.path("https://genome.med.nyu.edu/external/clinpathlab",
                     params$pactName))
}

gb$css_code <- base::readLines(
    "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts/consensus.css")
