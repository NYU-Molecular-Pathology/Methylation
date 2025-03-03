#!/usr/bin/env Rscript
## Script name: tableMaker.R
## Purpose: Source global parameters for research report tables and values
## Date Created: May 17, 2022
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2024

gb <- globalenv(); assign("gb", gb)

setDirectory <- function(foldr) {
    bsDir = paste("cd", foldr)
    mm2 = crayon::white$bgRed("Location Not Found:", foldr)
    if (dir.exists(foldr)) {
        system(bsDir); setwd(foldr); assign("runDir", foldr)
    } else{ warning(mm2) }
}


checkNeeded <- function(pkgList) {
    installedPkgs <- rownames(installed.packages())
    neededPkgs <- pkgList[!pkgList %in% installedPkgs]
    return(neededPkgs)
}


if (length(checkNeeded("librarian")) > 0) {
    install.packages(
        "librarian",
        dependencies = T,
        verbose = T,
        Ncpus = 4,
        quiet = T
    )
}

if (length(checkNeeded("pak")) > 0) {
    install.packages("pak", dependencies = T, ask = F, type = "binary")
}


# Load/install Packages
pkgs <- c(
    'knitr',
    'readr',
    'ggplot2',
    'devtools',
    "BiocManager",
    'parallel',
    'iterators',
    'doParallel',
    "foreach",
    "dplyr",
    "DT",
    "kableExtra",
    "foreach",
    "magrittr",
    "conumee",
    "enrichplot",
    "limma",
    "stringr",
    "tidyverse",
    "tidyr",
    "DOSE",
    "ggnewscale",
    "enrichplot",
    "ComplexHeatmap",
    "assertr",
    "Biobase",
    "irlba",
    "grid",
    "RColorBrewer",
    "Cairo",
    "gridExtra"
)

pkgs_needed <- checkNeeded(pkgs)

if (length(pkgs_needed) > 0) {
    pak::pkg_install(pkgs_needed, ask = F)
}


require("gridExtra")
require("Cairo")
require("minfi")
require("ComplexHeatmap")

EndDiv <- function(){
    cat("\n\n")
    cat(":::")
    cat("\n\n")
}

makeDt <- function(targets, gb = NULL, showHide = FALSE) {
    dtOpts <- list(
        columnDefs = list(list(className = 'dt-center', targets = "_all")),
        scrollX = T, scrollY = T, info = F, pageLength = 20, autoWidth = F,
        rownames = F, lengthChange = F, searchable = T
    )

    if (!is.null(gb)){
        colFilter <- c(gb$col_samNames, gb$col_samTypes, gb$col_shapes, "SentrixID_Pos")
        if (!is.null(gb$col_arrayType)){
            colFilter <- c(colFilter, gb$col_arrayType)
        }
        if (!is.null(gb$batch_col)){
            colFilter <- c(colFilter, gb$batch_col)
        }
        if (!is.null(gb$col_samGrp)){
            colFilter <- c(colFilter, gb$col_samGrp)
        }
        if (!is.null(gb$col_other)){
            colFilter <- c(colFilter, gb$col_other)
        }
        colFilter <- unique(colFilter)
        targets <- targets[, colFilter]
    }

    theDt <- DT::datatable(
        targets,
        selection = "single",
        autoHideNavigation = F,
        options = dtOpts,
        class = 'white-space: nowrap'
    )

    if(showHide){
        gb$MakeHideButton("SampleData")
        return(theDt)
    }else{
        return(theDt)
    }

}

smallTab <- function(dtObj) {
    dtTable <- knitr::kable(dtObj, row.names = F, "html")
    dtTable <- kableExtra::kable_styling(
        dtTable,
        bootstrap_options = c("striped", "condensed"),
        full_width = F,
        position = "left"
    )

    dtTable <- kableExtra::column_spec(dtTable, 1:ncol(dtObj), width = "4cm")
    return(dtTable)
}

# Adds colors to csv targets file to maintain same color scheme between plots
FixNullNaVars <- function(targets, varColumns) {
    if(all(varColumns %in% dimnames(targets)[[2]])==F){
        message("The column names you provided are not found in targets columns: ", paste(varColumns, collapse = ", "))
        stopifnot(all(varColumns %in% dimnames(targets)[[2]]))
    }
    for (variable in varColumns) {
        if (any(is.null(targets[, variable]))) {
            targets[is.null(targets[, variable]), variable] <- "NONE"
        }
        if (any(is.na(targets[, variable]))) {
            targets[is.na(targets[, variable]), variable] <- "NONE"
        }
        if (any(targets[, variable] == "")) {
            toFix <- targets[, variable] == ""
            targets[toFix, variable] <- "NONE"
        }
    }
    return(targets)
}


GetDefaultColors <- function(col_vect = NULL) {
    if (!is.null(col_vect)) {
        return(col_vect)
    }
    color_vect <- c(
        "#FF0733", # Vivid Red
        "#F76e1e", # Flame
        "#FFD300", # Cyber Yellow
        "#33CC33", # Lime Green
        "#004CFF", # Blue Ribbon
        "#FF00FF", # Magenta
        "#66EFFF", # Cyan
        "#8A2BE2", # Blue Violet
        "#5A5A5A", # Davy's Grey
        "#111789", # Ultra Blue
        "#C8FF00", # Electric Lime
        "#795a34", # Coffee
        "#FF93C7", # Pink Sherbet
        "#898A00", # Olive Drab
        "#858567", # Sage
        "#02AD24", # Green
        "#F8EF8C", # Banana Mania
        "#019965", # Green Sheen
        "#FF6D6D", # Bittersweet
        "#B1CC71", # Olive Green
        "#DC5E93", # Raspberry Pink
        "#ADD8E6", # Light Blue
        "#1D8932", # Forest Green
        "#766C95", # Mountbatten Pink
        "#FDBC7C", # Peach
        "#A8D9D0", # Light Cyan
        "#FFFF00", # Yellow
        "#73FFCC", # Aquamarine
        "#9B4513", # Saddle Brown
        "#91FF57", # Screamin' Green
        "#7F0000", # Maroon
        "#2B330B", # Pine Tree
        "#009FFF", # Vivid Sky Blue
        "#C6123E", # Ruby Red
        "#FFB79F", # Light Salmon Pink
        "#1F9698", # Blue-Green
        "#1C4D42", # Charleston Green
        "#4A6D40", # Fern Green
        "#FF1493", # Deep Pink
        "#BBCBB4", # Laurel Green
        "#FFA200"  # Orange
    )
    return(color_vect)
}


CheckColorCount <- function(varColumns, targets, col_vect) {
    for (vc in varColumns) {
        if (length(unique(targets[, vc])) > length(col_vect)) {
            warning(
                'Column "', vc, '" has more than 32 variables!',
                '\nYou need to input more colors in colorTargets(col_vect) or reduce number of variables!'
            )
        }
    }
}


GetColorVariables <- function(targets, varColumns, col_vect){
    CheckColorCount(varColumns, targets, col_vect)
    dat <- targets[, varColumns] # varColumns
    anno_df <- data.frame(dat)
    vars2Color <- as.list(lapply(dat, unique))
    colorValues <- lapply(vars2Color, function(x) {x = (col_vect)[1:(length(x))]})
    for (x in 1:length(vars2Color)) {
        for (varNum in 1:length(vars2Color[x])) {
            names(colorValues[x][[1]]) = c(vars2Color[x][[1]])
        }
    }
    return(colorValues)
}


FlipColorVector <- function(targets, colorColNames, col_vect){
    if (length(colorColNames) > 1) {
        for (varN in 1:length(colorColNames)) {
            if ((varN %% 2) == 0) {
                col2Change <- colorColNames[varN]
                oldColors2 <- unique(targets[, col2Change])
                col_vect <- rev(col_vect)
                newColAssign <- col_vect[1:length(oldColors2)]
                for (oldColN in 1:length(oldColors2)) {
                    oldCol <- oldColors2[oldColN]
                    newCol <- newColAssign[oldColN]
                    colorSelec <- which(targets[, col2Change] == oldCol)
                    targets[colorSelec, col2Change] <- newCol
                }
            }
        }
    }
    return(targets)

}

colorTargets <- function(targets, varColumns = c("Type","Origin"), col_vect = NULL){
    col_vect <- gb$GetDefaultColors(col_vect)
    targets <- gb$FixNullNaVars(targets, varColumns)

    hasType <- any("Type" %in% varColumns)
    if (length(varColumns) <= 1 & hasType == F) {
        targets$Type <- targets[, varColumns[1]]
        varColumns <- c(varColumns, "Type")
    }
    if (length(unique(varColumns)) == 1) {
        targets$NewCol <- targets[, varColumns[1]]
        varColumns <- c(varColumns, "NewCol")
    }
    colorColNames <- unlist(lapply(varColumns, paste0, "_color"))
    colorValues <- GetColorVariables(targets, varColumns, col_vect)
    targets$color <- NULL

    for (colorCol in colorColNames) {
        targets[, colorCol] <- NA
    }
    for (colNam in varColumns) {
        newColumnId <- paste0(colNam, "_color")
        varNameColVals <- names(colorValues[colNam][[1]])
        blankVals <- varNameColVals==""
        if(any(blankVals)){
            varNameColVals[blankVals] <- "NONE"
            names(colorValues[colNam][[1]])[names(colorValues[colNam][[1]])==""] <- "NONE"
        }

        for (samNam in names(colorValues[colNam][[1]])) {
            varToColor <- targets[, colNam] == samNam
            colorHex <- paste0(colorValues[[colNam]][[samNam]])
            targets[, newColumnId][varToColor] <- colorHex
        }
    }

    return(targets)
}


getColors <- function(samTypes) {
    library('RColorBrewer')
    qualCols = brewer.pal.info[brewer.pal.info$category == 'qual', ]
    col_vector = unlist(mapply(brewer.pal, qualCols$maxcolors, rownames(qualCols)))
    col_vector <- unique(col_vector)
    myColors = col_vector[1:length(unique(samTypes))]
    names(myColors) <- unique(samTypes)
    return(myColors)
}


StandardizeHeaders <- function(targets, samNames, sentrixs){
    # Create any missing header columns to standardize names
    targets[,"Barcode"] <- targets[,"Sentrix_ID"] <- targets[,"SentrixID_Pos"] <- NA
    targets[,"Sample_ID"] <- targets[,"Sample_Name"] <- targets[,"Sample_Group"] <- NA
    targets[,"Sample_ID"] <- targets[,"Sample_Name"] <- samNames
    targets[,"Barcode"] <- targets[,"Sentrix_ID"] <- targets[,"SentrixID_Pos"] <- sentrixs
    return(targets)
}


StripSheetSpaces <- function(samSh, samsheet, gb){
    colnames(samSh) <- gsub(pattern = " ", replacement = "_", colnames(samSh))
    samSh <- samSh %>% dplyr::mutate_all(stringr::str_replace_all, " ", "-")
    if(any(is.na(samSh))){
        samSh[is.na(samSh)] <- gb$blank_keywd
    }
    write.csv(samSh, samsheet, quote = F, row.names = F)
    targets <- read.csv(samsheet)
    if (class(targets) != "data.frame") {targets <- as.data.frame(targets)}
    return(targets)
}


ValidateColumns <- function(targets, gb) {
    if (is.null(gb$col_samTypes)) {
        if (is.null(gb$col_samGrp)) {
            targets$Unknown_Samples <- "Unknown"
            gb$col_samGrp <- gb$col_samTypes <- "Unknown_Samples"
        } else{
            gb$col_samTypes <- gb$col_samGrp
        }
    }
    if (gb$col_samNames %in% colnames(targets) == F) {
        stop("\nYour input col_samNames", ' "', gb$col_samNames, '" ',
             "was not found in the colnames(targets) available:\n",
             paste(colnames(targets), collapse = "\n"))
    }
    if (gb$col_samTypes %in% colnames(targets) == F){
        targets[, gb$col_samTypes] <- "Unknown"
    }
    return(targets)
}


ReadSheetType <- function(inputFi){
    if (stringr::str_detect(inputFi, ".xlsx")) {
        samSh <- readxl::read_excel(inputFi)
        samSh <- samSh %>% dplyr::mutate_all(stringr::str_replace_all, ",", "")
    } else{
        samSh <- read.csv(inputFi, strip.white = T)
    }
    return(samSh)
}


ValidateSentrix <- function(targets, gb) {
    if (is.null(gb$col_sentrix) | any(gb$col_sentrix %in% colnames(targets))) {
        message("Creating Sentrix ID column")
        gb$col_sentrix <- "SentrixID_Pos"
    }
    gb$GetCsvSheet(gb$needFi, gb$samsheet, gb$token, idatPath = gb$idatPath)
    targets <- gb$SetKeyColumns(
        targets,
        gb$col_samTypes,
        gb$col_samNames,
        gb$col_other,
        gb$col_shapes,
        gb$col_Grouping
    )
    targets <- FillMissingData(targets)
    return(targets)
}


MsgCaptureOut <- function(objDF){
    message(paste0(capture.output(objDF), collapse="\n"))
}


ValidateSampleIDs <- function(targets){
    if(any(duplicated(targets$Sample_ID))) {
        warning("Duplicated sample IDs will be dropped!")
        MsgCaptureOut(targets[duplicated(targets$Sample_ID), ])
        targets <- targets[!duplicated(targets$Sample_ID), ]
        row.names(targets) <- 1:nrow(targets)
    }
    return(targets)
}


sanitizeSheet <- function(inputFi, samsheet, gb) {
    library("magrittr")
    library("dplyr")
    samSh <- ReadSheetType(inputFi)
    targets <- StripSheetSpaces(samSh, samsheet, gb)
    targets <- ValidateColumns(targets, gb)
    targets <- ValidateSentrix(targets, gb)
    targets <- StandardizeHeaders(targets,
                                  samNames = targets[, gb$col_samNames],
                                  sentrixs = targets[, gb$col_sentrix]
                                  )
    targets <- ValidateSampleIDs(targets)
    write.csv(targets, file = samsheet, quote=F, row.names=F)
    return(targets)
}


#animation::ani.options(autobrowse = FALSE)
options(width=1200)
library("magrittr")
library("dplyr")
require("tidyverse")
require("plotly")
require("ggplot2")
library("DT")

printMissing <- function(theMissing, gb) {
    cat("#### Poor Quality or Dropped Samples:\n\n")
    if (!is.null(theMissing) & nrow(theMissing) > 0) {
        selCols <-
            c("Sample_Name", gb$col_sentrix, gb$col_samTypes, gb$col_samGrp, "Sample_Group")
        return(gb$smallTab(theMissing[, selCols]))
    } else{
        cat("None\n\n")
    }
}

loadHtmlTag <- function(){
    require("tidyverse")
    library("plotly")
    require("plotly")
    require("ggplot2")
    library("ggplot2")
    library("DT")
    htmltools::tagList(DT::datatable(
        cars, options = list(pageLength = 15, width = "100%")))
    par(ask = F)
    devAskNewPage(ask = F)
    doParallel::registerDoParallel(cores=2)
    return(htmltools::tagList(plotly::ggplotly(ggplot2::ggplot())))
}


sourceParams <- function(X = c("Params_input.R", "Params_output.R")) {
    paramFiles <- unlist(lapply(X, function(X) {file.path(".", "params", X)}))
    stopifnot(file.exists(paramFiles[1])|file.exists(paramFiles[2]))
    invisible(lapply(paramFiles, source))
}


GetCsvSheet <- function(needFi, samsheet, token, idatPath = NULL, outputFi = "samplesheet_og.csv"){
    if(is.null(idatPath)){
        idatPath <- file.path(getwd(), "idats")
    }
    # Using "pullRedcap_manual.R"
    rds <- gb$readInfo(inputSheet = samsheet) # inputSheet can be xlsx or csv
    stopifnot(length(rds) > 1)
    valid_rd <- stringr::str_detect(rds,"RD-")
    if(any(!valid_rd)){
        message("Some samples do not have valid RD-numbers and will be removed:")
        toDrop <- which(!valid_rd)
        message(capture.output(rds[toDrop]), collapse="\n")
        rds <- rds[-toDrop]
    }
    if(gb$needFi==T) {
        gb$grabRDCopyIdat(rd_numbers=rds, token, copyIdats=T, outputFi=outputFi)
        gb$MoveIdats()
    }else{
        result <- gb$search.redcap(rd_numbers=rds, token)
        result <- result[!is.na(result$barcode_and_row_column),]
        samplesheet_ID = as.data.frame(stringr::str_split_fixed(result[,"barcode_and_row_column"],"_",2))
        gb$makeSampleSheet(result, samplesheet_ID, bn = NULL, outputFi=outputFi)
    }
}


SubSetGroup <- function(strPatt, samPairs){
    idx <- which(grepl(samPairs, pattern = strPatt, ignore.case=T))
    return(samPairs[idx])
}


CheckTargetIdats <- function(targets) {
    stopifnot(gb$col_samNames %in% colnames(targets) == T)
    stopifnot(gb$col_samTypes %in% colnames(targets) == T)
    if (gb$needFi == T) {
        rds <-gb$readInfo(inputSheet = gb$samsheet) # inputSheet can be xlsx or csv
        stopifnot(length(rds) > 1 & stringr::str_detect(rds[1], "RD-"))
        gb$grabRDCopyIdat(rd_numbers = rds, gb$token,copyIdats = T)
        gb$MoveIdats()
    } else{
        targets <- as.data.frame(read.csv("samplesheet.csv"))
    }
    return(targets)
}


GetSamFreqTab <- function(targets, varCol1, varCol2 = NULL) {
    varColHead <- c("Sample Type", "Freq")
    sampleColumn <- targets[, varCol1]
    t1 <- setNames(as.data.frame(table(sampleColumn), row.names = NULL), varColHead)
    if (!is.null(varCol2)) {
        sampleColumn <- targets[, varCol2]
        t1 <- c(t1, setNames(as.data.frame(table(sampleColumn), row.names = NULL), varColHead))
    }
    return(list(t1))
}


MessageBatches <- function(targets, col_batchEffect) {
    if (!is.null(col_batchEffect)) {
        cat(paste0('#### Batch Effect Column Name: "', col_batchEffect, '" \n\n'))
        cat("#### Different Batches:\n")
        cat(paste(unique(targets[, col_batchEffect]), collapse = " & "))
        cat("\n\n")
    }else{
        cat("#### No Batch Effect Correction Performed\n\n")
    }
}

MessageArrayMix <- function(targets, col_arrayType) {
    if (!is.null(col_arrayType)) {
        cat(paste0('#### Arrays Type Column Name: "', col_arrayType, '" \n\n'))
        cat("#### Different Arrays: ")
        cat(paste(unique(targets[, col_arrayType]), collapse = " & "))
        cat("\n\n")
    }else{
        cat("No array probes merged, all microarrays are the same version\n\n")
    }
}


MessageBatchMix <- function(targets, gb){
    cat("\n\n")
    cat('# Data Quality and Analysis Parameters\n\n<hr style="border:4px solid darkgreen">\n\n')
    knitr::opts_chunk$set(out.width='50%')
    cat(paste("#### Samples Batch Corrected:", gb$batchEffect, "\n\n"))
    gb$MessageBatches(targets, gb$col_batchEffect)
    cat(paste("#### Sample 450k Probes Merged:", gb$mergeProbes, "\n\n"))
    gb$MessageArrayMix(targets, gb$col_arrayType)
    return(gb$printMissing(theMissing, gb))
}


PrintSamTypes <- function(targets, shCol) {
    if (!is.null(shCol)) {
        samF <- setNames(as.data.frame(table(targets[, shCol]), NULL), c(paste("Type:", shCol), "---Total"))
        return(samF)
    }
}


MsgSamFreq <- function(targets, gb){
    plotLi <- gb$PrintSamTypes(targets, gb$col_samTypes)
    if(gb$col_samGrp != gb$col_samTypes){
        plotLi <- list("Sample Types" = plotLi, "Sample Groups" = gb$PrintSamTypes(targets, gb$col_samGrp))
    }
    if (!is.null(gb$col_arrayType)){
        plotLi <- list(plotLi, "Array Types"= gb$PrintSamTypes(targets, gb$col_arrayType))
    }
    return(
        knitr::kable(plotLi, format = "html", table.attr = 'class="myTable"') %>%
            kableExtra::kable_styling(full_width = F, bootstrap_options = "bordered", position = "left")
    )
}


MsgNoobInfo <- function(getNoob) {
    if (getNoob == T) {
        cat(
            "<p>Noob (normal-exponential out-of-band) is a background correction method with dye-bias normalization for Illumina Infinium methylation arrays.</p><br />
<p>Functional normalization (FunNorm) is a between-array normalization method for the Illumina Infinium HumanMethylation450 platform. It removes unwanted variation by regressing out variability explained by the control probes present on the array.</p>"
        )
    }
}


if(Sys.info()[['sysname']]!="Darwin") {
    if (!require("rprofile")){remotes::install_github("csgillespie/rprofile", upgrade="never")}
    rprofile::set_startup_options(show.signif.stars = FALSE, useFancyQuotes = FALSE, Ncpus = parallel::detectCores()-2)
    Sys.setenv("R_PROFILE"=file.path(Sys.getenv("HOME"), "Rprofile.site"))
    magickPath <- paste0(system("which convert", intern=T))
    Sys.setenv(IMAGEMAGICK_V6_HOME=magickPath)
}


GetColorShape <- function(var1Col, var2Col){
    if(is.null(var2Col)){
        var2Col <- "None"
    }
    return(paste0("Color Labels = ", var1Col, " & ", "Shapes = ", var2Col))
}


ShowAnyMissed <- function(gb){
    cat("\n\n")
    cat("#### Samples Removed from Analysis with Missing or Duplicate idat files:\n\n")
    oldTargFile <- file.path(gb$runDir, "csv", "oldTargs.csv")
    oldTargs <- sanitizeSheet(gb$inputFi, oldTargFile, gb)
    oldTargs <- oldTargs[!c(oldTargs[,gb$col_samNames] %in% targets[,gb$col_samNames]),]
    if(nrow(oldTargs)>0){
        return(gb$smallTab(oldTargs))
    }else{
        return(cat("NONE\n\n"))
    }
}


PrintTableLists <- function(ta1, ta2, ...) {
    cat("\n\n")
    print(
        knitr::kable(
            list(ta1, ta2, ...), format = "html", table.attr = 'class="myTable"') %>%
            kableExtra::kable_styling(full_width = F)
    )
    cat("\n\n")
}


PrintPathHeader <-function(GroupSubset){
    cat("\n\n")
    cat("# Pathway Analysis")
    cat("\n\n")
    cat('<button class="btn btn-primary" data-toggle="collapse" data-target="#pathways">Show/Hide</button>')
    cat("\n\n")
    cat('::: {#pathways .collapse style="width:100%; height:auto; margin: 0px;"}')
    cat("\n\n")
    cat(paste("Cluster Condition: Tumor Type is", GroupSubset, "vs All\n\n"))
    cat('<p class="ncgnote"></p>')
    cat("\n\n")
    cat("## Cluster Profiler {.tabset}")
    cat("\n\n")
}


GetCatHeader <- function(tsneHead = NULL, tbset = T, preFix = "Unsupervised") {
    knitr::opts_chunk$set(out.width = '100%')
    cat('\n\n')
    if (is.null(tsneHead)) {
        tsneHead <- paste("All Groups", preFix , "TSNE (No Shapes/Symbols)")
    }
    if (tbset == T) {
        cat(paste('##', tsneHead, '{.tabset}', '\n\n'))
    } else{
        cat(paste('#', tsneHead, '\n\n'))
    }
}


MakeHideButton <- function(button_name){
    cat("\n\n")
    cat(paste0('<button class="btn btn-primary" data-toggle="collapse" data-target="#',
               button_name,'">Show/Hide</button>'))
    cat("\n\n")
    cat(paste0('::: {#', button_name,' .collapse style="width:100%; height:auto; margin: 0px;"}'))
    cat("\n\n")
}

MessageBatchMix <- function(targets, gb){
    cat("\n\n")
    cat('# Data Quality and Analysis Parameters\n\n')
    MakeHideButton("BatchMix")
    knitr::opts_chunk$set(out.width='50%')
    cat(paste("#### Samples Batch Corrected:", gb$batchEffect, "\n\n"))
    gb$MessageBatches(targets, gb$col_batchEffect)
    cat(paste("#### Sample 450k Probes Merged:", gb$mergeProbes, "\n\n"))
    gb$MessageArrayMix(targets, gb$col_arrayType)
    return(gb$printMissing(theMissing, gb))
}

MsgNullVar <- function(gb, colTxt, strTxt, varTxt) {
    if(is.null(gb[[colTxt]])){
        errMsg <- paste0(
            "\n>> The variable ", colTxt,' cannot be NULL', "\n",
            'Assign "', colTxt, '" to the name of column(s) in your targets/samplesheet dataframe that contains:\n',
            strTxt, '. For example, gb$', colTxt, ' <- ', varTxt
        )
        stop(errMsg)
    }
}

ShowKnitProgress <- function() {
    library("knitr")
    library("knitrProgressBar")
    library("progressr")
    knitr::opts_knit$set(progress = TRUE, verbose = TRUE)
    knitr::opts_hooks$set(
        log = function(options) {
            parsed <- parse(text = options$code, keep.source = FALSE)
            newCode <- sapply(as.character(parsed), function(x) {
                return(c('writeLog(sprintf("%s - format(Sys.time(), "%Y-%m-%d %H:%M:%S")), sep = NULL)',
                         sprintf("writeLog('%s')", x), x)
                )
            })
            options$code = as.vector(newCode)
            return(options)
        }
    )

    knitr::knit_hooks$set(time_it = local({
        now <- NULL
        function(before, options) {
            if (before) {
                # record the current time before each chunk
                now <<- Sys.time()
            } else {
                # calculate the time difference after a chunk
                res <- difftime(Sys.time(), now)
                # return a character string to show the time
                paste("Time for this code chunk to run:", res)
            }
        }
    }))

    options(
        knitr.progress.fun = function(total, labels) {
            p <- progressr::progressor(100, on_exit = FALSE)
            list(
                update = function(i) p(sprintf("chunk: %s", labels[i])),
                done = function() p(type = "finish")
            )
        }
    )
}

CatShapeColor <- function(var1Col = NULL, var2Col = NULL, tbset = TRUE, preFix = "Unsupervised") {
    knitr::opts_chunk$set(out.width = '100%')
    cat('\n\n')
    if (is.null(var2Col)) {
        var2Col <- "None"
    }
    if (tbset) {
        tsneHead <- paste0("Color Labels = ", var1Col, " & ", "Shapes = ", var2Col)
        cat(paste('##', tsneHead, '{.tabset}', '\n\n'))
    } else {
        tsneHead <- paste("All Groups", preFix, "TSNE (No Shapes/Symbols)")
        cat(paste('#', tsneHead, '\n\n'))
    }
}
