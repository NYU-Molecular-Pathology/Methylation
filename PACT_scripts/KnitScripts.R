gb <- globalenv(); assign("gb", gb)

# Librarian shelf loads or installs package from CRAN, BioConductor, & GitHub
if(!require("librarian")){install.packages("librarian", dependencies=T, verbose=T, Ncpus = 4, quiet=T)}
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
librarian::shelf(pkgs, ask=F)

library("kableExtra")

CheckMntDirs <- function(critialMnts, params) {
    outDir <- file.path(params$workDir, paste0(params$pactName,"_consensus"))
    failMount <-
        lapply(critialMnts, function(mnt) {
            ifelse(!dir.exists(mnt), return(T), return(F))
        })
    if (any(failMount == T)) {
        toFix <- paste(critialMnts[which(failMount == T)])
        cat("PATH does not exist, ensure drive is mounted:",
            toFix, sep = "\n")
        message(paste(
            "\nTry executing:\nsystem('diskutil mountdisk",toFix,"')"
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

StopMissingFile <- function(fiNam, isDir = F){
    if(isDir==F){
    if(!file.exists(fiNam)){
        stop(fiNam, " does not exist in the working directory!")
    }}else{
        if(!dir.exists(fiNam)){
            stop(fiNam, " does not exist as a directory! Check directory path name.")
        }
    }
}


CheckFiExist <- function(pactName, philipsFtp) {
    methSheet <- paste0(pactName, "_MethylMatch.xlsx")
    qcTsv <- paste0(pactName, "-QC.tsv")
    descrip <- paste0(pactName, "_desc.csv")
    samsheet <- list.files('.', "-SampleSheet.csv", T)[1]
    stopifnot(file.exists(samsheet) & file.exists(methSheet))
    StopMissingFile(samsheet)
    StopMissingFile(methSheet)
    StopMissingFile(qcTsv)
    StopMissingFile(descrip)
    StopMissingFile(philipsFtp, T)
}


GetMethDf <- function(pactName) {
    methSheet <- paste0(pactName, "_MethylMatch.xlsx")
    methData <- as.data.frame(readxl::read_excel(methSheet))
    methData <- methData[methData$report_complete == "YES", ]
    return(methData)
}

GetSamList <- function(pactName) {
    colFltr <- c("Test_Number","Specimen_ID","Tumor_Content")
    samsheet <- list.files('.', "-SampleSheet.csv", T)[1]
    message("Reading file: ", samsheet)
    samList <- read.csv(samsheet, skip=19)[, colFltr]
    isNGS <- !is.na(stringr::str_extract(samList$Test_Number, 'NGS'))
    toKeep <- samList$Test_Number!="0" & isNGS
    if(any(toKeep==F)){
        message(crayon::bgRed("The following cases have no NGS number and are being excluded:"),"\n")
        message(paste0(utils::capture.output(samList[!toKeep,]), collapse = "\n"))
        samList <- samList[toKeep,]
    }
    stopifnot(nrow(samList)>2)
    return(samList)
}

GrabSamples <- function(samList){
    sam <- unique(samList$Test_Number)
    samples <- sam[sam!=0 & !is.na(sam) & sam!="0"]
    return(samples)
}

ReadQcFile <- function(pactName){
    qcTsv <- paste0(pactName, "-QC.tsv")
    qcData <- read.delim(file=qcTsv)
    row.names(qcData) <- qcData$Row
    return(qcData)
}

makePdfTab <- function(pdfFi, outDir) {
    pdfFi <- file.path(outDir, "cnvpng", paste0(pdfFi, ".png"))
    cat("<span style='color: red;'>",
        "**Note**: CNV Amplifications are listed first,",
        "followed other Philips CNV abberations in table below this image\n\n"
    )
    cat(paste0("![CNV Plot](", pdfFi, "){width=100%}"))
    cat("\n\n")
}

makeRdTab <- function(rdNumb, sam) {
    cnvFi <- file.path(".","methCNV", paste0(rdNumb, "_cnv.png"))
    altTxt <- paste("![Methylation",rdNumb,sam,"CNV Plot](")
    cat(paste0(altTxt,cnvFi,"){width=120%}"))
    cat("\n\n")
}

makeSpecimenTab <- function(objDat) {
    newTa <- knitr::kable(objDat, row.names = F, "html")
    newTa <- kableExtra::kable_styling(
        newTa, bootstrap_options = c("condensed"), full_width = T, position = "left")
    newTa <- kableExtra::column_spec(newTa, 1:5, width = "2cm")
    print(newTa)
    cat("\n\n")
}

makeQCTab <- function(objDat){
    qcHeader <- c(
        "Sample", "Mapped Reads QC", "Deduplicated Reads QC",
        "Targets <50 Coverage", "Avg & Med Coverage", "Overlapped HOMO SNP"
    )
    colnames(objDat) <- qcHeader
    newTa <- knitr::kable(objDat, row.names = F, "html")
    newTa <- kableExtra::kable_styling(
        newTa, bootstrap_options = c("condensed"), full_width = F, position = "left"
    )
    newTa <- kableExtra::column_spec(newTa, 1:5, width = "5cm")
    print(newTa)
    cat("\n\n")
}

makeDefaultDt <- function(objDat) {
    dtOpts <- list(
        scrollX = T, scrollY=T, info = F, autoWidth = F, pageLength = 100,
        rownames=T, lengthChange = T, searchable = T
    )
    dtTab <- htmltools::tagList(DT::datatable(
        objDat, style="bootstrap", rownames = F,
        options = dtOpts, height = "120%", width = "120%"
    ))
    print(dtTab)
    cat("\n\n")
}

CheckTabName <- function(tabNam) {
    if (grepl("Indels|INDEL", tabNam)) {
        cat(
            "<span style='color: red;'>",
            "**Note** INDEL/Frameshift calls have the following filters:",
            "**Tumor freq >= 5%**, **Normal freq <2%**,",
            " and **Tumor Depth >50**</span> \n\n"
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

MakeRegularTab <- function(tabNam, objDat) {
    if (!is.null(objDat) & tabNam != "Methylation") {
        if (tabNam == "QC") {
            makeQCTab(objDat)
        } else{
            CheckTabName(tabNam)
            makeDefaultDt(objDat)
        }
    } else{
        cat("\n\nNo additonal results for this sample yet\n\n")
        cat("\n\n")
    }
}

# Prints out the table tab headings and tab information based on object parameters provided
makeDT <- function(tabNam, objDat, pdfFi = NULL, rdNumb = NULL, sam = NULL, outDir=NULL){
    MakeTabColor(tabNam)
    if (!is.null(pdfFi)) {
        makePdfTab(pdfFi, outDir)
    }
    if (!is.null(rdNumb)) {
        makeRdTab(rdNumb, sam)
    }
    if (stringr::str_detect(tabNam, "Info") == T) {
        makeSpecimenTab(objDat)
    } else{
        MakeRegularTab(tabNam, objDat)
    }
}

# Parses the philips CNV abberations export from Philips dump csv file
parsePhilipsCn <- function(cnvInfo, sam, cnvTab){
    cnvInfo$Test_Case <- sam
    cnvInfo$Variant <- "CNV"
    cnvInfo$Mutation.Type <- cnvInfo$AberrationType
    cnvInfo$Other <- cnvInfo$Chrom
    cnvInfo$Comments <- substr(cnvInfo$CopyNumber,1,3)
    cnvInfo$In.NYU <- cnvInfo$In.Philips <- NA
    cnvInfo$IGV <- "Copy Number"
    newTab <- cnvInfo[,names(cnvTab)]
    newTab <- rbind(cnvTab,newTab)
    cnvTab <- newTab
    dropMut <- cnvTab$Mutation.Type!="Hemizygous Loss"
    cnvTab <- cnvTab[dropMut,]
    cnvTab <- cnvTab[cnvTab$Mutation.Type!="",]
    return(cnvTab[order(cnvTab$Mutation.Type),])
}

parseCNV <- function(outPath, sam, cnvTab){
    cnvFi = file.path(outPath, sam, "aberration_cnv.csv")
    cnvInfo <- as.data.frame(readr::read_csv(
        cnvFi, col_types = readr::cols(.default = readr::col_character())
    ))
    if (nrow(cnvInfo) > 0) {
        cnvTab <- parsePhilipsCn(cnvInfo, sam, cnvTab)
    }
    return(cnvTab)
}

getDumpFiles <- function(outPath, sam,
                         cnvTab,
                         philipsFtp = "/Volumes/molecular/Molecular/Philips_SFTP") {
    samZip <- paste0(sam, ".zip")
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
    outPath <- file.path(getwd(), "zipfiles")
    if (!dir.exists(outPath)) {dir.create(outPath)}
    cnvTab<- getDumpFiles(outPath, sam , cnvTab)
    drpCol <- which(stringr::str_detect(colnames(cnvTab), "In."))
    cnvTab <- cnvTab[,-drpCol]
    return(cnvTab)
}

checkTumorPdf <- function(samList, outDir){
    pngOutDir <- file.path(outDir,"cnvpng") # output copy of cnvPNG files
    if(!dir.exists(pngOutDir)){
        dir.create(pngOutDir)
    }
    pdfList <- list.files(outDir, pattern = "*.pdf", recursive = T, full.names=T)
    if(length(pdfList)==0){
        message("No PDFs found in current directory:\n",getwd())
        message(
            "Make sure this Runs Facet PDFs are availible in the directory:\n",
            "/Volumes/molecular/Molecular/REDCap/cnv_facets\n",
            "or copy facetpdf files to working directory")
        stopifnot(length(pdfList)!=0)
    }
    toDrop <- samList$Tumor_Content != 0 & samList$Specimen_ID != "SC"
    tumors <- samList[toDrop, ] # drop controls/normals
    row.names(tumors)<- 1:nrow(tumors) # re-number rows
    if(length(pdfList)!=nrow(tumors)){
        message(
            "Number of sample rows ", nrow(tumors),
            " does not equal length of pdf files: ", length(pdfList),
            "\nCheck if any are missing:\n"
        )
        message("Tumors:")
        print(tumors$Specimen_ID)
        message("Files:")
        print(pdfList)
    }
    return(tumors)
}

# Converts FACETS PDF cnv facets to PNG format renamed as NGS name----------------------------------
convert.plots <- function(tumors, pdfList) {
    stopifnot(class(tumors)=="data.frame")
    ngsOrder = base::which(sapply(tumors$Specimen_ID, grepl, pdfList), arr.ind = T)[, "row"]
    invisible(lapply(X=1:length(ngsOrder), FUN=function(X){
        pngName <- paste0(tumors$Test_Number[X], ".png")
        pdfFile <- pdfList[ngsOrder[X]]
        tryCatch(
            expr = {suppressWarnings(pdftools::pdf_convert(
                pdfFile, filenames = file.path("cnvpng", pngName),dpi = 300))
            },
            error = function(e) {
                message("There was an error:\n", e,
                        "\nTry running:\nbrew install fontconfig --universal")
            }
        )}))
}

GetMethCnv <- function(params, methDir){
    methData <- gb$GetMethDf(params$pactName)
    outDir <- file.path(params$workDir, paste0(params$pactName,"_consensus"))
    pngDir <- file.path(outDir,"methCNV") # output copy of methylation png files
    methSamples <- paste(methData$record_id,"cnv.png",sep = "_")
    methSamples <- methSamples[methData$report_complete=="YES"]
    if(!dir.exists(pngDir)){dir.create(pngDir)}
    for (mSam in methSamples) {
        methPath <- file.path(methDir, mSam)
        methOut <- file.path(pngDir, mSam)
        if (!file.exists(methOut)) {
            message("Now trying to copy cnv file to...\n", methOut)
            if (file.exists(methPath)) {
                tryCatch(
                fs::file_copy(path = methPath, new_path = file.path(pngDir, mSam)),
                error= function(e){
                    warning("Methylation png failed to copy:\n", "From: ", methPath, "\nTo: ", methOut)
                    message(e)
                }
                )
            } else{
                warning("PNG file not found in Z-drive path:\n", methPath)
            }
        } else{
            message("File already copied to working directory: ", mSam)
        }
    }
}

# Checks if the facets pdfs have been converted to png
CopyPdfsPngs <- function(params) {
    outDir <- file.path(params$workDir, paste0(params$pactName,"_consensus"))
    samList <- gb$GetSamList(params$pactName)
    pdfDir <- file.path(outDir,"FACETpdfs") # input facet pdf directory
    tumors <- checkTumorPdf(samList, outDir)
    pdfList <- list.files(pattern = "*.pdf", recursive = T)
    if (length(pdfList)>0) {
        convert.plots(tumors, pdfList)
        dir.create(pdfDir)
        try(fs::file_move(pdfList, pdfDir), silent = T)
    }
}

# Generates a new Sample Tabbed row in html ------------------------------------
makeNewTab <- function(sam, samList, qcData){
    currSam <- samList[samList$Test_Number==sam,"Specimen_ID"]
    cat(paste0("\n\n# **", sam, "** {.tabset}","\n\n"))
    cat(paste0("(",currSam[1],")\n\n"))
    currQC <- grepl(paste(currSam, collapse="|"), qcData$Sample)
    qcTab <- qcData[currQC,2:ncol(qcData)]
    if (nrow(qcTab) == 0) {
        warning(sam, " is missing from qcTsv file")
    } else{
        makeDT("QC", qcTab)
    }
}

# Generates a new Sample Tabbed row in html  -------------------------------
makeMethTab <- function(sam, methCn, methData) {
    methCn <- if(nrow(methCn) == 0){NULL}
    if (any(sam %in% methData$Test_Number) == TRUE) {
        currRD <- methData$record_id[methData$Test_Number == sam]
        if (length(currRD) > 1) {
            message(
                sam, " has more than one methylation RD-number: ", paste(currRD, sep=", ")
            )
            for (rd in 1:length(currRD)) {
                rdTab <- paste0("Methylation", rd)
                makeDT(rdTab, methCn, pdfFi = NULL, currRD[rd], sam)
            }
        } else{
            makeDT("Methylation", methCn, NULL, currRD, sam)
        }
    } else{
        cat("## **No Methylation** \n\n")
        cat("\n\nNo additional data table for this sample tab\n\n")
    }
}

# Filters out Philips Abberations Tab -----------------------------
filterAbberations <- function(snvOut){
    snvOut$tumor_freq <- readr::parse_number(snvOut$tumor_freq)
    snvOut <- snvOut[!is.na(snvOut$tumor_freq),]
    snvOut <- snvOut[snvOut$tumor_freq>=5,]
    fltrStr <- stringr::str_detect(
        snvOut$AberrationType, pattern = "frame|delet|insert|indel")
    snvOut <- snvOut[fltrStr, ]
    drpFlt <- as.numeric(sub("%","",snvOut$normal_freq)) < 2.0
    snvOut <- snvOut[drpFlt,]
    drpFlt <- as.numeric(sub("%","",snvOut$tumor_freq)) >= 5.0
    snvOut <- snvOut[drpFlt,]
    drpFlt <- snvOut$tumor_dp >= 50
    snvOut <- snvOut[drpFlt,]
    return(snvOut)
}

printSnvs <- function(snvCsv){
    snvCsvCol <- c(
        "HGNC_gene", "AberrationType", "Coordinate",
        "tumor_freq", "normal_freq", "tumor_dp", "normal_dp",
        "THERAPY_AVAILABILITY", "HGVSp_Short", "SomaticStatus"
    )
    snvOutDf <- as.data.frame(read.delim(snvCsv, sep = ","))[, snvCsvCol]
    snvOut <- filterAbberations(snvOutDf)
    makeDT("Philips Filtered Indels", snvOut)
    cat("\n\n")
}

printSpecInfo <- function(samCsv) {
    samOut <- as.data.frame(read.delim(samCsv, sep = ","))[, c(1, 3:13)]
    makeDT("Philips Specimen Info", samOut)
    cat("\n\n")
}

printDiagInfo <- function(diagCsv){
    diagCsvCols <- c(
        "Diagnostic Order", "Orderer Name", "Orderer NPI",
        "Patient Name", "Gender", "Diagnosis Code", "Diagnosis",
        "Diagnosis for interpretation code",
        "Diagnosis for interpretation", "Tumor Stage"
    )
    samOut <- as.data.frame(read.csv(diagCsv, header = F))
    filterCol <- paste0(samOut[1, ]) %in% diagCsvCols
    samOut <- samOut[2, filterCol]
    colnames(samOut) <- diagCsvCols
    rownames(samOut) <- NULL
    makeDT("Philips Diagnostic Info", samOut)
    cat("\n\n")
}

# Makes Abberations Tab ------------------------------
makeAbTab <- function(sam,philipsFtp="/Volumes/molecular/Molecular/Philips_SFTP") {
    dumpDir <- file.path(getwd(), "zipfiles", sam)
    snvCsv <- file.path(dumpDir, "aberration_snv.csv")
    samCsv <- file.path(dumpDir, "specimen.csv")
    diagCsv <- file.path(dumpDir, "diagnosticorder.csv")
    if (file.exists(snvCsv) & file.exists(samCsv)) {
        printSnvs(snvCsv)
        printSpecInfo(samCsv)
        printDiagInfo(diagCsv)
    } else{
        warning("Data dump directory not found: ",
                file.path(philipsFtp, paste0(sam, ".zip")))
    }
}

makeBlankRow <- function(sam, snvDt) {
    nonMutant <- data.frame(matrix(ncol = ncol(snvDt), nrow = 1))
    colnames(nonMutant) <- colnames(snvDt)
    nonMutant[1, ] <- "None"
    nonMutant$Test_Case <- sam
    makeDT("In-House FrameShifts/INDEL", objDat = nonMutant)
}

LoopSampleTabs <-function(params){
    pactName <- params$pactName
    methData <- gb$GetMethDf(params$pactName)
    qcData <- gb$ReadQcFile(pactName)
    samList <- gb$GetSamList(pactName)
    samples <- gb$GrabSamples(samList)
    snvDt <- read.csv(paste0(pactName, "_desc.csv"))
    outDir <- file.path(params$workDir, paste0(params$pactName,"_consensus"))
    for (sam in samples) {
        gb$makeNewTab(sam, samList, qcData)
        if (sam %in% snvDt$Test_Case) {
            samRows <- snvDt$Test_Case == sam
            snvTab <- snvDt[samRows & snvDt$Variant == "SNV",]
            makeDT("In-House FrameShifts/INDEL", objDat = snvTab)
            cnvTab <- snvDt[samRows & snvDt$Variant == "CNV",]
            cnvTab <- checkDataDump(sam, cnvTab)
            makeDT("CNV", cnvTab, pdfFi = sam, outDir=outDir)
            methCn <- snvDt[samRows & snvDt$Variant == "Methylation",]
            makeMethTab(sam, methCn, methData)
        } else{
            message(sam, " is missing from your Description input File")
            makeBlankRow(sam, snvDt)
        }
        makeAbTab(sam)
    }
}

loadHtmlTag <- function(){
    require("tidyverse")
    library("plotly")
    require("plotly")
    require("ggplot2")
    library("ggplot2")
    library("DT")
    htmltools::tagList(DT::datatable(cars, options=list(pageLength = 100, width="100%")))
    par(ask=F); devAskNewPage(ask = F)
    #doParallel::registerDoParallel(cores=2)
    return(htmltools::tagList(plotly::ggplotly(ggplot2::ggplot())))
}

animation::ani.options(autobrowse = FALSE); options(width = 1600)
knitr::opts_chunk$set(echo = FALSE)
if(!require("xfun")){install.packages("xfun", ask=F)}
