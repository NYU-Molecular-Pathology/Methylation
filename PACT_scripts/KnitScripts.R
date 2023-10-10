gb <- globalenv(); assign("gb", gb)

CheckBrewLatex <- function() {
    # Checks if Brew is installed
    tryCatch(
        system("which brew", intern = T),
        warning = function(ww) {
            system('/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"')
        }
    )
    # Checks if Latex is installed
    tryCatch(
        system("which latex", intern = T),
        warning = function(ww) {
            system("brew install latex")
        }
    )
    
    tryCatch(
        system("which pdflatex", intern = T),
        warning = function(ww) {
            system("brew install pdflatex")
        }
    )
    
    try(system("which pdflatex", intern = T), T)
    try(system("which latex", intern = T), T)
}


# Function to check if a command exists
command_exists <- function(cmd) {
    exit_status <- system(paste("which", cmd, ">/dev/null 2>&1"))
    return (exit_status == 0)
}

# Function to check and install fontconfig via Homebrew
check_brew_fontconfig <- function() {
    brew_list <- system("brew list", intern = TRUE)
    if ("fontconfig" %in% brew_list) {
        message("fontconfig is already installed.")
    } else {
        message("fontconfig is not installed. Installing...")
        system("brew install fontconfig")
    }
}

# Function to check and install fontconfig via Homebrew
check_fontconfig <- function() {
    if (command_exists("brew")) {
        check_brew_fontconfig()
    } else {
        message("Homebrew is not installed on this system. Installing...")
        system(
            '/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"'
        )
        if (command_exists("brew")) {
            check_brew_fontconfig()
        } else {
            message("Failed to install Homebrew. Exiting.")
        }
    }
}

if (Sys.info()[['sysname']]=="Darwin") {
    invisible(check_fontconfig())
    invisible(CheckBrewLatex())
}

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

librarian::shelf(pkgs, ask=F, update_all = F, quiet = F, dependencies=T)
librarian::shelf(pdfToolsPkgs, ask=F, update_all = F, quiet = F, dependencies=T)

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


CheckMethPaths <- function(methData){
    for(i in 1:length(methData$`Report Path`)){
          currPath <- methData$`Report Path`[i]
          currSplit <- stringr::str_split_fixed(currPath, "/",11)[1,]
          if(stringr::str_detect(currSplit[10], "MGDM")==F) {
              next
          }
          runYear <- stringr::str_split_fixed(currSplit[10],"-", 2)[1,1]
          runYear <- paste0("20", runYear)
          currSplit[9] <- runYear
          newPath <- paste(currSplit, collapse = "/")
          methData[i, "Report Path"] <- newPath
      }
    
    checkPaths <- stringr::str_replace_all(methData$`Report Path`,
                                           "smb://shares-cifs.nyumc.org/apps/acc_pathology", "/Volumes")
    checkPaths <- checkPaths[checkPaths != ""]
    anyPathsFalse <- file.exists(checkPaths) == F
    
    if(any(anyPathsFalse)){
        message("Some paths need to be fixed 'Report Path' column of in ", methSheet, ":" )
        cat("Paths don't exist:\n")
        cat(paste(checkPaths[anyPathsFalse], collapse="\n"))
    }
    stopifnot(all(file.exists(checkPaths)))
    return(methData)
}


GetMethDf <- function(pactName) {
    methSheet <- paste0(pactName, "_MethylMatch.xlsx")
    methData <- as.data.frame(readxl::read_excel(methSheet))
    methData <- methData[methData$report_complete == "YES", ]
    if(any(is.na(methData$Test_Number))){
      samCsv <- dir(getwd(), "-SampleSheet.csv", full.names = T)
      wsData <- as.data.frame(read.csv(samCsv, skip = 19))
      toFix <- which(is.na(methData$Test_Number))
      qcData <- gb$ReadQcFile(pactName)
      for(xRow in toFix){
        nSam <- methData$accession_number[xRow]
        ngsRow <- which(stringr::str_detect(wsData$Specimen_ID, pattern = nSam))[[1]]
        methData$Test_Number[xRow] <- wsData$Test_Number[ngsRow]
      }
      
    }
    
    if(nrow(methData) > 0){
        methData <- CheckMethPaths(methData)
    }
    
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
        message(crayon::bgRed("The following cases have no NGS number and are being excluded:"), "\n")
        message(paste(utils::capture.output(as.data.frame(samList[!toKeep,])), collapse = "\n"))
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


makePdfTab <- function(pdfFi, cnvTab, outDir) {
    pdfFi <- file.path(outDir, "cnvpng", paste0(pdfFi, ".png"))
    cat(
        "<span style='color: red;'>",
        "**Note**: CNV Amplifications are listed first,",
        "followed other Philips CNV abberations in table below this image\n\n"
    )
    cat(paste0("![In-house Facets Plot (Above)](", pdfFi, "){width=100%}"))
    cat("\n\n")
    if (nrow(cnvTab) > 0) {
        newTa <- knitr::kable(cnvTab, row.names = F, "html", align="c")
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
    cnvFi <- file.path(".","methCNV", paste0(rdNumb, "_cnv.png"))
    altTxt <- paste("![Methylation",rdNumb,sam,"CNV Plot](")
    cat(paste0(altTxt,cnvFi,"){height=400px}"))
    cat("\n\n")
}

makeSpecimenTab <- function(objDat) {
    newTa <- knitr::kable(objDat, row.names = F, "html")
    newTa <- kableExtra::kable_styling(
        newTa, bootstrap_options = c("condensed"), full_width = T, position = "left")
    if(ncol(objDat) > 3){
    newTa <- kableExtra::column_spec(newTa, 1:5, width = "2cm")    
    }else{
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
sorted_objDat <- objDat[order(-rows_with_yes),]
dtOpts <- list(
        scrollX = T, scrollY=T, info = F, autoWidth = F, pageLength = 25,
        lengthChange = T, searchable = T
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
for (col in c("Strelka", "LoFreqSomatic", "Mutect")) {
    datTab <- datTab %>%
      formatStyle(col, backgroundColor = styleEqual("YES", "#64a463"), target = 'cell')
}

dtTab <- htmltools::tagList(datTab)

print(dtTab)
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
            "INDEL/Frameshift calls have the following filters:",
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

MakeRegularTab <- function(tabNam, objDat) {
    if (!is.null(objDat) & stringr::str_detect(tabNam, "Methylation") == F) {
        if (tabNam == "QC") {
            makeQCTab(objDat)
        } else{
            CheckTabName(tabNam)
            makeDefaultDt(objDat)
        }
    } else{
        cat("\n\nNo additonal results for this case yet\n\n")
        cat("\n\n")
    }
}


# Prints out the table tab headings and tab information based on object parameters provided
makeDT <- function(tabNam, objDat, pdfFi = NULL, rdNumb = NULL, sam = NULL, outDir=NULL){
    MakeTabColor(tabNam)
    if (!is.null(pdfFi)) {
        return(makePdfTab(pdfFi, cnvTab=objDat, outDir))
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


makeDT <- function(tabNam, objDat, pdfFi = NULL, rdNumb = NULL, sam = NULL, outDir=NULL){
    MakeTabColor(tabNam)
    if (!is.null(pdfFi)) {
        return(makePdfTab(pdfFi, cnvTab=objDat, outDir))
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


MakeHStab <- function(sam, hsDat, samList, outDir){
    ngsRows <- samList$Test_Number == sam
    ts_number <- samList$Specimen_ID[samList$Tumor_Content!=0 & ngsRows]
    hotspot_tsv <- file.path(outDir, 'hotspots', paste0(ts_number, "_Hotspots.tsv"))
    hsMain <- hsDat[hsDat$Test_Number == sam,]
    makeDT(tabNam = "Hotspots", objDat = hsMain, NULL, NULL, sam)
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
    if(!file.exists(cnvFi)){
      warning("File does not exist:\n", cnvFi, "\nSkipping Philips...")
      cnvTab$Gene <- "No Philips Data Dump CSV file availible, check ISPM for CNV"
      return(cnvTab)
    }
    cnvInfo <- as.data.frame(readr::read_csv(
        cnvFi, col_types = readr::cols(.default = readr::col_character())
    ))
    if (nrow(cnvInfo) > 0) {
        cnvTab <- parsePhilipsCn(cnvInfo, sam, cnvTab)
    }
    return(cnvTab)
}


getDumpFiles <- function(
    outPath, sam, cnvTab, philipsFtp = "/Volumes/molecular/Molecular/Philips_SFTP")
{
    samZip <- paste0(sam, ".zip")
    dumpDir <- file.path(philipsFtp, samZip)
    destDir <- file.path(outPath, paste0(basename(samZip)))
    if (file.exists(dumpDir) & !file.exists(destDir)) {
        try(fs::file_copy(path = dumpDir, new_path = destDir), silent = T)
        unzip(zipfile = destDir, exdir = outPath)
    }
    cnvTab <- parseCNV(outPath, sam, cnvTab)
    return(cnvTab)
}

# Appends any additional cnv abberations from data dump
checkDataDump <- function(sam, cnvTab) {
    if(nrow(cnvTab) == 0){
        cnvTab <- as.data.frame(matrix(nrow=1, ncol=length(cnvTab), dimnames=list(1, names(cnvTab))))
        cnvTab$Test_Case <- sam
    }
    outPath <- file.path(getwd(), "zipfiles")
    if (!dir.exists(outPath)) {dir.create(outPath)}
    cnvTab <- getDumpFiles(outPath, sam , cnvTab)
    drpCol <- which(stringr::str_detect(colnames(cnvTab), "In.|IGV"))
    cnvTab <- cnvTab[,-drpCol]
    names(cnvTab)[names(cnvTab) == 'Comments'] <- 'Copy Number'
    names(cnvTab)[names(cnvTab) == 'Other'] <- 'Chromosome'
    cnvTab$`Copy Number` <- as.numeric(cnvTab$`Copy Number`)
    cnvTab <- cnvTab[!is.na(cnvTab$Test_Case),]
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
        outputFile <- file.path("cnvpng", pngName)
        if(!file.exists(outputFile)){
            tryCatch(
                expr = {suppressWarnings(pdftools::pdf_convert(pdfFile, filenames = outputFile, dpi = 300))},
                error = function(e) {
                    message(e,"\nTry running:\nbrew install fontconfig --universal")
                    system("brew update --auto-update")
                })
            }}
        ))
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

# Checks if the facets pdfs have been converted to png ------------------------------------
CopyPdfsPngs <- function(params) {
    outDir <- file.path(params$workDir, paste0(params$pactName,"_consensus"))
    samList <- gb$GetSamList(params$pactName)
    pdfDir <- file.path(outDir,"FACETpdfs") # input facet pdf directory
    tumors <- checkTumorPdf(samList, outDir)
    pdfList <- list.files(pattern = "*.pdf", recursive = T)
    if (length(pdfList)>0) {
        convert.plots(tumors, pdfList)
        if(!dir.exists(pdfDir)){dir.create(pdfDir)}
        try(fs::file_move(pdfList, pdfDir), silent = T)
    }
}


CheckHotspots <- function(params){
    outDir <- file.path(params$workDir, paste0(params$pactName,"_consensus"))
    hotspot_tsv <- file.path(outDir, paste0(params$pactName,"_Hotspots.tsv"))
    if(file.exists(hotspot_tsv)){
        message("Hotspots exist")
    }else{
        message("Hotspots Missing:\n", hotspot_tsv)
    }
}


# Generates the link to the BAM file based on PACT run ------------------------------------
makeBamLink <- function(sam, pactID){
    samsheet <- list.files(".", "-SampleSheet.csv", T)[1]
    samCSVfi <- as.data.frame(read.csv(samsheet, skip = 19))
    allNames <- samCSVfi[,1]
    rowNGS <- which(samCSVfi$Test_Number == sam)
    bamFiNam <- paste0(allNames[rowNGS], ".dd.ra.rc.bam")
    bamFiUrl <- file.path("https://genome.med.nyu.edu/external/clinpathlab", pactID, bamFiNam)
    cat(paste0("```\nIGV Load from URL: ", bamFiUrl, "\n```\n\n"))
}


# Generates a new Sample Tabbed row in html ------------------------------------------------
makeNewTab <- function(sam, samList, qcData, pactID){
    currSam <- samList[samList$Test_Number==sam,"Specimen_ID"]
    cat(' <div class="boxed"> ')
    cat(paste0("\n\n# **", sam, "** {.tabset}","\n\n"))
    cat(paste0("(", currSam[1], ")\n\n"))
    makeBamLink(sam, pactID)
    currQC <- grepl(paste(currSam, collapse="|"), qcData$Sample)
    qcTab <- qcData[currQC,2:ncol(qcData)]
    if (nrow(qcTab) == 0) {
        warning(sam, " is missing from qcTsv file")
    } else{
        makeDT("QC", qcTab)
    }
}


GetMethMatchInfo <- function(methData, theRd){
    rdDtRows <- methData[methData$record_id == theRd, ]
    cat("\n\n")
    newTa <- knitr::kable(rdDtRows, row.names = F, "html")
    newTa <- kableExtra::kable_styling(newTa, bootstrap_options = c("bordered","condensed"), full_width = F, position = "left")
    newTa <- kableExtra::column_spec(newTa, 1:5, width = "5cm")
    cat("#### Methylation Full Report Link:")
    cat("\n\n")
    txtLink <- paste0("[", theRd, "](", rdDtRows$`Report Path`, ")")
    cat(txtLink)
    cat("\n\n")
    pathLink <- file.path(stringr::str_split_fixed(dirname(rdDtRows$`Report Path`), "acc_pathology", 2)[1, 2], 
                          basename(rdDtRows$`Report Path`))
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
        if (length(currRD) > 1) {
            message(sam, " has more than one RD-number:\n", paste(currRD, sep = " ", collapse=" "))
            for (rd in 1:length(currRD)) {
                rdTab <- paste0("Methylation_", rd)
                makeDT(rdTab, methCn, pdfFi = NULL, currRD[rd], sam)
                GetMethMatchInfo(methData, currRD[rd])
            }
        } else{
            makeDT("Methylation", methCn, NULL, currRD, sam)
            GetMethMatchInfo(methData, currRD)
        }
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

# Renders error when csv file fails to be parsed ------------------------------------------------------------
PrintParseErr <- function(csvPath, tabTxt){
    eTxt = "\n\n## **Failed to read Philips CSV File"
    snvErr = paste(eTxt, tabTxt, "**\n`", csvPath,"`\n\n")
    return(cat(snvErr))
}

# Makes Abberations Tab ------------------------------------------------------------
makeAbTab <-  function(sam, philipsFtp="/Volumes/molecular/Molecular/Philips_SFTP") {
    dumpDir <- file.path(getwd(), "zipfiles", sam)
    snvCsv <- file.path(dumpDir, "aberration_snv.csv")
    samCsv <- file.path(dumpDir, "specimen.csv")
    diagCsv <- file.path(dumpDir, "diagnosticorder.csv")
    if (file.exists(snvCsv) | file.exists(samCsv) | file.exists(diagCsv)) {
        tryCatch(printSnvs(snvCsv), error=function(e){PrintParseErr(snvCsv, "SNV")})
        tryCatch(printSpecInfo(samCsv), error=function(e){PrintParseErr(samCsv, "Sample")})
        tryCatch(printDiagInfo(diagCsv), error=function(e){PrintParseErr(diagCsv, "Diagnostics")})
    } else{
        zipFiN <- file.path(philipsFtp, paste0(sam, ".zip"))
        cat("\n\n## **No Philips Data Dump**\n\n")
        cat("Data dump not found:\n")
        cat(paste0(zipFiN, "\n\n"))
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


GrabHotspots <- function(params){
    outDir <- file.path(params$workDir, paste0(params$pactName,"_consensus"))
    hsTsv <- file.path(outDir, paste0(params$pactName, "_Hotspots.tsv"))
    if(file.exists(hsTsv)){
        hsDat <- as.data.frame(read.csv(hsTsv, sep='\t'))
    }else{
        message("No file found here:", hsTsv)
        hsDat <- NULL
    }
    return(hsDat)
}


LoopSampleTabs <- function(params){
    pactName <- params$pactName
    methData <- gb$GetMethDf(params$pactName)
    qcData <- gb$ReadQcFile(pactName)
    samList <- gb$GetSamList(pactName)
    samples <- gb$GrabSamples(samList)
    hsDat <- gb$GrabHotspots(params)
    snvDt <- read.csv(paste0(pactName, "_desc.csv"))
    outDir <- file.path(params$workDir, paste0(params$pactName,"_consensus"))
    for (sam in samples) {
      makeNewTab(sam, samList, qcData, params$pactName)
      if (!(sam %in% snvDt$Test_Case)) {
        message(sam, " is missing from your Description input File:\n",
                paste0(pactName, "_desc.csv"))
        message("Adding Sample as a blank row")
        snvDt <- makeBlankRow(sam, snvDt) 
      }
      samRows <- snvDt$Test_Case == sam
      snvTab <- snvDt[samRows & snvDt$Variant == "SNV",]
      makeDT("In-House FrameShifts/INDEL", objDat = snvTab)
      cnvSam <- snvDt[samRows & snvDt$Variant == "CNV",]
      cnvTab <- checkDataDump(sam, cnvSam)
      makeDT("CNV", cnvTab, pdfFi = sam, outDir=outDir)
      methCn <- snvDt[samRows & snvDt$Variant == "Methylation",]
      makeMethTab(sam, methCn, methData)
      makeAbTab(sam)
      if(!is.null(hsDat)){
          MakeHStab(sam, hsDat, samList, outDir)
      }
      cat("\n\n")
      cat(' </div> ')
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
