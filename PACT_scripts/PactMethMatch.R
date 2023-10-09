#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library("base"); gb <- globalenv(); assign("gb", gb)
dsh="\n================"
dsh2="================\n"

# Main arguments input in comandline
token <- args[1]
inputSheet <- args[2]

# Displays the Input args -----
message(dsh,"Parameters input",dsh2)
message("token: ",token)
message("inputSheet: ", inputSheet,"\n")

stopifnot(!is.na(token))
stopifnot(!is.na(inputSheet))

readFlag <- endsWith(inputSheet,".csv")==T
stopifnot(rlang::is_bool(readFlag))

# REDcap Heading Fields -----
flds = c("record_id","b_number","tm_number","accession_number","block","diagnosis",
         "organ","tissue_comments","run_number", "nyu_mrn")

# Load redcapAPI Package -----
if(!requireNamespace("redcapAPI")){
    params = list('nutterb/redcapAPI', dependencies = T, upgrade = "always", type = "source")
    do.call(devtools::install_github, c(params))
}

current_version <- utils::packageVersion("redcapAPI")
required_version <- "2.7.4"

if (is.na(current_version) || utils::compareVersion(as.character(current_version), required_version) < 0) {
    install.packages("redcapAPI", ask=FALSE, update=TRUE, dependencies=TRUE)
}

supM <- function(sobj){return(suppressMessages(suppressWarnings(sobj)))}

# FUN: Checks if z-drive is accessible to the Rscript
checkMounts <- function(){
    molecDrive = "/Volumes/molecular/MOLECULAR LAB ONLY"
    zDrive = "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular"
    failMount <- ifelse(dir.exists(molecDrive),T,F)
    if(failMount!=T){
        cat("PATH does not exist, ensure path is mounted:")
        cat(crayon::white$bgRed$bold(molecDrive))
        cat("You must mount the network Z-drive path:\n")
        cat(crayon::white$bgRed$bold(zDrive),"\n")
        stopifnot(!any(failMount==T))
    } else {message("\n",crayon::bgGreen("Z-drive path is accessible"),"\n")}
}


# Functions to load packages and get redcap info -----
loadPacks <- function(){
    pkgs = c("data.table", "foreach", "openxlsx", "jsonlite", "RCurl", "readxl", "stringr", "webshot",
             "tidyverse", "crayon", "plotly", "htmlwidgets", "tinytex")
    rlis = getOption("repos")
    rlis["CRAN"] = "http://cran.us.r-project.org"
    options(repos = rlis)
    invisible(lapply(pkgs, function(pk){
        if(supM(!requireNamespace(pk))){
            install.packages(pk, dependencies = T, verbose = T, repos = "http://cran.us.r-project.org", type = "both")
        }}))
    if(!requireNamespace("systemfonts")){devtools::install_github('r-lib/systemfonts')}
    if(!requireNamespace("redcapAPI")){install.packages("redcapAPI", dependencies = T, type="both",ask=F)}
    if(!requireNamespace("remotes")){install.packages("remotes", dependencies=T)}
    #if(!requireNamespace("chromote")){remotes::install_github("rstudio/chromote", upgrade ="never")}
    #if(!requireNamespace("webshot2")){remotes::install_github("rstudio/webshot2")}
    supM(library("redcapAPI"))
    #supM(library("chromote"))
    #try(tinytex::tlmgr_install(),silent=T)
    #try(tinytex:::install_prebuilt(),silent=T)
}

# API Call functions -----
grabAllRecords <- function(flds, rcon){
    message("Pulling REDCap data...")
    params = list(rcon, fields=flds, labels=F, dates = F, survey = F, dag = F, factors=F, form_complete_auto=F)
    dbCols <- do.call(redcapAPI::exportRecordsTyped, c(params))
    return(as.data.frame(dbCols))
}

# Database search function -----
searchDb <- function(queryList, db){
    v2f <- paste(queryList, collapse="|"); i=NULL
    res <- NULL
    for (idx in 1:length(queryList)) {
        item <- queryList[idx]
        ngsNumber = paste(names(item)[1])
        for (i in colnames(db)) {
            ngsMatch <- which(grepl(item, db[,i]))
            if(length(ngsMatch) > 0) {
                message("Match found for ",item, " (", db[ngsMatch,i], ") for ",
                        ngsNumber," in:\n", '"', i,'"', " column")
                dbMatch <- db[ngsMatch,]
                dbMatch$Test_Number <- names(item)
                if(is.null(res)){
                    res <- dbMatch
                }else{
                    res <- rbind(res, dbMatch)
                }
            }
        }
    }
    return(res)
}

# Filters list of possible files in the directory for worksheet
filterFiles <- function(potentialFi) {
    wbFiles <- stringr::str_which(basename(potentialFi), pattern = "xlsm")
    if (length(wbFiles) == 0) {
        message("\nNo .xlsm worksheet found. Checking .xlsx files...\n")
        wbFiles <- stringr::str_which(basename(potentialFi), pattern = "book")
    }
    potentialFi <- potentialFi[wbFiles]
    potentialFi <- potentialFi[!stringr::str_detect(potentialFi, "\\$")]
    return(potentialFi)
}

getAltPath <- function(inputFi){
    message("\n", crayon::bgRed("This PACT run worksheet was not found:"), "\n", inputFi)
    message(dsh, "Checking other files in PACT folder: ", basename(dirname(inputFi)), "================")
    potentialFi <- list.files(dirname(inputFi), full.names=T)
    print(basename(potentialFi))
    altFi <- filterFiles(potentialFi)
    if(length(altFi)>1){altFi <- altFi[1]}
    message(crayon::bgGreen("Using this workbook instead:"), basename(altFi), "\n")
    stopifnot(file.exists(altFi)==T)
    return(altFi)
}

# Returns the path to /Volumes/molecular/MOLECULAR LAB ONLY/NYU PACT Patient Data/Workbook/YYYY/PACT-YY-##"
getPactFolder <- function(inputSheet){
    drive <- file.path("", "Volumes", "molecular", "MOLECULAR LAB ONLY")
    folder <- file.path("NYU PACT Patient Data", "Workbook")
    runyr <- stringr::str_split_fixed(inputSheet,"-",3)[,2]
    runFolder <- file.path(drive, folder, paste0("20", runyr), inputSheet)
    print(runFolder); stopifnot(dir.exists(runFolder)==T); return(runFolder)
}

# Returns path to run worksheet if path is valid
getFilePath <- function(inputSheet){
    runFolder <- getPactFolder(inputSheet)
    inputFi <- file.path(runFolder, paste0(inputSheet, ".xlsm"))
    if(file.exists(inputFi)) {return(inputFi)} else{return(getAltPath(inputFi))}
}

# Parses the input file for the "PhilipsExport" tab
parseWorksheet <- function(inputFi){
    sheet2Read <- "PhilipsExport"
    message("Reading the file:\n", inputFi)
    shNames <- readxl::excel_sheets(inputFi)
    message("Excel sheet names:")
    message(paste(shNames, collapse="\n"))
    stopifnot(!is.null(shNames) & length(shNames) > 2)
    sh <- which(grepl(sheet2Read, shNames, ignore.case = T))[1]
    pactShCol <- c("Tumor Specimen ID", "Normal Specimen ID", "Tumor DNA/RNA Number", "MRN", "Test Number")
    vals2find <-  suppressMessages(as.data.frame(readxl::read_excel(
        inputFi, sheet=shNames[sh], skip=3, col_types ="text")[,pactShCol]))
    vals2find <- vals2find[!is.na(vals2find[,1]),]
    return(vals2find)
}

# Import csv Worksheet -----
getCaseValues <- function(inputSheet, readFlag) {
    isSamSheet <- stringr::str_detect(inputSheet, "-SampleSheet") == T
    if (readFlag) {
        if (isSamSheet == T) {
            vals2find <- utils::read.csv(inputSheet, skip = 19)[, c(6, 7, 9)]
            vals2find <- vals2find[!grepl("H20|SERACARE|HAPMAP", vals2find[, 2]),]
        } else{
            vals2find <- as.data.frame(read.csv(inputSheet))
            if (ncol(vals2find) > 1) {
                vals2find <- unlist(lapply(X = 1:ncol(vals2find), FUN = function(X) {vals2find[, X]}))
            }
        }
        vals2find <- as.data.frame(vals2find)
        return(vals2find)
    }
    if (stringr::str_detect(inputSheet, .Platform$file.sep) == T) {
        vals2find <- parseWorksheet(inputSheet)
    } else{
        inputFi <- getFilePath(inputSheet)
        vals2find <- parseWorksheet(inputFi)
    }
    return(vals2find)
}

genQuery <- function(dbCol,vals2find){
    currCol <- vals2find[, dbCol]
    toKeep <- which(currCol != 0 & !is.na(currCol) & currCol != "")
    q1 <- currCol[toKeep]
    names(q1) <- vals2find$`Test Number`[toKeep]
    return(q1)
}

# Get Methylation and Molecular Samples list ----
queryCases <- function(vals2find, db) {
    library("foreach")
    i=NULL
    queryList <- foreach::foreach(i=1:ncol(vals2find), .combine="c", .packages="foreach") foreach:::`%do%` {genQuery(i,vals2find)}
    tsTb <- stringr::str_detect(queryList, "TS|TB|TC")
    theTScases <- queryList[tsTb]
    for (x in 1:length(theTScases)) {
        y <- stringr::str_split_fixed(theTScases[x], "-", 3)[1:2]
        z <- paste(y,  collapse="-")
        theTScases[x] <- z
    }
    queryList <- c(queryList, theTScases)
    queryList <- queryList[!duplicated(queryList)]
    methQuery <- searchDb(queryList, db)
    return(unique(methQuery))
}

addOutputLinks <- function(output){
    winpath = "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular/Molecular/MethylationClassifier/"
    yearSplit <- stringr::str_split_fixed(output$run_number,"-",2)[,1]
    yearSplit <- stringr::str_remove_all(yearSplit, "MC")
    yearSplit <- yearSplit[!is.na(yearSplit)]
    yearPath <- lapply(yearSplit, function(yr) {
        rnum <- NULL
        if(nchar(yr)>2){rnum <- substring(yr, 3)}else{rnum <- yr}
        if(nchar(yr)>0){paste0("20",rnum)}else{rnum}
    })
    output$report_complete <- ifelse(!is.na(output$run_number), "YES", "NOT_YET_RUN")
    output$'Report Link' <- paste0(winpath, yearPath,"/",output$run_number,"/",output$record_id,".html")
    output$'Report Link'[is.na(output$run_number)] <- ""
    output$'Report Path'<-output$'Report Link'
    return(output)
}


FillMissingNGS <- function(output, vals2find){
    skipCols = F
    rows2fill <- which(is.na(output$Test_Number))
    outputRows <- output$accession_number[rows2fill]
    matchedTs <- which(grepl(paste(outputRows, collapse = "|"),
                             vals2find$`Tumor Specimen ID`, ignore.case = T))
    if(length(matchedTs)>0){
        foundTs <- vals2find$`Tumor Specimen ID`[matchedTs]
        ngsFound <- vals2find$`Test Number`[matchedTs]
        if(length(foundTs)==1){
            skipCols = T
        }
        foundTs <- stringr::str_split_fixed(foundTs, "-", 3)[,1:2]
        if(length(foundTs)>0){
            if(skipCols==F){
                foundTs <- paste(foundTs[,1], foundTs[,2], sep = "-")
            }else{
                foundTs <- paste(foundTs[1], foundTs[2], sep = "-")
            }
        }
        matchFound <- which(grepl(paste(foundTs, collapse = "|"),
                                  output$accession_number, ignore.case = T))
        output[matchFound,"Test_Number"] <- ngsFound
        rows2fill <- which(is.na(output$Test_Number))
    }
    if(length(rows2fill)>0){
        warning("Some samples still missing NGS Numbers:\n",
                paste(output$record_id[rows2fill], collapse = "\n"))
    }
    return(output)
}


modifyOutput <- function(output, vals2find) {
    #output$Test_Number <- NA
    if (length(vals2find$`Test Number`) == 0) {
        vals2find$`Test Number` <- ""
    }
    NGSmissing <- F
    if(all(output$Test_Number %in% vals2find$`Test Number`)){
        message("All NGS Found")
    }else{
        message("Not all NGS do not have methylation")
    }
    # for (i in 1:nrow(output)) {
    #   theVal = NA
    #   for (var in 1:ncol(vals2find)) {
    #     pat <- vals2find[, var]
    #     pat <- unique(pat)
    #     currRow <- paste(output[i,])
    #     theMatch <- which(pat %in% currRow[currRow != "0" & currRow != "NA"])
    #     if (length(theMatch) > 1) {
    #       warning(paste("Multiple records found matching pattern:",
    #                     paste(theMatch, collapse = " "),"\nPattern:",
    #                     paste(pat, collapse = " ")))
    #     }
    #     if (length(theMatch) > 0) {
    #       theVal <- vals2find$`Test Number`[theMatch]
    #     }
    #     output$Test_Number[i] <- theVal
    #   }
    #   if (is.na(theVal)) {
    #     warning(paste(output$record_id[i], "is missing NGS Number!"))
    #     NGSmissing <- T
    #   }
    # }
    if(NGSmissing==T){
        output <- FillMissingNGS(output, vals2find)
    }
    output <- addOutputLinks(output)
    return(output)
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
        message("Some paths need to be fixed 'Report Path' column of in the methylation sheet!" )
        cat("Paths don't exist:\n")
        cat(paste(checkPaths[anyPathsFalse], collapse="\n"))
    }

    return(methData)
}


addExcelLink <- function(output, fiLn, wb, runId){
    x <- c(output$'Report Link'[fiLn])
    names(x) <- paste0(output$record_id[fiLn],".html")
    class(x) <- "hyperlink"
    colNum <- which(colnames(output)=="Report Link")

    openxlsx::writeData(wb, sheet = runId, x = x, startCol = colNum, startRow = fiLn + 1)
}

createXlFile <- function(runId, output){
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, runId)
    openxlsx::writeData(wb, sheet = runId, x = output)
    output <- CheckMethPaths(methData = output)
    for (fiLn in 1:length(output$'Report Link')) {
        if(output$'Report Link'[fiLn]!=''){
            addExcelLink(output, fiLn, wb, runId)
        }
    }
    outFi <-file.path(fs::path_home(),"Desktop", paste0(runId,"_MethylMatch.xlsx"))
    openxlsx::saveWorkbook(wb, outFi, overwrite = T)
    return(outFi)
}

postData <- function(rcon, record){
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
    res <- RCurl::postForm(rcon$url,token=rcon$token,content='record',format='json',type='flat', data=datarecord, returnContent = 'nothing', returnFormat = 'csv')
    cat(res)
}

emailFile <- function(runId, outFi, rcon){
    record = data.frame(record_id = runId, run_number = runId)
    postData(rcon, record)
    isDone <- redcapAPI::exportRecordsTyped(rcon, factors=F, records=record$record_id, fields=c("record_id", "other_file"))
    if (length(isDone$other_file) == 0) {
        #rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", gb$token)
        body <- list(
            token = rcon$token,
            content = 'file',
            action = 'import',
            record = runId,
            field = "other_file",
            file = httr::upload_file(outFi),
            returnFormat = 'csv'
        )
        res <-
            tryCatch(
                httr::POST(url = rcon$url, body = body, config = rcon$config),
                error = function(cond){
                    list(status_code = "200")
                }
            )
        if(res$status_code=="200"){
            message("REDCap file upload successful: ", outFi)
        }else{
            message("REDCap file upload failed: ", outFi)
        }
    }
    record$comments <- "pact_sample_list_email"
    postData(rcon, record)
    message(dsh,"\nEmail Notification Created",dsh,"\n")
}

# Main workflow to get dataframe
getOuputData <- function(token, flds, inputSheet, readFlag){
    apiUrl = "https://redcap.nyumc.org/apps/redcap/api/"
    rcon <- redcapAPI::redcapConnection(apiUrl, token)
    vals2find <- getCaseValues(inputSheet, readFlag)
    db <- grabAllRecords(flds, rcon)
    if(nrow(db)==0){
        message("REDCap API connection failed!\n",
                "Check the Database for non-ASCII characters and verify the API Token: ", token)
        stopifnot(nrow(db)>0)
    }
    if(class(vals2find)!="data.frame"){
        vals2find <- as.data.frame(vals2find)
    }
    if(nrow(vals2find)!=0){
        output <- queryCases(vals2find, db)
    }
    if(nrow(output)==0){
        warning("No Methylation Cases on this Pact run, generating blank file")
        output[1,] <- "NONE"
    }else{
        rownames(output) <- 1:nrow(output)
    }
    if(readFlag==T){
        output <- modifyOutput(output, vals2find)
        write.csv(output, file="meth_sample_data.csv", quote=F, row.names=F)
        return(output)
    }
    runId <- ifelse(readFlag == T, substr(inputSheet, 1, nchar(inputSheet) - 4), inputSheet)
    output <- modifyOutput(output, vals2find)
    toDrop <- is.na(output$Test_Number)
    if(any(toDrop)){
        hasNGS <- which(grepl("NGS", output$tm_number))
        message("NGS Number found in 'tm_number' field of REDCap!")
        if(length(hasNGS) > 0){
            output$Test_Number[hasNGS] <- output$tm_number[hasNGS]
        }
    }
    outFi <- createXlFile(runId, output)
    emailFile(runId, outFi, rcon)
    return(output)
}

# FUN: Sets your directory and sources the helper functions
sourceFuns2 <- function(workingPath = NULL) {
    mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/"
    script.list <- c("R/SetRunParams.R","R/CopyInputs.R","PACT_scripts/generateCNV.R", "Research/cnvFunctions.R")
    if (is.null(workingPath)) {workingPath = getwd()}
    scripts <- paste0(mainHub, script.list)
    invisible(lapply(scripts, function(i){
        suppressPackageStartupMessages(devtools::source_url(i))}))
    supM(library("sest"))
    supM(library("mnp.v11b6"))
    #supM(require("plotly"))
    #require("htmlwidgets")
    gb$setDirectory(workingPath)
    return(gb$defineParams())
}

msgRDs <- function(rds, token){
    message("\nRD-numbers with idats:\n")
    message(paste(rds, collapse="\n"))
    assign("rds", rds)
    message(dsh, crayon::bgMagenta("Starting CNV PNG Creation"), dsh2)
    ApiToken <- token
    assign("ApiToken", ApiToken)
}


grabRDs1 <- function(rd_numbers, token){
    result <- gb$search.redcap(rd_numbers, token)
    result <- result[!is.na(result$barcode_and_row_column),]
    samplesheet_ID = as.data.frame(stringr::str_split_fixed(result[,"barcode_and_row_column"],"_",2))
    gb$writeFromRedcap(result, samplesheet_ID) # writes API export as minfi dataframe sheet
    gb$get.idats()  # copies idat files from return to current directory
}

msgCreated <- function(mySentrix){
    pngFiles <- paste0(file.path(fs::path_home(),"Desktop",mySentrix[, 1]),"_cnv.png")
    cnvMade <- file.exists(pngFiles)
    if(any(cnvMade==F)){
        message("The following failed to be created:")
        print(pngFiles[!cnvMade])
        message("Try running again or check GitHub troubleshooting")
    }
    if(any(cnvMade==T)){
        message("The following were created successfully:")
        print(pngFiles[cnvMade])
    }
}


loopCNV <- function(mySentrix, asPNG){
    for (sam in rownames(mySentrix)) {
        sampleName <- mySentrix[sam, 1]
        fn = file.path(fs::path_home(), "Desktop", paste0(sampleName, "_cnv.png"))
        if (file.exists(fn)) {
            message("\nFile already exists, skipping:", fn, "\n")
        } else{
            sentrix <- mySentrix[sam, "SentrixID_Pos"]
            message("\nGetting RGset for ", sentrix, "\n")
            RGsetEpic <- gb$grabRGset(getwd(), sentrix)
            tryCatch(
                expr = {gb$gen.cnv.png2(RGsetEpic, sampleName, asPNG)},
                error = function(e) {
                    erTxt <-paste0("An error occured with ", mySentrix[sam, 1]," png creation:")
                    message(crayon::bgRed(erTxt), "\n", e)
                    message(crayon::bgGreen("Trying next sample"))}
            )
        }
    }
}


CheckIfPngExists <- function(
        rds,
        outFolder = "/Volumes/molecular/Molecular/MethylationClassifier/CNV_PNG") {
    outpng <- paste0(rds, "_cnv.png")
    outFiles <- file.path(outFolder, outpng)
    finished <- file.exists(outFiles)
    if (any(finished)) {
        message(crayon::bgGreen(
            "The following samples are completed and will be skipped:"), "\n",
            paste(capture.output(outFiles[finished]), collapse = '\n')
        )
        rds <- rds[!finished]
    }
    return(rds)
}


gb$SaveConumeePACT <-  function(x, sampleImg, doXY=F){
    chrAll <- paste0("chr", 1:22)
    if(doXY==T){
        chrAll <- "all"
    }
    message("Saving file to:\n", sampleImg)
    png(filename = sampleImg, width = 1820, height = 1040, res=150)
    conumee::CNV.genomeplot(x, chr = chrAll)
    invisible(dev.off())
}


gb$SaveCNVplotsPACT <- function(samplename_data, sentrix.ids, i, idatPath = NULL, chrNum = NULL, doXY = F) {
    if (is.null(idatPath)) {
        idatPath <- getwd()
    }
    samName <- samplename_data[i]
    sampleEpic <- sentrix.ids[i]
    sampleImg <- file.path(fs::path_home(), "Desktop", paste0(samName, "_cnv.png"))
    pathEpic <- file.path(idatPath, sampleEpic)
    RGsetEpic <- minfi::read.metharray(pathEpic, verbose = T, force = T)
    if(RGsetEpic@annotation[['array']] == "IlluminaHumanMethylation450k"){
        library("mnp.v11b4")
        MsetEpic <- mnp.v11b4::MNPpreprocessIllumina(RGsetEpic, bg.correct = T, normalize = "controls")
    }else{
        library("mnp.v11b6")
        MsetEpic <- mnp.v11b6::MNPpreprocessIllumina(RGsetEpic, bg.correct = T, normalize = "controls")
    }
    x <- gb$customCNV(MsetEpic, samName, NULL)
    slot(x, 'detail', check = FALSE) <- NULL
    invisible(format(object.size(x), units = 'auto'))
    gb$SaveConumeePACT(x, sampleImg, F)
}


LoopSavePlainCNV3 <- function(targets) {
    samplename_data <- as.character(targets[,1])
    sentrix.ids <- as.character(targets$SentrixID_Pos)
    for (i in 1:length(sentrix.ids)) {
        gb$SaveCNVplotsPACT(samplename_data, sentrix.ids, i)
    }
}


TryCnvMaker <- function(myDt) {
    tryCatch(
        expr = {
            gb$makeCNV(myDt)
        },
        error = function(e) {
            message("The following error occured:\n", e)
            message("\n\nTry checking the troubleshooting section on GitHub:\n")
            message("https://github.com/NYU-Molecular-Pathology/Methylation/PACT_scripts/README.md\n"
            )
        },
        finally = {
            gb$copyOutputPng()
        }
    )
}


makeCNV <- function(myDt, asPNG = T) {
    mySentrix <- myDt[myDt[, "SentrixID_Pos"] %like% "_R0", ]
    if (nrow(mySentrix) > 0) {
        #loopCNV(mySentrix, asPNG)
        gb$LoopSavePlainCNV3(myDt)
    } else{
        message("The RD-number(s) do not have idat files in REDCap:/n")
        print(myDt)
    }
    msgCreated(mySentrix)
    while (!is.null(dev.list())) {dev.off()}
}


CheckIfPngExists <- function(rds, outFolder = "/Volumes/molecular/Molecular/MethylationClassifier/CNV_PNG") {
    outpng <- paste0(rds, "_cnv.png")
    outFiles <- file.path(outFolder, outpng)
    finished <- file.exists(outFiles)
    if (any(finished)) {
        message(
            crayon::bgGreen("The following samples are completed and will be skipped:"),
            "\n",
            paste(capture.output(outFiles[finished]), collapse = '\n')
        )
        rds <- rds[!finished]
    }
    return(rds)
}


TryCnvMaker <- function(myDt) {
    tryCatch(
        expr = {
            gb$makeCNV(myDt)
        },
        error = function(e) {
            message("The following error occured:\n", e)
            message("\n\nTry checking the troubleshooting section on GitHub:\n")
            message(
                "https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/PACT_scripts/README.md\n"
            )
        },
        finally = {
            gb$copyOutputPng()
        }
    )
}

GetSampleList <- function(rds, sampleSheet="samplesheet.csv"){
    myDt <- read.csv(sampleSheet)
    toDrop <- myDt$Sample_Name %in% rds
    myDt <- myDt[toDrop,]
    rownames(myDt) <- 1:length(rds)
    return(myDt)
}

GreenMsg <- function(strMsg){
    message(crayon::bgGreen(strMsg))
}

QueCnvMaker <- function(output, token) {
    rds <- output$record_id[output$report_complete == "YES"]
    if (all(!is.null(rds)) == T & all(!is.na(rds)) == T & length(rds) > 0) {
        msgRDs(rds, token)
        rds <- CheckIfPngExists(rds)
        if (length(rds) > 0) {
            sourceFuns2()
            grabRDs1(rds, token)
            myDt <- GetSampleList(rds)
            TryCnvMaker(myDt)
            try(gb$copyOutputPng(), silent=T)
        } else{
            GreenMsg("No CNV png images to generate.  Check the output directory.")
        }
    } else{
        GreenMsg("The PACT run has no cases with methylation completed.")
        GreenMsg("No CNV png images will generate.")
    }
}

# Search REDCap Worksheets for MRN Match for output -------------------------------------
loadPacks()
checkMounts()
output <- getOuputData(token, flds, inputSheet, readFlag)

# CNV PNG Creation -------------------------------------
if(output[1,1]!="NONE"){
    QueCnvMaker(output, token)
}
