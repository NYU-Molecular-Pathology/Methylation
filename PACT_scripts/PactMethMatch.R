#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library("base"); gb <- globalenv(); assign("gb", gb)
dsh="\n================"
dsh2="================\n"

# Main arguments input in comandline
#token<-args[1]; inputSheet<-args[2]
token<-args[1]; inputSheet<-args[2]

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
if(suppressWarnings(!require("redcapAPI"))){
    params=list('nutterb/redcapAPI', dependencies=T, upgrade="always", type="source")
    do.call(devtools::install_github,c(params))
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
    pkgs = c("data.table", "foreach", "openxlsx","jsonlite","RCurl","readxl","stringr","webshot","tidyverse","crayon","plotly","htmlwidgets","tinytex")
    rlis = getOption("repos")
    rlis["CRAN"] = "http://cran.us.r-project.org"
    options(repos = rlis)
    invisible(lapply(pkgs, function(pk){
        if(suppressWarnings(!require(pk, character.only=T))){
            install.packages(pk,dependencies=T, verbose=T, repos="http://cran.us.r-project.org", type="both")
        }}))
    if(!require("systemfonts")){devtools::install_github('r-lib/systemfonts')}
    if(!require("redcapAPI")){install.packages("redcapAPI", dependencies = T, type="both",ask=F)}
    if(!require("remotes")){install.packages("remotes", dependencies=T)}
    if(!require("chromote")){remotes::install_github("rstudio/chromote", upgrade ="never")}
    if(!require("webshot2")){remotes::install_github("rstudio/webshot2")}
    library("redcapAPI")
    library("chromote")
    #try(tinytex::tlmgr_install(),silent=T)
    #try(tinytex:::install_prebuilt(),silent=T)
}

# API Call functions -----
grabAllRecords <- function(flds, rcon){
    params = list(rcon, fields=flds, labels=F, dates = F, survey = F, dag = F, factors=F, form_complete_auto=F)
    dbCols <- do.call(redcapAPI::exportRecords, c(params))
    return(as.data.frame(dbCols))
}

# Database search function -----
searchDb <- function(vals, db){
    v2f <- paste(vals, collapse="|"); i=NULL
    return(foreach::foreach(i=1:ncol(db), .combine='rbind') %do% {db[grepl(v2f,db[,i]),]})
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
    inputFi <- file.path(runFolder,paste0(inputSheet,".xlsm"))
    if(file.exists(inputFi)) {return(inputFi)} else{return(getAltPath(inputFi))}
}

# Parses the input file for the "PhilipsExport" tab
parseWorksheet <- function(inputFi){
    shNames <- readxl::excel_sheets(inputFi)
    message("Excel sheet names:"); print(shNames)
    sheet2Read <- "PhilipsExport"
    sh <- which(grepl(sheet2Read,shNames,ignore.case=T))[1]
    pactShCol <- c("Tumor Specimen ID","Tumor DNA/RNA Number", "MRN", "Test Number")
    vals2find <-  as.data.frame(readxl::read_excel(inputFi, sheet=shNames[sh], skip=3, col_types ="text")[,pactShCol])
    vals2find <- vals2find[!is.na(vals2find[,1]),]
    return(vals2find)
}

# Import csv Worksheet -----
getCaseValues <- function(inputSheet, readFlag){
    if(readFlag){
        vals2find <- utils::read.csv(inputSheet, skip=19)[,c(6,7,9)]
        vals2find <- vals2find[!grepl("H20|SERACARE|HAPMAP", vals2find[,2]),]
        return(vals2find)
    }else{
        inputFi <- getFilePath(inputSheet)
        vals2find <- parseWorksheet(inputFi)
        return(vals2find)
    }
}

genQuery <- function(dbCol,vals2find){
    currCol <- vals2find[, dbCol]
    q1 <- currCol[currCol != 0 & !is.na(currCol) & currCol != ""]
    return(unique(q1))
}

# Get Methylation and Molecular Samples list ----
queryCases <- function(vals2find, db) {
    i=NULL
    queryList <- foreach::foreach(i=1:ncol(vals2find), .combine="c") %do% {genQuery(i,vals2find)}
    methQuery <- searchDb(queryList, db)
    return(unique(methQuery))
}

addOutputLinks <- function(output){
    winpath = "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular/Molecular/MethylationClassifier/"
    yearPath <- lapply(stringr::str_split_fixed(output$run_number,"-",2)[,1], function(yr) {
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

modifyOutput <- function(output,vals2find){
    output$Test_Number <- NA
    for(i in 1:nrow(output)) {
        theVal = NA
        for (var in 1:ncol(vals2find)) {
            pat <- vals2find[, var]
            currRow <- paste(output[i, ])
            theMatch <- which(pat %in% currRow[currRow != "0" & currRow != "NA"])
            if (length(theMatch) > 0) {
                theVal <- vals2find$`Test Number`[theMatch]
            }
            }
        output$Test_Number[i] <- theVal
        }
    output <- addOutputLinks(output)
    return(output)
}

addExcelLink <- function(output,fiLn,wb,runId){
    x <- c(output$'Report Link'[fiLn])
    names(x) <- paste0(output$record_id[fiLn],".html")
    class(x) <- "hyperlink"
    colNum <- which(colnames(output)=="Report Link")
    openxlsx::writeData(wb, sheet=runId, x=x, startCol=colNum, startRow=fiLn+1)
}

createXlFile <- function(runId,output){
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, runId)
    openxlsx::writeData(wb, sheet=runId, x=output)
    for (fiLn in 1:length(output$'Report Link')) {
        if(output$'Report Link'[fiLn]!=''){addExcelLink(output,fiLn,wb,runId)}
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
    suppressMessages(redcapAPI::importFiles(rcon,file= outFi, runId, field="other_file", repeat_instance=1))
    record$comments <- "pact_sample_list_email"
    postData(rcon, record)
    message(dsh,"\nEmail Notification Created",dsh,"\n")
}

# Main workflow to get dataframe
getOuputData <- function(token, flds, inputSheet, readFlag){
    apiUrl="https://redcap.nyumc.org/apps/redcap/api/"
    rcon <- redcapAPI::redcapConnection(apiUrl, token)
    vals2find <- getCaseValues(inputSheet, readFlag)
    db <- grabAllRecords(flds, rcon)
    output <- queryCases(vals2find, db)
    runId <- ifelse(readFlag == T, substr(inputSheet,1,nchar(inputSheet)-4), inputSheet)
    output <- modifyOutput(output, vals2find)
    outFi <- createXlFile(runId, output)
    emailFile(runId, outFi, rcon)
    return(output)
}

# FUN: Sets your directory and sources the helper functions
sourceFuns2 <- function(workingPath = NULL) {
    mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/"
    script.list <- c("SetRunParams.R","CopyInputs.R","PACT_scripts/generateCNV.R")
    if (is.null(workingPath)) {workingPath = getwd()}
    scripts <- paste0(mainHub, script.list)
    invisible(lapply(scripts, function(i){devtools::source_url(i)}))
    library("sest")
    library("mnp.v11b6")
    require("plotly")
    require("htmlwidgets")
    gb$setDirectory(workingPath)
    return(gb$defineParams())
}

startCNVmaker <- function(output, token) {
    rds <- output$record_id[output$report_complete == "YES"]
    if (length(rds) > 0) {
        assign("rds", rds)
        message(dsh, crayon::bgMagenta("Starting CNV PNG Creation"),dsh2)
        sourceFuns2()
        ApiToken = gb$token
        tryCatch(
            expr = {gb$save.png.files(gb$rds, token)},
            error = function(e) {
                message("The following error occured:\n",e)
                message("\n\nTry checking the troubleshooting section on GitHub:\n")
                message("https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/PACT_scripts/README.md\n")
            },
                 finally={gb$copyOutputPng()}
        )
    } else{message(crayon::bgGreen("This PACT run has no cases with methylation complete.  No CNV png images to generate."))}
}

# Search REDCap Worksheets for MRN Match for output -------------------------------------
loadPacks()
checkMounts()
output <- getOuputData(token, flds, inputSheet, readFlag)

# CNV PNG Creation -------------------------------------
startCNVmaker(output, token)
