#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library("base"); gb <- globalenv(); assign("gb", gb)
dsh="\n================"

# Main arguments input in comandline
#token<-args[1]; inputSheet<-args[2]
token<-args[1]; inputSheet<-args[2]

# Displays the Input args -----
message(dsh,"\nParameters input",dsh)
message("token: ",token)
message("inputSheet: ", inputSheet,"\n")

stopifnot(!is.na(token))
stopifnot(!is.na(inputSheet))

readFlag <- endsWith(inputSheet,".csv")==T
stopifnot(rlang::is_bool(readFlag))

# REDcap Heading Fields -----
flds = c("record_id","b_number","tm_number","accession_number","block","diagnosis","organ","tissue_comments","run_number")

# Load redcapAPI Package -----
if(suppressWarnings(!require("redcapAPI"))){
    params=list('nutterb/redcapAPI', dependencies=T, upgrade="always", type="source")
    do.call(devtools::install_github,c(params))
}

# Functions to load packages and get redcap info -----
loadPacks <- function(){
    pkgs = c("data.table", "foreach", "openxlsx","jsonlite","RCurl","readxl","stringr")
    rlis = getOption("repos")
    rlis["CRAN"] = "http://cran.us.r-project.org"
    options(repos = rlis)
    invisible(lapply(pkgs, function(pk){
        if(suppressWarnings(!require(pk, character.only=T))){
            install.packages(pk,dependencies=T, verbose=T, repos="http://cran.us.r-project.org", type="both")
        }}))
    if(!require("webshot2")){remotes::install_github("rstudio/webshot2")}
}

# API Call functions -----
grabAllRecords <- function(flds, rcon){
    params = list(rcon, fields=flds, labels=F, dates = F, survey = F, dag = F, factors=F, form_complete_auto=F)
    dbColumns <- do.call(redcapAPI::exportRecords, c(params))
    return(as.data.frame(dbColumns))
}

# Database search function -----
searchDb <- function(vals, db){
    v2f <- paste(vals, collapse="|"); i=NULL
    return(foreach::foreach(i=1:ncol(db), .combine='rbind') %do% {db[grepl(v2f,db[,i]),]})
}

getFilePath <- function(inputSheet, normFlag=F){
    drive = file.path("", "Volumes", "molecular", "MOLECULAR LAB ONLY")
    folder <- file.path("NYU PACT Patient Data", "Workbook")
    runyr <- stringr::str_split_fixed(inputSheet,"-",3)[,2]
    if(normFlag==F){
    return(file.path(drive, folder, paste0("20", runyr),inputSheet,paste0(inputSheet,".xlsx")))
    }else{
        return(file.path(drive, folder, paste0("20", runyr),inputSheet))
        }
}

parseWorksheet <- function(inputFi){
    shNames <- readxl::excel_sheets(inputFi)
    sh <- which(grepl("SampleSheet",shNames,ignore.case=T))[1]
    vals2find <-  as.data.frame(readxl::read_excel(inputFi, sheet=shNames[sh], skip=19, col_types ="text")[,c(6,7,9)])
    vals2find <- vals2find[!grepl("H20|SERACARE|HAPMAP", vals2find[,2]),]
    vals2find <- vals2find[!is.na(vals2find[,1]),]
    return(vals2find)
}

# Import csv Worksheet -----
getCaseValues <- function(inputSheet,readFlag){
    if(readFlag){
        vals2find <- utils::read.csv(inputSheet, skip=19)[,c(6,7,9)]
        vals2find <- vals2find[!grepl("H20|SERACARE|HAPMAP", vals2find[,2]),]
        return(vals2find)
    }else{
        inputFi <- getFilePath(inputSheet)
        if(file.exists(inputFi)){
            vals2find <- parseWorksheet(inputFi)
            return(vals2find)
        }else{
            message(dsh,"The PACT run worksheet was not found:\n", inputFi, dsh)
            inputFi <- getFilePath(inputSheet, T)
            allFi <- list.files(path=inputFi, pattern="*.xlsx")
            allFi <- allFi[!grepl( "~$", allFi, fixed = T)]
            allFi <- allFi[!grepl( "export", allFi, fixed = F,ignore.case=T)]
            allFi <- allFi[grepl( "Book", allFi, fixed = F,ignore.case=T)]
            if(length(allFi)>1){allFi <- allFi[1]}
            inputFi <- file.path(inputFi,allFi)
            vals2find <- parseWorksheet(inputFi)
            return(vals2find)
            }
    }
}

# Get Methylation and Molecular Samples list ----
queryCases <- function(vals2find, db) {
    query1 <- vals2find[, 1][vals2find[, 1] != 0 & vals2find[, 1] != ""]
    methResA <- searchDb(query1, db)
    query2 <- vals2find[, 2][vals2find[, 2] != 0 & vals2find[, 2] != ""]
    message("\n~~~ Queried Cases:\n")
    print(query1)
    print(query2)
    methResB <- searchDb(query2, db)
    return(unique(rbind(methResA, methResB)))
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
                theVal <- vals2find$Test_Number[theMatch]
            }
        }
        output$Test_Number[i] <- theVal
    }
    output <- addOutputLinks(output)
    return(output)
}

createXlFile <- function(runId,output){
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, runId)
    openxlsx::writeData(wb, sheet=runId, x= output)
    for (fiLn in 1:length(output$'Report Link')) {
        if(output$'Report Link'[fiLn]!=''){
            x <- c(output$'Report Link'[fiLn])
            names(x) <- paste0(output$record_id[fiLn],".html")
            class(x) <- "hyperlink"
            openxlsx::writeData(wb, sheet = runId, x = x, startCol = 11, startRow=fiLn+1)
        }
    }
    outFi <-file.path("~","Desktop", paste0(runId,"_MethylMatch.xlsx"))
    openxlsx::saveWorkbook(wb, outFi, overwrite = T)
    return(outFi)
}

postData <- function(rcon, record){
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
    res <- RCurl::postForm(rcon$url,token=rcon$token,content='record',format='json',type='flat',
                           data=datarecord, returnContent = 'ids', returnFormat = 'csv')
    cat(res,sep="\n")
}

emailFile <- function(runId, outFi, rcon){
    record = data.frame(record_id = runId, run_number = runId)
    postData(rcon, record)
    suppressMessages(redcapAPI::importFiles(rcon,file= outFi, runId, field="other_file", repeat_instance=1))
    record$comments <- "pact_sample_list_email"
    postData(rcon, record)
    message(dsh,"\nEmail Notification Created",dsh)
}

loadPacks()
rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", token)
vals2find <- getCaseValues(inputSheet,readFlag)
db <- grabAllRecords(flds, rcon)
output <- queryCases(vals2find, db)
runId <- ifelse(readFlag==T, paste0(head(read.csv(inputSheet))[3,2]),
                paste0(head(suppressMessages(readxl::read_excel(getFilePath(inputSheet), sheet=7)))[3,2]))
output <- modifyOutput(output,vals2find)
outFi <- createXlFile(runId,output)
emailFile(runId, outFi, rcon)

#db[db$record_id %in% output$record_id,]
rds <- output$record_id
assign("rds", rds)

supM <- function(sobj){return(suppressMessages(suppressWarnings(sobj)))}

# FUN: Sets your directory and sources the helper functions for cnv
sourceFuns2 <- function(workingPath = NULL) {
    mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/"
    script.list <- c("SetRunParams.R","CopyInputs.R","PACT_scripts/generateCNV.R")
    if(!require("redcapAPI")){install.packages("redcapAPI", dependencies = T, type="both",ask=F)
        library("redcapAPI")}
    if (is.null(workingPath)) {workingPath = getwd()}
    scripts <- paste0(mainHub, script.list)
    invisible(lapply(scripts, function(i){devtools::source_url(i)}))
    gb$setDirectory(workingPath)
    return(gb$defineParams())
}

library("conumee");library("sest");library("mnp.v11b6")
require("plotly");require("htmlwidgets")
sourceFuns2()
ApiToken = gb$token
gb$save.png.files(gb$rds, token)
