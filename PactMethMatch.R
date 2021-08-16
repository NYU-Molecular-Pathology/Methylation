#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library("base"); gb <- globalenv(); assign("gb", gb)
dsh="\n================"

# Main arguments input in comandline
#token<-args[1]; inputSheet<-args[2]
token<-args[2]; inputSheet<-args[1]

# Displays the Input args -----
message(dsh,"\nParameters input",dsh)
message("token: ",token)
readFlag <- endsWith(inputSheet,".csv")==T
library("rlang")
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
}

# API Call functions -----
grabAllRecords <- function(token, flds, rcon){
    params = list(rcon, fields=flds, labels=F, dates = F, survey = F, dag = F, factors=F, form_complete_auto=F)
    dbColumns <- do.call(redcapAPI::exportRecords, c(params))
    return(as.data.frame(dbColumns))
}

# Database search function -----
searchDb <- function(vals, db){
    v2f <- paste(vals, collapse="|"); i=NULL
    return(foreach::foreach(i=1:ncol(db), .combine='rbind') %do% {db[grepl(v2f,db[,i]),]})
}

# Import csv Worksheet -----
inputFi = NULL
if(readFlag){
    vals2find <- utils::read.csv(inputSheet, skip=19)[,c(6,7,9)]
    vals2find <- vals2find[!grepl("H20|SERACARE|HAPMAP", vals2find[,2]),]
}else{
    drive = file.path("", "Volumes", "molecular", "MOLECULAR LAB ONLY")
    folder <- file.path("NYU PACT Patient Data", "Workbook")
    runyr <- stringr::str_split_fixed(inputSheet,"-",3)[,2]
    inputFi <- file.path(drive, folder, paste0("20", runyr),inputSheet,paste0(inputSheet,".xlsx"))
    if(file.exists(inputFi)){
        shNames <- readxl::excel_sheets(inputFi)
        sh <- which(grepl("SampleSheet",shNames,ignore.case=T))[1]
        vals2find <-  as.data.frame(readxl::read_excel(inputFi, sheet=shNames[sh], skip=19, col_types ="text")[,c(6,7,9)])
        vals2find <- vals2find[!grepl("H20|SERACARE|HAPMAP", vals2find[,2]),]
        vals2find <- vals2find[!is.na(vals2find[,1]),]
    }else{message("The PACT run worksheet was not found:\n",inputFi)}
}
loadPacks()

# Get Methylation and Molecular Samples list ----
rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", token)
db <- grabAllRecords(token, flds, rcon)
query1 <- vals2find[,1][vals2find[,1]!=0 & vals2find[,1] !=""]
message("\n~~~ Queried Cases:\n")
print(query1)
methResA <- searchDb(query1, db)
query2 <- vals2find[,2][vals2find[,2]!=0 & vals2find[,2] !=""]
print(query2)
methResB <- searchDb(query2, db)
output <- unique(rbind(methResA, methResB))

if(readFlag){
    runId <- paste0(head(read.csv(inputSheet))[3,2])
}else{
    runId <- paste0(head(readxl::read_excel(inputFi, sheet=7))[3,2])
}
message("======",runId,"======")
yearPat <- stringr::str_split_fixed(output$run_number,"-",2)[,1]
yearPath <- lapply(yearPat, function(yr) {
    rnum <- NULL
    if(nchar(yr)>2){rnum <- substring(yr, 3)}else{rnum <- yr}
    if(nchar(yr)>0){paste0("20",rnum)}else{rnum}
})
winpath = "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular/Molecular/MethylationClassifier/"
ngsNumb <- vals2find[,c(1,2,3)]

v2f <- lapply(X=c(1,2,4), FUN=function(X){
    newval<- output[,X]
    newval<- newval[!is.na(newval)&newval!=""&newval!="0"]
    c(paste(newval,collapse="|"))
})
output$Test_Number <- NA
for(i in 1:nrow(output)){
    theVal = NA
    for (var in 1:ncol(vals2find)) {
        pat<-vals2find[,var]
        currRow <- paste(output[i,])
        theMatch <- which(pat %in% currRow[currRow!="0"&currRow!="NA"])
        if(length(theMatch)>0){
            theVal <- vals2find$Test_Number[theMatch]
        }
    }
    output$Test_Number[i] <- theVal
}

output$report_complete <- ifelse(!is.na(output$run_number), "YES", "NOT_YET_RUN")
output$'Report Link' <- paste0(winpath, yearPath,"/",output$run_number,"/",output$record_id,".html")
output$'Report Link'[is.na(output$run_number)] <- ""
output$'Report Path'<-output$'Report Link'
outFi <- paste0(runId,"_MethylMatch.xlsx")

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, runId)
openxlsx::writeData(wb, sheet=runId, x= output)
for (fiLn in 1:length(output$'Report Link')) {
    if(output$'Report Link'[fiLn]!=''){
        x <- c(output$'Report Link'[fiLn])
        names(x) <- paste0(output$record_id[fiLn],".html")
        class(x) <- "hyperlink"
        writeData(wb, sheet = runId, x = x, startCol = 11, startRow=fiLn+1)
    }
}
openxlsx::saveWorkbook(wb, file.path("~/Desktop",outFi), overwrite = TRUE)

record = data.frame(record_id = runId, run_number = runId)
datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
res <- RCurl::postForm(rcon$url,token=rcon$token,content='record',format='json',
                       type='flat',data=datarecord, returnContent = 'count', returnFormat = 'csv')
cat(res,sep="\n")
redcapAPI::importFiles(rcon,file= file.path("~/Desktop",outFi), runId, field="other_file", repeat_instance=1)

record$comments <- "pact_sample_list_email"
datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
res<-RCurl::postForm(rcon$url,token=rcon$token,content='record',format='json',type='flat',
                     data=datarecord,returnContent = 'ids', returnFormat = 'csv')
cat(res)
message("\n====Email Notification Created====")
