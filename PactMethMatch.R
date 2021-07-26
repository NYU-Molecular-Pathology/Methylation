#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library("base"); gb <- globalenv(); assign("gb", gb)
token=NULL;inputSheet=NULL
# Main Parameters trailing commandline
token<-args[1]; inputSheet<-args[2]
# Check Input Params -----
message("====Parameters input====\ntoken: ",token,"\n","Worksheet: ", inputSheet,"\n"); stopifnot(!is.null(token))

# REDcap Heading Fields -----
flds = c("record_id","b_number","tm_number","accession_number","block","diagnosis","organ","tissue_comments","run_number")

# Load redcapAPI Package -----
if(suppressWarnings(!require("redcapAPI"))){
    params=list('nutterb/redcapAPI', dependencies=T, upgrade="always", type="source")
    do.call(devtools::install_github,c(params))
}

# Functions to load packages and get redcap info -----
loadPacks <- function(){
    pkgs = c("data.table", "foreach", "xlsx","jsonlite","RCurl")
    rlis = getOption("repos")
    rlis["CRAN"] = "http://cran.us.r-project.org"
    options(repos = rlis)
   invisible(lapply(pkgs, function(pk){
        if(suppressWarnings(!require(pk, character.only=T))){

        install.packages(pk,dependencies=T, verbose=T, repos = "http://cran.us.r-project.org", type="source")}}))
}

# API Call functions -----
grabAllRecords <- function(token, flds, rcon){
    params = list(rcon, fields=flds, labels=F, dates = F, survey = F, dag = F, factors=F, form_complete_auto=F)
    dbColumns <- do.call(redcapAPI::exportRecords, c(params))
    return(as.data.frame(dbColumns))
}

# Database search function -----
searchDb <- function(vals, db){
    vals2find <- paste(vals, collapse="|"); i=NULL
    return(foreach::foreach(i=1:ncol(db), .combine='rbind') %do% {db[grepl(vals2find,db[,i]),]})
}

# Import csv Worksheet -----
vals2find <- utils::read.csv(inputSheet, skip=19)[,c(6,7)]
vals2find <- vals2find[!grepl("H20|SERACARE|HAPMAP", vals2find[,2]),]

loadPacks()

# Get Methylation and Molecular Samples list ----
rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", token)
db <- grabAllRecords(token, flds, rcon)
query1 <- vals2find[,1][vals2find[,1]!=0 & vals2find[,1] !=""]
print(query1)
methResA <- searchDb(strtrim(query1, 10), db)
query2 <- vals2find[,2][vals2find[,2]!=0 & vals2find[,2] !=""]
print(query2)
methResB <- searchDb(query2, db)
output <- unique(rbind(methResA, methResB))
runId <- paste0(head(read.csv(inputSheet))[3,2])
winpath = "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular/MOLECULAR/NYU-METHYLATION/Results/"
output$report_complete <- ifelse(!is.na(output$run_number), "YES", "NOT_YET_RUN")
output$reportPath <- paste0(winpath,format(Sys.Date(),"%Y"),"/",output$run_number,"/",output$record_id,".html")
output$reportPath[is.na(output$run_number)] <- "NOT_YET_RUN"
outFi <- paste0(runId,"_MethylMatch.xlsx")
xlsx::write.xlsx(output, file.path("~/Desktop",outFi), sheetName = runId, row.names = F)

record = data.frame(record_id = runId, run_number = runId)
datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
RCurl::postForm(rcon$url,token=rcon$token,content='record',format='json',type='flat',data=datarecord)
redcapAPI::importFiles(rcon,file= file.path("~/Desktop",outFi), runId, field="other_file", repeat_instance=1)

record = data.frame(record_id = runId, comments="pact_sample_list_email")
datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
RCurl::postForm(rcon$url,token=rcon$token,content='record',format='json',type='flat',data=datarecord)
message("====Email Notification Created====")
