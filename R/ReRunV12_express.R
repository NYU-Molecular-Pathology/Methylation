#!/usr/bin/env Rscript
library("base"); gb <- globalenv(); assign("gb", gb)
args <- commandArgs(T)
if(!require("devtools", warn.conflicts = F)) {
    install.packages("devtools", quiet = T)
}

try(args[1] -> token)
try(args[2] -> samSheetIn)
try(args[3] -> generateSarcReport)

runID <- "MR23-rerun_v12"

if(!length(generateSarcReport)>0){
    generateSarcReport <- F
}else{
    sarcRepFlag <- as.character(generateSarcReport)
    generateSarcReport <- ifelse(sarcRepFlag=="T", T, F)
}
gb$token <- gb$ApiToken <- token

# Default Variables ---------------------------------------------------------------
runLocal=T; selectRDs = NULL; baseFolder = NULL; redcapUp <- redcapUpload <- F
reportMd <- file.path(fs::path_home(),"report.Rmd")
assign("token", token, envir = gb); assign("ApiToken", token, envir = gb)

ghLink <- "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R"
script <- "Report-Scripts/SourceLoadGitHub.R"
devtools::source_url(file.path(ghLink, script))

gb$SetBaseFolder(token, baseFolder, runID)
if(file.exists("samplesheet.csv")){file.remove("samplesheet.csv")}

gb$MakeLocalSampleSheet(runID, token, samSheetIn)

selectRDs <- gb$AssignArgs(runID, baseFolder, token=token, selectRDs, redcapUp, gb)
gb$StartCustomRun(redcapUp=redcapUp)

if(generateSarcReport==T) {
  targets <- as.data.frame(read.csv("samplesheet.csv"))
  gb$MakeSarcomaReport(targets = targets)
}

file.list <- dir(getwd(), pattern = ".html", full.names = T)
print(as.data.frame(file.list))

try(args[1] -> token)
gb$ForceUploadToRedcap(file.list, token, F)
)
