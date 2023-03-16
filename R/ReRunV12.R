#!/usr/bin/env Rscript
library("base"); gb <- globalenv(); assign("gb", gb)
if(!require("devtools")){install.packages("devtools", quiet=T)}

runID <- "MR23-rerun_v12"
generateSarcReport <- F
gb$token <- gb$ApiToken <- token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"

# Default Variables ---------------------------------------------------------------
runLocal=T; selectRDs = NULL; baseFolder = NULL; redcapUp <- redcapUpload <- F
reportMd <- file.path(fs::path_home(),"report.Rmd")
assign("token", token, envir = gb); assign("ApiToken", token, envir = gb)

ghLink <- "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R"
script <- "Report-Scripts/SourceLoadGitHub.R"
devtools::source_url(file.path(ghLink, script))

gb$SetBaseFolder(token, baseFolder, runID)
if(file.exists("samplesheet.csv")){file.remove("samplesheet.csv")}

gb$MakeLocalSampleSheet(runID, token)

selectRDs <- gb$AssignArgs(runID, baseFolder, token, selectRDs, redcapUp, gb)
gb$StartCustomRun(redcapUp=redcapUp)

if(generateSarcReport==T) {
  targets <- as.data.frame(read.csv("samplesheet.csv"))
  gb$MakeSarcomaReport(targets = targets)
}

file.list <- dir(getwd(), pattern = ".html", full.names = T)
print(as.data.frame(file.list))
gb$ForceUploadToRedcap(file.list, F)
