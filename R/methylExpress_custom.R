#!/usr/bin/env Rscript
library("base"); gb <- globalenv(); assign("gb", gb)
if(!require("devtools")){install.packages("devtools", quiet = T)}

# Input Parameters --------------------------------------------------
rd_numbers = c("RD-17-140") # Default NULL, a list of any individual RD-numbers to run in REDCap
rd_csv = "/Users/NAME/Desktop/my_rdnumbers.csv" # path to a .csv file with no header listing RD-numbers in column 1

# Run Flags ------------------------------------------------------------------------------------------
runSarcoma   <- F # To Generate Sarcoma Classifier? If FALSE, Brain Classifier report will output
makeNewSheet <- T # To Generate a new samplesheet or use existing Samplesheet in current directory
forcedUpload <- T # To upload an overwrite existing html report in REDCap with html in current directory

# Default Variables --------------------------------------------------
runFolder <- getwd() # Can be modified to any other directory to copy idats and run from
runID <- "MR23-rerun_v12" # Anything with MR##-, as a new MGDM will trigger an email notification
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" -> token -> gb$token -> gb$ApiToken # Replace with REDCap API token

runLocal=T; selectRDs = NULL; baseFolder = NULL; redcapUp  <- F; setwd(runFolder)
ghLink <- "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R"
devtools::source_url(file.path(ghLink, "Report-Scripts/SourceLoadGitHub.R"))

if(makeNewSheet==T){gb$MakeLocalSampleSheet(runID, token, rd_csv, rd_numbers)}

if(generateSarcReport == T) {
    gb$MakeSarcomaReport()
} else{
    gb$StartCustomRun(redcapUp = redcapUp)
}

if(forcedUpload ==T) {
    gb$ForceUploadToRedcap(file.list = dir(runFolder, pattern = ".html", full.names = T), token, F)
}
