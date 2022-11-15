if(!require(compiler)){install.packages("compiler", dependencies=T)}
library(compiler)
enableJIT(3)

mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/"
script.list = c("LoadInstallPackages.R", "SetRunParams.R", "MakeSampleSheet.R", "CopyInputs.R", "CopyOutput.R", "pipelineHelper.R", "CustomRuns.R")
scripts = paste0(mainHub, script.list)
lapply(scripts, function(i){message("Sourcing: ", i); devtools::source_url(i)})

mainPage = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/Report-Scripts/"
scripts <- c("ClassTables.R", "MLH1_Functions.R", "PipeLineU.R", "RedcapOutput.R", "TsneFunctions.R", "cnvggplotly.R")
scripts <- paste0(mainPage, scripts)
lapply(scripts, function(i) {devtools::source_url(i); return(message("Sourced: ", i))})
