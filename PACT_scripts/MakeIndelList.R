library("base")
args <- commandArgs(TRUE)

# Parameters Input trailing commandline -------------------------------------------------------
args[1] -> pactRunName

if(!require("stringr")){install.packages("stringr", dependencies = T)}
if(!require("dplyr")){install.packages("dplyr", dependencies = T)}
library("dplyr")

monthDay <- gsub("0(\\d)", "\\1", format(Sys.Date(), "%m_%d_"))
thisYear <- format(Sys.Date(), "%Y.csv")

thisFile <- paste0("export_mytable_", paste0(monthDay, thisYear))

inputFile <- dir(path = file.path(fs::path_home(), "Downloads"), pattern = thisFile, full.names = T)

if(length(inputFile)==0){
    message("File not found: ", thisFile, "\n" ,"\nTrying other files in downloads folder:")
    inputFile <- sort(dir(
        path = file.path(fs::path_home(), "Downloads"), pattern = "export_mytable", full.names = T),
        decreasing = T)[1]
    message(inputFile)
}

variantsData <-as.data.frame(read.csv(inputFile, strip.white = T))

varColumns <- c("Test_Number", "Gene.refGene", "IGV", "ExonicFunc.refGene")
exonicFilter <- stringr::str_detect(variantsData$ExonicFunc.refGene, "frame|delet|insert")
indelsList <- variantsData[exonicFilter, varColumns]
indelsList$Gene.refGene <- stringr::str_replace_all(indelsList$Gene.refGene, ",", " ")
indelsList$IGV <- stringr::str_replace_all(indelsList$IGV, ".png", "")
indelsList$nyu <- "Yes"
indelsList$philips <- ""
indelsList$Variant <- "SNV"
indelsList <- indelsList %>% dplyr::distinct(IGV, .keep_all =T)

for(ngs in unique(variantsData$Test_Number)){
    if(!(ngs %in% unique(indelsList$Test_Number))){
        newSNVRow <- c(ngs, "No Indels or FrameShifts called in-house", "", "", "", "", "SNV")
        indelsList[nrow(indelsList) + 1,] = newSNVRow
    }
    newCnvRow <- c(ngs, "", "", "", "NA", "", "CNV")
    indelsList[nrow(indelsList) + 1,] = newCnvRow
}

varsToCheck <- data.frame(
    "Test_Case" = indelsList$Test_Number,
    "Gene" = indelsList$Gene.refGene,
    "Mutation Type" = indelsList$ExonicFunc.refGene,
    "Other" = indelsList$IGV,
    "In NYU" = indelsList$nyu,
    "In Philips" = indelsList$philips,
    "IGV" = '',
    "Comments" = '',
    "Variant" = indelsList$Variant
)

outPutFile <- file.path(fs::path_home(), "Desktop", paste0(pactRunName, "_desc.csv"))

write.csv(varsToCheck, file=outPutFile, quote=F, row.names=F)
message("Output file saved: '", outPutFile, "'\n\n", "Copy this file into your consensus Rmd directory")
