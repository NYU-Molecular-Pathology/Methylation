library("base")
args <- commandArgs(TRUE)

# Parameters Input trailing commandline -------------------------------------------------------
args[1] -> pactRunName
args[2] -> concensusDir

if(is.na(concensusDir)){
    concensusDir <- paste0("/Volumes/CBioinformatics/jonathan/pact/consensus/", pactRunName, "_consensus")
}

molecOut <- "/Volumes/molecular/MOLECULAR LAB ONLY/NYU PACT Patient Data/Results/Bioinformatics"
csvPath <- file.path(molecOut, format(Sys.Date(), "%Y"), pactRunName)

stopifnot(dir.exists(csvPath))

if(!require("stringr", warn.conflicts = F)){install.packages("stringr", dependencies = T)}
if(!require("dplyr", warn.conflicts = F)){install.packages("dplyr", dependencies = T)}
library("dplyr", mask.ok = T)

#monthDay <- gsub("0(\\d)", "\\1", format(Sys.Date(), "%m_%d_"))
#thisYear <- format(Sys.Date(), "%Y.csv")
#thisFile <- paste0("export_mytable_", paste0(monthDay, thisYear))
#inputFile <- dir(path = file.path(fs::path_home(), "Downloads"), pattern = thisFile, full.names = T)
#path = file.path(fs::path_home(), "Downloads")

thisFile <- "somatic_variants.csv"
inputFile <- inputFile <- dir(path = csvPath, pattern = thisFile, full.names = T)

if(length(inputFile)==0){
    message("File not found: ", thisFile, "\n" , "Ensure the file was copied here:")
    message(inputFile)
}

variantsData <-as.data.frame(read.csv(inputFile, strip.white = T))

varColumns <- c("Test_Number", "Gene.refGene", "Variant", "ExonicFunc.refGene")
exonicFilter <- stringr::str_detect(variantsData$ExonicFunc.refGene, "frame|delet|insert")
indelsList <- variantsData[exonicFilter, varColumns]
indelsList$Gene.refGene <- stringr::str_replace_all(indelsList$Gene.refGene, ",", " ")

positions <- stringr::str_split_fixed(indelsList$Variant, ":", 3)[,1:2]
indelsList$Position <- paste(positions[,1], positions[,2], sep="_")
indelsList$nyu <- "Yes"
indelsList$philips <- ""
indelsList$Variant <- "SNV"
indelsList <- indelsList %>% dplyr::distinct(Position, .keep_all =T)

for(ngs in unique(variantsData$Test_Number)){
    if(!(ngs %in% unique(indelsList$Test_Number))){
        newSNVRow <- c(ngs, "No Indels or FrameShifts called in-house", "", "", "", "", "SNV")
        indelsList[nrow(indelsList) + 1,] = newSNVRow
    }
    newCnvRow <- c(ngs, "", "", "", "NA", "", "", "", "CNV")
    indelsList[nrow(indelsList) + 1,] = newCnvRow
}

varsToCheck <- data.frame(
    "Test_Case" = indelsList$Test_Number,
    "Gene" = indelsList$Gene.refGene,
    "Mutation Type" = indelsList$ExonicFunc.refGene,
    "Other" = indelsList$Position,
    "In NYU" = indelsList$nyu,
    "In Philips" = indelsList$philips,
    "IGV" = '',
    "Comments" = '',
    "Variant" = indelsList$Variant
)

#outPutFile <- file.path(fs::path_home(), "Desktop", paste0(pactRunName, "_desc.csv"))
outPutFile <- file.path(concensusDir, paste0(pactRunName, "_desc.csv"))
write.csv(varsToCheck, file=outPutFile, quote=F, row.names=F)
message("Output file saved: '", outPutFile, "'")
