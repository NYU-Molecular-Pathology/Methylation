library("dplyr")
inputFile <- sort(dir(path="~/Downloads", pattern="export_mytable",full.names=T),decreasing=T)[1]
variantsData <- as.data.frame(read.csv(inputFile, strip.white = T))
varColumns <- c("Test_Number", "Gene.refGene", "IGV", "ExonicFunc.refGene")
exonicFilter <- stringr::str_detect(variantsData$ExonicFunc.refGene, "frame|delet|insert")
indelsList <- variantsData[exonicFilter, varColumns]
indelsList$Gene.refGene <- stringr::str_replace_all(indelsList$Gene.refGene, ",", " ")
indelsList$IGV <- stringr::str_replace_all(indelsList$IGV, ".png", "")
indelsList$nyu <- "Yes"
indelsList$philips <- ""
indelsList$Variant <- "SNV"
indelsList <- indelsList %>% dplyr::distinct(IGV, .keep_all =T)
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
write.csv(varsToCheck, file="~/Desktop/Pact_desc.csv", quote=F, row.names=F)
