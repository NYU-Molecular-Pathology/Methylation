writeTargets <- function(targets, csvFi="samplesheet.csv"){
  write.csv(targets, file, quote=F, row.names=F)
}


setKeyColumns <- function(targets, col_samTypes, col_samNames, col_other){
  targets$Type <- targets[ , col_samTypes] # Creates any new "Type" column
  targets$Sample_Name <- targets$Sample_ID <- targets[ , col_samNames] # generates Sample_ID column if doesn't exist
  targets$Sample_Group <- targets[,col_other]
}
