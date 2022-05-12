writeTargets <- function(targets, csvFi="samplesheet.csv"){
  write.csv(targets, file, quote=F, row.names=F)
}


setKeyColumns <- function(targets, col_samTypes, col_samNames, col_other){
  targets$Type <- targets[ , col_samTypes] # Creates any new "Type" column
  targets$Sample_Name <- targets$Sample_ID <- targets[ , col_samNames] # generates Sample_ID column if doesn't exist
  targets$Sample_Group <- targets[,col_other]
}

checkSamNam <- function(samNames, targets){
  if(anyDuplicated(samNames) > 0){
  warning("Sample IDs contain duplicates, assigning unique names")
  samNames <- targets$Sample_ID = make.unique(samNames, sep="_")
  return(samNames) }else{return(samNames)}
}

getTargCsv <- function(csvFi = "samplesheet.csv") {
  gc(verbose = F)
  targets <- read.csv(csvFi, strip.white = T, row.names = NULL)
  if (class(targets) != "data.frame") {targets <- as.data.frame(targets)}
  return(targets)
}
