gb <- globalenv(); assign("gb", gb)

dfTargets <- function(targets){
  if(class(targets)!="data.frame"){targets<- as.data.frame(targets)}
  return(targets)
  }

writeTargets<- function(targets, csvFi="samplesheet.csv"){
  targets <- dfTargets(targets)
  write.csv(targets, file=csvFi, quote=F, row.names=F)
}

CheckSamNames <- function(samNames, targets){
  if(anyDuplicated(samNames) > 0){
      warning("Sample IDs contain duplicates, assigning unique names")
      samNames <- targets$Sample_ID = make.unique(samNames, sep="_")
      targets$Sample_Name <- targets$Sample_ID <- samNames
    }
  return(targets)
}

SetKeyColumns <- function(targets, col_samTypes, col_samNames, col_other, col_shapes, sam.grp.type=NULL) {
    targets <- dfTargets(targets)
    targets$Type <- targets[, col_samTypes] # Creates any new "Type" column
    targets$Sample_Name <- targets$Sample_ID <- targets[, col_samNames] # generates Sample_ID column if doesn't exist
    targets$Other_Group <- targets$Sample_Group <- targets[, col_other]
    targets$Sym_Shape <- targets[,col_shapes]
    targets <- CheckSamNames(targets$Sample_ID, targets)
    if (!is.null(sam.grp.type)) {
      targets$Sample_Group <- sam.grp.type
    }
    return(targets)
  }

getTargCsv <- function(csvFi = "samplesheet.csv") {
  gc(verbose = F)
  targets <- read.csv(csvFi, strip.white = T, row.names = NULL)
  targets <- dfTargets(targets)
  return(targets)
}

checkIdats<- function(samsheet, token, needFi=F){
    if(needFi==T) {
      rds <- gb$readInfo(inputSheet = samsheet) # inputSheet can be xlsx or csv
      stopifnot(length(rds)>1 & stringr::str_detect(rds[1],"RD-"))
      gb$grabRDCopyIdat(rd_numbers=rds, token, copyIdats=T) 
    }
}

setKnitDir <- function(runDir) {
  syscmd <- paste("cd", runDir)
  system(syscmd); setwd(runDir)
  knitr::opts_knit$set(root.dir = runDir)
}

MoveIdats <- function(pathName=NULL){
    if(is.null(pathName)){pathName<- file.path(getwd(),"idats")}
    if(!dir.exists(pathName)){dir.create(pathName)}
    idatFiles <- dir(path=getwd(),pattern = ".idat", full.names = T)
    if(length(idatFiles)>0){
    gb$supM(lapply(idatFiles, function(x){file.rename(from = x, to = file.path(pathName, basename(x)))}))    
    }
}

CheckOriginal <- function(ogSheet) {
    if (!file.exists(ogSheet)) {
        rds <- gb$readInfo(inputSheet = gb$samsheet)
        result <- gb$search.redcap(rds, token)
        result <- result[!is.na(result$barcode_and_row_column),]
        samplesheet_ID = as.data.frame(stringr::str_split_fixed(result[, "barcode_and_row_column"], "_", 2))
        gb$makeSampleSheet(result, samplesheet_ID, bn = NULL, outputFi = ogSheet)
    }
}
