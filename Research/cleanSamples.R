gb <- globalenv(); assign("gb", gb)

dfTargets <- function(targets){
  if(class(targets)!="data.frame"){targets<- as.data.frame(targets)}
  return(targets)
  }

writeTargets<- function(targets, csvFi="samplesheet.csv"){
  targets <- dfTargets(targets)
  write.csv(targets, file=csvFi, quote=F, row.names=F)
}


CheckSamNames <-  function(samNames, targets){
  if(anyDuplicated(samNames) > 0){
      warning("Sample IDs contain duplicates, assigning unique names\n")
      message(paste0(capture.output(samNames[which(duplicated(samNames))]), collapse="\n"))

      samNames <- targets$Sample_ID <- make.unique(samNames, sep="_")
      targets$Sample_Name <- targets$Sample_ID <- samNames
    }
  return(targets)
}


SetKeyColumns <- function(targets, col_samTypes, col_samNames, col_other, col_shapes, col_Grouping=NULL) {
    targets <- dfTargets(targets)
    if(is.null(col_samTypes)){
        targets$Type <- "All Types"
    }else{
        targets$Type <- targets[, col_samTypes] # Creates any new "Type" column
    }
    if(is.null(col_samNames)){
        targets$Sample_Name <- targets$Sample_ID <-  targets[, 1]
    }else{
        targets$Sample_Name <- targets$Sample_ID <- targets[, col_samNames]
    }
    if(is.null(col_other)){
        targets$Other_Group <- targets$Sample_Group <- "All Groups"
    }else{
        targets$Other_Group <- targets$Sample_Group <- targets[, col_other]
    }
    if(is.null(col_shapes)){
        targets$Sym_Shape <- "All Samples"
    }else{
        targets$Sym_Shape <- targets[,col_shapes]
    }
    targets <- CheckSamNames(targets$Sample_ID, targets)
    if (is.null(col_Grouping)) {
      targets$Sample_Group <- "All Samples"
    }else{
      targets$Sample_Group <- col_Grouping  
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

FillMissingData <- function(targets, col_samNames="Sample_Name", originalFi="samplesheet_og.csv"){
  newTarg <- read.csv(originalFi, strip.white=T, row.names=NULL)
  targets <- merge(newTarg,targets, by=col_samNames, all=F, suffixes = c("",".xyzq"))
  dupeDrop <- grepl(".xyzq", colnames(targets))==F
  targets <- targets[, dupeDrop]
  extraCol <- c("DNA_Number", "RunID", "MP_num", "Date")
  toDrop <- !colnames(targets) %in% extraCol
  targets <- targets[, toDrop]
  write.csv(targets, file="samplesheet.csv", quote=F, row.names=F)
  targets <- read.csv("samplesheet.csv", strip.white=T, row.names=NULL)
  if(class(targets)!="data.frame"){targets <- as.data.frame(targets)}
  return(targets)
}


FixNullNaVars <- function(targets, varColumns) {
  if(all(varColumns %in% dimnames(targets)[[2]])==F){
    message("The column names you provided are not found in targets columns: ", paste(varColumns, collapse = ", "))
    stopifnot(all(varColumns %in% dimnames(targets)[[2]]))
    }
  for (variable in varColumns) {
    if (any(is.null(targets[, variable]))) {
      targets[is.null(targets[, variable]), variable] <- "NONE"
    }
    if (any(is.na(targets[, variable]))) {
      targets[is.na(targets[, variable]), variable] <- "NONE"
    }
    if (any(targets[, variable] == "")) {
        toFix <- targets[, variable] == ""
        targets[toFix, variable] <- "NONE"
    }
  }
  return(targets)
}


FixBaseName <- function(targets, runDir = NULL, col_sentrix) {
    if (is.null(runDir)) {runDir <- file.path(getwd(), "idats")}
    if (class(targets) != "data.frame") {targets <- as.data.frame(targets)}
    stopifnot(col_sentrix %in% colnames(targets))
    senCol <- min(which(grepl(col_sentrix, colnames(targets)) == T))
    targets$Basename <- file.path(runDir, targets[, senCol])
    return(targets)
}


ModifyTargetColumns <- function(targets, gb) {
    if (gb$needFi == T) {
        gb$GetCsvSheet(gb$needFi, gb$samsheet, gb$token, idatPath = gb$idatPath)
        targets <- gb$SetKeyColumns(
            targets,
            gb$col_samTypes,
            gb$col_samNames,
            gb$col_other,
            gb$col_shapes,
            gb$col_Grouping
        )
        targets <- FillMissingData(targets)
    } else{
        targets <- gb$SetKeyColumns(
            targets,
            gb$col_samTypes,
            gb$col_samNames,
            gb$col_other,
            gb$col_shapes,
            gb$col_Grouping
        )
    }
    
    if(is.null(gb$selectedVars)){
        gb$selectedVars <- "Type"
    }
    targets <- gb$FixNullNaVars(targets, selectedVars)
    targets <- gb$colorTargets(targets, varColumns = gb$selectedVars)
    if (is.null(gb$col_sentrix)) {
        gb$col_sentrix <- "SentrixID_Pos"
    }
    targets <-
        gb$FixBaseName(targets, runDir = gb$idatPath, gb$col_sentrix)
    return(targets)
}


FilterArrayKind <- function(targets, array_column, arrayToDrop = "450k"){
    toDrop <- targets[, array_column] != "450k"
    targets <- targets[toDrop,]
    row.names(targets) <- 1:nrow(targets)
    write.csv(targets, "samplesheet.csv", row.names = F, quote = F)
    return(targets)
}


DropMissingIdats <- function(targets, gb){
    idatReal <- file.exists(file.path(gb$idatPath, paste0(targets$SentrixID_Pos, "_Grn.idat")))
    targets <- targets[idatReal, ]
    rownames(targets) <- 1:nrow(targets)
    idatReal <- file.exists(file.path(gb$idatPath, paste0(targets$SentrixID_Pos, "_Red.idat")))
    targets <- targets[idatReal, ]
    rownames(targets) <- 1:nrow(targets)
    if(length(unique(targets$ArrayType)) > 1){
      message("Multiple Array types on this sample list!  EPIC and 450K")
      targets <- FilterArrayKind(targets, gb$col_arrayType)
    }
    if(!is.null(gb$col_sentrix)){
      targets[, gb$col_sentrix] <- targets$SentrixID_Pos
    }
    write.csv(targets, gb$samsheet, quote=F, row.names=F)
    targets <- as.data.frame(read.csv(gb$samsheet))
    return(targets)
}


GenerateTargets <- function(gb){
  targets <- gb$sanitizeSheet(gb$inputFi, gb$samsheet, gb) # Drop commas in xlsx sheet
  targets <- gb$ModifyTargetColumns(targets, gb)
  targets <- gb$GetArrayTypes(targets, arrayColumn = gb$col_arrayType, idatPath = gb$idatPath)
  targets <- gb$DropMissingIdats(targets, gb)
  targets <- gb$colorTargets(targets, gb$selectedVars)
  write.csv(targets, gb$samsheet, quote=F, row.names=F)
  targets <- as.data.frame(read.csv(gb$samsheet))
  return(targets)
}

