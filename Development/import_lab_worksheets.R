#!/usr/bin/env Rscript
## ---------------------------
## Script name: import_lab_worksheets.R
## Purpose: Read and Import all Excel Worksheets into REDCap
## Author: Jonathan Serrano
## Created: April 20, 2021
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------

apiToken = "8XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" # REDCap API token
ws.type = NULL # Default is NULL "DNA" "RNA" "PACT" "FUSION" will upload if NULL
doAll = T # Set to TRUE to import all worksheets from this month
allMths <- T # Set to TRUE to import all worksheets from the input year
td = Sys.Date() #Sys.Date() # today's date is default can be a string date "2022-01-26"
#td = "2022-01-26"

formals(library)$quietly <- T
formals(library)$warn.conflicts <- F
formals(require)$warn.conflicts <- F
formals(install.packages)$dependencies <- T
formals(install.packages)$verbose <- T
formals(install.packages)$ask <- F
options("install.packages.compile.from.source" = "never")
options("install.packages.compile.from.source" = "never")

rqPk <- function(pk) {
    if (!requireNamespace(pk, quietly = TRUE)) {
        install.packages(pk, warn.conflicts = F, ask = F, dependencies = T, Ncpus = 5)
    }
    library(pk, warn.conflicts = F, mask.ok = T, quietly = T, character.only = T)
}

check.packages <- function(pkgs = NULL) {
  if (is.null(pkgs)) {
    pkgs <- c(
      "readxl", "stringi", "readr", "parallel", "RCurl", "foreach",
      "redcapAPI", "doParallel", "tidyr", "chron", "crayon", "compiler",
      "purrr", "tidyverse", "cronR", "plyr", "doSNOW", "doParallel",
      "magrittr", "dplyr", "lubridate", "pbmcapply"
    )
  }
  rqPk("purrr")
  purrr::walk(pkgs, rqPk)
}

check.packages()

Sys.setenv(R_ENABLE_JIT = TRUE)
compiler::enableJIT(3)
rqPk("devtools")
rqPk("pacman")

if (!requireNamespace("Rmpi", quietly = TRUE)) {
  install.packages("Rmpi", type = "source")
}

doParallel::registerDoParallel(parallel::detectCores() - 2)


mountMsg <- function(msg) {
  if (msg == 1) {
    cat(crayon:::bgBlue("\nChecking files... "))
  }
  if (msg == 2) {
    cat(crayon::red("\nDirectory not found, stopping execution\n"))
    cat(
      crayon::white$bgRed$bold(
        "Mount the network drive or edit your input path to worksheets folder:"
      )
    )
    cat("\n'/molecular/MOLECULAR LAB ONLY/'\n")
  }
}

Check.File.Paths <- function(ws.type, ws) {
  cat(paste0("\n", ws.type, " INPUT PATH:"))
  message(crayon::red$bgBlack$bold(ws), "\n")
  cat(crayon::yellow$bgBlack$bold(paste("Does your", ws.type, "input path exist?")))
  folder.found <- dir.exists(ws)
  cat(" ", folder.found)
  if (folder.found) {
    mountMsg(1)
  } else {
    mountMsg(2)
  }
  stopifnot(folder.found)
}

removePattern <- function(fileList, file.patterns = NULL) {
  if (is.null(file.patterns)) {
    file.patterns = c(
      "~$",
      ".xml",
      "TRAIN",
      "Qubit",
      "test",
      "List",
      "SAM",
      "orig",
      "Tumor",
      "INUSE",
      "TEMPLATE",
      "combined",
      "TBDELETED",
      "dsDNA",
      "copy",
      "review",
      "Backup",
      "exported",
      "Exported",
      "-Sample-List",
      "SampleList",
      "NOTUSED",
      "old",
      "orginal",
      "PM",
      "pm",
      "comp",
      "Tranining",
      "ROUND",
      "TISSUE",
      "BLOOD",
      "VALIDATION",
      "BS",
      "Tranining"
    )
 }

valid.files <- as.character(parallel::mclapply(fileList, function(i) {
      keepFile <- parallel::mclapply(file.patterns, function(j) {
        opn = stringi::stri_opts_fixed(case_insensitive = T, overlap = F)
        ! stringi::stri_detect_coll(
          i,
          j,
          strength = 1L,
          case_level = F,
          opts_fixed = opn
        )
      })
      if (suppressWarnings(all(keepFile))) {
        return(i)
      }
    }))

  allFiles <- unique(valid.files[valid.files != "NULL"])
  return(allFiles)
}

set.headers <- function(ws.type) {
  naHeader <-
    switch(
      ws.type,
      "RNA" = c(
        "sticker",
        "record_id",
        "rna_number",
        "source",
        "tech",
        "date",
        "ffpe_block_id_2",
        "tumor_percent",
        "rna_comments",
        "date2"
      ),
      "PACT" = c(
        "Sample_ID",
        "Sample_Name",
        "Paired_Normal",
        "I7_Index_ID",
        "index",
        "record_id",
        "b_number",
        "epic_id_pact",
        "test_number_pact",
        "run_number_pact_id",
        "Sequencer_ID",
        "Chip_ID",
        "pact_run_number",
        "tumor_percentage",
        "notes_pact",
        "Description",
        "GenomeFolder"
      ),
      "FUSION" = c("archer_id", "record_id", "rna_number", "archer_lane"),
      "METH" = c(
        "rd_number",
        "record_id",
        "b_number",
        "conc",
        "barcode_and_row_column"
      ),
      c(
        "sticker",
        "record_id",
        "b_number",
        "source",
        "tech",
        "date",
        "ffpe_block_id",
        "tumor_percent",
        "comments_dna",
        "date2"
      )
    )
  return(naHeader)
}

GetSheetPath <- function(ws.type, td) {
    stopifnot(!is.null(ws.type))
    message("td is: ", td, "\n", "Year: ", paste(format(td, "%Y")))
    drive = file.path("", "Volumes", "molecular", "MOLECULAR LAB ONLY")
    folder <- switch(
        ws.type,
        "PACT" = file.path("NYU PACT Patient Data", "Workbook"),
        "FUSION" = file.path("NYU FUSION SEQer", paste(format(td, "%Y"), "reports")),
        "METH" = file.path(drive, "NYU-METHYLATION/WORKSHEETS", paste(format(td, "%Y"))),
        "MYLO" = file.path(
            drive,
            "NYU-MyeloSeqer",
            "Worksheets-Wetlab",
            paste(format(td, "%Y"))
        ),
        file.path(
            "Restricted-Molecular LAB staff ONLY",
            "PROCESSING&EXTRACTION",
            ws.type
        )
    )

    finalDir <- switch(
        ws.type,
        "FUSION" = file.path(drive, folder),
        "PACT" = file.path(drive, folder, format(td, "%Y")),
        "METH" = folder,
        "MYLO" = folder,
        file.path(drive, folder, format(td, "%Y"))
    )

    if (!dir.exists(finalDir)) {
        if (as.numeric(format(td, "%Y")) == 2022) {
            finalDir <- file.path(finalDir, "2022")
        }
    }

    if (!dir.exists(finalDir)) {
        message("Directory not found: ", finalDir)
        stopifnot(dir.exists(finalDir))
    }

    if (ws.type == "METH") {
        return(finalDir)
    }

    if (ws.type == "PACT") {
        if (as.numeric(format(td, "%Y")) < 2022) {
            dayMonth <- paste(format(td, "%m"), format(td, "%B"))
            finalDir <- file.path(finalDir , dayMonth)
        }
        return(finalDir)
    }

    monthDir = paste0("(?i)", format(td, "%B"))
    allDir <- dir(finalDir)
    subDir <- allDir[stringr::str_detect(allDir, regex(monthDir))]
    if(length(subDir) == 0) {
      return(finalDir)
    }
    return(file.path(finalDir, subDir))

}

CleanFusion <- function(fileList, fiNum, mgfs) {
  keywds <- c(".r_and_d_results.txt" = "", "/var/www/analysis/" = "")
  qsFix <- c("/" = " ", "_S" = " S", "-Q2" = " Q2")
  dataRead <- read.table(
    fileList[fiNum],
    sep = "\t",
    header = F,
    row.names = NULL,
    stringsAsFactors = F
  )
  finam <- as.list(dataRead[1, ][-1])
  finam <- stringr::str_replace_all(finam, keywds)
  finam <- stringr::str_split_fixed(stringr::str_replace_all(finam, qsFix), " ", 4)
  colnames(finam) <- set.headers("FUSION")
  finam <- as.data.frame(finam)
  finam$fs_batch <- mgfs[fiNum]
  finam$file_name <- fileList[fiNum]
  return(finam)
}

ReadFusionWs <- function(fileList, worksheet) {
    cntrls <- c("NTC-H2O", "SERASEQ-125ng", "NC-IVS0035")
    cntrls <- paste(cntrls, collapse = "|")
    spltSheets <- stringr::str_replace_all(fileList, worksheet, "")
    mgfs <- stringr::str_split_fixed(spltSheets, "/", 3)[, 2]
    dat.io <- foreach::foreach(
      fiNum = 1:length(fileList),
      .combine = "rbind",
      .packages = c("readxl", "foreach", "doParallel")
    ) %dopar% {
      return(CleanFusion(fileList, fiNum, mgfs))
    }
    drpSam <- !grepl(cntrls, dat.io[, 2], ignore.case = T)
    dat.io <- dat.io[drpSam,]
    recordIDs <- paste0(dat.io[, "record_id"])
    rnaNumber <- paste0(dat.io[, "rna_number"])
    recordIDs <- paste(recordIDs, rnaNumber, sep = "_")
    dat.io[, "record_id"] <- make.unique(recordIDs, sep = "_")
    return(dat.io)
}

msgTime <- function(modTime, accessTime) {
  recently <- modTime <= accessTime + 1
  if (any(recently)) {
    mm1 <-
      crayon::green(paste("File last imported/accessed:", accessTime, "hours ago\n\n"))
    mm2 <-
      crayon::blue(paste("File last modified:", modTime, "hours ago\n", mm1))
    cat(mm2)
    return(recently)
  } else {
    return(F)
  }
}

minusTi <- function(rf, fiNum) {
  return(chron::hours(Sys.time()) -
           chron::hours(switch(
             fiNum, "1" = file.mtime(rf), "2" = file.info(rf)$atime
           )))
}

get.mod.time <- function(fileList) {
  if (is.null(fileList)) {
    return(NULL)
  }
  thisMon <- chron::days(Sys.time()) == chron::days(file.mtime(fileList))
  #difftime(Sys.Date(),file.mtime(fileList))<0
  if (any(thisMon)) {
    recentFi <- fileList[thisMon] #files modified today
    if (length(recentFi) > 0) {
      todays <- file.info(recentFi)$atime - file.mtime(recentFi) < 2
      #difftime(file.info(recentFi)$atime,file.mtime(recentFi))<24
      if (any(todays)) {
        recentFi <- recentFi[todays]
      } else {
        return(NULL)
      }
      fiTime <- msgTime(minusTi(recentFi, 1), minusTi(recentFi, 2))
      if (any(fiTime)) {
        return(recentFi[fiTime])
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}


grab.today <- function(fileList) {
  recentMod <- (Sys.time() - file.mtime(fileList))
  units(recentMod) <- "days"
  if (any(recentMod <= 4)) {
    fileList <- fileList[recentMod <= 4]
    lastAccessed <- Sys.time() - file.info(fileList)$atime
    units(lastAccessed) <- "hours"
    if (any(lastAccessed >= 1)) {
      fileList <- fileList[lastAccessed >= 1]
      return(fileList)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}


GetFusionPath <- function(worksheet, td, doAll=F) {
  run.fld <- list.dirs(path = worksheet, full.names = T, recursive = F)
  toKeep <- stringr::str_detect(run.fld, "MGFS")
  run.fld <- run.fld[toKeep]

  if(doAll==F){
    run.fld <- run.fld[which.max(file.mtime(run.fld))]
  }

  frl <- list.dirs(run.fld, full.names = T, recursive = F)
  run.fld <- file.path(frl, "summaries")
  run.fis <- dir(run.fld, full.names = T)
  fileList <- run.fis[grepl("Summary.r_and_d_results.txt", run.fis)]

  if(doAll == F){
    fileList <- grab.today(fileList)
  }

  if (length(fileList) == 0) {
    return(NULL)
  }
  print(fileList)
  return(fileList)
}


GetPactPath <- function(worksheet, fend) {
    allSheetDirs <- list.dirs(worksheet)
    if(length(allSheetDirs)==0){
        fileList <- dir(worksheet, pattern = fend, full.names = T, recursive = T)
    }else{
        fileList <- dir(allSheetDirs, pattern = fend, full.names = T, recursive = T)
    }
    file.patterns = c("VAL", "PM", "dist", "Test", "-TumorLibraries", "-TIPSTEST", "initial")
    fileList <- removePattern(fileList, file.patterns) #[-1]
    fileList <- removePattern(fileList)
    return(fileList)
}


# Select only files created today (read only new samples)
GrabRecentFiles <- function(ws.type, td = NULL, doAll = F) {
    if (class(td) != "Date") {
        td <- as.Date(td)
        doAll <- T
    }

    worksheet <- GetSheetPath(ws.type, td)
    fend = ifelse(grepl("PACT", ws.type) == F, "*.xlsm", "*.xlsm")

    Check.File.Paths(ws.type, worksheet)

    if(ws.type == "FUSION"){
      fileList <- GetFusionPath(worksheet, td, doAll)
      return(fileList)
    }

    allSheetDirs <- list.dirs(worksheet)
    fileList <- dir(allSheetDirs, pattern = fend, full.names = T, recursive = T)
    fileList <- removePattern (fileList)

    if (ws.type == "PACT") {
        fileList <- GetPactPath(worksheet, fend)
    }

    if (ws.type == "MYLO") {
        file.patterns = c("dnr", "-dis")
        fileList <- removePattern(fileList, file.patterns) #[-1]
    }

    if (doAll == T) {
        return(fileList)
    }

    recentMod <- grab.today(fileList)

    if (!is.null(recentMod)) {
        return(recentMod)
    }

    recentFi <- get.mod.time(recentMod)
    return(recentFi)
}

# Filter non-OLD, broken, or duplicate files
filter_files <- function(fileList){
  nameFlags <- c('Old', 'OLD', 'old', 'NOTUSED', 'initial', 'Copy', 'copy',
                 'philip', 'Philip', '\\$', '\\~')
  nameFlags <- paste(nameFlags, collapse = "|")
  fileList <- unique(fileList[str_detect(fileList, nameFlags, negate = T)])
  toKeep <- endsWith(fileList, ".xlsm")
  fileList <- fileList[toKeep]
  return(fileList)
}


# Get sheet name
get_sheet_name <- function(fileList, type){
  switch(
    type,
    "METH" = rep(3, length(fileList)),
    "PACT" = rep("PhilipsExport", length(fileList)),
    "MYLO" = rep("DNA", length(fileList)),
    rep("INFO_ENTRY", length(fileList))
  )
}


# Process METH files
process_meth_files <- function(fileList, shNam, xl.range, naHeader) {
  map_dfr(
    1:length(fileList),
      ~read_excel(fileList[.], shNam[.], xl.range, naHeader, "text") %>%
      filter(rd_number != 0 & rd_number != "MP-20-##") %>%
      select(-4)
  )
}


process_mylo_fi <- function(fiPath){
    col2keep <- c(1, 2, 3, 4, 10)
    dataRead <- readxl::read_excel(path = fiPath, sheet = "DNA", skip = 5, col_types = "text")[,col2keep]
    dataRead2 <- readxl::read_excel(path = fiPath, sheet = "RNA", skip = 5, col_types = "text")[,col2keep]
    dataRead$file_name <- fiPath
    dataRead2$file_name <- fiPath
    allDat <- dplyr::bind_rows(dataRead, dataRead2)
    colnames(allDat) <- c("mylo_barcode", "accession_number", "b_number", "tumor_percent_mylo", "mylo_chip", "file_name", "rna_number")
    allDat <- allDat[!is.na(allDat[,2]),]
    return(as.data.frame(allDat))
}


ParseReadSheet <- function(dataRead, currFi) {
  if(any(stringr::str_detect(colnames(dataRead), "Sticker|sticker"))) {
      dataRead <- dataRead[, 2:ncol(dataRead)]
      dataRead$file_name <- currFi
      return(as.data.frame(dataRead))
  }

  if("Patient Name" %in% colnames(dataRead)){
      runIDs <- as.data.frame(dataRead[,"run_number_pact_id"])
      dataRead <- dataRead[, 1:25]
      dataRead$file_name <- currFi
      dataRead$run_number_pact_id <- paste(runIDs[1,1])
      return(as.data.frame(dataRead))
  }

  if("I7_Index_ID" %in% colnames(dataRead)){
      dataRead <- dataRead[, 1:19]
      dataRead$file_name <- currFi
      return(as.data.frame(dataRead))
  }
}


# Process excel file
process_excel_file <- function(currFi, shNam, xl.range){
    fileRead <- R.utils::withTimeout(read_excel(currFi, sheet = shNam[1], col_types = "text", skip = 0), timeout=10, onTimeout="warning")
    if(shNam[1] == "PhilipsExport"){
        runNumb <- R.utils::withTimeout(
            read_excel(currFi, sheet = 1, col_types = "text", skip = 5), timeout=10,
            onTimeout="warning") %>% as.data.frame
        runNumRow <- which(stringr::str_detect(runNumb[,6],"Run ID"))[1]
        runId <- paste(runNumb[runNumRow, 9])
        fileRead <- as.data.frame(fileRead)
        fileRead$run_number_pact_id <- runId
    }
    return(fileRead)
}


# Message on error
process_error <- function(currFi, shNam, xl.range){
    cat(crayon::bgYellow(paste0("File was missing the sheet named: ", shNam)), " reading Sheet1 instead")

    dataRead <- tryCatch(

        expr = R.utils::withTimeout(
            readxl::read_excel(
                path = currFi, sheet = 1, range = xl.range, col_types = "text"),
            timeout=10, onTimeout="error"),

        error = function(e){
            return(R.utils::withTimeout(
            readxl::read_excel(path = currFi, sheet = 1, skip = 5, col_types = "text"), timeout=10, onTimeout="silent"))
        }
    )

    if(nrow(dataRead)==0){
        return(R.utils::withTimeout(
            readxl::read_excel(path = currFi, sheet = 1, skip = 5, col_types = "text"),
            timeout=10, onTimeout="error")
            )
    }

    return(dataRead)
}


ProcessMyloDoc <- function(fileList) {

  allReadData <- pbmcapply::pbmclapply(X = 1:length(fileList), function(ws){
    fiPath = fileList[ws]
    dataRead <- process_mylo_fi(fiPath)
    if(is.null(dataRead)){
      return(NULL)
    }
    return(dataRead)
  })

  files.read <- unique(as.data.frame(dplyr::bind_rows(allReadData)))

  return(files.read)
}



# Process docs
process_docs <- function(fileList, shNam, xl.range){

    allReadData <- pbmcapply::pbmclapply(X = 1:length(fileList), function(ws){
            currFi = fileList[ws]
            dataRead <- tryCatch(
                expr = process_excel_file(currFi, shNam[1], xl.range),
                error = function(e) {return(process_error(currFi, shNam[ws], xl.range))}
                )
            if(is.null(dataRead)){
              return(NULL)
            }
            dataRead <- dataRead[!is.na(dataRead[,2]),]

        if (nrow(dataRead) > 0 & ncol(dataRead) >= 15) {
            dataRead <- ParseReadSheet(dataRead, currFi)
            dataRead <- as.data.frame(dataRead)
            if(class(dataRead) == "data.frame" & nrow(dataRead) > 4){
                return(dataRead)
            }
        }
            return(NULL)
    })
    toDrop <- unlist(lapply(allReadData, is.null))
    allReadData <- allReadData[!toDrop]
    files.read <- unique(as.data.frame(dplyr::bind_rows(allReadData)))
  return(files.read)
}


ParseWorkbook <- function(fileList, xl.range, ws.type) {
  naHeader <- set.headers(toupper(ws.type))
  if(length(fileList)==0){
    return(NULL)
  }
  message(paste0(capture.output(fileList), collapse="\n"))
  shNam <- get_sheet_name(fileList, type = ws.type)

  if (ws.type == "MYLO"){
    files.read <- ProcessMyloDoc(fileList)
    return(files.read)
  }

  if (ws.type == "METH") {
    files.read <- process_meth_files(fileList, shNam, xl.range, naHeader)
    return(files.read)
  }

  files.read <- process_docs(fileList, shNam, xl.range)

  if (ws.type == "PACT") {
      #shNam <- stringr::str_remove_all(basename(fileList), ".xlsm")
      #files.read2 <- process_docs(fileList, shNam, xl.range)
  }

  return(files.read)
}


# Function to replace specific strings in data
replace_in_data <- function(data, replacements){
  data[is.na(data)] <- ""
  toupper(stringr::str_replace_all(trimws(data), replacements))
}

# Function to correct typos in the data
autocorrect_typos <- function(data){
  autocorrect <- c(" " = "-", "NGS2" = "NGS-2", "RD2" = "RD-2", "RD1" = "RD-1", "NGS1" = "NGS-1")
  replace_in_data(data, autocorrect)
}

# Function to remove illegal characters
remove_illegals <- function(data){
  illegals <- c("/" = " ", "\\+" = "", "&" = " ", "#" = "num", "\\\\"  = " ", "<" = "", ">" = "", "\\*"="")
  make.unique(as.character(replace_in_data(paste0(data), illegals)), sep = "_")
}

GetColumnPatterns <- function(ws.type){
  dnaRnaPatt <- c("Accession", "DNA", "RNA", "Container", "Source", "Product", "by", "Block",
               "Tumor", "Comment", "Test", "file_name")
  methPat <- ""
  pactPat <- c('Patient', 'Gender', 'Assignee', 'MRN', 'Specimen', 'DNA', 'RNA',
               'Test\\sNumber', 'Diagnosis', 'Epic', 'Percentage', 'Type', 'Site', "file_name", "run_number_pact_id")

  pattern <- switch(
    ws.type,
    "METH" = methPat,
    "PACT" = pactPat,
    dnaRnaPatt
  )

  return(paste(pattern, collapse = "|"))
}

GetNewColNames <- function(ws.type){
  dnaColumns <- c(
    "accession_number", "b_number", "container", "sample_source", "products",
    "processing_tech", "block_number", "tumor_percent", "comments", "tests", "file_name"
  )
  rnaColumns <- c(
    "accession_number", "rna_number", "container", "sample_source", "products",
    "processing_tech", "block_number", "tumor_percent", "comments", "file_name", "tests"
  )

  methColumns <- ""
  defaultCols <- ""
  pactColumns <- c("philips_patient_name", "philips_gender", "philips_assignee",
                   "philips_mrn", "philips_tumor_specimen_id", "philips_normal_specimen_id",
                   "philips_test_number", "philips_diagnosis", "philips_epic_order",
                   "philips_diagnosis_interp", "philips_tumor_dna_num",
                   "philips_normal_dna_num", "philips_tumor_percent", "philips_tumor_type",
                   "philips_normal_type", "philips_tumor_site", "philips_normal_site", "file_name", "run_number_pact_id")

  newColNames <- switch(
    ws.type,
    "DNA" = dnaColumns,
    "RNA" = rnaColumns,
    "METH" = methColumns,
    "PACT" = pactColumns,
    defaultCols
  )

  return(newColNames)
}


CleanOutput <- function(files.read, ws.type) {
  # Corrections and removals
  spaceFix <- c("\\," = "", "\\." =";", "  " = " ")
  cntrls <- c("NTC" = "", "SC" = "", "NC" = "")
  allCols <- colnames(files.read)

  if(ws.type!="MYLO" & ws.type!="METH"){
      pattern <- GetColumnPatterns(ws.type)
      dat.io <- files.read[, allCols[stringr::str_detect(allCols, pattern)]]
      newColnames <- GetNewColNames(ws.type)
      if(ncol(dat.io) != length(newColnames)){
          newColnames <- newColnames[-10]
      }
      if(ncol(dat.io) != length(newColnames)){
          newColnames <- c("accession_number", "b_number", "sample_source", "processing_tech",
                           "block_number", "tumor_percent", "comments", "file_name")
      }
      colnames(dat.io) <- newColnames
  }else{
    dat.io <- files.read
  }

  for(n in 1:ncol(dat.io)){
    dat.io[, n] <- replace_in_data(dat.io[, n], spaceFix)
  }

  dat.io[, 1] <- autocorrect_typos(dat.io[, 1])
  dat.io <- dat.io[dat.io[, 1] != 0, ]
  dat.io[, 1] <- remove_illegals(dat.io[, 1])

  return(as.data.frame(dat.io))
}


Excel2csv <- function(fileList, xl.range, ws.type) {
    message("\nReading Worksheets...")
    if (length(fileList) == 0) {
      warning("File not readable, no files to import")
      return(NULL)
    }
    fileList <- filter_files(fileList)
    files.read <- ParseWorkbook(fileList, xl.range, ws.type)
    dat.io <- CleanOutput(files.read, ws.type)
    dat.io <- as.data.frame(dat.io)
    return(dat.io)
}


pullRedcap <- function(dat.io, rcon, col2swap) {
  bqCurr <- redcapAPI::exportRecordsTyped(
      rcon,
      fields = c("record_id", col2swap),
      survey = F,
      form_complete_auto = F,
      dag=F
    )
  rownames(bqCurr) <- bqCurr[, 1]
  return(bqCurr[dat.io[, 1], ])
}


CombineBQs <- function(dat.io, bqCurr, rcrd, col2swap) {
  theRecord <- bqCurr[rcrd, col2swap]

  if (is.na(theRecord)) {
    return(dat.io[rcrd, 2])
  }

  alreadyIn <- grepl(dat.io[rcrd, 2], bqCurr[rcrd, 2])

  if (alreadyIn & is.na(bqCurr[rcrd, 2])) {
    return(dat.io[rcrd, 2])
  }

  if (!is.na(dat.io[rcrd, 2]) & !alreadyIn) {
    return(paste(bqCurr[rcrd, 2], dat.io[rcrd, 2]))
  }

  return(bqCurr[rcrd, 2])
}


DropAllNAvals <- function(dat.io) {
  idx <- which(dat.io == "NA" | is.na(dat.io), arr.ind = TRUE)
  if (length(idx) > 0) {
    dat.io[idx] <- ""
  }
  return(dat.io)
}


CheckBQDupes <- function(dat.io, rcon, col2swap){
    bqCurr <- pullRedcap(dat.io, rcon, col2swap)
    if (!is.null(bqCurr) & nrow(bqCurr) > 0) {
        bqSwap <- foreach::foreach(rcd = 1:nrow(bqCurr), .combine = "rbind") %dopar% {
            CombineBQs(dat.io, bqCurr, rcd, col2swap)
        }
        bqSwap <- gsub("NA ", "", paste(bqSwap))
        dat.io[, col2swap] <- bqSwap
    }
    return(dat.io)
}


Check4Dupes <- function(rcon, dat.io, block = NULL) {
    if (!is.null(dat.io) | nrow(dat.io) == 0) {
        return(dat.io)
    }

    dat.io <- DropAllNAvals(dat.io)

    if(!any(c("b_number","dna_number", "rna_number") %in% colnames(dat.io))){
        return(dat.io)
    }

    if("b_number" %in% colnames(dat.io)){
        dat.io <- CheckBQDupes(dat.io, rcon, "b_number")
    }

    if("rna_number" %in% colnames(dat.io)){
        dat.io <- CheckBQDupes(dat.io, rcon, "rna_number")
    }

    return(as.data.frame(dat.io))
}

MsgImport <- function(dat.io = NULL, ws.type = NULL, rez = NULL) {
  if (!is.null(rez)){
    message("Imported Records:", paste(rez, sep = "\n"))
  }
    if(ws.type == "METH"){
        output = file.path(fs::path_home(),"Desktop", paste0(Sys.Date(),"_meth_file.csv"))
        readr::write_csv(dat.io, output, na = "", quote = "none")
    }else{
  if (!is.null(ws.type)) {
    io_name <- if (!is.null(dat.io) && length(dat.io) > 0) {
      wbNam <- paste(stringr::str_split(dat.io$file_name[1], "/")[[1]][6:7], collapse = "_")
      extraNm <- stringr::str_replace_all(wbNam, " ", "_")
      csvFi = paste(format(Sys.Date(), "%b%d"), ws.type, extraNm, "sheets.csv", sep = "_")
      output = file.path(fs::path_home(),"Desktop", csvFi)
      readr::write_csv(dat.io, output, na = "", quote = "none")
      paste0("Saved CSV: ", output)
      } else {
        paste0(ws.type, ": Import is up-to-date")
      }
    message(crayon::bgGreen(io_name))
  }}
}


pushRedcap <- function(raw_data_csv, rcon) {
  message("\nRecords importing...\n")
  redcapAPI::importRecords(rcon, raw_data_csv)
}


CleanAndSaveData <- function(dat.io, rcon){
    rownames(dat.io) <- 1:nrow(dat.io)
    dat.io <- Check4Dupes(rcon, dat.io)
    dat.io <- DropAllNAvals(dat.io)
    MsgImport(dat.io, ws.type = ws.type)
    return(dat.io)
}


MakeMyloRecords <- function(raw_data_csv) {
    raw_data_csv$record_id <- raw_data_csv$accession_number
    rnaVals <- raw_data_csv$rna_number!=""
    dnaVals <- raw_data_csv$b_number!=""
    raw_data_csv$record_id[rnaVals] <-
      paste(raw_data_csv$record_id[rnaVals], raw_data_csv$rna_number[rnaVals], sep="_")
    raw_data_csv$record_id[dnaVals] <-
      paste(raw_data_csv$record_id[dnaVals], raw_data_csv$b_number[dnaVals], sep="_")
    return(raw_data_csv)
}


MakePactRecords <- function(raw_data_csv) {
    recordID <- paste(raw_data_csv$philips_tumor_specimen_id,
                      raw_data_csv$philips_tumor_dna_num, sep = "_")
    raw_data_csv$record_id <- recordID
    runNames <- paste(stringr::str_split_fixed(basename(raw_data_csv$file_name), ";", 2)[,1])
    raw_data_csv$pact_run_number <- runNames
    return(raw_data_csv)
}


csv2redcap <- function(dat.io, apiToken, ws.type) {
  rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", apiToken)

  if (is.null(dat.io)) {return(NULL)}

  raw_data_csv <- CleanAndSaveData(dat.io, rcon)

  if (ws.type == "METH") {
    raw_data_csv$record_id <- paste(raw_data_csv$record_id, raw_data_csv$b_number, sep = "_")
  }

  if (ws.type == "MYLO") {
    raw_data_csv <- MakeMyloRecords(raw_data_csv)
  }

  if (ws.type == "PACT") {
    raw_data_csv <- MakePactRecords(raw_data_csv)
  }

  if (ws.type == "DNA") {
    raw_data_csv$record_id <-
      paste(raw_data_csv$accession_number, raw_data_csv$b_number, sep = "_")
  }

  if (ws.type == "RNA") {
    raw_data_csv$record_id <-
      paste(raw_data_csv$accession_number, raw_data_csv$rna_number, sep = "_")
  }

  dupes <- duplicated(raw_data_csv$record_id)

  if(any(dupes)){
    raw_data_csv <- raw_data_csv[!dupes,]
  }

  if (ws.type != "METH") {
    raw_data_csv$file_name <- stringr::str_replace_all(raw_data_csv$file_name, ";", ".")
  }
  pushRedcap(raw_data_csv, rcon)
}

FsGrab <- function(ws.type, apiToken, td, doAll) {
  fileList <- GrabRecentFiles(ws.type, td, doAll)
  dat.io <- ReadFusionWs(fileList, worksheet)
  csv2redcap(dat.io, apiToken, ws.type)
  return(NULL)
}


AutomateFunCalls <- function(ws.type = "DNA", apiToken = NULL, td = NULL, doAll = F) {

    message("Checking ", ws.type)

    xl.range <- switch(
        ws.type,
        "PACT" = "W7:AM41",
        "MYLO" = "B7:E42",
        "METH" = "A2:E97",
        "FUSION" = FsGrab(ws.type, apiToken, td, doAll),
        "A2:J86"
    )

    dat.io <- fileList <- NULL
    if(doAll==F){
        td <- as.Date(td)
    }

    if (!is.null(xl.range)) {
        fileList <- GrabRecentFiles(ws.type, td, doAll)
    }

    if (!is.null(fileList) & length(fileList) > 0) {
        dat.io <- Excel2csv(fileList, xl.range, ws.type)
    }

    gc(verbose = F)

    if (!is.null(dat.io) & !is.null(apiToken) & length(dat.io) > 0) {
        csv2redcap(as.data.frame(dat.io), apiToken, ws.type)
    }

    MsgImport(NULL, ws.type, NULL)
    }


SheetsIntoRedcap <- function(apiToken = NULL, ws.type = NULL, doAll = F, td) {
  if (is.null(apiToken)) {
    stop("apiToken is NULL")
  }
  if (is.null(ws.type)) {
    ws.type <- c("DNA", "RNA", "PACT", "FUSION", "MYLO", "METH")
  }
  lapply(ws.type, AutomateFunCalls, apiToken = apiToken, td = td, doAll = doAll)
  #AutomateFunCalls("PACT", apiToken, td, doAll)
}


GetAllMonths <- function(yr) {
  monthLi <- sprintf("%02d", 1:12)
  thisYear = as.character(format(as.Date(Sys.Date()), "%Y"))
  if(yr == thisYear){
      thisMonth = as.numeric(format(Sys.Date(), "%m"))
      monthLi <- monthLi[1:thisMonth]
  }
  return(monthLi)
}

DoWholeYear <- function(yr = "2021", ws.type = NULL, apiToken) {
    if(ws.type == "PACT" | ws.type == "FUSION" | ws.type == "METH" | ws.type == "MYLO"){
        SheetsIntoRedcap(apiToken, ws.type, doAll=T, td=paste(yr, "01", "01", sep = "-"))
    }else{
        monthLi <- GetAllMonths(yr)
        for (mnth in monthLi) {
            td <- paste(yr, mnth, "01", sep = "-")
            message("Processing Month: ", paste(lubridate::month(td, T, F)),"\n")
            SheetsIntoRedcap(apiToken, ws.type, doAll=T, td)
        }
    }
}

if (allMths == T) {
  if (is.null(ws.type)) {
    allTypes <- c("DNA", "RNA", "PACT", "FUSION", "MYLO", "METH")
  }
  yr <- as.character(format(as.Date(td), "%Y"))
  if (length(allTypes) > 1) {
    for (nSheet in allTypes) {
      ws.type <- nSheet
      assign("ws.type", nSheet, envir = .GlobalEnv)
      DoWholeYear(yr, ws.type, apiToken)
    }
  } else{
    DoWholeYear(yr, ws.type, apiToken)
  }
} else{
  SheetsIntoRedcap(apiToken, ws.type, doAll, td)
}

# # cronR Job start --------------------------------------------------
# if(!require("cronR")){install.packages("cronR")};library("cronR")
# scriptPath <- "/Volumes/CBioinformatics/jonathan/Rprojects/Methylation-scripts/import_lab_worksheets.R"
# cmd <- cronR::cron_rscript(scriptPath)
# cronR::cron_add(cmd, frequency = "*/5 * * * *", days_of_week=c(1:6))
# cronR::cron_njobs(); cronR::cron_ls()
