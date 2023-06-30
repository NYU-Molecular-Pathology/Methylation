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
    install.packages(
      pk,
      warn.conflicts = FALSE,
      ask = FALSE,
      dependencies = TRUE,
      Ncpus = 5
    )
  }
  library(pk, warn.conflicts = FALSE, mask.ok = TRUE, quietly = TRUE, character.only = TRUE)
}

check.packages <- function(pkgs = NULL) {
  if (is.null(pkgs)) {
    pkgs <- c(
      "readxl", "stringi", "readr", "parallel", "RCurl", "foreach",
      "redcapAPI", "doParallel", "tidyr", "chron", "crayon", "compiler",
      "purrr", "tidyverse", "cronR", "plyr", "doSNOW", "doParallel",
      "doMPI", "magrittr", "dplyr"
    )
  }

  # Ensure 'purrr' is installed before calling purrr::walk
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

removePattern <- function(file.list, file.patterns = NULL) {
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
      "old",
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
      "SampleList"
    )
  }

  valid.files <- as.character(parallel::mclapply(file.list, function(i) {
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
  return(valid.files[valid.files != "NULL"])
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

get.sheet.path <- function(ws.type, td) {
  stopifnot(!is.null(ws.type))
  message("td is: ", td, " ", paste(format(td, "%Y")))
  drive = file.path("", "Volumes", "molecular", "MOLECULAR LAB ONLY")
  folder <- switch(
    ws.type,
    "PACT" = file.path("NYU PACT Patient Data", "Workbook"),
    "FUSION" = file.path("NYU FUSION SEQer", paste(format(td, "%Y"), "reports")),
    "METH" = file.path(drive, "NYU-METHYLATION/WORKSHEETS", paste(format(td, "%Y"))),
    file.path(
      "Restricted-Molecular LAB staff ONLY",
      "PROCESSING&EXTRACTION",
      ws.type
    )
  )
  finalDir <- switch(
    ws.type,
    "FUSION" = return(file.path(drive, folder)),
    "PACT" = return(file.path(drive, folder, format(td, "%Y"))),
    "METH" = return(folder),
    return(file.path(drive, folder, format(td, "%Y")))
  )
  if (ws.type == "PACT") {
    if (as.numeric(format(td, "%Y")) < 2022) {
      return(file.path(finalDir , format(td, "%B")))
    } else{
      return(file.path(finalDir , paste(format(td, "%m"), format(td, "%B"))))
    }
  } else{
    return(file.path(finalDir, format(td, "%B")))
  }
}

CleanFusion <- function(file.list, fiNum, mgfs) {
  keywds <- c(".r_and_d_results.txt" = "", "/var/www/analysis/" = "")
  qsFix <- c("/" = " ", "_S" = " S", "-Q2" = " Q2")
  dataRead <- read.table(
    file.list[fiNum],
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
  finam$file_name <- file.list[fiNum]
  return(finam)
}

ReadFusionWs <- function(file.list, worksheet) {
  cntrls <- c("NTC-H2O", "SERASEQ-125ng", "NC-IVS0035")
  mgfs <-
    stringr::str_split_fixed((stringr::str_replace_all(file.list, worksheet, "")), "/", 3)[, 2]
  dat.io <-
    foreach::foreach(
      fiNum = 1:length(file.list),
      .combine = "rbind",
      .packages = c("readxl", "foreach", "doParallel")
    ) %dopar% {
      return(CleanFusion(file.list, fiNum, mgfs))
    }
  drpSam <- !grepl(paste(cntrls, collapse = "|"), dat.io[, 2], ignore.case = T)
  return(dat.io[drpSam, ])
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

get.mod.time <- function(file.list) {
  if (is.null(file.list)) {
    return(NULL)
  }
  thisMon <- chron::days(Sys.time()) == chron::days(file.mtime(file.list))
  #difftime(Sys.Date(),file.mtime(file.list))<0
  if (any(thisMon)) {
    recentFi <- file.list[thisMon] #files modified today
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


gtd <- function(file.list) {
  recentMod <- (Sys.time() - file.mtime(file.list))
  units(recentMod) <- "days"
  if (any(recentMod <= 4)) {
    file.list <- file.list[recentMod <= 4]
    lastAccessed <- Sys.time() - file.info(file.list)$atime
    units(lastAccessed) <- "hours"
    if (any(lastAccessed >= 1)) {
      file.list <- file.list[lastAccessed >= 1]
      return(file.list)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}


GetFusionPath <- function(worksheet, td, doAll=F) {
  run.fld <- list.dirs(path = worksheet,
                       full.names = T,
                       recursive = F)
  if(doAll==F){
    run.fld <- run.fld[which.max(file.mtime(run.fld))]
  }
  frl <- list.dirs(run.fld, full.names = T, recursive = F)
  run.fld <- file.path(frl, "summaries")
  run.fis <- dir(run.fld, full.names = T)
  file.list <- run.fis[grepl("Summary.r_and_d_results.txt", run.fis)]

  if(doAll == F){
    file.list <- grab.today(file.list)
  }

  if (length(file.list) == 0) {
    return(NULL)
  }
  print(file.list)
  dat.io <- ReadFusionWs(file.list, worksheet)
  dat.io[, "record_id"] <- make.unique(paste0(dat.io[, "record_id"]), sep = "_")
  return(dat.io)
}


# Select only files created today (read only new samples)
grab.recent <- function(ws.type, td = NULL, doAll = F) {
  if (class(td) != "Date") {
    td <- as.Date(td)
    doAll <- T
  }
  worksheet <- get.sheet.path(ws.type, td)
  fend = ifelse(grepl("PACT", ws.type) == F, "*.xlsm", "*.xlsm")
  Check.File.Paths(ws.type, worksheet)
  if (ws.type == "FUSION") {
    dat.io <- GetFusionPath(worksheet, td, doAll)
    return(dat.io)
  }

  if (ws.type == "PACT") {
    file.list <- list.dirs(ws)
    file.list <- dir(file.list, pattern = fend, full.names = T, recursive = T)
    file.patterns = c("VAL", "PM", "dist", "Test", "-TumorLibraries", "-TIPSTEST")
    file.list <- removePattern(file.list, file.patterns)[-1]
    #xlFiles <- file.path(file.list, paste0(basename(file.list), ".xlsm"))
    #file.list<- grab.today(xlFiles[file.exists(xlFiles)])
    return(file.list)
    }

  file.list <- dir(worksheet, full.names = T, recursive = T, include.dirs = T)
  file.list <- removePattern(file.list)
  if (doAll == T) {
    return(file.list)
  }
  recentMod <- grab.today(file.list)
  if (is.null(recentMod)) {
    return(recentMod)
  } else {
    recentFi <- gmt(recentMod)
    return(recentFi)
  }
}

# Filter non-OLD, broken, or duplicate files
filter_files <- function(file.list){
  nameFlags <- c('Old', 'OLD', 'old', 'NOTUSED', 'initial', 'Copy', 'copy',
                 'philip', 'Philip', '\\$', '\\~')
  nameFlags <- paste(nameFlags, collapse = "|")
  file.list <- unique(file.list[str_detect(file.list, nameFlags, negate = T)])
  toKeep <- endsWith(file.list, ".xlsm")
  file.list <- file.list[toKeep]
  return(file.list)
}

# Get sheet name
get_sheet_name <- function(file.list, type){
  switch(
    type,
    "METH" = rep(3, length(file.list)),
    "PACT" = rep("PhilipsExport", length(file.list)),
    rep("INFO_ENTRY", length(file.list))
  )
}

# Process excel file
process_excel_file <- function(file, shNam, xl.range){
  return(read_excel(file, sheet = shNam[1], col_types = "text", skip = 0))
}

# Message on error
process_error <- function(file, shNam, xl.range){
  cat(bgYellow(paste0("File was missing the sheet named: ", shNam)), " reading Sheet1 instead")
  return(read_excel(file, 1, xl.range, TRUE, "text"))
}

# Process METH files
process_meth_files <- function(file.list, shNam, xl.range, naHeader) {
  map_dfr(
    1:length(file.list),
      ~read_excel(file.list[.], shNam[.], xl.range, naHeader, "text") %>%
      filter(rd_number != 0 & rd_number != "MP-20-##") %>%
      select(-4)
  )
}

# Process docs
process_docs <- function(file.list, shNam, xl.range){
  alldataRead <- foreach(ws = 1:length(file.list)) %dopar% {

    dataRead <- tryCatch(
      expr = process_excel_file(file.list[ws], shNam[1], xl.range),
      error = process_error(file.list[ws], shNam[1], xl.range)
    )

    toKeep <- !is.na(dataRead[,2])
    dataRead <- dataRead[toKeep,]

    if (nrow(dataRead) > 0 & ncol(dataRead) >= 19) {

      if(any(stringr::str_detect(colnames(dataRead), "Sticker|sticker"))) {
        dataRead <- dataRead[, 2:ncol(dataRead)]
        dataRead$fileName <- file.list[ws]
        return(as.data.frame(dataRead))
      }

      if("Patient Name" %in% colnames(dataRead)){
        dataRead <- dataRead[,1:25]
        dataRead$fileName <- file.list[ws]
        return(as.data.frame(dataRead))
      }
    }
  }

  boundData <- as.data.frame(dplyr::bind_rows(alldataRead))
  return(boundData)
}


parse.ws <- function(file.list, xl.range, ws.type) {
  naHeader <- set.headers(toupper(ws.type))
  file.list <- filter_files(file.list)
  if(length(file.list)==0){
    return(NULL)
  }
  message(paste0(capture.output(file.list), collapse="\n"))
  shNam <- get_sheet_name(file.list, type = ws.type)
  if (ws.type == "METH") {
    files.read <- process_meth_files(file.list, shNam, xl.range, naHeader)
  }else{
    files.read <- process_docs(file.list, shNam, xl.range)
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
               "Tumor", "Comment", "Test", "fileName")
  methPat <- ""
  pactPat <- c('Patient', 'Gender', 'Assignee', 'MRN', 'Specimen', 'DNA', 'RNA',
               'Test\\sNumber', 'Diagnosis', 'Epic', 'Percentage', 'Type', 'Site', "fileName")

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
                   "philips_normal_type", "philips_tumor_site", "philips_normal_site", "file_name")

  newColNames <- switch(
    ws.type,
    "DNA" = dnaColumns,
    "RNA" = rnaColumns,
    "METH" = methPattern,
    "PACT" = pactColumns,
    defaltPatt
  )

  return(newColNames)
}


CleanOutput <- function(files.read, ws.type) {
  # Corrections and removals
  spaceFix <- c("\\," = "", "\\." =";", "  " = " ")
  cntrls <- c("NTC" = "", "SC" = "", "NC" = "")
  allCols <- colnames(files.read)

  pattern <- GetColumnPatterns(ws.type)
  dat.io <- files.read[, allCols[stringr::str_detect(allCols, pattern)]]
  colnames(dat.io) <- GetNewColNames(ws.type)

  for(n in 1:ncol(dat.io)){
    dat.io[, n] <- replace_in_data(dat.io[, n], spaceFix)
  }

  dat.io[, 1] <- autocorrect_typos(dat.io[, 1])
  dat.io <- dat.io[dat.io[, 1] != 0, ]
  dat.io[, 1] <- remove_illegals(dat.io[, 1])

  return(as.data.frame(dat.io))
}

xL2csv <- function(file.list, xl.range, ws.type) {
  message("\nReading Worksheets...")

  if (length(file.list) > 0) {
    files.read <- parse.ws(file.list, xl.range, ws.type)
    dat.io <- CleanOutput(files.read, ws.type)
    return(dat.io)
    } else {
      warning("File not readable, no files to import")
    }
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

  bqCurr[rcrd, 2]
}

DropAllNAvals <- function(dat.io) {
  idx <- which(dat.io == "NA" | is.na(dat.io), arr.ind = TRUE)
  if (length(idx) > 0) {
    dat.io[idx] <- ""
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

  col2swap <- ifelse(any(grepl("b_number", colnames(dat.io))), "b_number", "rna_number")
  bqCurr <- pullRedcap(dat.io, rcon, col2swap)

  if (!is.null(bqCurr) & nrow(bqCurr) > 0) {
    bqSwap <- foreach::foreach(rcd = 1:nrow(bqCurr), .combine = "rbind") %dopar% {
      CombineBQs(dat.io, bqCurr, rcd, col2swap)
    }
    bqSwap <- gsub("NA ", "", paste(bqSwap))
    dat.io[, col2swap] <- bqSwap
  }

  return(as.data.frame(dat.io))
}

MsgImport <- function(dat.io = NULL, ws.type = NULL, rez = NULL) {
  if (!is.null(rez)){
    message("Imported Records:", paste(rez, sep = "\n"))
  }
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
  }
}


pushRedcap <- function(raw_data_csv, rcon) {
  message("\nRecords importing...\n")
  redcapAPI::importRecords(rcon, raw_data_csv)
}

CleanAndSaveData <- function(dat.io){
  rownames(dat.io) <- 1:nrow(dat.io)
  dat.io <- Check4Dupes(rcon, dat.io)
  dat.io <- DropAllNAvals(dat.io)
  MsgImport(dat.io, ws.type = ws.type)
  return(dat.io)
}

csv2redcap <- function(dat.io, apiToken, ws.type) {
  rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", apiToken)

  if (is.null(dat.io)) {return(NULL)}

  raw_data_csv <- CleanAndSaveData(dat.io)

  if (ws.type == "METH") {
    raw_data_csv <- raw_data_csv[, 1:3]
  }

  if (ws.type == "PACT") {
    recordID <- paste(raw_data_csv$philips_tumor_specimen_id,
                      raw_data_csv$philips_tumor_dna_num, sep = "_")
    raw_data_csv$record_id <- recordID
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

  pushRedcap(raw_data_csv, rcon)
}

FsGrab <- function(ws.type, apiToken, td, doAll) {
  dat.io <- grab.recent(ws.type, td, doAll)
  csv2redcap(dat.io, apiToken, ws.type)
  return(NULL)
}


AutomateFunCalls <- function(ws.type = "DNA", apiToken = NULL,
                             td = NULL, doAll = F) {
    xl.range <- switch(
      ws.type,
      "PACT" = "W7:AM41",
      "MYLO" = "B7:E42",
      "METH" = "A2:E97",
      "FUSION" = FsGrab(ws.type, apiToken, td, doAll),
      "A2:J86"
    )
    if (!is.null(xl.range)) {
      file.list <- grab.recent(ws.type, td, doAll)
      if (!is.null(file.list) & length(file.list) > 0) {
        dat.io <- xL2csv(file.list, xl.range, ws.type)
        gc(verbose = F)
        if (!is.null(dat.io) & !is.null(apiToken) & length(dat.io) > 0) {
          csv2redcap(dat.io, apiToken, ws.type)
        }
      }
      MsgImport(NULL, ws.type, NULL)
    }
}


SheetsIntoRedcap <- function(apiToken = NULL, ws.type = NULL, doAll = F, td) {
  if (is.null(apiToken)) {
    stop("apiToken is NULL")
  }
  if (is.null(ws.type)) {
    ws.type <- c("DNA", "RNA", "PACT", "FUSION")
  }
  lapply(ws.type, AutomateFunCalls, apiToken = apiToken, td = td, doAll = doAll)
  #AutomateFunCalls("PACT", apiToken, td, doAll)
}

DoWholeYear <- function(yr = "2021", ws.type, apiToken) {
  if(ws.type == "PACT" | ws.type == "FUSION"){
    td <- paste(yr, "01", "01", sep = "-")
    return(
      SheetsIntoRedcap(apiToken, ws.type, doAll=T, td)
      )
  }

  monthLi <- sprintf("%02d", 1:12)
  thisYear = as.character(format(as.Date(Sys.Date()), "%Y"))
  # Filter up to current month if current year
  if(yr == thisYear){
    thisMonth = as.numeric(format(Sys.Date(), "%m"))
    monthLi <- monthLi[1:thisMonth]
  }
  for (mnth in monthLi) {
    td <- paste(yr, mnth, "01", sep = "-")
    SheetsIntoRedcap(apiToken, ws.type, doAll=T, td)
  }
}

gmt <- compiler::cmpfun(get.mod.time)
grab.today <- compiler::cmpfun(gtd)

if (allMths == T) {
  yr <- as.character(format(as.Date(td), "%Y"))
  DoWholeYear(yr, ws.type, apiToken)
} else{
  SheetsIntoRedcap(apiToken, ws.type, doAll, td)
}


# # cronR Job start --------------------------------------------------
# if(!require("cronR")){install.packages("cronR")};library("cronR")
# scriptPath <- "/Volumes/CBioinformatics/jonathan/Rprojects/Methylation-scripts/import_lab_worksheets.R"
# cmd <- cronR::cron_rscript(scriptPath)
# cronR::cron_add(cmd, frequency = "*/5 * * * *", days_of_week=c(1:6))
# cronR::cron_njobs(); cronR::cron_ls()
