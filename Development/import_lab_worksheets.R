#!/usr/bin/env Rscript
## ---------------------------
## Script name: import_lab_worksheets.R
## Purpose: Read and Import all Excel Worksheets into REDCap
## Author: Jonathan Serrano
## Created: April 20, 2021
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------

api.tkn = "8XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" # REDCap API token
NULL -> ws.type     # Default is NULL "DNA" "RNA" "PACT" "FUSION" will upload if NULL
T ->    doAll       # Set TRUE to import all worksheets from current month of td
F ->    allMths     # Set TRUE to import all worksheets from the input year
td <-   Sys.Date()  # default Sys.Date() td (today's date) can be a string ie "2022-01-26"

formals(library)$quietly <- T
formals(library)$warn.conflicts <- F
formals(require)$warn.conflicts <- F
formals(install.packages)$dependencies <- T
formals(install.packages)$verbose <- T
formals(install.packages)$ask <- F
options("install.packages.compile.from.source" = "never")

rqPk <- function(pk) {if (!require(pk, character.only = T)) {
    install.packages(pk, warn.conflicts = F, ask = F, dependencies = T, Ncpus = 5)
    } else{library(pk, warn.conflicts = F, mask.ok = T, character.only = T, quietly = T)}
}

check.packages <- function(pkgs = NULL) {
  if (is.null(pkgs)) {
    pkgs <- c(
      "readxl",
      "stringi",
      "readr",
      "parallel",
      "RCurl",
      "foreach",
      "redcapAPI",
      "doParallel",
      "tidyr",
      "chron",
      "crayon",
      "compiler",
      "purrr",
      "tidyverse",
      "cronR",
      "plyr",
      "tidyverse"
    )
  }
    if (!requireNamespace("easypackages", quietly = T)) {install.packages("easypackages", dependencies = T, verbose = T, upgrade = "always")}
    tryCatch(
        expr = {easypackages::libraries(pkgs)},
        warning = function(cond) {
            missed <- !pkgs %in% rownames(installed.packages())
            if (any(missed)) {easypackages::install_packages(pkgs[missed])}
            }
        )
}


check.packages()
Sys.setenv(R_ENABLE_JIT = T)
compiler::enableJIT(3)
rqPk("devtools")
rqPk("pacman")
if(!requireNamespace("Rmpi", quietly = T)){install.packages("Rmpi",type="source")}
do.call(pacman::p_load, list("doSNOW", "doParallel","doMPI","magrittr","dplyr"))

mountMsg <- function(msg) {
  if (msg == 1) {cat(crayon:::bgBlue("\nChecking files... "))}
  if (msg == 2) {
    cat(crayon::red("\nDirectory not found, stopping execution\n"))
    cat(crayon::white$bgRed$bold("Mount the network drive or edit your input path to worksheets folder:"))
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


removePattern <- function(file.list,file.patterns=NULL) {
  if(is.null(file.patterns)){
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
      ),
      "METH"= c("rd_number","record_id", "b_number", "conc", "barcode_and_row_column")
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
    "METH" = file.path(drive,"NYU-METHYLATION/WORKSHEETS",paste(format(td, "%Y"))),
    file.path("Restricted-Molecular LAB staff ONLY", "PROCESSING&EXTRACTION", ws.type)
    )
    finalDir <- switch(ws.type,
         "FUSION" = return(file.path(drive, folder)),
         "PACT" = return(file.path(drive, folder, format(td, "%Y"))),
         "METH" = return(folder),
         return(file.path(drive, folder, format(td, "%Y")))
    )
    if(ws.type=="PACT"){
        if(as.numeric(format(td, "%Y")) < 2022){
            return(file.path(finalDir , format(td,"%B")))
        }else{return(file.path(finalDir , paste(format(td,"%m"), format(td,"%B"))))}
    }else{return(file.path(finalDir, format(td,"%B")))}

}

clean.fusion <- function(run.fis,ws,mgfs) {
  keywds <- c(".r_and_d_results.txt" = "",
              "/var/www/analysis/" = "")
  qsFix <- c("/" = " ", "_S" = " S", "-Q2" = " Q2")
  dataRead <- read.table(
      run.fis[ws],
      sep = "\t",
      header = F,
      row.names = NULL,
      stringsAsFactors = F
  )
  finam <- as.list(dataRead[1,][-1])
  finam <- stringr::str_replace_all(finam,keywds)
  finam <- stringr::str_split_fixed(stringr::str_replace_all(finam,qsFix)," ",4)
  colnames(finam) <- set.headers("FUSION")
  finam <- as.data.frame(finam)
  finam$fs_batch <- mgfs[ws]
  return(finam)
}

read.fusion.ws <- function(run.fis, ws) {
    cntrls <- c("NTC-H2O", "SERASEQ-125ng", "NC-IVS0035")
    mgfs <-
        stringr::str_split_fixed((stringr::str_replace_all(run.fis, ws, "")), "/", 3)[, 2]
    dat.io <-
        foreach::foreach(
            ws = 1:length(run.fis),
            .combine = "rbind",
            .packages = c("readxl", "foreach", "doParallel")
        ) %dopar% {
            return(clean.fusion(run.fis, ws, mgfs))
        }
    drpSam <- !grepl(paste(cntrls, collapse = "|"), dat.io[, 2], ignore.case = T)
    return(dat.io[drpSam, c(2, 3, 1, 4, 5)])
}

msgTime <- function(modTime, accessTime) {
  recently <- modTime <= accessTime + 1
  if (any(recently)) {
    mm1 <- crayon::green(paste("File last imported/accessed:", accessTime, "hours ago\n\n"))
    mm2 <- crayon::blue(paste("File last modified:", modTime, "hours ago\n", mm1))
    cat(mm2)
    return(recently)} else {return(F)}
}

minusTi <- function(rf, fiNum) {
  return(
    chron::hours(Sys.time()) -
      chron::hours(switch(fiNum, "1" = file.mtime(rf), "2" = file.info(rf)$atime))
    )
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
      #chron::hours(file.info(recentFi)$atime)-chron::hours(file.mtime(recentFi))<2
      if (any(todays)) {
        recentFi <- recentFi[todays]
      } else {
        return(NULL)}
      fiTime <- msgTime(minusTi(recentFi, 1), minusTi(recentFi, 2))
      if (any(fiTime)) {
        print(recentFi[fiTime])
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

gmt <- compiler::cmpfun(get.mod.time)

gtd <- function(file.list) {
  #return(file.list)
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

grab.today <- compiler::cmpfun(gtd)

grab.fusion.ws <- function(ws, td) {
  run.fld <- list.dirs(path = ws, full.names = T, recursive = F)
  run.fld <- run.fld[which.max(file.mtime(run.fld))]
  frl <- list.dirs(run.fld, full.names = T, recursive = F)
  run.fld <- file.path(frl, "summaries")
  run.fis <- dir(run.fld, full.names = T)
  file.list <-
    run.fis[grepl("Summary.r_and_d_results.txt", run.fis)]
  file.list <- grab.today(file.list)
  #file.list <-gmt(run.fis)
  if (length(file.list) > 0) {
    print(file.list)
    dat.io <- read.fusion.ws(file.list, ws)
    dat.io[, 1] <- make.unique(paste0(dat.io[, 1]), sep = "_")
    return(dat.io)
  }else {return(NULL)}
}

# Select only files created today (read only new samples)
grab.recent <- function(ws.type, td = NULL, doAll = F) {
    if (class(td) != "Date") {
        td <- as.Date(td)
        doAll <- T
        }
  ws <- get.sheet.path(ws.type, td)
  fend = ifelse(grepl("PACT", ws.type) == F, "*.xlsm", "*.xlsx")

  Check.File.Paths(ws.type, ws)
  if (ws.type == "FUSION") {return(grab.fusion.ws(ws, td))
  } else
  if (ws.type == "PACT") {
    file.list <- list.dirs(ws)
    file.patterns = c("VAL", "PM","dist","Test","-TumorLibraries","-TIPSTEST")
    file.list <-removePattern(file.list,file.patterns)[-1]
    xlFiles <-file.path(file.list, paste0(basename(file.list),".xlsx"))
    file.list<- grab.today(xlFiles[file.exists(xlFiles)])
    return(file.list)
  } else {
    file.list <- dir(ws, pattern = fend, full.names = T)
  }
  if (doAll == T) {return(removePattern(file.list))}
  file.list <- removePattern(file.list)
  recentMod <- grab.today(file.list)
  if (is.null(recentMod)) {return(recentMod)} else {
    recentFi <- gmt(recentMod)
    return(recentFi)
    }
}

parse.ws <- function(file.list, naHeader, xl.range, ws.type) {
  shNam <- switch(
    ws.type,
    "METH" = rep(3, length(file.list)),
    "PACT" = substr(basename(file.list), 1, nchar(basename(file.list)) - 5),
    rep("INFO_ENTRY", length(file.list))
    )
  ws = NULL
  if(ws.type=="METH"){
    dataRead <- 1:length(file.list) %>%
      purrr::map_df(function(i){
      return((readxl::read_excel(file.list[i], shNam[i], xl.range, naHeader, "text") %>%
                  dplyr::filter(rd_number != 0 & rd_number !="MP-20-##"))[,-4])})
    return(dataRead)
    }
  docsRead <-
      foreach::foreach(ws = 1:length(file.list), .combine = "rbind", .packages = c("readxl","foreach","doParallel")) %dopar% {
              dataRead <- tryCatch(
          expr = {readxl::read_excel(file.list[ws], shNam[ws], xl.range, naHeader, "text")},
          error = function(cond) {
            msg <-crayon::bgYellow(paste0("File was missing the sheet named: ", shNam[ws]))
            cat(msg); cat(" reading Sheet1 instead")
            return(readxl::read_excel(file.list[ws], sheet = 1, xl.range, naHeader, "text"))
          }
        )
        doc <- as.data.frame(dataRead[!is.na(dataRead[, 2]),])
        if(nrow(doc)>0){
            doc$fileName <- paste0(file.list[ws])
            doc
        }
      }
  return(docsRead)
}

# This function filters columns for input and auto-correct typos in the data
clean.output <- function(files.read, ws.type) {
  autocorrect <-
    c(
      " " = "-",
      "NGS2" = "NGS-2",
      "RD2" = "RD-2",
      "RD1" = "RD-1",
      "NGS1" = "NGS-1"
    )
  illegals <-
    c(
      "/" = " ",
      "\\+" = "",
      "&" = " ",
      "#" = "",
      "\\\\"  = " ",
      "<" = "",
      ">" = ""
    )
  spaceFix <- c("\\," = "", "  " =" ")
  cntrls <- c("NTC" = "", "SC" = "", "NC" = "")
  if (grepl("PACT", toupper(ws.type)) == T) {
    files.read[, 6] <- stringr::str_replace_all(trimws(files.read[, 6]), cntrls)
    dat.io <-as.data.frame(files.read[c(files.read[, 6] != ""), c(6:10, 13:15,ncol(files.read))])
  }
  if(ws.type=="METH"){
    dat.io <- as.data.frame(files.read[!is.na(files.read[, 2]), c(2,3,1,4,ncol(files.read))])
  } else {
    dat.io <- as.data.frame(files.read[!is.na(files.read[, 2]), c(2, 3, 7:9,ncol(files.read))])
    dat.io[, 5] <- toupper(stringr::str_replace_all(trimws(dat.io[, 5]), spaceFix))
  }
  dat.io[, 1] <- toupper(stringr::str_replace_all(trimws(dat.io[, 1]), autocorrect))
  dat.io[, 1] <- make.unique(as.character(stringr::str_replace_all(paste0(dat.io[, 1]), illegals)), sep = "_")
  dat.io[, 4] <- paste0(stringr::str_replace_all(paste0(dat.io[, 4]), illegals)) #tumor_percent number
  return(as.data.frame(dat.io))
}

xL2csv <- function(file.list, xl.range, ws.type) {
  message("\nReading Worksheets...")
  if (length(file.list) > 0) {
    print(file.list)
    files.read <- parse.ws(file.list, set.headers(toupper(ws.type)), xl.range, ws.type)
    return(clean.output(files.read, ws.type))
  } else {
    warning("File not readable, no files to import")
  }
}

pullRedcap <- function(dat.io, rcon, col2swap) {
  fl <- c("record_id", col2swap)
  bqCurr <-
    redcapAPI::exportRecords(
      rcon,
      fields = fl,
      survey = F,
      form_complete_auto = F
    )
  bqCurr<-as.data.frame(bqCurr)
  rownames(bqCurr) <- bqCurr[, 1]
  bqCurr <- bqCurr[dat.io[, 1], ]
  return(bqCurr)
}

combBQ <- local(function(dat.io, bqCurr, rcrd, col2swap) {
  datSwap = NULL
  theRecord <- paste0(bqCurr[rcrd, col2swap])
  if (is.na(theRecord)) {
    datSwap = paste0(dat.io[rcrd, 2]) # if current B/Q is empty paste new
  } else {
    alreadyIn <- grepl(paste0(dat.io[rcrd, 2]), paste0(bqCurr[rcrd, 2]))
    if (alreadyIn == T & is.na(bqCurr[rcrd, 2])) {
      datSwap = paste0(dat.io[rcrd, 2]) # check if B/Q is already in record
    } else {
      if (!is.na(dat.io[rcrd, 2]) & alreadyIn == F) {
        datSwap = paste(bqCurr[rcrd, 2], dat.io[rcrd, 2])
      } else {
        datSwap = paste0(bqCurr[rcrd, 2])
      } # if the new B/Q is not <NA>
    }
  }
  return(datSwap)
})

DropAllNAvals <- function(dat.io) {
  if (length(names(table(dat.io == "NA"))) > 1) {
    idx <- which(dat.io == "NA", arr.ind=TRUE)
    if(length(idx)>0){dat.io[idx] <- ""}
    idx <- which(is.na(dat.io), arr.ind=TRUE)
    if(length(idx)>0){dat.io[idx] <- ""}
  }
  #if (length(names(table(is.na(dat.io)))) > 1) {dat.io[is.na(dat.io)] <- ""}
  return(dat.io)
}

Check4Dupes <- function(rcon, dat.io, block = NULL) {
  #if(block){col2swap = "block"
  if (!is.null(dat.io) & nrow(dat.io) > 0) {
    dat.io <- DropAllNAvals(dat.io)
    col2swap <- ifelse(any(grepl("b_number", colnames(dat.io)) == T), "b_number", "rna_number")
    bqCurr <- pullRedcap(dat.io, rcon, col2swap)
    if (!is.null(bqCurr) & nrow(bqCurr) > 0) {
      rcd = NULL
      bqSwap <- foreach::foreach(rcd = 1:nrow(bqCurr), .combine = "rbind",.packages = c("readxl","foreach","doParallel")) %dopar% {
        combBQ(dat.io, bqCurr, rcd, col2swap)
      }
      bqSwap <- paste(bqSwap)
      bqSwap <-gsub("NA ","",bqSwap)
      dat.io[, 2] <- bqSwap
    }
  }
  return(as.data.frame(dat.io))
}

MsgIoImp <- function(dat.io = NULL, ws.type = NULL, rez = NULL) {
  if (!is.null(rez)){
    message("Imported Records:")
    cat(rez, sep = "\n")
  }
  if (is.null(dat.io) & !is.null(ws.type)) {
    cat(crayon::bgGreen(paste(ws.type, ": Import is up-to-date\n")))
  }
  if (length(dat.io) > 0 & !is.null(ws.type)) {
    io_name <- paste(format(Sys.Date(), "%b%d"), ws.type, "sheets.csv", sep = "_")
    message("\nSaving CSV: ", io_name, "\n")
    dat.io <- DropAllNAvals(dat.io)
    readr::write_csv(dat.io, file.path("~/Desktop",io_name), na = "", quote="none")
  }
}

pushRedcap <- function(dat, rcon) {
  result <- RCurl::postForm(
    uri = rcon$url,
    token = rcon$token,
    content = 'record',
    format = 'json',
    overwriteBehavior = 'normal',
    forceAutoNumber = 'false',
    data = dat,
    returnContent = 'ids',
    type = 'flat',
    returnFormat = 'csv'
  )
  MsgIoImp(NULL, NULL, rez = result)
}

csv2redcap <- function(dat.io, api.tkn, ws.type) {
  ur = "https://redcap.nyumc.org/apps/redcap/api/"
  rcon <- redcapAPI::redcapConnection(ur, api.tkn)
  if (is.null(dat.io)) {return(NULL)}
  dat.io <- Check4Dupes(rcon, dat.io)
  MsgIoImp(dat.io, ws.type, NULL)
  dat.io <- dat.io[,!names(dat.io) %in% "fileName"]
  dat.io <- DropAllNAvals(dat.io)
  if(ws.type=="METH"){
    dat.io <- dat.io[,1:3]
  }
  if(ws.type=="PACT"){
    colnames(dat.io) <- c("record_id", "paired_normal_pact", "b_number", "epic_id_pact", "test_number_pact")
  }
  pushRedcap(jsonlite::toJSON(dat.io), rcon)
}

FsGrab <- function(wst, tk) {
  dat.io <- grab.recent(wst, td) #Sys.Date()
  csv2redcap(dat.io, tk, wst)
  return(NULL)
}

AutomateFunCalls <- function(ws.type = "DNA", api.tkn = NULL, td = NULL, doAll = F) {
  xl.range <- switch(
    ws.type,
    "PACT" = "W7:AM41",
    "MYLO" = "B7:E42",
    "METH" = "A2:E97",
    "FUSION" = FsGrab(ws.type, api.tkn),
    "A2:J86"
    )
  if (!is.null(xl.range)) {
    file.list <- grab.recent(ws.type, td, doAll)
    if (!is.null(file.list) & length(file.list) > 0) {
      dat.io <- xL2csv(file.list, xl.range, ws.type)
      if (!is.null(dat.io) & !is.null(api.tkn) & length(dat.io) > 0){
        csv2redcap(dat.io, api.tkn, ws.type)}}
  }
  MsgIoImp(NULL, ws.type, NULL)
}

SheetsIntoRedcap <- function(api.tkn = NULL,ws.type = NULL,doAll = F, td) {
    if (is.null(api.tkn)) {warning("apitoken was NULL")}
    doParallel::registerDoParallel(parallel::detectCores() - 2)
    if (is.null(ws.type)) {ws.type <- c("DNA", "RNA")}
    AutomateFunCalls(ws.type[1], api.tkn, td, doAll)
    if (!is.na(ws.type[2])) {
      AutomateFunCalls(ws.type[2], api.tkn, td, doAll)
      AutomateFunCalls("PACT", api.tkn, td, doAll)
      AutomateFunCalls("FUSION", api.tkn, td, doAll)
      AutomateFunCalls("METH", api.tkn, td, doAll)
    }
}

DoWholeYear <- function(yr = "2021"){
    mnLi <- c(paste0("0",1:9), "11", "12")
    for (mnth in mnLi) {
        td <- paste(yr, mnth, "01", sep="-")
        assign("td", td)
        SheetsIntoRedcap(api.tkn, ws.type, doAll, td)
    }
}


if(allMths==T){
    td <- as.character(format(as.Date(td),"%Y"))
    DoWholeYear(td)
}else{
    SheetsIntoRedcap(api.tkn, ws.type, doAll, td)
}


# # cronR Job start --------------------------------------------------
# if(!require("cronR")){install.packages("cronR")};library("cronR")
# scriptPath <- "/Volumes/CBioinformatics/jonathan/Rprojects/Methylation-scripts/import_lab_worksheets.R"
# cmd <- cronR::cron_rscript(scriptPath)
# cronR::cron_add(cmd, frequency = "*/5 * * * *", days_of_week=c(1:6))
# cronR::cron_njobs(); cronR::cron_ls()
