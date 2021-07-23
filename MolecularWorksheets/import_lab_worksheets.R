## Read Excel Worksheets
## Jonathan Serrano
## April 20, 2021
Sys.setenv(R_ENABLE_JIT = T); compiler::enableJIT(3)

api.tkn = "#####moleclularDbPath#######" # REDCap API token
ws.type = NULL # Default is NULL "DNA" "RNA" "PACT" "FUSION" will upload if NULL
doAll = F # Set to TRUE to import all worksheets from this month

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
      "cronR"
    )
  }
  if (!require("easypackages")) {
    install.packages(
      "easypackages",
      dependencies = T,
      verbose = T,
      upgrade = "always"
    )
  }
  tryCatch(
    expr = {
      easypackages::libraries(pkgs)
    },
    warning = function(cond) {
      missed <- !pkgs %in% rownames(installed.packages())
      if (any(missed)) {
        easypackages::install_packages(pkgs[missed])
      }
    }
  )
}

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
  message(ws)
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

removePattern <- function(file.list) {
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

get.sheet.path <- function(ws.type, td) {stopifnot(!is.null(ws.type))
  drive = file.path("", "Volumes", "molecular", "MOLECULAR LAB ONLY")
  folder <- switch(
    ws.type,
    "PACT" = file.path("NYU PACT Patient Data", "Workbook"),
    "FUSION" = file.path("NYU FUSION SEQer", paste(format(td, "%Y"), "reports")),
    "METH" = file.path(drive,"NYU-METHYLATION/WORKSHEETS",paste(format(td, "%Y"))),
    file.path("Restricted-Molecular LAB staff ONLY", "PROCESSING&EXTRACTION", ws.type)
  )
  switch(ws.type,
         "FUSION" = return(file.path(drive, folder)),
         "PACT" = return(file.path(drive, folder, format(td, "%Y"))),
         "METH" = return(folder),
         return(file.path(drive, folder, format(td, "%Y"),format(td,"%B")))
  )
}

clean.fusion <- function(run.fis,ws,mgfs) {
  keywds <- c(".r_and_d_results.txt" = "",
              "/var/www/analysis/" = "")
  qsFix <- c("/" = " ", "_S" = " S", "-Q2" = " Q2")
  dataRead <-
    read.table(
      run.fis[ws],
      sep = "\t",
      header = F,
      row.names =NULL,
      stringsAsFactors = F
    )
  finam <- as.list(dataRead[1,][-1])

  finam <- stringr::str_replace_all(finam,keywds)
  finam <-
    stringr::str_split_fixed(stringr::str_replace_all(finam,qsFix)," ",4)
  colnames(finam) <- set.headers("FUSION")
  finam <- as.data.frame(finam)
  finam$fs_batch <- mgfs[ws]
  return(finam)
}

read.fusion.ws <- function(run.fis,ws) {
  cntrls <- c("NTC-H2O","SERASEQ-125ng","NC-IVS0035")
  mgfs <-
    stringr::str_split_fixed((stringr::str_replace_all(run.fis,ws,"")), "/", 3)[, 2]
  dat.io <-
    foreach::foreach(ws = 1:length(run.fis), .combine = "rbind") %dopar% {
      return(clean.fusion(run.fis, ws, mgfs))
    }
  drpSam <-
    !grepl(paste(cntrls, collapse = "|"), dat.io[, 2], ignore.case = T)
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
  #difftime(Sys.Date(),file.mtime(file.list)) <=1
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
  run.fld <- list.dirs(path = ws,
                       full.names = T,
                       recursive = F)
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

  } else {
    return(NULL)
  }
}

# Select only files created today (read only new samples)
grab.recent <- function(ws.type, td = NULL, doAll = F) {
  if (class(td) != "Date") {td = as.Date(td);doAll <- T}
  ws <- get.sheet.path(ws.type, td)
  fend = ifelse(grepl("PACT", ws.type) == F, "*.xlsm", "*.xlsx")
  Check.File.Paths(ws.type, ws)
  if (ws.type == "FUSION") {return(grab.fusion.ws(ws, td))
  } else
  if (ws.type == "PACT") {file.list <- dir(ws, pattern = fend, full.names = T)
  } else {file.list <- dir(ws, pattern = fend, full.names = T)}
  if (doAll == T) {return(removePattern(file.list))}
  file.list <- removePattern(file.list)
  recentMod <- grab.today(file.list)
  if (is.null(recentMod)) {
    return(recentMod)
    } else {
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
    dataRead <- 1:length(file.list) %>% purrr::map_df(function(i){
      return((readxl::read_excel(file.list[i], shNam[i], xl.range, naHeader, "text") %>% dplyr::filter(
        rd_number != 0 & rd_number !="MP-20-##"))[,-4])
      })
    return(dataRead)
    }
  foreach::foreach(ws = 1:length(file.list), .combine = "rbind", .packages = "readxl") %dopar% {
    dataRead <- tryCatch(
      expr = {
        readxl::read_excel(file.list[ws], shNam[ws], xl.range, naHeader, "text")
      },
      error = function(cond) {
        msg <-crayon::bgYellow(paste0("File was missing the sheet named: ", shNam[ws]))
        cat(msg); cat(" reading Sheet1 instead")
        return(readxl::read_excel(file.list[ws], sheet = 1, xl.range, naHeader, "text"))
      }
    )
    as.data.frame(dataRead[!is.na(dataRead[, 2]),])
  }
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
    dat.io <-as.data.frame(files.read[c(files.read[, 6] != ""), c(6:10, 13:15)])
  }
  if(ws.type=="METH"){
    dat.io <- as.data.frame(files.read[!is.na(files.read[, 2]), c(2,3,1,4)])
  } else {
    dat.io <- as.data.frame(files.read[!is.na(files.read[, 2]), c(2, 3, 7:9)])
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
    as.data.frame(redcapAPI::exportRecords(
      rcon,
      fields = fl,
      survey = F,
      form_complete_auto = F
    ))
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

dropNAvals <- function(dat.io) {
  if (length(names(table(dat.io == "NA"))) > 1) {
    dat.io[dat.io == "NA"] <- ""
  }
  if (length(names(table(is.na(dat.io)))) > 1) {
    dat.io[is.na(dat.io)] <- ""
  }
  return(dat.io)
}

checkDupes <- function(rcon, dat.io, block = NULL) {
  #if(block){col2swap = "block"
  if (!is.null(dat.io) & nrow(dat.io) > 0) {
    dat.io <- dropNAvals(dat.io)
    col2swap <- ifelse(any(grepl("b_number", colnames(dat.io)) == T), "b_number", "rna_number")
    bqCurr <- pullRedcap(dat.io, rcon, col2swap)
    if (!is.null(bqCurr) & nrow(bqCurr) > 0) {
      rcd = NULL
      bqSwap <- foreach::foreach(rcd = 1:nrow(bqCurr), .combine = "rbind") %dopar% {combBQ(dat.io, bqCurr, rcd, col2swap)}
      dat.io[, 2] <- paste(bqSwap)
    }
  }
  return(as.data.frame(dat.io))
}

msg.io <- function(dat.io = NULL, ws.type = NULL, rez = NULL) {
  if (!is.null(rez)){
    message("Imported Records:")
    cat(rez, sep = "\n")
  }
  if (is.null(dat.io) & !is.null(ws.type)) {
    cat(crayon::bgGreen(paste(ws.type, ": Import is up-to-date\n")))
  }
  if (length(dat.io) > 0 & !is.null(ws.type)) {
    io_name <-
      paste(format(Sys.Date(), "%b%d"), ws.type, "sheets.csv", sep = "_")
    message("\nSaving CSV: ", io_name, "\n")
    readr::write_csv(dat.io, io_name, na = "")
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
  msg.io(NULL, NULL, rez = result)
}

csv2redcap <- function(dat.io, api.tkn, ws.type) {
  ur = "https://redcap.nyumc.org/apps/redcap/api/"
  rcon <- redcapAPI::redcapConnection(ur, api.tkn)
  if (is.null(dat.io)) {return(NULL)}
  dat.io <- checkDupes(rcon, dat.io)
  msg.io(dat.io, ws.type, NULL)
  pushRedcap(jsonlite::toJSON(dat.io), rcon)
}

fsGrab <- function(wst, tk) {
  dat.io <- grab.recent(wst, Sys.Date())
  csv2redcap(dat.io, tk, wst)
  return(NULL)
}

automate.calls <- function(ws.type = "DNA", api.tkn = NULL, td = NULL, doAll = F) {
  xl.range <- switch(
    ws.type,
    "PACT" = "W7:AM41",
    "MYLO" = "B7:E42",
    "METH" = "A2:E97",
    "FUSION" = fsGrab(ws.type, api.tkn),
    "A2:J86"
  )
  if (!is.null(xl.range)) {
    file.list <- grab.recent(ws.type, Sys.Date(), doAll)
    if (!is.null(file.list) & length(file.list) > 0) {
      dat.io <- xL2csv(file.list, xl.range, ws.type)
      if (!is.null(dat.io) & !is.null(api.tkn) & length(dat.io) > 0) {
        csv2redcap(dat.io, api.tkn, ws.type)
      }
    }
  }
  msg.io(NULL, ws.type, NULL)
}

sheets.to.redcap <-
  function(api.tkn = NULL,ws.type = NULL,doAll = F,td = Sys.Date()) {
    check.packages()
    if (is.null(api.tkn)) {warning("apitoken was NULL")}
    doParallel::registerDoParallel(parallel::detectCores() - 2)
    if (is.null(ws.type)) {ws.type <- c("DNA", "RNA")}
    automate.calls(ws.type[1], api.tkn, td, doAll)
    if (!is.na(ws.type[2])) {
      automate.calls(ws.type[2], api.tkn, td, doAll)
      automate.calls("PACT", api.tkn, td, doAll)
      automate.calls("FUSION", api.tkn, td, doAll)
      automate.calls("METH", api.tkn, td, doAll)
    }
  }

sheets.to.redcap(api.tkn,ws.type,doAll)

# cronR Job start --------------------------------------------------
#if(!require("cronR")){install.packages("cronR")};library("cronR")
# scriptPath <- "/path/TO/DIR/import_lab_worksheets.R"
# cmd <- cronR::cron_rscript(scriptPath)
# cronR::cron_add(cmd, frequency = "*/5 * * * *", days_of_week=c(1:6))
# cronR::cron_njobs(); cronR::cron_ls()

