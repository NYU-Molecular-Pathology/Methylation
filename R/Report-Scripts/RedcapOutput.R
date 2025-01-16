#!/usr/bin/env Rscript
## Script name: RedcapOutput.R
## Purpose: source of global functions for methylation REDCap database data import
## Date Created: August 11, 2022
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025

library(verbose = F, warn.conflicts = F, quietly = T, package = "jsonlite")
library(verbose = F, warn.conflicts = F, quietly = T, package = "crayon")
library(verbose = F, warn.conflicts = F, quietly = T, package = "RCurl")
library(verbose = F, warn.conflicts = F, quietly = T, package = "redcapAPI")

makePost <- function(dfNewRed, params){

    apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
    tk <- params$token
    data <- dfNewRed[1,]
    datarecord = jsonlite::toJSON((as.list(dfNewRed[1,])), auto_unbox = T)
    message("~~",crayon::bgBlue("Record Uploaded:"),"\n", datarecord)
    RCurl::postForm(
        apiLink, token = tk, content = 'record', format = 'csv', type = 'flat',
        data = datarecord, returnFormat = 'csv', overwriteBehavior = 'normal'
    )

    rcon <- redcapAPI::redcapConnection(apiLink, tk)
    res <- redcapAPI::importRecords(rcon, data, "normal", "ids", logfile = "REDCapImportLog.txt")
    message("REDCap Response:\n", res)
}

supM <- function(objTing){return(suppressMessages(suppressWarnings(objTing)))}
gb <- globalenv(); assign("gb", gb)

writeRedcapPred <- function(run_id = NULL, dfNewRed) {
    stopifnot(length(run_id) > 0 & !is.na(run_id) & !is.null(run_id))
    redDir <- file.path(fs::path_home(), "Desktop", run_id)
    if (!dir.exists(redDir)) {dir.create(redDir, recursive = T)}

    redcsv <- file.path(redDir, paste0(run_id, "_Redcap.csv"))
    if (file.exists(redcsv)) {
        dfRedcap = read.csv(redcsv, header = T, row.names = NULL)
        dfRedcap <- as.data.frame(dfRedcap, row.names = NULL)
        redDF <- rbind(dfRedcap, dfNewRed)
    } else{
        redDF <- dfNewRed
    }
    row.names(redDF) = NULL
    write.csv(redDF, redcsv, row.names = F)
}

SetDesktopOutput <- function(run_id){
    redfolder <- file.path(fs::path_home(),"Desktop", run_id)
    redcsv <- file.path(redfolder, paste0(run_id,"_v11_Redcap.csv"))
    redcapDF <- data.frame(
        record_id = character(),
        b_number = character(),
        barcode_and_row_column = character(),
        array_type = character(),
        classifier_sex = character(),
        classifier_score = character(),
        classifier_value = character(),
        subgroup = character(),
        subgroup_score = character(),
        mgmt_status = character(),
        mlh1_status = character(),
        mlh1_pos_loci = character(),
        second_tech = character(),
        primary_tech = character(),
        run_number = character(),
        tm_number = character()
    )
    if (!dir.exists(redfolder)) {dir.create(redfolder)}
    if (!file.exists(redcsv)) {write.csv(x = redcapDF, file = redcsv, row.names = F)}
}

FormatSuppInfo <- function(suppinfo){
    names(suppinfo) <- c(
        "ID",
        "Batch ID",
        "Tech",
        "B-number",
        "MP-TN#",
        "Sentrix ID",
        "Array type",
        "Material type",
        "Gender"
    )
    suppinfo <- as.data.frame(suppinfo)
    colnames(suppinfo) <- ""
    return(suppinfo)
}


DebugDataFrame <- function(e, gb) {
    eMsg <- crayon::bgRed("Potentially missing variable(s) in REDCap dataframe:")
    message(e,"\n",eMsg,"\n")
    fixNull <- function(obj) {
        if (is.null(obj) |
            length(obj) == 0) {
            return("NONE or Missing")
        } else{
            return(obj)
        }
    }

    gb$is450k <- gb$RGset@annotation[["array"]] != "IlluminaHumanMethylationEPIC"
    array_opt1 <- ifelse(gb$is450k, yes = "450k", no = "EPIC")
    array_opt <- ifelse(gb$RGset@annotation[["array"]] == "IlluminaHumanMethylationEPICv2", yes = "EPICV2", no = array_opt1)

    if (array_opt == "EPICV2") {
        familia <- gb$outList["family", "predicted"]
        fscore <- gb$outList["family", "maxscore"]
        subfam <- gb$outList["subclass", "predicted"]
        subScore <- gb$outList["subclass", "maxscore"]
        mgmtStat1 <- as.data.frame(gb$mgmtValues)

    } else{
        out <- fixNull(gb$outList$out)
        familia <- fixNull(gb$outList$out_class_family$`Methylation Family`[1])
        fscore <- fixNull(gb$outList$out_class_family$`Class Score`[1])
        subfam <- fixNull(gb$out$`Methylation Subgroup`[1])
        subScore <- fixNull(gb$out$`Subgroup Score`[1])
        mgmtStat1 <- fixNull(gb$mgmtValues$mgmtVal)
    }

    mlh_status <- gb$mlh1Pred$theValue$m.reslt
    mlh_total <- gb$mlh1Pred$theValue$MLH1.pos.loci
    mlh1_status <- fixNull(paste0(mlh_status))
    mlh1_pos_loci <- fixNull(paste0(mlh_total))
    run_id <- fixNull(paste(gb$dat$run_id))

    message("sampleID: ", gb$dat$sampleID)
    message("paste(dat$bnumber): ", paste(gb$dat$bnumber))
    message("colnames(RGset): ", colnames(gb$RGset))
    message("is450k: ", array_opt)
    message("sex: ", tolower(gb$msetDat$sex))
    message("familia: ", familia)
    message("fscore: ", fscore)
    message("subfam: ", subfam)
    message("subScore: ", subScore)
    message("mgmtStat1$Status: ", paste0(mgmtStat1$Status))
    message("mlh1_status: ", mlh1_status)
    message("mlh1_pos_loci: ", mlh1_pos_loci)
    message("second_tech: ", paste(gb$dat$tech2))
    message("primary_tech: ", paste(gb$dat$tech))
    message("run_id: ", run_id)
    message("mlh1_pos_loci: ", paste(gb$dat$mp_number))
}


GetRedcapDF <- function(gb) {
    gb$is450k <- gb$RGset@annotation[["array"]] != "IlluminaHumanMethylationEPIC"
    array_opt1 <- ifelse(gb$is450k, yes = "450k", no = "EPIC")
    array_opt <- ifelse(gb$RGset@annotation[["array"]] == "IlluminaHumanMethylationEPICv2", yes = "EPICV2", no = array_opt1)

    familia <- gb$outList["family", "predicted"]
    fscore <- gb$outList["family", "maxscore"]
    subfam <- gb$outList["subclass", "predicted"]
    subScore <- gb$outList["subclass", "maxscore"]
    mgmtStat1 <- as.data.frame(gb$mgmtValues)
    mlh_status <- gb$mlh1Pred$theValue$m.reslt
    mlh_total <- gb$mlh1Pred$theValue$MLH1.pos.loci

    dfNewRed <- data.frame(
        record_id = paste0(gb$dat$sampleID),
        b_number = paste0(gb$dat$bnumber),
        barcode_and_row_column = paste0(colnames(gb$RGset)),
        array_type = array_opt,
        classifier_sex = paste0(tolower(gb$msetDat$sex)),
        classifier_score = gsub(",","", familia),
        classifier_value = gsub(",","", fscore),
        subgroup = gsub(",","", subfam),
        subgroup_score = gsub(",","", subScore),
        mgmt_status = paste0(mgmtStat1$Status),
        mlh1_status = paste0(mlh_status),
        mlh1_pos_loci = paste0(mlh_total),
        second_tech = paste0(gb$dat$tech2),
        primary_tech = paste0(gb$dat$tech),
        run_number =  paste0(gb$dat$run_id),
        tm_number = paste0(gb$dat$mp_number)
    )
    stopifnot(nrow(dfNewRed) > 0 & !is.null(dfNewRed$run_number))
    return(dfNewRed)
}


TryREDCap <- function(gb) {
    tryCatch(
        gb$writeRedcapPred(gb$dat$run_id, dfNewRed = gb$GetRedcapDF(gb)),
        error = function(e) {
            gb$DebugDataFrame(e, gb)
            stop("REDCap csv saving failed!")
        }
    )
}

CheckScoreCsv <- function(targets){
    deskDir <- file.path(fs::path_home(), "Desktop", targets$RunID[1])
    if (!dir.exists(deskDir)) {dir.create(deskDir)}
}
