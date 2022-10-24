makePost <- function(dfNewRed, params){
    library(verbose=F, warn.conflicts = F, quietly = T, package= "jsonlite")
    library(verbose=F, warn.conflicts = F, quietly = T, package= "crayon")
    library(verbose=F, warn.conflicts = F, quietly = T, package= "RCurl")
    ur = "https://redcap.nyumc.org/apps/redcap/api/"
    tk <- params$token
    datarecord = jsonlite::toJSON((as.list(dfNewRed[1,])), auto_unbox=T)
    message("~~",crayon::bgBlue("Record Uploaded:"),"\n", datarecord)
    RCurl::postForm(
        ur, token = tk, content = 'record', format = 'csv', type = 'flat',
        data = datarecord, returnFormat = 'csv', overwriteBehavior='normal'
    )
}
supM <- function(objTing){return(suppressMessages(suppressWarnings(objTing)))}
gb <- globalenv(); assign("gb", gb)

writeRedcapPred <- function(run_id, dfNewRed) {
    redfolder <- file.path("~","Desktop", run_id)
    redcsv <- list.files(path=redfolder, pattern="_Redcap.csv", full.names=T)[1]
    if (!is.na(redcsv)) {
        dfRedcap = read.csv(redcsv, header = T, row.names = NULL)
        dfRedcap <- as.data.frame(dfRedcap, row.names = NULL)
        redDF <- rbind(dfRedcap, dfNewRed)
        row.names(redDF) = NULL
        write.csv(x = redDF, file = redcsv, row.names = F)
    } else{
        message("Saving the REDCap csv data failed for sample:", sampleID,
                "second upload attempt will be made at end of run."
        )
    }
}

SetDesktopOutput <- function(run_id){
    redfolder <- file.path("~","Desktop", run_id)
    redcsv <- file.path(redfolder, paste0(run_id,"_Redcap.csv"))
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

LoadMnpData <- function(is450k){
    message("Loading MNP Data...")
if(is450k==T){
    load("/Volumes/CBioinformatics/Methylation/Methylation_classifier_v11b4/mnp.v11b4/data/rfpred.v11b4.RData")
}else{
    load("/Volumes/CBioinformatics/Methylation/in_house/mnp.v116/mnp.v11b6/data/rfpred.v11b6.RData")
}
}

DebugDataFrame <-function(e, gb){
    eMsg <- crayon::bgRed("Potentially missing variable(s) in REDCap dataframe:")
    message(e,"\n",eMsg,"\n")
    fixNull <- function(obj) {if (is.null(obj)|length(obj)==0) {return("NONE or Missing")} else{return(obj)}}
    out <- fixNull(gb$outList$out)
    familia <- fixNull(gb$outList$out_class_family$`Methylation Family`[1])
    fscore <- fixNull(gb$outList$out_class_family$`Class Score`[1])
    subfam <- fixNull(gb$out$`Methylation Subgroup`[1])
    subScore <- fixNull(gb$out$`Subgroup Score`[1])
    mgmtStat <- fixNull(gb$mgmtValues$mgmtVal)
    mlh1_status <- fixNull(paste0(gb$mlhP[[1]][1]))
    mlh1_pos_loci <- fixNull(paste0(gb$mlhP[[1]][2]))
    run_id <- fixNull(paste(gb$dat$run_id))
    message("sampleID: ", gb$sampleID)
    message("paste(dat$bnumber): ", paste(gb$dat$bnumber))
    message("colnames(RGset): ", colnames(gb$RGset))
    message("is450k: ", paste0(ifelse(gb$is450k, yes = "450k", no = "EPIC")))
    message("sex: ", tolower(gb$msetDat$sex))
    message("familia: ", familia)
    message("fscore: ", fscore)
    message("subfam: ", subfam)
    message("subScore: ", subScore)
    message("mgmtStat$Status: ", mgmtStat$Status)
    message("mlh1_status: ", mlh1_status)
    message("mlh1_pos_loci: ", mlh1_pos_loci)
    message("second_tech: ", paste(gb$dat$tech2))
    message("primary_tech: ", paste(gb$dat$tech))
    message("run_id: ", run_id)
    message("mlh1_pos_loci: ", paste(gb$dat$mp_number))
}

GetRedcapDF <- function(gb) {
    out <- gb$outList$out
    familia <- gb$outList$out_class_family$`Methylation Family`[1]
    fscore <- gb$outList$out_class_family$`Class Score`[1]
    subfam <- gb$out$`Methylation Subgroup`[1]
    subScore <- gb$out$`Subgroup Score`[1]
    mgmtStat <- gb$mgmtValues$mgmtVal
    
    dfNewRed <- data.frame(
        record_id = gb$sampleID,
        b_number = paste(gb$dat$bnumber),
        barcode_and_row_column = colnames(gb$RGset),
        array_type = ifelse(gb$is450k, yes = "450k", no = "EPIC"),
        classifier_sex = tolower(gb$msetDat$sex),
        classifier_score = familia,
        classifier_value = fscore,
        subgroup = subfam,
        subgroup_score = subScore,
        mgmt_status = mgmtStat$Status,
        mlh1_status = paste0(gb$mlhP[[1]][1]),
        mlh1_pos_loci = paste0(gb$mlhP[[1]][2]),
        second_tech = paste(gb$dat$tech2),
        primary_tech = paste(gb$dat$tech),
        run_number =  paste(gb$dat$run_id),
        tm_number = paste(gb$dat$mp_number)
    )
    
    return(dfNewRed)
}

TryREDCap <- function(gb){
    gb$supM(tryCatch(
    gb$writeRedcapPred(gb$dat$run_id, dfNewRed = gb$GetRedcapDF(gb)),
    error = function(e) {
        gb$DebugDataFrame(e, gb)
    }
))

}

CheckScoreCsv <- function(targets){
    deskDir <- file.path("~", "Desktop", targets$RunID[1])
    if(!dir.exists(deskDir)){dir.create(deskDir)}
    #scoreFile <- file.path(deskDir, paste0(targets$RunID[1], "_v12.csv"))
    #if(!file.exists(scoreFile)){write.csv(v12df, file = scoreFile, row.names=F)}
}
