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
        message("REDCap Upload failed for sample:", sampleID,
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
if(is450k==T){
    load("/Volumes/CBioinformatics/Methylation/Methylation_classifier_v11b4/mnp.v11b4/data/rfpred.v11b4.RData")
}else{
    load("/Volumes/CBioinformatics/Methylation/in_house/mnp.v116/mnp.v11b6/data/rfpred.v11b6.RData")
}    
}
