#!/usr/bin/env Rscript
library("base"); gb <- globalenv(); assign("gb", gb)
if(!require("optparse")){install.packages("optparse", dependencies = T, quiet = F)}
library(optparse)

option_list <- list(
    optparse::make_option(c("-t", "--token"),
                          default = NULL,
                          help = "REDCap API token (default %default)",),
    optparse::make_option(c("-m", "--molecular"),
                          default = NULL,
                          help = "Single-spaced list of molecular TS or TM (default %default)")
)

parser <- optparse::OptionParser(option_list=option_list)
arguments <- optparse::parse_args(parser, positional_arguments = T)
opt <- arguments$options
args <- arguments$args

token <- opt$token
myMolecular <- strsplit(opt$molecular, " ")[[1]]

dsh="\n================"
dsh2="================\n"

# Displays the Input args -----
message(dsh,"Parameters input",dsh2)
message("token: ",token)
message("Args input:")
print (myMolecular)
message("")

stopifnot(!is.na(token))
stopifnot(!is.na(myMolecular))

# REDcap Heading Fields -----
flds = c("record_id","b_number","tm_number","accession_number","block","diagnosis",
         "organ","tissue_comments","run_number", "nyu_mrn")

# Load redcapAPI Package -----
if(suppressWarnings(!require("redcapAPI"))){
    params=list('nutterb/redcapAPI', dependencies=T, upgrade="always", type="source")
    do.call(devtools::install_github,c(params))
}

supM <- function(sobj){return(suppressMessages(suppressWarnings(sobj)))}
supPk <- function(sobj){return(suppressPackageStartupMessages(suppressWarnings(sobj)))}

# FUN: Checks if z-drive is accessible to the Rscript
checkMounts <- function(){
    molecDrive = "/Volumes/molecular/MOLECULAR LAB ONLY"
    zDrive = "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular"
    failMount <- ifelse(dir.exists(molecDrive),T,F)
    if(failMount!=T){
        cat("PATH does not exist, ensure path is mounted:")
        cat(crayon::white$bgRed$bold(molecDrive))
        cat("You must mount the network Z-drive path:\n")
        cat(crayon::white$bgRed$bold(zDrive),"\n")
        stopifnot(!any(failMount==T))
    } else {message("\n",crayon::bgGreen("Z-drive path is accessible"),"\n")}
}

# Functions to load packages and get redcap info -----
loadPacks <- function(){
    pkgs = c("data.table", "foreach", "openxlsx","jsonlite","RCurl","readxl","stringr",
             "webshot","tidyverse","crayon","plotly","htmlwidgets","tinytex")
    rlis = getOption("repos")
    rlis["CRAN"] = "http://cran.us.r-project.org"
    options(repos = rlis)
    invisible(lapply(pkgs, function(pk){
        if(suppressPackageStartupMessages(suppressWarnings(!require(pk, character.only=T)))){
            install.packages(pk,dependencies=T, verbose=T, repos="http://cran.us.r-project.org",
                             type="both")
        }}))
    if(!require("systemfonts")){devtools::install_github('r-lib/systemfonts')}
    if(!require("redcapAPI")){install.packages("redcapAPI", dependencies = T, type="both",ask=F)}
    if(!require("remotes")){install.packages("remotes", dependencies=T)}
    if(!require("chromote")){remotes::install_github("rstudio/chromote", upgrade ="never")}
    if(suppressPackageStartupMessages(suppressWarnings(!require("webshot2")))){
        remotes::install_github("rstudio/webshot2")
        }
    library("redcapAPI")
    library("chromote")
    #try(tinytex::tlmgr_install(),silent=T)
    #try(tinytex:::install_prebuilt(),silent=T)
}

# API Call functions -----
grabAllRecords <- function(flds, rcon){
    params = list(rcon, fields=flds, labels=F, dates = F, survey = F, dag = F, factors=F, form_complete_auto=F)
    dbCols <- do.call(redcapAPI::exportRecords, c(params))
    return(as.data.frame(dbCols))
}

# Database search function -----
searchDb <- function(vals, db){
    v2f <- paste(vals, collapse="|"); i=NULL
    return(foreach::foreach(i=1:ncol(db), .combine='rbind') %do% {db[grepl(v2f,db[,i]),]})
}

genQuery <- function(dbCol,vals2find){
    currCol <- vals2find[, dbCol]
    q1 <- currCol[currCol != 0 & !is.na(currCol) & currCol != ""]
    return(unique(q1))
}

# Get Methylation and Molecular Samples list ----
queryCases <- function(vals2find, db) {
    i=NULL
    if(!is.null(ncol(vals2find))){
    queryList <- foreach::foreach(i=1:ncol(vals2find), .combine="c") %do% {genQuery(i,vals2find)}
    methQuery <- searchDb(queryList, db)
    return(unique(methQuery))
    }else{print(vals2find)}

}

addOutputLinks <- function(output){
    winpath = "smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular/Molecular/MethylationClassifier/"
    yearPath <- lapply(stringr::str_split_fixed(output$run_number,"-",2)[,1], function(yr) {
        rnum <- NULL
        if(nchar(yr)>2){rnum <- substring(yr, 3)}else{rnum <- yr}
        if(nchar(yr)>0){paste0("20",rnum)}else{rnum}
    })
    output$report_complete <- ifelse(!is.na(output$run_number), "YES", "NOT_YET_RUN")
    output$'Report Link' <- paste0(winpath, yearPath,"/",output$run_number,"/",output$record_id,".html")
    output$'Report Link'[is.na(output$run_number)] <- ""
    output$'Report Path'<-output$'Report Link'
    return(output)
}

modifyOutput <- function(output,vals2find){
    output$Test_Number <- NA
    for(i in 1:nrow(output)) {
        theVal = NA
        for (var in 1:ncol(vals2find)) {
            pat <- vals2find[, var]
            currRow <- paste(output[i, ])
            theMatch <- which(pat %in% currRow[currRow != "0" & currRow != "NA"])
            if (length(theMatch) > 0) {
                theVal <- vals2find$`Test Number`[theMatch]
                }
        }
        if(!is.null(theVal)){
        output$Test_Number[i] <- theVal
        }
    }
    output <- addOutputLinks(output)
    return(output)
}

# Main workflow to get dataframe
getOuputData <- function(token, flds, inputSheet){
    apiUrl="https://redcap.nyumc.org/apps/redcap/api/"
    rcon <- redcapAPI::redcapConnection(apiUrl, token)
    vals2find <- data.frame(molecular= c(inputSheet))
    db <- grabAllRecords(flds, rcon)
    if(nrow(vals2find)!=0){
        output <- queryCases(vals2find, db)
    if (nrow(output) != 0) {
        output <- modifyOutput(output, vals2find)
        return(output)
    } else{
        return(NULL)
    }
    }else{
        return(NULL)
    }
}

# FUN: Sets your directory and sources the helper functions
sourceFuns2 <- function(workingPath = NULL) {
    mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/"
    script.list <- c("SetRunParams.R","CopyInputs.R","PACT_scripts/generateCNV.R")
    if (is.null(workingPath)) {workingPath = getwd()}
    scripts <- paste0(mainHub, script.list)
    invisible(lapply(scripts, function(i){supPk(devtools::source_url(i))}))
    #supPk(library("sest"))
    #supPk(library("mnp.v11b6"))
    supPk(require("plotly"))
    supPk(require("htmlwidgets"))
    gb$setDirectory(workingPath)
    return(gb$defineParams())
}

msgRDs <- function(rds,token){
    message("\nRD-numbers with idats:")

    assign("rds", rds)
rdPrnt <- paste(rds, collapse=" ")
    message(rdPrnt)
    message(dsh, crayon::bgMagenta("Starting CNV PNG Creation"),dsh2)
    sourceFuns2()
    ApiToken = token
    assign("ApiToken", ApiToken)
}

grabRDs <- function(rd_numbers, token){
    result <- gb$search.redcap(rd_numbers, token)
    result <- result[!is.na(result$barcode_and_row_column),]
    samplesheet_ID = as.data.frame(stringr::str_split_fixed(result[,"barcode_and_row_column"],"_",2))
    gb$writeFromRedcap(result, samplesheet_ID) # writes API export as minfi dataframe sheet
    gb$get.idats()  # copies idat files from return to current directory
}

msgCreated <- function(mySentrix){
    pngFiles <- paste0(file.path("~","Desktop",mySentrix[, 1]),"_cnv.png")
    cnvMade <- file.exists(pngFiles)
    if(any(cnvMade==F)){
        message("The following failed to be created:")
        print(pngFiles[!cnvMade])
        message("Try running again or check GitHub troubleshooting")
    }
    if(any(cnvMade==T)){
        message("The following were created successfully:")
        print(pngFiles[cnvMade])
    }
}

loopCNV <- function(mySentrix, asPNG){
    for (sam in rownames(mySentrix)) {
        sampleName <- mySentrix[sam, 1]
        fn = file.path("~", "Desktop", paste0(sampleName, "_cnv.png"))
        if (file.exists(fn)) {
            message("\nFile already exists, skipping:", fn, "\n")
        } else{
            sentrix <- mySentrix[sam, "SentrixID_Pos"]
            message("\nGetting RGset for ", sentrix, "\n")
            RGsetEpic <- gb$grabRGset(getwd(), sentrix)
            tryCatch(
                expr = {gb$gen.cnv.png2(RGsetEpic, sampleName, asPNG)},
                error = function(e) {
                    erTxt <-paste0("An error occured with ", mySentrix[sam, 1]," png creation:")
                    message(crayon::bgRed(erTxt), "\n", e)
                    message(crayon::bgGreen("Trying next sample"))}
            )
        }
    }
}

makeCNV <- function(myDt, asPNG = T) {
    mySentrix <- myDt[myDt[, "SentrixID_Pos"] %like% "_R0", ]
    if (nrow(mySentrix) > 0) {
        loopCNV(mySentrix, asPNG)
    } else{
        message("The RD-number(s) do not have idat files in REDCap:/n")
        print(myDt)
    }
    msgCreated(mySentrix)
    while (!is.null(dev.list())) {dev.off()}
}

startCNVmaker <- function(output, token) {
    rds <- output$record_id[output$report_complete == "YES"]
    if (length(rds) > 0) {
        msgRDs(rds, token)
        grabRDs(rds, token)
        myDt <- read.csv("samplesheet.csv")
        tryCatch(
            expr = {gb$makeCNV(myDt)},
            error = function(e) {
                message("The following error occured:\n", e)
                message("\n\nTry checking the troubleshooting section on GitHub:\n")
                message("https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/PACT_scripts/README.md\n")
            },
            finally={gb$copyOutputPng()}
        )
    } else{message(crayon::bgGreen("This run has no cases with methylation complete.
                                   No CNV png images to generate."))}
}

# Search REDCap Worksheets for MRN Match for output -------------------------------------
loadPacks()
checkMounts()
output <- getOuputData(token, flds, myMolecular)

# CNV PNG Creation -------------------------------------
if (!is.null(output)) {
    try(write.csv(output, file = "~/Desktop/methQueryResults.csv"))
    startCNVmaker(output, token)
} else{
    message("\nNo matching methylation idat files found for input:")
    print(myMolecular)
}
