#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library("base"); gb <- globalenv(); assign("gb", gb)
dsh="\n================"

# Main arguments input in comandline
#token<-args[1]; inputSheet<-args[2]
token<-args[1]; inputSheet<-args[2]

# Displays the Input args -----
message(dsh,"\nParameters input",dsh)
message("token: ",token)
message("inputSheet: ", inputSheet,"\n")

stopifnot(!is.na(token))
stopifnot(!is.na(inputSheet))

readFlag <- endsWith(inputSheet,".csv")==T
stopifnot(rlang::is_bool(readFlag))

# REDcap Heading Fields -----
flds = c("record_id","b_number","tm_number","accession_number","block","diagnosis",
         "organ","tissue_comments","run_number", "nyu_mrn")

# Load redcapAPI Package -----
if(suppressWarnings(!require("redcapAPI"))){
    params=list('nutterb/redcapAPI', dependencies=T, upgrade="always", type="source")
    do.call(devtools::install_github,c(params))
}

# Functions to load packages and get redcap info -----
loadPacks <- function(){
    pkgs = c("data.table", "foreach", "openxlsx","jsonlite","RCurl","readxl","stringr")
    rlis = getOption("repos")
    rlis["CRAN"] = "http://cran.us.r-project.org"
    options(repos = rlis)
    invisible(lapply(pkgs, function(pk){
        if(suppressWarnings(!require(pk, character.only=T))){
            install.packages(pk,dependencies=T, verbose=T, repos="http://cran.us.r-project.org", type="both")
        }}))
    if(!require("webshot2")){remotes::install_github("rstudio/webshot2")}
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

getFilePath <- function(inputSheet){
    drive = file.path("", "Volumes", "molecular", "MOLECULAR LAB ONLY")
    folder <- file.path("NYU PACT Patient Data", "Workbook")
    runyr <- stringr::str_split_fixed(inputSheet,"-",3)[,2]
    inputFi <- file.path(drive, folder, paste0("20", runyr),inputSheet,paste0(inputSheet,".xlsm"))
    if(file.exists(inputFi)) {
        return(inputFi)
    } else{
        inputFi <- file.path(drive, folder, paste0("20", runyr), inputSheet)
        message(dsh, "The PACT run worksheet was not found:\n", inputFi, dsh)
        allFi <- list.files(path=inputFi, pattern="*.xlsm")
        allFi <- allFi[!grepl( "~$", allFi, fixed = T)]
        allFi <- allFi[!grepl( "export", allFi, fixed = F,ignore.case=T)]
        allFi <- allFi[grepl( "Book", allFi, fixed = F,ignore.case=T)]
        if(length(allFi)>1){allFi <- allFi[1]}
        inputFi <- file.path(inputFi,allFi)
        message(dsh,"Using this workbook instead:\n", inputFi, dsh)
        return(inputFi)
        }
}

parseWorksheet <- function(inputFi){
    shNames <- readxl::excel_sheets(inputFi)
    sheet2Read <- "PhilipsExport"
    sh <- which(grepl(sheet2Read,shNames,ignore.case=T))[1]
    pactShCol <- c("Tumor Specimen ID","Tumor DNA/RNA Number", "MRN", "Test Number")
    vals2find <-  as.data.frame(readxl::read_excel(inputFi, sheet=shNames[sh], skip=3, col_types ="text")[,pactShCol])
    vals2find <- vals2find[!is.na(vals2find[,1]),]
    return(vals2find)
}

# Import csv Worksheet -----
getCaseValues <- function(inputSheet,readFlag){
    if(readFlag){
        vals2find <- utils::read.csv(inputSheet, skip=19)[,c(6,7,9)]
        vals2find <- vals2find[!grepl("H20|SERACARE|HAPMAP", vals2find[,2]),]
        return(vals2find)
    }else{
        inputFi <- getFilePath(inputSheet)
        vals2find <- parseWorksheet(inputFi)
        return(vals2find)
    }
}

genQuery <- function(dbCol,vals2find){
    currCol <- vals2find[, dbCol]
    q1 <- currCol[currCol != 0 & !is.na(currCol) & currCol != ""]
    return(unique(q1))
}

# Get Methylation and Molecular Samples list ----
queryCases <- function(vals2find, db) {
    i=NULL
    queryList <- foreach::foreach(i=1:ncol(vals2find), .combine="c") %do% {genQuery(i,vals2find)}
    methQuery <- searchDb(queryList, db)
    return(unique(methQuery))
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
        
        output$Test_Number[i] <- theVal
    }
    output <- addOutputLinks(output)
    return(output)
}

createXlFile <- function(runId,output){
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, runId)
    openxlsx::writeData(wb, sheet=runId, x= output)
    for (fiLn in 1:length(output$'Report Link')) {
        if(output$'Report Link'[fiLn]!=''){
            x <- c(output$'Report Link'[fiLn])
            names(x) <- paste0(output$record_id[fiLn],".html")
            class(x) <- "hyperlink"
            colNum <- which(colnames(output)=="Report Link")
            openxlsx::writeData(wb, sheet = runId, x = x, startCol = colNum, startRow=fiLn+1)
        }
    }
    outFi <-file.path("~","Desktop", paste0(runId,"_MethylMatch.xlsx"))
    openxlsx::saveWorkbook(wb, outFi, overwrite = T)
    return(outFi)
}

postData <- function(rcon, record){
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
    res <- RCurl::postForm(rcon$url,token=rcon$token,content='record',format='json',type='flat',
                           data=datarecord, returnContent = 'ids', returnFormat = 'csv')
    cat(res,sep="\n")
}

emailFile <- function(runId, outFi, rcon){
    record = data.frame(record_id = runId, run_number = runId)
    postData(rcon, record)
    suppressMessages(redcapAPI::importFiles(rcon,file= outFi, runId, field="other_file", repeat_instance=1))
    record$comments <- "pact_sample_list_email"
    postData(rcon, record)
    message(dsh,"\nEmail Notification Created",dsh)
}

# Searching REDCap Worksheets for MRN Match -------------------------------------

loadPacks()
rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", token)
vals2find <- getCaseValues(inputSheet, readFlag)
db <- grabAllRecords(flds, rcon)
output <- queryCases(vals2find, db)
theExcel <- getFilePath(inputSheet)

runId <- ifelse(readFlag == T, paste0(head(read.csv(inputSheet))[3, 2]),
                paste0(head(suppressMessages(
                    readxl::read_excel(theExcel, sheet = 7)
                ))[3, 2]))
output <- modifyOutput(output,vals2find)
outFi <- createXlFile(runId,output)
emailFile(runId, outFi, rcon)

# CNV PNG Creation -------------------------------------

rds <- output$record_id[output$report_complete=="YES"]
assign("rds", rds)

supM <- function(sobj){return(suppressMessages(suppressWarnings(sobj)))}

# FUN: Sets your directory and sources the helper functions
sourceFuns2 <- function(workingPath = NULL) {
    mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/"
    script.list <- c("SetRunParams.R","CopyInputs.R","PACT_scripts/generateCNV.R")
    if(!require("redcapAPI")){install.packages("redcapAPI", dependencies = T, type="both",ask=F)
        library("redcapAPI")}
    if (is.null(workingPath)) {workingPath = getwd()}
    scripts <- paste0(mainHub, script.list)
    invisible(lapply(scripts, function(i){devtools::source_url(i)}))
    gb$setDirectory(workingPath)
    return(gb$defineParams())
}

message("\nStarting CNV PNG Creation\n")

sourceFuns2()

gen.cnv.png <- function(RGsetEpic, sampleName) {
    RGset=RGsetEpic
    imgName <- paste(sampleName, "cnv.png", sep="_")
    savePath <- file.path("~","Desktop")
    fn=file.path(savePath,imgName)
    tempPathFi <- file.path(savePath,"temp.html")
    if(file.exists(fn)){
        message("File already exists, skipping:", fn)
    }else{
        Mset <- mnp.v11b6::MNPpreprocessIllumina(RGsetEpic)
        Mset@annotation=c(array="IlluminaHumanMethylationEPIC", annotation="ilm10b4.hg19")
        FFPE <- mnp.v11b6::MNPgetFFPE(RGsetEpic)
        Mset_ba <- mnp.v11b6::MNPbatchadjust(Mset, FFPE)
        detP <- minfi::detectionP(RGsetEpic)
        bs <- minfi::getBeta(Mset)
        sexEstimate <-as.data.frame(signif(sest::get.proportion_table(bs, detP), digits = 2))
        yest <- as.double(sexEstimate$`p.Y:(-18,-5]`) >= 0.75
        yest1 <- as.double(sexEstimate$`Y:(0,0.1]`) >= 0.12
        sex <- ifelse((yest == TRUE && yest1 == TRUE), "male", "female")
        message("Generating ", sampleName, " cnv plot...")
        xx <- mnp.v11b6::MNPcnv(Mset,sex = sex,main = sampleID)
        thePlot<-supM(mnp.v11b6::MNPcnvggplotly(xx, getTables = F))
        p<-supM(plotly::ggplotly(thePlot))
        supM(htmlwidgets::saveWidget(widget=plotly::as.widget(p), file=tempPathFi))
        supM(webshot2::webshot(url=tempPathFi, file = fn, cliprect = "viewport", vwidth = 1152, vheight = 672))
        #dev.off()
        message("File saved:\n",imgName,"\n")
    }
}


grabRGset <- function(runPath, sentrix){
    barcode = stringr::str_split_fixed(sentrix, "_",2)[1]
    RGsetEpic <- minfi::read.metharray(file.path(runPath, sentrix), verbose = T, force = T)
    aEpic=c(array="IlluminaHumanMethylationEPIC", annotation="ilm10b4.hg19")
    a450k=c(array="IlluminaHumanMethylation450k", annotation="ilmn12.hg19")
    if (barcode >= as.numeric("204220033000")) {RGsetEpic@annotation=aEpic}
    if (RGsetEpic@annotation['array']=="IlluminaHumanMethylation450k") {RGsetEpic@annotation=a450k}
    return(RGsetEpic)
}

copyOutputPng <- function(){
    unlink("~/Desktop/temp.html")
    the.cnvs <- dir(path="~/Desktop",pattern="_cnv.png", full.names=T)
    outFolder <- "/Volumes/molecular/Molecular/MethylationClassifier/CNV_PNG"
    savePath <- file.path(outFolder, basename(the.cnvs))
    message("Copying png files to Molecular folder:\n")
    message(outFolder)
    try(fs::file_copy(path=the.cnvs,new_path=savePath),silent = T)
    if(any(!file.exists(savePath))){
        message("The following failed to copy from the desktop:\n")
        print(basename(savePath[!file.exists(savePath)]))
    }
    # while (!is.null(dev.list()))  dev.off()
}


library("conumee");library("sest");library("mnp.v11b6")
require("plotly");require("htmlwidgets")
ApiToken = gb$token
gb$save.png.files(gb$rds, token)
