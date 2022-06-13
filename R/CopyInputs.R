#!/usr/bin/env Rscript
gb <- globalenv(); assign("gb", gb)
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
cpInLnk = "https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/CopyInputs.R"
defaultDir = "/Volumes/CBioinformatics/Methylation/Clinical_Runs"

msgFunName <- function(pthLnk, funNam){
message("\nExecuting function: ", crayon::black$bgYellow(funNam), " from RScript in:\n", pthLnk,"\n")
}

msgParams <- function(...){message("\nParams passed: ", crayon::bgGreen(paste(..., sep = ",")), "\n")}

# Sets the methylation run directory named by the new run name
setRunDir <- function(runID=NULL, workFolder=NULL){
    msgFunName(cpInLnk, "setRunDir")
    msgParams("runID", "workFolder")
    msgParams(runID, workFolder)

    if(!is.null(workFolder)){gb$methDir <- workFolder}
    if(is.null(gb$methDir)){gb$methDir<- gb$defaultDir}
    workFolder <- gb$ckNull(workFolder, gb$methDir, deparse(substitute(workFolder,env=gb)))
    if(is.null(workFolder)){
        workFolder <- gb$defaultDir
        assign("workFolder", gb$defaultDir)
    }
    newRun <- file.path(workFolder, runID)
    if(stringr::str_detect(newRun, "NULL")==T){
         newRun <- file.path(gb$methDir, runID)
    }

    assign("newRunPath", newRun)
    if(grepl("TEST",runID)){
        if(dir.exists(newRun)){unlink(newRun, T, T)}
        dir.create(newRun)
        try(unlink(file.path("~/Desktop",runID), T, T),silent = T)
    }

    message(newRun,"\n", crayon::bgGreen("newRunPath: "),newRunPath)

    if(!dir.exists(newRun)){
        #base::dir.create(newRun, recursive=T)
        cmdDir = paste('mkdir -p -m 777', newRun)
        system(cmdDir)
        cat("Creating folder: ",newRun)
        gb$setDirectory(newRun)
    } else {gb$setDirectory(newRun)}
    return(newRun)
}

# FUN: Returns a list of idat files given an idat drive location -
getAllFiles <- function(idatDir, csvNam=NULL) {
    msgFunName(cpInLnk, "getAllFiles")
    msgParams("idatDir","csvNam")

    if(!is.null(csvNam)){ssheet=read.csv(csvNam,strip.white=T)
    barcode=as.vector(ssheet$Sentrix_ID)} else {
        ssheet=read.csv(csvNam,strip.white=T)
        barcode=as.vector(ssheet$Sentrix_ID)}
    allFi=NULL
    senP=ssheet$SentrixID_Pos
    for (i in 1:length(idatDir)){
        dirNames = file.path(idatDir[i], barcode)
        green.files <- file.path(dirNames, paste0(senP, "_Grn.idat"))
        red.files <- file.path(dirNames, paste0(senP, "_Red.idat"))
        allFi <- c(allFi, green.files,red.files)
    }
    return(allFi)
}

# FUN: Copies .idat files to your current directory using sample sheet
copyBaseIdats <- function(allFi) {
    msgFunName(cpInLnk, "copyBaseIdats")
    cat(crayon::white$bgCyan("Copying idats to current directory..."),"\n")
    fs::file_copy(allFi, file.path(getwd()),overwrite=T)
    idcs = basename(allFi)
    idatsCopied <- idcs[idcs != ""]
    success = file.exists(idatsCopied)
    all(success)
    message(".idat files that failed to copy:")
    if (all(success)) {cat("none","\n")} else {print(idatsCopied[!success])}
}

# Helper FUN called during copying idats to notify if a network mount is not found
warnMount <- function(idat.dir){
    msgFunName(cpInLnk, "warnMount")
    cat(crayon::bgRed("Directory not found, ensure the idat folder location is accessible:"),idat.dir,sep="\n")
}

# FUN: Returns a list of idat files that exist on Molecular and Snuderl lab drives -
get.idats <-function(csvNam = "samplesheet.csv"){
    msgFunName(cpInLnk, "get.idats")
    rsch.idat <- gb$rsch.idat;clin.idat <- gb$clin.idat
    if(!dir.exists(rsch.idat)){warnMount(rsch.idat)}; if(!dir.exists(clin.idat)){warnMount(clin.idat)}
    stopifnot(dir.exists(rsch.idat)|dir.exists(clin.idat))
    if (file.exists(csvNam)) {
        allFi <- getAllFiles(idatDir = c(rsch.idat, clin.idat), csvNam = csvNam)
        allFi = allFi[file.exists(allFi)]
        if (length(allFi) > 0) {
            message("Files found: "); print(allFi)
            cur.idat <- dir(pattern = "*.idat$")
            bcds <- paste0(basename(allFi))
            if (all(bcds %in% cur.idat)) {message(".idat files already copied")}
            if (!all(bcds %in% cur.idat)) {copyBaseIdats(allFi[!(bcds %in% cur.idat)])} #length(cur.idat) < length(allFi)
        } else {cat(crayon::bgRed("No .idat files found! Check worksheet for barcode if found in the folder path:\n"))
               cat(crayon::bgRed("/Volumes/molecular/Molecular/iScan/"))
                filesFound = F
                stopifnot(filesFound==T)
               }
    } else {message(paste("Cannot find your sheet named:", csvNam))
           filesFound = F
                stopifnot(filesFound==T)
           }
}

# FUN: Copies samplesheet to Desktop folder
moveSampleSheet <- function(methDir, runID=NULL) {
    msgFunName(cpInLnk, "moveSampleSheet")
    if (is.null(runID)){runID=paste0(basename(getwd()))}
    baseFolder=paste0("~/Desktop/",runID,"/")
    if(!dir.exists(baseFolder)){dir.create(baseFolder)}
    currDir=paste0(methDir, "/",runID)
    endDir = paste0(baseFolder, runID, "_samplesheet.csv")
    fs::file_copy(path=paste0(currDir,"/samplesheet.csv"),new_path=baseFolder,overwrite=T)
    file.rename(from=paste0(baseFolder,"samplesheet.csv"), to=endDir)
}

#  Copy idats and Worksheets creation
writeFromRedcap <- function(df, samplesheet_ID, bn = NULL) {
     msgFunName(cpInLnk, "writeFromRedcap")
    if (is.null(bn)) {bn = file.path(getwd(), df$barcode_and_row_column)}
    message("~~~Writing from redcap samplesheet.csv:")
    names(df)
    samplesheet_csv = data.frame(
        Sample_Name = df[, "record_id"],
        DNA_Number = df[,"b_number"],
        Sentrix_ID = samplesheet_ID[, 1],
        Sentrix_Position = samplesheet_ID[, 2],
        SentrixID_Pos = df[, "barcode_and_row_column"],
        Basename = paste0(bn),
        RunID = df$run_number,
        MP_num = df$accession_number,
        tech = df$primary_tech,
        tech2 = df$second_tech,
        Date = df$arrived
    )
    print(samplesheet_csv)
    write.csv(
        samplesheet_csv,
        file = "samplesheet.csv",
        quote = F,
        row.names = F
        )
}

#' FUN: Returns dataframe of redcap search using default worksheet header and fields
search.redcap <- function(rd_numbers, token=NULL, flds=NULL) {
    msgFunName(cpInLnk, "search.redcap")
    if(!require("redcapAPI")){install.packages("redcapAPI", dependencies = T, type="both",ask=F);library("redcapAPI")}
    if(is.null(token)){message("You must provide an ApiToken!")}
    stopifnot(!is.null(token))
    rcon <- redcapAPI::redcapConnection(gb$apiLink, token)
    if (is.null(flds)){
        flds = c("record_id","b_number","primary_tech","second_tech","run_number",
                 "barcode_and_row_column","accession_number","arrived")
        }
    result <- redcapAPI::exportRecords(
        rcon,records = rd_numbers,fields = flds,dag = F,factors = F,
        labels = F,dates = F, form_complete_auto = F,format = 'csv')
    return(as.data.frame(result))
}

# FUN: Copies .idat files to your directory and saves samplesheet.csv
get.rd.info <- function(rd_numbers=NULL, token=NULL, sh_name=NULL){
    msgFunName(cpInLnk, "get.rd.info")
    if (is.null(rd_numbers)){message("Input RD-numbers using get.rd.info(rd_numbers)")}
    if (is.null(sh_name)) {sh_name = "samplesheet.csv"}
    if(is.null(token)){message("You must provide an ApiToken!")};stopifnot(!is.null(token))
    result <- search.redcap(rd_numbers, token, NULL)
    samplesheet_ID = as.data.frame(stringr::str_split_fixed(result[,"barcode_and_row_column"],"_",2))
    writeFromRedcap(result, samplesheet_ID) # writes API export as minfi dataframe sheet
    get.idats(csvNam = sh_name)  # copies idat files from return to current directory

    return(result)
}
