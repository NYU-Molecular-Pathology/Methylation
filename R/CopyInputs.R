#!/usr/bin/env Rscript
gb <- globalenv(); assign("gb", gb)
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"
cpInLnk = "https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/CopyInputs.R"
defaultDir = "/Volumes/CBioinformatics/Methylation/Clinical_Runs"

msgFunName <- function(pthLnk, funNam){
message("\nExecuting function: ", crayon::black$bgYellow(funNam), " from RScript in:\n", pthLnk,"\n")
}

msgParams <- function(...){
    message("\nParam passed: ", crayon::bgGreen(paste(deparse(substitute(...)), "=", ...)), "\n")
}

# Changes the working directory using the system CD command
setDirectory <- function(foldr) {
    msgFunName(cpInLnk,"setDirectory")
    bsDir = paste("cd", foldr)
    mm2 = crayon::white$bgRed("Location Not Found:", foldr)

    if (dir.exists(foldr)) {
        system(bsDir)
        setwd(foldr)
        assign("workDir", foldr)
        }else{warning(mm2)}
}

CreateRunDir <- function(newRun) {
    if (endsWith(newRun, "/")) {
        newRun <- substr(newRun, 1, nchar(newRun) - 1)
    }
    message(crayon::bgGreen("New Run Path:"), "\n", newRun)
    if (!dir.exists(newRun)) {
        dir.create(newRun, recursive = T)
        Sys.chmod(newRun, "0777", use_umask = TRUE)
    }
    if(Sys.info()[['sysname']]=="Darwin"){
    cmd <- paste("chmod -R 770", newRun)
    system(cmd)
    cmd1 <-
        paste0(
            'chmod -R +ai ',
            "'group:NYUMC\\shared-rsc-CBioinfo-Meth-Clinicalruns",
            " allow list,add_file,search,add_subdirectory,delete_child,readattr,writeattr,readextattr,writeextattr,readsecurity,file_inherit,directory_inherit'",
            " '",
            newRun,
            "'","/*"
        )
    try(system(cmd1),silent=T)
    setDirectory(newRun)
    cmd2 <-
        paste0(
            'chgrp ',
            file.path(
                "'NYUMC",
                "shared-rsc-CBioinfo-Meth-Clinicalruns'",
                fsep = "\\"
            ),
            " '",
            newRun,
            "'"
        )
    try(system(cmd2),silent=T)

    cmd3 <-
        paste0(
            'chmod -R +ai ',
            file.path(
                "'group:NYUMC",
                "shared-rsc-CBioinfo-Meth-Clinicalruns",
                fsep = "\\"
            ),
            " allow list,add_file,search,delete,add_subdirectory,readattr,writeattr,readextattr,writeextattr,readsecurity,file_inherit,directory_inherit'",
            " '",
            newRun,
            "'","/*"
        )
    try(system(cmd3),silent=T)
    }

    setDirectory(newRun)
    return(newRun)
}

# Sets the methylation run directory named by the new run name
setRunDir <- function(runID=NULL, workFolder=NULL){
    msgFunName(cpInLnk, "setRunDir")
    msgParams(runID); msgParams(workFolder)
    workFolder <- ifelse(is.null(workFolder),gb$defaultDir, workFolder)
    if(is.null(runID)){
        runID <- paste0(basename(getwd()))
    }
    newRun <- file.path(workFolder, runID)
    if(grepl("TEST",runID)){
        if(dir.exists(newRun)){unlink(newRun, T, T)}
        dir.create(newRun)
        try(unlink(file.path(fs::path_home(),"Desktop",runID), T, T),silent = T)
    }
    methylPath <- CreateRunDir(newRun)
    return(methylPath)
}

# FUN: Returns a list of idat files given an idat drive location -
getAllFiles <- function(idatDir, csvNam=NULL) {
    msgFunName(cpInLnk, "getAllFiles")
    msgParams(idatDir);msgParams(csvNam)
    if (!is.null(csvNam)) {
        ssheet = read.csv(csvNam, strip.white = T)
        barcode = as.vector(ssheet$Sentrix_ID)
    } else {
        ssheet = read.csv("samplesheet.csv", strip.white = T)
        barcode = as.vector(ssheet$Sentrix_ID)
    }
    allFi = NULL
    senP = ssheet$SentrixID_Pos
    for (i in 1:length(idatDir)){
        dirNames = file.path(idatDir[i], barcode)
        green.files <- file.path(dirNames, paste0(senP, "_Grn.idat"))
        red.files <- file.path(dirNames, paste0(senP, "_Red.idat"))
        allFi <- c(allFi, green.files, red.files)
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
    if(!dir.exists(idat.dir)){
        message(crayon::bgRed("Directory not found, ensure the idat path is accessible:"), "\n", idat.dir)
        stopifnot(dir.exists(idat.dir))
        }
}

DataFrameMessage <- function(dat){
    return(message(paste0(capture.output(as.data.frame(dat)), collapse="\n")))
}


# FUN: Returns a list of idat files that exist on Molecular and Snuderl lab drives -
get.idats <- function(csvNam = "samplesheet.csv", runDir=NULL){
    msgFunName(cpInLnk, "get.idats")
    rsch.idat <- gb$rsch.idat
    clin.idat <- gb$clin.idat
    extr.idat <- file.path(gb$rsch.idat,"External")
    warnMount(rsch.idat)
    warnMount(clin.idat)
    if (is.null(runDir)){ runDir <- getwd()}
    if (file.exists(csvNam)) {
        allFi <- getAllFiles(idatDir = c(rsch.idat, clin.idat), csvNam = csvNam)
        allFi = allFi[file.exists(allFi)]
        ssheet = read.csv(csvNam, strip.white = T)
        if (length(allFi) > 0) {
            message("Files found: ")
            DataFrameMessage(allFi)
            message("Checking if idats exist in run directory...")
            cur.idat <- basename(dir(path=runDir, pattern = "*.idat$", recursive = F))
            bcds <- paste0(basename(allFi))
            if(nrow(ssheet) * 2 != length(allFi)){
                message("Still missing some idats! Checking External Folder: ", extr.idat)
                barcode = as.vector(ssheet$SentrixID_Pos)
                #allToFind <- unique(basename(getAllFiles(idatDir = c(rsch.idat, clin.idat), csvNam = csvNam)))
                toSearch <- paste(barcode, collapse = "|")
                stillMissing <- stringr::str_detect(allFi, pattern=toSearch, negate=T)
                otherIdat <- dir(extr.idat, pattern = ".idat", full.names = T, recursive = T)
                toBeFound <- barcode[stillMissing]
                foundIdat <- stringr::str_detect(otherIdat, pattern=toBeFound)
                if(any(foundIdat)==T){
                    message("Found extra idats in External folder: ", extr.idat)

                }else{
                    message("Still missing idat files not in External folder: ", extr.idat)

                }
            }

            if (all(bcds %in% cur.idat)) {
                message(".idat files already copied")
            }else{
                copyBaseIdats(allFi[!(bcds %in% cur.idat)])
            }
        } else {
            warning(crayon::bgRed("No .idat files found!"))
            message("Check worksheet for barcode if found in the folder path:")
            message(rsch.idat, "\nor\n", clin.idat)
            stopifnot(length(allFi) > 0)
        }
    } else {
        message("Cannot find your sheet named:", csvNam)
        stopifnot(file.exists(csvNam))
    }
}


# FUN: Copies samplesheet to Desktop folder
moveSampleSheet <- function(methDir, runID = NULL) {
    msgFunName(cpInLnk, "moveSampleSheet")
    if (is.null(runID)) {
        runID = paste0(basename(getwd()))
        message("Setting runID=", runID)
    }
    deskDir = file.path(fs::path_home(),"Desktop", runID)
    if (!dir.exists(deskDir)) {
        dir.create(deskDir)
    }
    if(is.null(methDir)){methDir <- "/Volumes/CBioinformatics/Methylation/Clinical_Runs"}
    currDir = file.path(methDir, runID)
    outFile = paste0(runID, "_samplesheet.csv")
    message("Copying from currDir:\n",currDir)
    message("To:\n",file.path(deskDir, "samplesheet.csv") )
    fs::file_copy(
        path = file.path(currDir, "samplesheet.csv"),
        new_path = file.path(deskDir, "samplesheet.csv"),
        overwrite = T
    )
    file.rename(from = file.path(deskDir, "samplesheet.csv"),
                to = file.path(deskDir, outFile))
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
        rcon, records = rd_numbers, fields = flds, dag = F, factors = F,
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

CopyFromDmn <- function(runID, runYear){
    kerberos <- Sys.info()[["user"]]
    cmdStart <- "rsync --protect-args --size-only -avzPe ssh"
    sshUser <- paste0(kerberos, "@dmn-0002:")
    runSheet <- paste0(runID,".xlsm")
    startPath <- file.path(
        "/mnt", kerberos, "molecular/MOLECULAR LAB ONLY/NYU-METHYLATION/WORKSHEETS", runYear, runSheet)
    destPath <- file.path(
        "/gpfs/home", kerberos, "molecpathlab/production/Methylation/Clinical_Runs", runID)
    cmdSync <- paste0(cmdStart, " '", sshUser, startPath, "' ", destPath)
    system(cmdSync)
}
