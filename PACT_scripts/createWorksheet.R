
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
    } else {message("Z-drive path is accessible")}
}

# Function to silently load library without conflict warnings
libLoad <- function(libName) {
    lib.opts <- list(package = libName, character.only = T, verbose = F, warn.conflicts = F)
    suppressWarnings(suppressPackageStartupMessages(do.call(library, c(lib.opts))))
    message(libName, " ...loaded successful")
}

# Functions to load or install missing required packages -----
loadPacks <- function(){
    if(suppressWarnings(!require("redcapAPI"))){
        params=list('nutterb/redcapAPI', dependencies=T, upgrade="always", type="source")
        do.call(devtools::install_github,c(params))
    } else{libLoad("redcapAPI")}
    pkgs = c("foreach","jsonlite","RCurl","readxl","stringr")
    rlis = getOption("repos");rlis["CRAN"] = "http://cran.us.r-project.org"
    options(repos = rlis)
    invisible(lapply(pkgs, function(pk){
        if(suppressWarnings(!require(pk, character.only=T))){
            install.packages(pk,dependencies=T, verbose=T, repos="http://cran.us.r-project.org", type="both")
            libLoad(pk)
        }else{libLoad(pk)}}))
}

# Returns Path to xlsx file -----
getExcelPath <- function(inputSheet, pathType=1){
    drive = file.path("", "Volumes", "molecular", "MOLECULAR LAB ONLY")
    folder = file.path("NYU PACT Patient Data", "Workbook")
    runyr <- stringr::str_split_fixed(inputSheet, "-", 3)[, 2]
    ending <- ifelse(pathType==1,".xlsm","_FinalExportedList.xlsx")
    return(file.path(drive, folder, paste0("20", runyr),inputSheet,
                     paste0(inputSheet, ending)))
}

# Removes and fixes newlines, commas, and blanks from samplesheet
sanitizeSheet <- function(sheetVals){
    mainSheet <- sheetVals[!is.na(sheetVals[,1]),]
    mainSheet$Tumor_Type <- gsub(" ", "-", mainSheet$Tumor_Type)
    for(i in 1:ncol(mainSheet)){
        mainSheet[,i] <- sapply(mainSheet[,i], function(x) { gsub("[\r\n]", "", x) })
        mainSheet[,i] <- sapply(mainSheet[,i], function(x) { gsub(",", "", x) })
        mainSheet[,i] <- sapply(mainSheet[,i], function(x) { gsub(" ", "", x) })
    }
    mainSheet$Paired_Normal[mainSheet$Paired_Normal==0|is.na(mainSheet$Paired_Normal)] <-""
    mainSheet$Tumor_Type[mainSheet$Tumor_Type==0] <- "NA"
    mainSheet[,1:16] <- sapply(mainSheet[,1:16], function(x){gsub("\\\\", "-", x) })
    controlNames <- "NTC_H20|SC_SERACARE|NC_HAPMAP"
    controlSamples <- grepl(pattern=controlNames, mainSheet$Sample_Name)
    if(table(controlSamples)[['TRUE']]!=3){
        warning("There are not 3 control samples, either NTC_H20, SC_SERACARE, or NC_HAPMAP is missing or added")
    }else{
        controlIndexes <- which(controlSamples==T)
        mainSheet[controlIndexes,'Paired_Normal'] <- ""
    }
    if(any(duplicated(mainSheet$Sample_ID))){
        warning("There are duplicated Sample_ID in the SampleSheet")
        print(mainSheet$Sample_ID[duplicated(mainSheet$Sample_ID)])
    }
    return(mainSheet)
}

# Reads the genders and outputs a .tsv file on desktop
getPhilipsGender <- function(mainSheet,inputFi, sh2){
    cnvSheet <- mainSheet[,1:15]
    runId <- cnvSheet[1,"Sample_Project"]
    philipVals <- as.data.frame(readxl::read_excel(inputFi,sheet = sh2,skip = 3,col_types = "text"))
    cnvSheet$Gender <- philipVals$Gender[match(cnvSheet$Test_Number, philipVals$`Test Number`)]
    cnvPath <- "/Volumes/molecular/Molecular/MethylationClassifier/CNV_PNG"
    write.table(cnvSheet,quote=F, sep='\t', file=file.path(cnvPath,paste0(runId,".tsv")),row.names=F)
}

# Parses xlsx file and writes as csv file -----
parseExcelFile <- function(inputFi){
    shNames <- readxl::excel_sheets(inputFi)
    sh <- which(grepl("SampleSheet", shNames, ignore.case = T))[1]
    sh2 <- which(grepl("Philips", shNames, ignore.case = T))[1]
    sheetHead <- as.data.frame(readxl::read_excel(
        inputFi,sheet = shNames[sh], na="", range="A1:B17", col_types = "text", col_names=F))
    sheetHead[is.na(sheetHead)] <- ""
    sheetHead <- rbind(sheetHead, c("",""), c("[Data]",""))
    sheetVals <-
        as.data.frame(readxl::read_excel(
            inputFi,
            sheet = shNames[sh],
            skip = 19,
            col_types = "text"
        ))
    mainSheet <- sanitizeSheet(sheetVals)
    getPhilipsGender(mainSheet,inputFi, sh2)
    outFile <- file.path("~","Desktop",paste(mainSheet[1,"Run_Number"],"SampleSheet.csv",sep="-"))
    write.table(sheetHead, sep=",", file=outFile, row.names=F, col.names=F, quote=F)
    controls <- which(grepl(pattern="NTC_H20|NC_HAPMAP|SC_SERACARE", mainSheet$Sample_Name))
    mainSheet$Paired_Normal[controls] <- ""
    suppressWarnings(write.table(mainSheet,sep=",", file=outFile, row.names=F, col.names=T, append=T,quote=F))
    return(c(runId=mainSheet[1,"Sample_Project"], outFile=outFile))
}

# Generate Email notification and attach csv file
emailNotify <- function(record,rcon){
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
    res<-RCurl::postForm(rcon$url, token=rcon$token, content='record',format='json',type='flat',
                         data = datarecord, returnContent = 'ids', returnFormat = 'csv')
    cat(res);message(dsh,"\nEmail Notification Created",dsh)
}

# Connect to REDCap and send email attachments of csv file ----
pushToRedcap <- function(runId,outFile,token){
    rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", token)
    record = data.frame(record_id = runId, pact_run_number = runId)
    datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
    res <- RCurl::postForm(rcon$url,token=rcon$token,content='record',format='json',
                           type='flat',data=datarecord, returnContent = 'count', returnFormat = 'csv')
    cat(res,sep="\n")
    redcapAPI::importFiles(rcon=rcon, file=outFile, record=runId, field="pact_csv_sheet", repeat_instance=1)
    record$pact_csv_email<-"pact_csv_email"
    emailNotify(record, rcon)
    #unlink(outFile)
}

# Gets dataframe and saves as CSV file -----
writeSampleSheet <- function(inputSheet, token){
    inputFi <- getExcelPath(inputSheet)
    if (file.exists(inputFi)) {
        outVals <- suppressMessages(parseExcelFile(inputFi))
        pushToRedcap(runId=outVals[[1]], outFile=outVals[[2]], token)
    } else {
        message(crayon::bgRed("The PACT run worksheet was not found:"),"\n", inputFi, dsh)
        inputFi <- gsub(paste0("/",basename(inputFi)),"",inputFi)
        potentialFi <- list.files(path=inputFi,full.names=T)
        message(crayon::bgRed("Checking the following files:"),"\n", dsh)
        if(length(potentialFi)>=1){
            print(potentialFi)
            wbFiles <- stringr::str_which(basename(potentialFi),pattern="xlsm")
            potentialFi <- potentialFi[wbFiles]
            potentialFi <- potentialFi[!stringr::str_detect(potentialFi,"\\$")]
        }
        if (file.exists(potentialFi[1])) {
            message("Now trying to read:\n",potentialFi[1])
            outVals <- suppressMessages(parseExcelFile(inputFi=potentialFi[1]))
            pushToRedcap(runId=outVals[[1]], outFile=outVals[[2]], token)
        }else{
            message(crayon::bgRed("The PACT run worksheet was not found:"),"\n", inputFi, dsh)
            stopifnot(file.exists(potentialFi[1]))
        }
    }
}

loadPacks()
checkMounts()
