cbVol = "/Volumes/CBioinformatics"
moVol = "/Volumes/molecular"
rsVol = "/Volumes/snudem01labspace"
rsch.idat = paste0(file.path(rsVol,"idats"))
clin.idat = paste0(file.path(moVol, "MOLECULAR/iScan"))

getAllFiles <- function(idatDir, csvNam="samplesheet.csv") {
    ssheet <- read.csv(csvNam, strip.white=T)
    senP=ssheet[,grepl("Pos|Sentrix",colnames(ssheet))]
    barcode = stringr::str_split_fixed(senP,"_",2)[,1]
    i <- NULL
    allFi <-foreach::foreach (i = 1:length(idatDir), .combine="c") %do% {
        folder.names = file.path(idatDir[i], barcode)
        green.files <-file.path(folder.names, paste0(senP, "_Grn.idat"))
        red.files <- file.path(folder.names, paste0(senP, "_Red.idat"))
        c(green.files,red.files)
    }
    return(allFi)
}

# FUN: Copies .idat files to your current directory using sample sheet
copyBaseIdats <- function(allFi) {
    message("Files copying...")
    fs::file_copy(allFi, file.path(getwd()))
    idcs = basename(allFi)
    idatsCopied <- idcs[idcs != ""]
    success = file.exists(idatsCopied)
    message(".idat files that failed to copy:")
    if (length(success) > 0) {print(idatsCopied[!success])} else{cat("none\n")}
}

get.idats <-function(csvNam = "samplesheet.csv", clin.idat, rsch.idat) {
    if (!dir.exists(rsch.idat)) {message("Path not found:");cat(rsch.idat)}
    if (!dir.exists(clin.idat)) {message("Path not found:");cat(clin.idat)}
    stopifnot(dir.exists(rsch.idat)|dir.exists(clin.idat))
    if (file.exists(csvNam)) {
        allFi <- getAllFiles(idatDir = c(rsch.idat, clin.idat), csvNam = csvNam)
        allFi = allFi[file.exists(allFi)]
        if (length(allFi) > 0) {
            message("Files found: "); print(allFi)
            cur.idat <- list.files(pattern = "*.idat")
            bcds <- paste0(basename(allFi))
            if (all(bcds %in% cur.idat)) {message(".idat files already copied")}
            if (length(cur.idat) < length(allFi)) {copyBaseIdats(allFi[!(bcds %in% cur.idat)])}
        } else{message("No .idat files found! Check worksheet and input folder")}
    } else {message(paste("Cannot find your sheet named:", csvNam))}
}
