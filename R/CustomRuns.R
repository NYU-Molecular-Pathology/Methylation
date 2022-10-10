gb <- globalenv(); assign("gb", gb)
reportMd <- "/Volumes/CBioinformatics/Methylation/report_v12.Rmd"
sarcRmdFile = "/Volumes/CBioinformatics/Methylation/SarcReport.Rmd"
GetLocalData <- function(rg){
    dat <- data.frame(
        sampleID = paste0(rg),
        bnumber = "NONE",
        senLi = paste0(rg),
        run_id = paste0(gb$runID),
        mp_number = "NONE",
        tech = "NONE",
        tech2 = "NONE",
        outFi = paste0(rg, ".html")
    )
    return(dat)
}

GetTargetData <- function(data) {
    runDt <- data.frame(
        sampleID = paste0(data[, "Sample_Name"]),
        bnumber = paste0(data[, "DNA_Number"]),
        senLi = paste0(data[, "SentrixID_Pos"]),
        run_id = paste0(data[, "RunID"]),
        mp_number = paste0(data[, "MP_num"]),
        tech = paste0("NONE"),
        tech2 = paste0("NONE"),
        outFi = paste0(data[, 1], ".html")
    )
    return(runDt)
}

KnitReportRmd <- function(dat, token, reportMd){
    message(paste0(capture.output(dat), collapse = "\n"))
    RGsetEpic <- suppressWarnings(gb$getRGset(getwd(), dat$senLi))
    RGset <- RGsetEpic[,1]
    rmarkdown::render(
        input = reportMd,
        output_format = "html_document",
        output_file = dat$outFi,
        output_dir = getwd(),
        quiet = F,
        output_options = c("self_contained = TRUE"),
        params = list(
            token = token,
            rundata = dat,
            RGsetEpic = RGsetEpic,
            knitDir = getwd()
        )
    )
}

loop_targets <- function(targets, reportMd="/Volumes/CBioinformatics/Methylation/report_v12.Rmd"){
    require(compiler)
    enableJIT(3)
    mainPage = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/Report-Scripts/"
    scripts <- c("ClassTables.R", "MLH1_Functions.R", "PipeLineU.R",
                 "RedcapOutput.R", "TsneFunctions.R", "cnvggplotly.R")
    scripts <- paste0(mainPage, scripts)
    lapply(scripts, function(i) {message("Sourcing: ", i);devtools::source_url(i)})
    gb$moveSampleSheet(gb$baseFolder, gb$runID)
    gb$TryLoadUniD()
    try(silent=T, unloadNamespace("sarc.v12b6"))
    gb$LoadMnpData(F)
    for(i in 1:nrow(targets)){
        message("Sample ", i, " of ", nrow(targets))
        sample = 1
        dat <- GetTargetData(data = targets[i,])
        message("dat:")
        message(paste(capture.output(dat), collapse="\n"))
        KnitReportRmd(dat, token = gb$ApiToken, reportMd)
    }
}

loop_local <- function(RGSet){
    reportMd <- "/Volumes/CBioinformatics/Methylation/report_v12.Rmd"
    for (rg in colnames(RGSet)) {
        thisSam <- RGSet[, rg]
        dat <- GetLocalData(rg)
        sentrix=dat$senLi
        RGsetEpic<-RGset<-thisSam
        rmarkdown::render(
            reportMd, "html_document", dat$outFi, getwd(), quiet = FALSE, output_options = c("self_contained = TRUE"),
            params = list(token = gb$ApiToken, rundata = dat, RGsetEpic = RGsetEpic, knitDir=getwd())
        )
    }
}

PromptInputCsv <- function(runID) {
    message('No idat files in current directory and no .xlsm file was found with the runID named "', runID,'"')
    message("Enter the full path name to a local csv file without quotes that has a list of RD-numbers and press return")
    message("The file should have no header and just list RD-numbers in Column A (Example: /Users/myName/Desktop/myFileList_rd_numbers.csv)")
    csvFilePath <- readline("Paste the full path to your csv file and hit return/Enter: ")
    csvFilePath <- as.character(csvFilePath)
    if(!file.exists(csvFilePath)){
        warning("The file ",csvFilePath," does not exist!")
        csvFilePath <- readline("Try pasting the full path again and hit return/Enter: ")
        csvFilePath <- as.character(csvFilePath)
    }
    rd_numbers <- read.csv(file = csvFilePath, header = F, colClasses="character")[,1]
    rd_numbers <- as.data.frame(rd_numbers)
    message("The following RD-numbers were entered:")
    print(rd_numbers$rd_numbers)
    return(rd_numbers$rd_numbers)
}

MakeLocalSampleSheet <- function(runID){
    idatScript <- "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/Research/pullRedcap_manual.R"
    rd_numbers <- PromptInputCsv(runID)
    stopifnot(length(rd_numbers)>=1 & length(rd_numbers)!=0 & stringr::str_detect(rd_numbers[1],"RD-"))
    message("Sourcing: ", idatScript)
    devtools::source_url(idatScript)
    gb$grabRDCopyIdat(rd_numbers, gb$token, copyIdats=T, outputFi="samplesheet.csv")
}

RunFromSamplesheet <- function(samSh="samplesheet.csv"){
    targets <- read.csv(file=file.path(getwd(),samSh), strip.white = T)
    targets[is.na(targets)] <- "NONE"
    write.csv(targets,samSh, quote = F,row.names = F)
    targets<- minfi::read.metharray.sheet(base=getwd(), pattern = samSh)
    dupedSamples <- anyDuplicated(targets$Basename)
    if(dupedSamples!=0){
        warning("Removing Samples with identical (duplicated) Sentrix IDs")
        message("Some basenames are duplicated and will be dropped:")
        tgBn <- targets$Basename
        print(targets[tgBn==tgBn[anyDuplicated(tgBn)],])
        targets <- targets[-dupedSamples,]
    }
    rownames(targets) <- targets$Sample_Name
    return(targets)
}

RunLocalIdats <- function(runID, token){
    if(!file.exists(file.path(getwd(), paste0(runID,".xlsm")))){
        idatFiles <- dir(path = getwd(), pattern = ".idat", full.names = T)
        if(!file.exists(file.path(getwd(), "samplesheet.csv")) & length(idatFiles)!=0){
            idatBase <- unique(substring(idatFiles, 1, nchar(idatFiles) - 9))
            RGSet <- minfi::read.metharray(basenames =idatBase, force=TRUE, verbose = T)
            loop_local(RGSet)
        }else{
            if(!file.exists(file.path(getwd(), "samplesheet.csv"))){
                MakeLocalSampleSheet(runID)
            }
            targets <- RunFromSamplesheet("samplesheet.csv")
            loop_targets(targets)
        }
    }else{
        gb$readSheetWrite(runID = runID)
        gb$moveSampleSheet(baseFolder, runID)
    }
}

# Work Directory Functions ----------------------------------------------------
CheckBaseDir <- function(baseFolder){
    if(is.null(baseFolder)){
        gb$baseDir <- gb$methDir <- gb$baseFolder <- "/Volumes/CBioinformatics/Methylation/Clinical_Runs"
    }else{gb$baseDir <- gb$methDir <- gb$baseFolder <- baseFolder}
    isDesktop <- stringr::str_detect(baseFolder, "Desktop")
    if(length(isDesktop)!=0)
    if(is.null(baseFolder) & isDesktop==T) {
        warning("Trying to run methylation from Desktop working directory is not allowed")
        message("Try setting baseFolder to '~/Documents/' instead")
        stopifnot(isDesktop == F)
    }
    return(gb$baseFolder)
}

# Sets the working folder directory
SetBaseFolder <- function(token, baseFolder, runID){
    baseFolder <- CheckBaseDir(baseFolder)
    methylPath <- gb$setRunDir(runID=gb$runID, workFolder = baseFolder)
    message("Working directory set to:\n", crayon::bgGreen(methylPath), "\n")
    gb$methDir <- gb$workFolder <- baseFolder
    gb$setVar("workFolder", baseFolder)
    gb$setVar("ApiToken", token) # assign the ApiToken & print params
    setwd(file.path(baseFolder, runID))
    return(baseFolder)
}

PrepareRun <- function(token, baseFolder=NULL, runID, runLocal=F, rdInput=F){
    if(runLocal==F){
        gb$checkMounts()
        gb$checkValidRun(runID)
        gb$SetBaseFolder(token, baseFolder, runID)
        gb$copyWorksheetFile(runID = runID) # copies the xlsm file
        gb$readSheetWrite(runID = runID) # reads xlsm and generates input .csv samplesheet
        gb$get.idats() # Copy idat files to current folder from molecular and snuderlabspace to cwd
        gb$moveSampleSheet(baseFolder, runID) #copies outputs temp to desktop for QC.Rmd
    } else{
        gb$SetBaseFolder(token, baseFolder, runID)
        if(rdInput==F){
            gb$RunLocalIdats(runID, token)
        }
    }
}



SarcomaReport <- function(RGset, sampleID, output_dir = getwd()) {
    output_file <- paste0(sampleID, "_sarc.html")
    RGset <- RGset[,1]
    rmarkdown::render(input= gb$sarcRmdFile, output_dir = output_dir, output_file = output_file)
}

PromptRDnumbers <- function(){
    message("Type your RD-numbers without quotes below")
    message("Separate with spaces & press Enter twice to exit")
    rd_numbers <- scan(what = " ")
    message("Values entered:\n", rd_numbers)
    rd_numbers <- rd_numbers[which(grepl(rd_numbers, pattern="^RD-"))]
    message("Valid RD-numbers entered:\n", paste0(capture.output(rd_numbers), collapse = '\n' ))
    if(length(rd_numbers)==0) {
        message("No valid RD-numbers entered! Try again.")
        PromptRDnumbers()
    }else(return(rd_numbers))
}

MakeBlankRun <- function(rd_numbers, token, outputFi="samplesheet.csv"){
    stopifnot(length(rd_numbers)>=1 & length(rd_numbers)!=0 & stringr::str_detect(rd_numbers[1],"RD-"))
    if(!exists("grabRDCopyIdat")){
        idatScript <- "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/Research/pullRedcap_manual.R"
        message("Sourcing: ", idatScript)
        devtools::source_url(idatScript)
    }
    gb$grabRDCopyIdat(rd_numbers, token, outputFi=outputFi)
    return(as.data.frame(read.csv(outputFi)))
}

MakeSarcomaReport <- function(worksheet = "samplesheet.csv", targets = NULL) {
    if (is.null(targets)) {
        if (file.exists(worksheet)) {
            targets <- as.data.frame(read.csv(worksheet))
        } else{
            rd_numbers <- PromptRDnumbers()
            targets <- MakeBlankRun(rd_numbers, gb$token, worksheet)
            }
    }
    for (samIdx in 1:nrow(targets)) {
        message("Running ", samIdx, " of ", nrow(targets))
        RGset <- suppressWarnings(minfi::read.metharray(targets$Basename[samIdx]))
        SarcomaReport(RGset, sampleID = targets[samIdx, 1])
    }
}

StartRun <- function(selectRDs = NULL, emailNotify = T, redcapUp = T) {
    gb$msgFunName(paste0(mainHub,"methylExpress.R"),"startRun")
    # Re-order sample report generation for priority
    if (!is.null(selectRDs)) {runOrder <- gb$reOrderRun(selectRDs) }else{runOrder <- NULL}
    gb$makeReports.v11b6(
        skipQC = F,            # Don't skip QC generation
        email = emailNotify,   # to email after Run complete
        cpReport = F,          # Flag to copy files to network drive
        selectSams = runOrder, # Prioritize specific RD-numbers
        redcapUp = redcapUp    # Flag to import files to REDCap
    )
}

