gb <- globalenv(); assign("gb", gb)

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

loop_targets <- function(targets){
    reportMd <- "/Users/serraj10/report_v12.Rmd"
    gb$moveSampleSheet(gb$baseFolder, gb$runID)
    for(i in 1:nrow(targets)){
        dat <- GetTargetData(data= targets[i,])
        print(dat)
        RGsetEpic <- suppressWarnings(gb$getRGset(getwd(), dat$senLi))
        rmarkdown::render(
            reportMd,
            "html_document",
            dat$outFi,
            getwd(),
            quiet = T,
            output_options = c("self_contained = TRUE"),
            params = list(token = gb$ApiToken, rundata = dat, RGsetEpic=RGsetEpic)
        )
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
            params = list(token = gb$ApiToken, rundata = dat, RGsetEpic=RGsetEpic)
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
    idatScript <-
        "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/Research/pullRedcap_manual.R"
    rd_numbers <- PromptInputCsv(runID)
    stopifnot(length(rd_numbers)>1 & stringr::str_detect(rd_numbers[1],"RD-"))
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
        if(length(idatFiles)>0){
            idatBase <- unique(substring(idatFiles, 1, nchar(idatFiles) - 9))
            RGSet<-minfi::read.metharray(basenames =idatBase, force=TRUE, verbose = T)
            loop_local(RGSet)
        }else{
            MakeLocalSampleSheet(runID)
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
    if(!is.null(baseFolder)){
        isDesktop <- stringr::str_detect(baseFolder, "Desktop")
    } else{isDesktop==F}
    if(isDesktop==T) {
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
}
