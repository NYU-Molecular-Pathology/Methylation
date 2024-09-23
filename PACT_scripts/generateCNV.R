gb <- globalenv(); assign("gb", gb)
supM <- function(sobj){return(suppressMessages(suppressWarnings(sobj)))}
instLin <- "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/LoadInstall_new.R"

require("minfi")
require("sest")
library(data.table)
library(dplyr)
library(tibble)

current_version <- utils::packageVersion("redcapAPI")
required_version <- "2.7.4"

if (is.na(current_version) || utils::compareVersion(as.character(current_version), required_version) < 0) {
    install.packages("redcapAPI", ask=FALSE, update=TRUE, dependencies=TRUE)
}

#  Copy idats and Worksheets creation
writeFromRedcap <- function(df, samplesheet_ID, bn = NULL) {
    if (is.null(bn)) {bn = file.path(getwd(), df$barcode_and_row_column)}
    message(crayon::bgCyan("~~~Writing from redcap samplesheet.csv using dataframe:"))

    df<- df[!is.na(df[, "barcode_and_row_column"]),]
    df<- df[!is.null(df[, "barcode_and_row_column"]),]
	    print(df)
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
    samplesheet_csv <- samplesheet_csv[!is.na(samplesheet_csv$SentrixID_Pos),]
    print(samplesheet_csv)
    write.csv(samplesheet_csv, file = "samplesheet.csv", quote = F,row.names = F)
}

search.redcap <- function(rd_numbers, ApiToken=NULL) {
    if(is.null(ApiToken)){message("You must provide an ApiToken!")}
    rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", ApiToken)
    flds = c("record_id","b_number","primary_tech","second_tech","run_number","barcode_and_row_column","accession_number","arrived")
    result <- redcapAPI::exportRecordsTyped(
	    rcon, records = rd_numbers, fields = flds, dag = F, factors = F, form_complete_auto = F, format = 'csv')
    result <- as.data.frame(result)
    return(result)
}

warnMount <- function(idat.dir){
    message("\n",crayon::bgRed("Directory not found, you may need this idat folder location mounted:"),"\n", idat.dir,"\n")
}

get.idats2<-function(csvNam = "samplesheet.csv"){
    rsch.idat <- gb$rsch.idat;clin.idat <- gb$clin.idat
    if(!dir.exists(rsch.idat)){warnMount(rsch.idat)}; if(!dir.exists(clin.idat)){warnMount(clin.idat)}
    stopifnot(dir.exists(rsch.idat)|dir.exists(clin.idat))
    if (file.exists(csvNam)) {
        allFi <- gb$getAllFiles(idatDir = c(rsch.idat, clin.idat), csvNam = csvNam)
        allFi = allFi[file.exists(allFi)]
        if (length(allFi) > 0) {
            message("Files found: "); print(allFi)
            cur.idat <- dir(pattern = "*.idat$")
            bcds <- paste0(basename(allFi))
            if (all(bcds %in% cur.idat)) {message(".idat files already copied")}
            if (!all(bcds %in% cur.idat)) {gb$copyBaseIdats(allFi[!(bcds %in% cur.idat)])}
        } else {message("No .idat files found! Check worksheet and input folder path")}
    } else {message(paste("Cannot find your sheet named:", csvNam))}
}

# FUN: Copies .idat files to your directory and saves samplesheet.csv
get.rd.info <- function(rd_numbers=NULL, token=NULL, sh_name=NULL){
    if (is.null(rd_numbers)){
	    message("No RD-numbers found, Input RD-numbers using get.rd.info(rd_numbers)")
	    return(NULL)
    }
	print(rd_numbers)
    if (is.null(sh_name)) {sh_name = "samplesheet.csv"}
    result <- gb$search.redcap(rd_numbers, token)
    samplesheet_ID = as.data.frame(stringr::str_split_fixed(result[, "barcode_and_row_column"], "_", 2))
    writeFromRedcap(result, samplesheet_ID) # writes API export as minfi dataframe sheet
    gb$get.idats2(csvNam = sh_name)  # copies idat files from return to current directory
    return(result)
}

grabSexEst <- function(bs, detP){
    sexEstimate <-as.data.frame(signif(sest::get.proportion_table(bs, detP), digits = 2))
    yest <- as.double(sexEstimate$`p.Y:(-18,-5]`) >= 0.75
    yest1 <- as.double(sexEstimate$`Y:(0,0.1]`) >= 0.12
    sex <- ifelse((yest == TRUE && yest1 == TRUE), "male", "female")
    return(sex)
}


getHtmlCN <- function(fn){
    newFi = paste0(substr(fn, 1, nchar(fn) - 3), "html")
    supM(htmlwidgets::saveWidget(plotly::as_widget(p), newFi))
    message("\nSaving html cnv ", newFi, "\n")
}

getPngFile <- function(fn, tempPathFi){
    message("\nSaving WebShot from temp.html...\n")
    webshot2::webshot(url = tempPathFi, file = fn,cliprect = "viewport", delay = 3.5,
        vwidth = 1560,
        vheight = 896
    )
}

getCnWebshot <- function(xx, fn, asPNG=T) {
    tempPathFi <- file.path("~", "Desktop", "temp.html")
    thePlot <- supM(mnp.v11b6::MNPcnvggplotly(xx, getTables = F))
    p <- supM(plotly::ggplotly(thePlot))
    supM(htmlwidgets::saveWidget(plotly::as_widget(p), tempPathFi))
    if (asPNG == F) {getHtmlCN(fn)} else{getPngFile(fn, tempPathFi)}
    try(file.remove(tempPathFi), silent = T)
    invisible(gc()) # garbage cleanup
}

calculateCnv <- function(RGsetEpic, sampleName) {
    Mset <- mnp.v11b6::MNPpreprocessIllumina(RGsetEpic)
    Mset@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19")
    FFPE <- mnp.v11b6::MNPgetFFPE(RGsetEpic)
    Mset_ba <- mnp.v11b6::MNPbatchadjust(Mset, FFPE)
    bs <- minfi::getBeta(Mset_ba) #Mset
    sex <- grabSexEst(bs, minfi::detectionP(RGsetEpic))
    xx <- mnp.v11b6::MNPcnv(Mset, sex = sex, main = sampleName)
    return(xx)
}

# gets the mset and converts mnp cnv analysis obj to png
gen.cnv.png2 <- function(RGsetEpic, sampleName, asPNG=T) {
    imgName <- paste(sampleName, "cnv.png", sep = "_")
    fn = file.path("~", "Desktop", imgName)
    if (file.exists(fn)) {
        message("\nFile already exists, skipping:", fn, "\n")
    } else{
        cnvMsg <- paste("Generating", sampleName, "cnv plot...")
        message("\n", crayon::bgYellow(cnvMsg), "\n")
        xx <- calculateCnv(RGsetEpic, sampleName)
        getCnWebshot(xx, fn,asPNG)
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

copyOutputPng <-
    function(outFolder = "/Volumes/molecular/Molecular/MethylationClassifier/CNV_PNG") {
        try(unlink("~/Desktop/temp.html"), silent = T)
        the.cnvs <- dir("~/Desktop", "_cnv.png", full.names = T) %>% file.info() %>%
            rownames_to_column() %>% filter(as.Date(ctime) == Sys.Date()) %>% pull(rowname)
        if (length(the.cnvs) > 0) {
            savePath <- file.path(outFolder, basename(the.cnvs))
            message("\nCopying png files to Molecular folder: ", outFolder, "\n")
            cat(the.cnvs, sep = "\n")
            try(fs::file_copy(path = the.cnvs, new_path = savePath), T)
            if (any(!file.exists(savePath))) {
                message("The following failed to copy from the desktop:\n")
                print(basename(savePath[!file.exists(savePath)]))
            }
        } else{message("No CNV files found on Desktop to copy")}
        # while (!is.null(dev.list()))  dev.off()
}

save.png.files <- function(rds, token, asPNG=T){
    get.rd.info(rd_numbers=rds, token=token,sh_name=NULL) # input your RD-numbers here rd_numbers = c("RD-21-21")
    myDt <- gb$search.redcap(rd_numbers = rds, token)
    mySentrix <- myDt[myDt[,2] %like% "_R0", ]
    if(nrow(mySentrix)>0){
    for (sam in rownames(mySentrix)) {
        sampleName<-mySentrix[sam,1]
        fn = file.path("~", "Desktop", paste0(sampleName, "_cnv.png"))
        if (file.exists(fn)) {
            message("\nFile already exists, skipping:", fn, "\n")
        } else{
            sentrix <- mySentrix[sam,2]
            message("\nGetting RGset for ",sentrix,"\n")
            RGsetEpic <- grabRGset(getwd(),sentrix)
            tryCatch(
                expr= {gen.cnv.png2(RGsetEpic, sampleName, asPNG)},
                error= function(e){
                    erTxt <- paste0("An error occured with ", mySentrix[sam,1]," png creation:")
                    message(crayon::bgRed(erTxt),"\n",e)
                    message(crayon::bgGreen("Trying next sample"))
                }
            )
        }
    }}else{
        message("The RD-number(s) provided do not have idat files listed in REDCap:/n")
        print(rds)
        }
    while (!is.null(dev.list()))  dev.off()
}
