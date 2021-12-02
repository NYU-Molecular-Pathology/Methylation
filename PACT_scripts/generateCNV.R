gb <- globalenv(); assign("gb", gb)
supM <- function(sobj){return(suppressMessages(suppressWarnings(sobj)))}

#  Copy idats and Worksheets creation
writeFromRedcap <- function(df, samplesheet_ID, bn = NULL) {
    if (is.null(bn)) {bn = file.path(getwd(), df$barcode_and_row_column)}
    message("~~~Writing from redcap samplesheet.csv:")
    names(df)
    df<- df[!is.na(df[, "barcode_and_row_column"]),]
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
    write.csv(
        samplesheet_csv,
        file = "samplesheet.csv",
        quote = F,
        row.names = F
    )
}

search.redcap <- function(rd_numbers, ApiToken=NULL) {
    if(is.null(ApiToken)){message("You must provide an ApiToken!")}
    rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", ApiToken)
    flds = c("record_id","b_number","primary_tech","second_tech","run_number","barcode_and_row_column","accession_number","arrived")
    result <- redcapAPI::exportRecords(rcon,records = rd_numbers,fields = flds,dag = F,factors = F,
                                       labels = F,dates = F, form_complete_auto = F,format = 'csv')
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
    if (is.null(rd_numbers)){message("Input RD-numbers using get.rd.info(rd_numbers)")}
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

getCnWebshot <- function(xx, fn){
    tempPathFi <- file.path("~","Desktop","temp.html")
    thePlot<-supM(mnp.v11b6::MNPcnvggplotly(xx, getTables = F))
    p<-plotly::ggplotly(thePlot)
    htmlwidgets::saveWidget(widget=plotly::as.widget(p), file=tempPathFi)
    message("\nSaving WebShot\n")
    webshot2::webshot(url=tempPathFi, file = fn, cliprect = "viewport", delay = 2.5, vwidth = 2340, vheight = 1344)
    dev.off()
    file.remove(tempPathFi)
    message("\nFile saved: ",fn,"\n")
}

gen.cnv.png2 <- function(RGsetEpic, sampleName) {
    RGset=RGsetEpic
    imgName <- paste(sampleName, "cnv.png", sep="_")
    fn = file.path("~","Desktop",imgName)
    if(file.exists(fn)){
        message("\nFile already exists, skipping:", fn,"\n")
    } else{
        Mset <- mnp.v11b6::MNPpreprocessIllumina(RGsetEpic)
        Mset@annotation=c(array="IlluminaHumanMethylationEPIC", annotation="ilm10b4.hg19")
	FFPE <- mnp.v11b6::MNPgetFFPE(RGset)
	Mset_ba <- mnp.v11b6::MNPbatchadjust(Mset, FFPE)
        FFPE <- mnp.v11b6::MNPgetFFPE(RGsetEpic)
        detP <- minfi::detectionP(RGsetEpic)
        bs <- minfi::getBeta(Mset_ba) #Mset
        sex = grabSexEst(bs, detP)
        message("\n~~~~~~~~~~~~~~~Generating ", sampleName, " cnv plot...\n")
        xx <- mnp.v11b6::MNPcnv(Mset,sex = sex,main = sampleID)
        getCnWebshot(xx, fn)
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

save.png.files <- function(rds, token){
    get.rd.info(rd_numbers=rds, token=token,sh_name=NULL) # input your RD-numbers here rd_numbers = c("RD-21-21")
    mySentrix <- gb$search.redcap(rd_numbers = rds, token)
    mySentrix <- mySentrix[!is.na(mySentrix$barcode_and_row_column),]
    mySentrix <- mySentrix[!is.null(mySentrix$barcode_and_row_column),]
    for (sam in rownames(mySentrix)) {
        sentrix <- mySentrix[sam,2]
        RGsetEpic <- grabRGset(getwd(),sentrix)
        tryCatch(
            expr= {gen.cnv.png2(RGsetEpic, sampleName=mySentrix[sam,1])},
            error= function(e){
                erTxt <- paste0("An error occured with ", mySentrix[sam,1]," png creation:")
                message(crayon::bgRed(erTxt),"\n",e,"\n")
                message(crayon::bgGreen("Trying next sample"))
            }
        )
    }
}

