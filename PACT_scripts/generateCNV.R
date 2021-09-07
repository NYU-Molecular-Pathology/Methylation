gb <- globalenv(); assign("gb", gb)

# FUN: Sets your directory and sources the helper functions
sourceFuns <- function(workingPath = NULL) {
    mainHub = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/"
    script.list <- c("SetRunParams.R","CopyInputs.R")
    if(!require("redcapAPI")){install.packages("redcapAPI", dependencies = T, type="both",ask=F)
        library("redcapAPI")}
    if (is.null(workingPath)) {workingPath = getwd()}
    scripts <- paste0(mainHub, script.list)
    invisible(lapply(scripts, function(i){devtools::source_url(i)}))
    gb$setDirectory(workingPath)
    return(gb$defineParams())
}

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

# FUN: Copies .idat files to your directory and saves samplesheet.csv
get.rd.info <- function(rd_numbers=NULL, token=NULL, sh_name=NULL){
    if (is.null(rd_numbers)){message("Input RD-numbers using get.rd.info(rd_numbers)")}
    if (is.null(sh_name)) {sh_name = "samplesheet.csv"}
    result <- gb$search.redcap(rd_numbers, token)
    samplesheet_ID = as.data.frame(stringr::str_split_fixed(result[, "barcode_and_row_column"], "_", 2))
    writeFromRedcap(result, samplesheet_ID) # writes API export as minfi dataframe sheet
    gb$get.idats(csvNam = sh_name)  # copies idat files from return to current directory
    return(result)
}

gen.cnv.png <- function(RGsetEpic, sampleName) {
    RGset=RGsetEpic
    imgName <- paste(sampleName, "cnv.png", sep="_")
    Mset <- mnp.v11b6::MNPpreprocessIllumina(RGsetEpic)
    Mset@annotation=c(array="IlluminaHumanMethylationEPIC", annotation="ilm10b4.hg19")
    FFPE <- mnp.v11b6::MNPgetFFPE(RGset)
    Mset_ba <- mnp.v11b6::MNPbatchadjust(Mset, FFPE)
    detP <- minfi::detectionP(RGset)
    bs <- minfi::getBeta(Mset)
    sexEstimate <-as.data.frame(signif(sest::get.proportion_table(bs, detP), digits = 2))
    yest <- as.double(sexEstimate$`p.Y:(-18,-5]`) >= 0.75
    yest1 <- as.double(sexEstimate$`Y:(0,0.1]`) >= 0.12
    sex <- ifelse((yest == TRUE && yest1 == TRUE), "male", "female")
    message("Generating ", sampleName, " cnv plot...")
    xx <- mnp.v11b6::MNPcnv(Mset,sex = sex,main = sampleID)
    thePlot<-mnp.v11b6::MNPcnvggplotly(xx, getTables = F)
    p<-plotly::ggplotly(thePlot)
    fn=file.path("~/Desktop",imgName)
    htmlwidgets::saveWidget(plotly::as.widget(p), "~/Desktop/temp.html")
    webshot2::webshot("~/Desktop/temp.html", file = fn, cliprect = "viewport", vwidth = 1152, vheight = 672)
    dev.off(); message("File saved:\n",imgName,"\n")
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

save.png.files <- function(rds, token){
    get.rd.info(rds, token=token) # input your RD-numbers here rd_numbers = c("RD-21-21")
    mySentrix <- gb$search.redcap(rd_numbers = rds, token)
    mySentrix <- mySentrix[!is.na(mySentrix$barcode_and_row_column),]
    
    for (sam in rownames(mySentrix)) {
        sentrix <- mySentrix[sam,2]
        RGsetEpic <- grabRGset(getwd(),sentrix)
        gen.cnv.png(RGsetEpic, sampleName=mySentrix[sam,1])
    }
    unlink("~/Desktop/temp.html")
    the.cnvs <- dir(path="~/Desktop",pattern="_cnv.png", full.names=T)
    print(the.cnvs)
}

library("conumee");library("sest");library("mnp.v11b6")
require(plotly)
require(htmlwidgets)
sourceFuns()
