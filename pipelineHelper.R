## Normal Methylation Run
## Jonathan Serrano Jan 14, 2020

#devtools::install_github("RGLab/Rtsne.multicore")}
#rstudioapi::executeCommand("toggleDocumentOutline")
rq <- function(pkg){return(!require(pkg, character.only = T))}
ld <- function(libName){return(library(libName, character.only = T))}
pk.in <- function(pkg){install.packages(pkg,dependencies=T,verbose=T,ask=F,suppressUpdates=F,type="binary")}
gh.in <- function(ghLnk){devtools::install_github(ghLnk, dependencies = T, type = "source")}
# FUN: Installs package from Bioconductor
bioInst <- function(pknm){if(rq(pknm)){BiocManager::install(pknm, update=T, ask=F)};ld(pknm)}

cranInstall <- function(pkg) {
  tryCatch(
    expr = {message("Checking Package...");cat(paste0(pkg, '\n'))
      if(length(find.package(pkg,lib.loc = NULL,quiet = T))==0){rq(pkg)}
    },
    error = function(cond) {
      message("Installing missing pacakge...")
      install.packages(pkg,dependencies=T,verbose=T,ask=F,suppressUpdates=F,type="source")
      ld(pkg)
    },warning = function(cond){message(cond)},finally = {ld(pkg)}
  )
}

notifyMissing <- function(varNam){
  message("\nThe variable ",varNam, " was not found, try assinging the value using:")
  mm1 <- paste0("setVar(\'",varNam, "\',\'value of variable\')")
  cat(mm1)
}
# Loads default packages or custom if input provided
loadPacks <- function(ezLibs = NULL, ghpack = NULL, bioPks = NULL) {
  pkgs <-
    c("devtools",
      "tidyverse",
      "remotes",
      "sjmisc",
      "stringi",
      "digest",
      "RCurl",
      "parallel"
    )
  if (is.null(ezLibs)) {
    ezLibs <-
      c("remotes",
        "tidyverse",
        "sjmisc",
        "stringi",
        "redcapAPI",
        "needs",
        "digest",
        "RCurl",
        "docstring",
        "minfi",
        'magick'
      )
  }
  if (is.null(ghpack)) {
    ghpack <-
      data.frame(
        docstring = "dasonk/docstring",
        rstudioapi = "rstudio/rstudioapi",
        easypackages = "jakesherman/easypackages",
        redcapAPI = "nutterb/redcapAPI"
      )
  }
  if (is.null(bioPks)) {bioPks <- c("minfi")}
  cranInstall("devtools")
  if (rq("easypackages")) {devtools::install_github("jakesherman/easypackages")}
  easypackages::packages("parallel", "doParallel", "foreach", "compiler", prompt=F)
  compiler::enableJIT(3); compiler::compilePKGS(enable = T);compiler::setCompilerOptions(suppressAll = T, optimize = 3)
  x=NULL
  # Installs packages from CRAN
  foreach::foreach(x = 1:length(pkgs)) %do% {cranInstall(pkgs[x])}
  # Installs packages from GitHub Repos
  foreach::foreach(x = 1:length(ghpack)) %do% {if (rq(names(ghpack[x]))) {gh.in(ghpack[[1, x]])}}
  foreach::foreach(x = 1:length(bioPks)) %do% {bioInst(bioPks[x])}
  fldx <- function(nmn = "foldAll") {return(rstudioapi::executeCommand(nmn))}
  assign(x = "fldx", value = fldx, envir = .GlobalEnv)
  easypackages::libraries(ezLibs)
}
installMgmt <- function(){
  url="https://github.com/badozor/mgmtstp27/raw/master/trunk/Rpackage/mgmtstp27_0.6-2.tar.gz"
  destfile <- file.path("~/Desktop", "mgmtstp27_0.6-2.tar.gz")
  download.file(url, destfile); message("~~Manually installing from desktop~~")
  install.packages(destfile,repos=NULL, type="source")
}

checkQCpkg <- function(){
  library("knitr")
  library("kableExtra")
  library("magick")
  library("webshot")
  library("plyr")
  library("beepr")
  library("ggplot2")
  library("knitr")
  library("reshape2")
  library("data.table")
  library("DT")
  library("plotly")
  library("MethylAid")
  library("minfi")
  library("scales")
  library("IlluminaHumanMethylation450kmanifest")
  library("IlluminaHumanMethylationEPICmanifest")
  library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
  library("Biobase")
  library("RColorBrewer")
  library("limma")
  library("ggfortify")
  library("Rtsne")
  library("qdapTools")
  library("gplots")
  library("readxl")
  library("stringr")
  library("ggrepel")
  library("Polychrome")
  library("tinytex")
  library("gridExtra") 
  library("rmarkdown")
  library("tinytex")
  library("BiocParallel")
  require("ggplot2")
  require("pals")
  require("scales")
  require("stringr")
  require("scales")
  library('grid')
  require('grid')
}



# Define Default input parameters
setRunType <- function(isMC){
  mv11="Install or update mnp.v11b6 package: use install.or.load()"
  sexEst="https://github.com/jungch/sest/raw/master/sest.tar"
  if (isMC){
    loadPacks()
    if (rq("BiocManager")){cranInstall("BiocManager")}
    # FUN: Installs package from Bioconductor
    bioInst <- function(pknm){if(rq(pknm)){BiocManager::install(pknm)}}
    # FUN: You can use easy libraries here by passings loadPacks(libs)
    ghInstall <- function(pkgName, ghLink) {
      if(rq(pkgName)){devtools::install_github(ghLink, dependencies = T, type = "source")}else{ld(pkgName)}
    }
    if (rq("mgmtstp27")){bioInst("IlluminaHumanMethylation450kmanifest")}
    if (rq("rmarkdown")){ghInstall("rmarkdown",'rstudio/rmarkdown')}
    if (rq("mgmtstp27")){installMgmt()}
    if (rq("mnp.v11b6")){message(mv11)}
    if(rq("sest")){install.packages(pkgs = sexEst,repos = NULL,type = "source", dependencies = T)}
    tryCatch(
      expr={checkQCpkg()},
      error=function(cond){
        message(cond)
      },
      warning=function(cond){
        message(cond)
      }
    )
    }else{loadPacks()}
}

checkScript <-function(mntPath){
  ifelse(file.exists(mntPath), source(mntPath), print(paste("Script not found at:", mntPath)))
}

pathMsg <- function(varName=NULL, foldr=NULL, wrong=NULL) {
  mm1 = "!!!This location does not exist:"
  mm2 = "\n\nAssign correct path using the fuction setVar():\n\n"
  mm3 = paste0("setVar(\"", varName, "\", 'yourCustomPathname')")
  mm4=paste(varName, "=")
  if (wrong) {
    cat("", sep="\n\n"); message(mm1)
    cat(foldr, sep="\n"); cat(mm2)
    cat("", sep="\n");message(mm3)
    } else{message(mm4);cat(foldr);cat("", sep="\n\n")}
}

checkMapping <- function(varName, foldr) {
  isReal = NULL
  isQcRmd <- stringr::str_detect(foldr, "_QC.Rmd")
  isTkn <- nchar(foldr)==32
  if (isTkn) {
    message("Token is:")
    cat(foldr, sep = "\n\n")
  } else{
    isReal <-
      ifelse(isQcRmd == T &
               length(isQcRmd) == 1,
             file.exists(foldr),
             dir.exists(foldr))
    if (!isReal) {
      pathMsg(varName, foldr, T)
    } else{
      pathMsg(varName, foldr, F)
    }
  }
}
setGlobalPaths <- function(dirPaths) {
  cat("", sep="\n\n")
  valNames <- names(dirPaths)
  for (vv in 1:length(dirPaths)){
    val <- paste0(dirPaths[1,vv])
    varnm <- paste0(valNames[vv])
    assign(varnm,val, envir=.GlobalEnv)
    checkMapping(varnm,val)
  }
}

# Assigns variable passed to global environment
setVar<-function(valueName,val){return(assign(valueName, val, envir=.GlobalEnv))}

assignVar <- function(varStr, assignedVal){
  tryCatch(
    expr = {if(!is.null(get(varStr))){message(varStr," =\n",assignedVal,"\n")}},
    error = {setVar(varStr,assignedVal)}
    )
}

ckNull <- function(nullVar, subVar, varName){
  if (is.null(nullVar)){
    message("\nRe-assigning variable:");cat(varName,sep = "\n")
    cat(paste("to", subVar,"\n"))
    setVar(as.character(varName),as.character(subVar))
    return(paste0(subVar))
    }else{return(paste0(nullVar))}
}

getDefaults <- function() {

  cbVol = "/Volumes/CBioinformatics"
  moVol = "/Volumes/molecular"

  defaultParams <- data.frame(
    mnp.pk.loc = paste0(file.path(cbVol, "Methylation/in_house/mnp.v116/mnp.v11b6")),
    ApiToken = "",
    methDir = paste0(file.path(cbVol, "jonathan/Rprojects/Clinical_Methylation")),
    clinDrv = paste0(file.path(moVol, "MOLECULAR LAB ONLY/NYU-METHYLATION")),
    rschOut = paste0(file.path(cbVol, "jonathan/Rprojects/Research_runs/")),
    clinOut = paste0(file.path(moVol, "MOLECULAR/MethylationClassifier")),
    rsch.idat = "/Volumes/snudem01labspace/idats",
    clin.idat = paste0(file.path(moVol, "MOLECULAR/iScan")),
    QC_file = paste0(system.file('Methyl_QC.Rmd', package = "mnp.v11b6")),
    baseDir = paste0(file.path(cbVol, "jonathan/Rprojects")),
    stringsAsFactors=F
  )
  return(defaultParams)
}

defineParams <-
  function(mnp.pk.loc = NULL,
           ApiToken = NULL,
           methDir = NULL,
           clinDrv = NULL,
           rschOut = NULL,
           clinOut = NULL,
           rsch.idat = NULL,
           clin.idat = NULL,
           QC_file = NULL,
           isMC = T,
           baseDir = NULL,
           runID = NULL) {

    defVars <- getDefaults()
    inVars <- list(mnp.pk.loc, ApiToken, methDir, clinDrv, rschOut, clinOut,
                   rsch.idat, clin.idat, QC_file, baseDir)
    i=1:length(inVars)
    invisible(lapply(
      i,function(x){if(!is.null(inVars[[x]])){setVar(names(defVars[x]), inVars[[x]])}}))
    setRunType(isMC)
    i=1:length(defVars)
    invisible(sapply(i,FUN=function(i){
      y <- defVars[i]
      assignVar(names((y)), paste0(defVars[,i]))}))

    if(!isMC) {methDir=rschOut}
    #setGlobalPaths(defVars)
}

setDirectory <- function(foldr) {
  bsDir = paste("cd", foldr)
  mm2 = paste("Location Not Found:", foldr)
  if (dir.exists(foldr)) {system(bsDir)
    setwd(foldr)
    } else{warning(mm2)}
}

getLastRunID <- function(methDir, returnAll=F) {
  prevMC <- list.dirs(methDir)
  prevMC <- as.data.frame(stringi::stri_split_fixed(prevMC, .Platform$file.sep))[7,]
  row.names(prevMC)="Run Numbers:"
  prevMC <- as.data.frame(prevMC[stringr::str_detect(prevMC, "MGDM")])
  names(prevMC)=NULL;newestRun = prevMC[length(prevMC)]
  if(returnAll){for(mcv in prevMC){cat(mcv,sep="\n")}}else{return(paste0(newestRun))}
}

printParams <-function(varsNames){
  cat(crayon::bgGreen("\nThese are your variables:\n"))
  invisible(foreach::foreach(i=1:length(varsNames))%do%{
    cat(paste(crayon::bgBlue("\nVariable:"),
              names(varsNames)[i],
              crayon::bgCyan("\nAssigned value:\n"),
              varsNames[1,i],"\n"))
  })
  message("\nUse ", crayon::bgBlue('gb$setVar("varNameString", varValue)'),
          " to change any variables\n")
}

copyWorksheetFile <- function(runID=NULL, runYear=NULL) {
  if (is.null(runYear)){runYear=paste0(format(Sys.Date(), "%Y"))}
  if (is.null(runID)){runID=paste0(basename(getwd()))}else(runID=runID)
  isMC = sjmisc::str_contains(runID, "MGDM") |sjmisc::str_contains(runID, "MC")
  if(!isMC) {mountLoc <- "/Volumes/snudem01lab/snudem01labspace/Methylation_Worksheets"}
  if(isMC) {mountLoc <- file.path(clinDrv,"WORKSHEETS")}
  message("\nCopying file from: ");cat(mountLoc)
  if (dir.exists(mountLoc)) {
    fileLoc <- file.path(mountLoc,runYear,paste0(runID,".xlsm"))
    if (file.exists(fileLoc)) {
      message("\n\nCopying worksheet from Worksheets Folder:")
      if(!file.exists(file.path(getwd(),basename(fileLoc)))){
      fs::file_copy(fileLoc, getwd())}
    }else{
      message("File not found:");cat(fileLoc)
      message("\nMake sure your path is correct, try print(gb$copyWorksheetFile)")
      }
  } else{
    message("Folder not found:\n");cat(mountLoc)
    message("\nEnsure your volume path is correct, check print(copyWorksheetFile)")
  }
}

setRunDir <- function(runID=NULL, workFolder=NULL){
  runID <- ckNull(runID, paste0(basename(getwd())), deparse(substitute(runID,env=.GlobalEnv)))
  workFolder <- ckNull(workFolder, methDir, deparse(substitute(workFolder,env=.GlobalEnv)))
  newRun <- file.path(workFolder, runID)
  if(!dir.exists(newRun)){
    dir.create(newRun);cat("creating folder: ",newRun)
    setDirectory(newRun)}else{setDirectory(newRun)}
  return(newRun)
}

# Returns Total Sample Count in the run
getTotalSamples <- function(){
  thisSh <- list.files(getwd(), "*.xlsm")
  temp <- stringi::stri_detect_fixed(thisSh, "~$")
  thisSh <- thisSh[!temp]
  if(length(thisSh)==0){print("No .xlsm sheet, defaulting to 16");return(16)}
  worksheet <- readxl::read_excel(thisSh[1], sheet=1, col_names=F, range="B4:B4",trim_ws=T)
  if (length(worksheet) == 0) {
    message("Samplesheet ", thisSh[1]," is invalid format, manually edit")
    fs::file_copy(
      path=file.path(baseDir, "Clinical_Methylation/methylation_run_TEMPLATE_new.xlsm"), new_path=getwd())
    } else {message("Total sample count found is: ", worksheet[1])}
  return(paste0(worksheet[1]))
}

writeSampleSheet <- function(df, samplesheet_ID, bn = NULL, sampleName, dnaNumber, Sentrix) {
  if (is.null(bn)) {bn = file.path(getwd(), df$Sentrix)}
  samplesheet_csv = data.frame(
    Sample_Name =  df[, sampleName],
    DNA_Number =df[, dnaNumber],
    Sentrix_ID = samplesheet_ID[, 1],
    Sentrix_Position = samplesheet_ID[, 2],
    SentrixID_Pos = df[, Sentrix],
    Basename = paste0(bn),
    RunID = df$Batch,
    MP_num = df$MP_number,
    tech = df$Tech,
    tech2 = df$Tech2,
    Date = df$Date
  )
  write.csv(samplesheet_csv, file = "samplesheet.csv", quote = F, row.names = F)
}

# Check if the QC File will be read
checkRunOutput <- function(runID) {
  location <- file.path("~/Desktop",runID,paste0(runID,"_Redcap.csv"))
  if (!file.exists(location)) {
    message("File not found: ")
    cat(location,sep="\n")
    message("QC Summary Table will not Knit")
  } else{fs::file_copy(path=location, new_path=getwd(), overwrite = T)}
}

# Read samplesheet
readSampleSheet <- function(runID=F, totalSam=F, wks=F) {
 file.list <- list.files(path=getwd(), "*.xlsm")
 temps <- stringi::stri_detect_fixed(file.list, "~$")
 file.list <- file.list[!temps]
 sampleSheet <- paste0(file.list[1])
 message(paste0("Reading worksheet named: ", sampleSheet))
 worksheet <- readxl::read_excel(sampleSheet, sheet=2, col_names=T, col_types="text", trim_ws=T)
 wsDate <-  as.data.frame(
   readxl::read_excel(sampleSheet, sheet=1, col_names=F, range="F4:F4", trim_ws=T))[1]
 names(wsDate)="Date"
 colnames(worksheet)
  worksheet$Date <- paste0(wsDate$Date[1])
 if (runID  ==  T) {return(worksheet$Project[1])}
 if (totalSam  ==  T){
   runNum <- readxl::read_excel(sampleSheet, sheet=1, col_names=F, range="B4")
   runNum <- as.numeric(runNum)
   return(runNum)
   }
 if (wks  ==  T) {return(worksheet)}
}

# Sends an email notification that the run is complete from redcap admin
launchEmailNotify <-function(runID){
  rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", ApiToken)
  ur=paste0(rcon$url);tk=rcon$token
  sam_id = paste0(runID,"_QC")
  com = "sample_qc"
  record= data.frame(record_id=sam_id,comments=com)
  datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T);print(datarecord)
  res<-RCurl::postForm(ur,token=tk,content='record',format='json',type='flat',data=datarecord)
  message("Run Emailed for QC");cat(res)
}

create.QC.record <- function(runID=NULL){
  rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", ApiToken)
  uri=paste0(rcon$url); tk=rcon$token
  if(is.null(runID)){runID<-paste0(basename(getwd()))}
  record = c(record_id=paste0(runID,"_QC"), run_number=runID)
  qcdata <- jsonlite::toJSON(list(as.list(record)), auto_unbox=T)
  rr<-RCurl::postForm(uri,token=tk,content='record',format='json', type='flat',data=qcdata)
  rr;message("Created QC Record"); print(qcdata)
}

# Imports the xlsm sheet 3 data
importRedcapStart <- function(nfldr){
  rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", ApiToken)
  uri=paste0(rcon$url); tk=rcon$token
  samSh <- list.files(path=getwd(), full.names=T, ".xlsm")
  sampleNumb <- getTotalSamples()
  sh_Dat <-as.data.frame(readxl::read_excel(samSh,sheet=3,range="A1:N97", col_types=c("text")))
  sh_Dat <- sh_Dat[,1:13]
  sampleNumb=as.integer(sampleNumb);sh_Dat = sh_Dat[1:sampleNumb,]
  runName=basename(getwd())
  filnm =paste0(sh_Dat$record_id, "_cnv.png")
  pathNam=file.path(nfldr,paste0(runName,"_CNVs"),filnm)
  rms=paste(c("control_","low_"),collapse ='|')
  pathNam <- stringr::str_replace(pathNam, rms, "")
  pathNam <- stringr::str_replace(pathNam, "//", "/")
  sh_Dat$cnv_file_path <- pathNam

  if (!is.null(sh_Dat)){
    for(n in 1:nrow(sh_Dat)){
      print(n);record=c(sh_Dat[n,])
      datarecord = jsonlite::toJSON(list(as.list(record)), auto_unbox=T);print(datarecord)
      res<-RCurl::postForm(uri,token=tk,content='record',format='json',type='flat',data=datarecord)
      message("Record Uploaded"); print(res)
      }
    create.QC.record(runName)
    }
}

readSheetWrite <- function(sampleNumb= NULL, runID = NULL) {
  if(is.null(sampleNumb)){sampleNumb<-getTotalSamples()}
  if(is.null(runID)){runID<-paste0(basename(getwd()))}
  current.run.Folder <- paste0(methDir, "/",runID)
  ww1="No tech name found: check df$Tech in samplesheet"
  ww2="No tech name found: check df$MP_number in samplesheet"
  if (!file.exists("samplesheet.csv")) {
    sampleNumb=as.integer(sampleNumb)
    worksheet=readSampleSheet(wks=T)
    hdrs = dimnames(worksheet)
    print(data.frame(HEADERS=hdrs[[2]]))
    sampleName=paste(hdrs[[2]][9])
    Sentrix=paste(hdrs[[2]][1])
    dnaNumber=paste(hdrs[[2]][10])
    mpnum=paste(hdrs[[2]][8])
    Var.names=c("~~~~~~~",sampleName,Sentrix, dnaNumber, mpnum);cat("\n\n")
    print(data.frame(Var.names))
    df=as.data.frame(worksheet)[1:sampleNumb, ]
    head(df)
    samplesheet_ID=as.data.frame(stringr::str_split_fixed(df[, Sentrix], "_", 2))
    print(samplesheet_ID)
    bn <- file.path(current.run.Folder, df[, Sentrix])
    message("Basename layout: ");message(bn[1])
    note=paste(df$Notes[1])
    df$Notes <- note
 if (is.null(df$Tech)){warning(ww1)
   df$Tech <- "NA"}
 if (is.null(df$MP_number)){warning(ww2)
   df$MP_number <- "none"}
    writeSampleSheet(df, samplesheet_ID=samplesheet_ID,bn=bn, sampleName, dnaNumber, Sentrix)
    } else{message("samplesheet.csv already exists in this directory!")}
}

# FUN: Returns a list of idat files given an idat drive location -
getAllFiles <- function(idatDir, csvNam=NULL) {
  if(!is.null(csvNam)){ssheet=read.csv(csvNam,strip.white=T)
    barcode=as.vector(ssheet$Sentrix_ID)} else{
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
  message("Files copying...")
  fs::file_copy(allFi, file.path(getwd()))
  idcs = basename(allFi)
  idatsCopied <- idcs[idcs != ""]
  success = file.exists(idatsCopied)
  message(".idat files that failed to copy:")
  if ((length(success) > 0)) {print(idatsCopied[!success])} else{cat("none\n")}
}

# FUN: Returns a list of idat files that exist on Molecular and Snuderl lab drives -
get.idats <-function(csvNam = "samplesheet.csv",
                     rsch.idat = "/Volumes/snudem01labspace/idats",
                     clin.idat = file.path("/Volumes/molecular", "MOLECULAR/iScan"))
{
  if (!dir.exists(rsch.idat)) {message("\nDirectory not found, make sure rsch.idat folder location is mounted:\n", rsch.idat)}
  if (!dir.exists(clin.idat)) {message("\nDirectory not found, make sure clin.idat folder location is mounted:\n", clin.idat)}
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

# Sets Parallel parameters for Mac Darwin system
checkRver <- function(){
  if (!nzchar(Sys.getenv("RSTUDIO_TERM")) && Sys.info()["sysname"]=="Darwin" && getRversion()>='4.1.0')
  {parallel:::setDefaultClusterOptions(setup_strategy="sequential")
    }else{warning("You do not have R 4.1 installed")}
}

# FUN: Check Packages & Environment -
checkPackages <- function(pathtoFile) {
  thisHist=file.path(getwd(),".Rhistory")
  print(pathtoFile)
  if (!file.exists(pathtoFile)){warning("You must re-map path to /Volumes/CBioinformatics")}
  if (file.exists(thisHist)){
    file.remove(thisHist)
    Sys.setenv(R_ENABLE_JIT=T)
    options(mc.cores=parallel::detectCores())
      }
  if (rq("doParallel")) {install.packages("doParallel", dependencies=T)}
  easypackages::libraries("parallel","foreach","doParallel")
  if (rq("pkgbuild")) {install.packages("pkgbuild", dependencies=T)}
  require(compiler)
  compiler::enableJIT(3);compiler::compilePKGS(enable=T)
  compiler::setCompilerOptions(suppressAll=T, optimize=3)
  checkRver()
}

# FUN: Copies samplesheet to Desktop folder
moveSampleSheet <- function(methDir, runID=NULL) {
  if (is.null(runID)){runID=paste0(basename(getwd()))}
  baseFolder=paste0("~/Desktop/",runID,"/")
  if(!dir.exists(baseFolder)){dir.create(baseFolder)}
  currDir=paste0(methDir, "/",runID); endDir = paste0(baseFolder, runID, "_samplesheet.csv")
  fs::file_copy(path=paste0(currDir,"/samplesheet.csv"),new_path=baseFolder,overwrite=T)
  file.rename(from=paste0(baseFolder,"samplesheet.csv"), to=endDir)
}
# FUN: Load Classifier Versions -
install.or.load <- function(pathtoFile=NULL, instNew=F, rmpkg=F) {
  pathtoFile<-ckNull(pathtoFile, mnp.pk.loc, deparse(substitute(pathtoFile,env=.GlobalEnv)))
  checkPackages(pathtoFile)
  if(rmpkg){remove.packages("mnp.v11b6")}
  if(instNew){
    install.packages(pathtoFile, repos=NULL, type="source")
    install.or.load(instNew=F)
    }else{ld("mnp.v11b6")}
}

# FUN: Generate CNV image
saveCNVhtml <- function(data) {
    # get RGSet for sample
    RGsetEpic <- read.metharray(file.path(getwd(),'idats',unlist(data[4])), verbose = F, force = T)
    # Save the CNV as a png file
    MsetEpic <- preprocessRaw(RGsetEpic)
    sample_id <- unlist(data[1])
    xx <- MNPcnv(MsetEpic,sex = unlist(data[6]),main = sample_id)
    hg <- MNPcnvggplotly(xx,getTables = F)
    htmlwidgets::saveWidget(hg, paste(sample_id, "_cnv.html", sep = ""))
}

# Save the CNV as a png file
generateCNVpng <- function(RGsetEpic, sampleName) {
  imgName <- paste(sampleName, "cnv.png", sep="_")
  MsetEpic <- minfi::preprocessRaw(RGsetEpic)
  png(filename=imgName,width=1820, height=1040)
  mnp.v11b6::MNPcnvplot(MsetEpic)
  dev.off()
}

testReport <- function(i = NULL, varsOnly = T) {
  if (is.null(i)) {i = 1}
  data = read.csv("samplesheet.csv", strip.white = T)
  tech = as.character(data$tech)
  tech2 = as.character(data$tech2)
  mp_num = as.character(data$MP_num)
  run_id = as.character(data$RunID)
  barcode = as.numeric(data$Sentrix_ID)
  samplename_data = as.character(data$Sample_Name)[i]
  tech = tech[i]
  tech2 = tech2[i]
  mp_num = mp_num[i]
  run_id = run_id[i]
  bnumber = as.character(data$DNA_Number)[i]
  barcode = barcode[i]
  runPath = getwd()
  senP = as.character(data$SentrixID_Pos)
  newTotal = paste("Running", i, "of", length(senP))
  thisrun <-
    paste(samplename_data, senP[i], tech, mp_num, run_id, bnumber, runPath, sep = ", ")
  print(thisrun)
  message(newTotal)
  cat(" ", sep = "\n\n")
  pathEpic = paste(runPath, data$SentrixID_Pos[i], sep = '/')
  currpath = getwd()

  RGsetE <- minfi::read.metharray(pathEpic, verbose = T, force = T)
  RGsetEpic <-checkBarcode(RGsetE, barcode)
  sampleID = paste(samplename_data, sep = "")
  RGset = RGsetEpic
  tsne = T
  sampleID = samplename_data
  run_id = run_id
  tech = tech
  tech2 = tech2
  bnumber = bnumber
  mp_number = mp_num
  FFPE = "FFPE"
  RGset <- RGset[, 1]

  if (varsOnly==T) {
    setVar("RGset", RGsetEpic)
    setVar("i", i)
    setVar("tsne", T)
    setVar("FFPE", "FFPE")
    setVar("sampleID", samplename_data)
    setVar("run_id", run_id)
    setVar("tech", tech)
    setVar("tech2", tech2)
    setVar("bnumber", bnumber)
    setVar("mp_number", mp_num)
  }else{
  system.time(
    mnp.v11b6::MNPreport(
      RGset = RGsetEpic,
      tsne = T,
      sampleID = sampleID,
      run_id = run_id,
      tech = tech,
      tech2 = tech2,
      bnumber = bnumber,
      mp_number = mp_num,
      FFPE = "FFPE"
    )
  )
  }
}

# New B5 or B4 Manifest
checkBarcode <- function(RGsetEpic, barcode) {
  aEpic=c(array="IlluminaHumanMethylationEPIC", annotation="ilm10b4.hg19")
  #bEpic=c(array="IlluminaHumanMethylationEPIC", annotation="ilm10b5.hg38")
  a450k=c(array="IlluminaHumanMethylation450k", annotation="ilmn12.hg19")
  if (barcode >= as.numeric("204220033000")) {
    RGsetEpic@annotation=aEpic;message("Barcode is new EPIC B5 manifest")
    }
  if (RGsetEpic@annotation['array']=="IlluminaHumanMethylation450k") {
    RGsetEpic@annotation=a450k; message("Barcode is 450k using B4 manifest")
  }else {message("Barcode is Default manifest")}
  return(RGsetEpic)
}
# Creates the html file report given variables
saveReportOutput <- function(samplename_data,sentrix_pos_list,tech,tech2,mp_num, run_id,bnumber,runPath,i, barcode,genCn=T) {
    newTotal <- paste("Now running sample", paste0("#",i), "of", length(sentrix_pos_list))
    thisrun <- paste(samplename_data, sentrix_pos_list[i], tech, mp_num, run_id, bnumber, runPath, sep=", ")
    print(thisrun);message(newTotal);cat(" ", sep = "\n\n")
    pathEpic <- file.path(runPath, sentrix_pos_list[i])
    RGsetE <- minfi::read.metharray(pathEpic, verbose = T, force = T)
    RGsetEpic <- checkBarcode(RGsetE, barcode)
    sampleID=paste0(samplename_data)
    if(genCn){generateCNVpng(RGsetEpic, paste0(sampleID))} # Saves the CNV as a png file if param = TRUE
  # Generates report html tsne
  system.time(mnp.v11b6::MNPreport(
    RGset=RGsetEpic, tsne=T, sampleID=sampleID, run_id=run_id, tech=tech, tech2=tech2,
    bnumber=bnumber,mp_number=mp_num,FFPE="FFPE"
    )
  )
  totfin <-paste("Finished #", i, "\n")
  message(totfin)
}

# QC REPORT maker
generateQCreport <- function(runID=NULL, qc=NULL) {
  runID<-ckNull(nullVar = runID, subVar=basename(getwd()), deparse(substitute(runID,env=.GlobalEnv)))
  QC_file <- system.file('Methyl_QC.Rmd', package = "mnp.v11b6")
  if (!file.exists(QC_file)){message("Check Working directory, QC_file.rmd not found")}
  fs::file_copy(QC_file, getwd(), overwrite = T)
  currentQC = list.files(getwd(),"*QC.Rmd", full.names=T)
  qcFile = paste0(runID,"_QC.html") # output file name
  rmarkdown::render(currentQC, output_file=file.path(dirname(currentQC), qcFile), params=list(runID=runID))
  currentQC <- stringr::str_replace_all(string=currentQC, ".Rmd", "_cache")
  unlink(currentQC, recursive=T) #clear cache
}
# FUN: Clean wd Returns Directory Junk Files and Copies redcap file
removeExtras <- function(runID) {
  folder.list = list.dirs(path=getwd(), full.names=T, recursive=T)
  temp.folder.list <- folder.list[stringi::stri_detect_fixed(folder.list, "_files")]
  thisReport <- paste0(runID, "_Redcap.csv")
  file.list<-list.files(path=paste0("~/Desktop/",runID,"/"), thisReport, full.names=T)
  fs::file_copy(file.list, getwd(), overwrite=T)
  print(as.data.frame(temp.folder.list))
  unlink(temp.folder.list, recursive=T)
  return(temp.folder.list)
}
# Helper archive function: creates a "previous" folder when reports already exists
create.prev.folder<-function(prevs,oldFi){
  dir.create(prevs)
  fs::file_copy(path=oldFi, new_path=prevs)
}

# FUN: Copies Reports to Z drive
copy.to.clinical <- function(clinOut, runID, runYear) {
  newFolder <- file.path(clinOut, runYear,runID); message(newFolder)
  if (!dir.exists(newFolder)) {dir.create(newFolder)}
  if (dir.exists(newFolder)) {
    oldFi=list.files(path=newFolder, full.names=T)
    prevs=file.path(newFolder,"previous")
    if (length(oldFi) > 0) {create.prev.folder(prevs,oldFi)}
    #file.list=list.files(path=getwd(), full.names=T, "*.html")
    message(paste0("Copying Reports to: ", newFolder))
    #newFolder <- rep(paste(newFolder,sep="/"), length(file.list))
    #fs::file_copy(file.list, newFolder, overwrite = T)
    sys_cmd <- paste("cp",paste(getwd(),"*.html",sep="/"), newFolder)
    system(sys_cmd)
  }
  #copy.cnv.files(newFolder, runID, runYear="2020")
}
# REDCap: API call & Upload
# uploads the redcap classifier values must convert to JSON first
importRedcapOutput <- function(rcon,samsheet=NULL) {
  rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", ApiToken)
  ur=paste0(rcon$url);tk=rcon$token
  if(is.null(samsheet)){
    samsheet=list.files(path=getwd(),full.names=T,"_Redcap.csv",recursive=F)}else{samsheet=samsheet}
  if (length(samsheet) < 1) {message("Redcap headers csv file not found")}
  if (length(samsheet) == 1) {

    data<-read.csv(samsheet)
    for (n in 1:nrow(data)) {

      datarecord = jsonlite::toJSON((as.list(data[n,])), auto_unbox=T);print(datarecord)
      message("Record Uploaded:")
      RCurl::postForm(
        ur,
        token = tk,
        content = 'record',
        format = 'csv',
        type = 'flat',
        data = datarecord
      )
      cat(redcapAPI::importRecords(rcon,data,"normal","ids",returnData = F))
    }
  }else{message("no redcap file found")}
}

# Copy Output Files
copy.cnv.files <- function(newFolder, runID, runYear=NULL) {
  if (is.null(runYear)){runYear=paste0(format(Sys.Date(), "%Y"))}
  cnv_folder <- file.path(newFolder,paste0(runID,"_CNVs/"))
  cnvNames <- list.files(path=getwd(), full.names=T, "*_cnv.png")
  if(length(cnvNames)>2){
    message(paste0("Copying PNG to: ", cnv_folder));print(as.data.frame(cnvNames))
    if (!dir.exists(cnv_folder)){dir.create(cnv_folder)}
    if (dir.exists(cnv_folder)){fs::file_copy(cnvNames,cnv_folder)}
  }
}


uploadCnPng <- function() {
  rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", ApiToken)
  samSh <- list.files(path=getwd(), full.names=T, ".xlsm")
  sampleNumb <- getTotalSamples()
  sh_Dat <-as.data.frame(readxl::read_excel(samSh,sheet=3,range="A1:N97", col_types=c("text")))
  sampleNumb=as.integer(sampleNumb)
  sh_Dat = sh_Dat[1:sampleNumb,]
  runName=basename(getwd())
  filnm =paste0(sh_Dat$tm_number,"_", sh_Dat$b_number, "_cnv.png")
  pathNam=file.path(getwd(),filnm)
  rms=paste(c("control_","low_"),collapse ='|')
  pathNam <- stringr::str_replace(pathNam, rms, "")
  pathNam <- stringr::str_replace(pathNam, "//", "/")
  records <- sh_Dat$record_id
  for (idx in 1:length(records)) {
    pth = pathNam[idx]
    recordName = paste0(records[idx])
    message("Importing Record:");cat(recordName, sep = "\n");cat(pth, sep = "\n")
    redcapAPI::importFiles(rcon, pth, recordName, field="methyl_cn", repeat_instance=1)
  }
}

# FUN: Copies Reports to Z drive
copy.to.outFolder <-
  function(clinDrv = NULL,
           runID,
           runYear = NULL) {
    if (is.null(runYear)) {
      runYear = paste0(format(Sys.Date(), "%Y"))
    }
    if (is.null(clinDrv)) {
      clinDrv = clinDrv
    }
    isMC = sjmisc::str_contains(runID, "MGDM") |
      sjmisc::str_contains(runID, "MC")
    rsrch = file.path("/Volumes/snudem01labspace/FINAL_PDF_Reports_Brain",
                      runID)
    clin = file.path(clinDrv, "Results", runYear, runID)
    runYear = ifelse(isMC, paste0("20", stringr::str_split_fixed(runID, "-", 2)[1]), runYear)
    newFolder = ifelse(isMC, clin, rsrch)
    message(newFolder)
    if (!dir.exists(newFolder)) {
      dir.create(newFolder)
    }
    oldFi = list.files(path = newFolder, full.names = T)
    prevs = file.path(newFolder, "previous")
    if (length(oldFi) > 0) {create.prev.folder(prevs, oldFi)} # saves any old files
    file.list = list.files(path = getwd(), "*.html", full.names = T)
    message("Copying Existing Reports to Folder")
    fs::file_copy(file.list, newFolder, overwrite = T)
    if (isMC) {
      cnList <- list.files(getwd(), "_cnv.png",recursive = F)
      hasCn <- length(cnList) > 2
      if (hasCn) {
        copy.cnv.files(newFolder, runID, runYear)
        uploadCnPng()
      }
      clinOut = file.path(stringr::str_split_fixed(clinDrv, " ", 2)[1],"MethylationClassifier")
      importRedcapStart(clinOut)
      copy.to.clinical(clinOut, runID, runYear)
    }
    return(file.list)
  }

# Creates QC record and uploads reports to redcap
uploadToRedcap <- function(file.list) {
  rcon <- redcapAPI::redcapConnection("https://redcap.nyumc.org/apps/redcap/api/", ApiToken)
  importRedcapOutput(rcon)
  #sla = .Platform$file.sep
  #slashes <- as.integer(stringr::str_count(file.list[1], sla) + 1)
  recordName <-stringr::str_replace_all(string = paste0(file.list), ".html", "")
  #runIDs <-stringr::str_split_fixed(string = recordName, "/", n = slashes)[, (slashes - 1)]
  runIDs <- rep(basename(getwd()), length(recordName))
  #records <-stringr::str_split_fixed(string = recordName, "/", n = slashes)[, slashes]
  records <- basename(recordName)
  for (idx in 1:length(records)) {
    pth = file.list[idx]
    recordName = paste0(records[idx])
    runID = paste0(runIDs[idx])
    message("\nImporting Record:")
    cat(recordName, sep = "\n")
    cat(pth, sep = "\n")
    data = data.frame(record_id = recordName, run_number = runID)
    redcapAPI::importRecords(
      rcon, data, overwriteBehavior = "normal",
      returnContent = "ids",
      returnData = T
    )
    redcapAPI::importFiles(
      rcon = rcon,
      file = pth,
      record = recordName,
      field = "classifier_pdf",
      repeat_instance = 1
    )
  }
}

do_report <-function(data = NULL) {
  require(rmarkdown)
  if(!is.null(data)){
  samplename_data = paste0(data[,1])
  sentrix_pos_list = (data[,5])
  tech = paste0(data[,9])
  tech2 = paste0(data[,10])
  mp_number = paste0(data[,8])
  run_id = paste0(data[,7])
  bnumber = paste0(data[,2])
  runPath = getwd()
  barcode = as.numeric(data[,3])
  pathEpic = file.path(runPath, sentrix_pos_list)
  RGsetEpic <- minfi::read.metharray(pathEpic, verbose = T, force = T)
  aEpic=c(array="IlluminaHumanMethylationEPIC", annotation="ilm10b4.hg19")
  a450k=c(array="IlluminaHumanMethylation450k", annotation="ilmn12.hg19")
  if (barcode >= as.numeric("204220033000")) {RGsetEpic@annotation=aEpic}
  if (RGsetEpic@annotation['array']=="IlluminaHumanMethylation450k") {RGsetEpic@annotation=a450k}
  RGset = RGsetEpic
  sampleID=paste0(samplename_data)
  FFPE = NULL
  output_dir = getwd()
  sample = 1
  output_file <- paste0(sampleID,".html")
  message(paste("Now running:", samplename_data, run_id, barcode, pathEpic))
  rmarkdown::render(
    #"/Volumes/CBioinformatics/Methylation/in_house/mnp.v116/mnp.v11b6/inst/report.Rmd",
                     #rmarkdown::render( "/Volumes/CBioinformatics/Methylation/in_house/mnp.v116/mnp.v11b6/inst/report.Rmd",
                     system.file("report.Rmd",package="mnp.v11b6"),
    #system.file('report.Rmd', package = "mnp.v11b6"),
    output_dir = output_dir,
    output_format = "html_document",
    output_file = output_file,
    clean = T
  )

  }else{message("your data is null")}
}


quickReport <- function(selectRDs, sh=NULL, getAll=F){
  if(is.null(selectRDs)){return(NULL)}
  if(is.null(sh)){sh<-"samplesheet.csv"}
  allRd <- as.data.frame(read.csv(sh))
  if(getAll==T){return(which(!(allRd[,1] %in% selectRDs)))}
  return(which(allRd[,1] %in% selectRDs))
}


# REPORT: Generates Html reports to cwd with samplesheet.csv
#' @param runPath The location of samplesheet.csv and idats
#' @param sheetName name of samplesheet if it is not "samplesheet.csv"
#' @param selectSams vector of integer index of specific samples to run from samplesheet
#' @param genCn also generate the CNV plot as a PNG file
makeReports.v11b6 <-

  function(runPath = NULL,
           sheetName = NULL,
           selectSams = NULL,
           genCn = F,
           skipQC = F,
           email = T,
           cpReport = F,
           redcapUp = T) {

  if (is.null(runPath)) {runPath=getwd()}
  if (is.null(sheetName)) {sheetName="samplesheet.csv"}
  data <- read.csv(sheetName, strip.white=T)
  runID<-paste0(as.character((data$RunID)[1]))
    if (is.null(selectSams)) {
      listedSams <- 1:length(as.character(data$SentrixID_Pos))
    } else {listedSams <- selectSams}
    for (i in listedSams) {
      message("\n-----------\nNow Running ",i, " of ", length(listedSams),"\n-----------")
      do_report(data=data[i, ])
      }
    message("RUN COMPLETE")
    checkRunOutput(runID)

  if(skipQC==F){create.QC.record();generateQCreport()}#Knits the QC RMD file
  if(cpReport==T){file.list <- copy.to.outFolder(clinDrv, runID)}
  if(redcapUp==T){file.list = list.files(path = getwd(), "*.html", full.names = T)
  uploadToRedcap(file.list)
  }
  if(email==T){launchEmailNotify(runID)}
}

automateRun <- function(runID=NULL){
  if (is.null(runID)){runID=paste0(basename(getwd()))}
  setRunDir(runID)
  readSheetWrite()
  get.idats("samplesheet.csv",rsch.idat,clin.idat)
  moveSampleSheet(methDir)
  install.or.load(instNew = F, rmpkg = F)
  makeReports.v11b6(skipQC = F, email=T)
}

startRun <- function(selectRDs=NULL){
  prioritySam <- quickReport(selectRDs) # returns index in samplelist get("quickReport")
  if(stringr::str_detect(runID, "MGDM")){
    #Clinical Run
    makeReports.v11b6(skipQC=T, email=F, cpReport=F, selectSams=prioritySam, redcapUp=F)
    if(!is.null(selectRDs)){
      remainingSam <- quickReport(selectRDs,getAll=T)
      makeReports.v11b6(skipQC=F, email=T, cpReport=T, selectSams=remainingSam, redcapUp=T)
    }
  }else{
    #Research Run
    makeReports.v11b6(skipQC=T, email=F, cpReport=T, selectSams = NULL, redcapUp=T)}
}

defineParams()
