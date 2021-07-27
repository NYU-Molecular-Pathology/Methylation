#!/usr/bin/env Rscript
options(Ncpus = 6)
# Cran Packages ----
cranPkgs <-
  c(
    'devtools',
    'remotes',
    'sjmisc',
    'stringi',
    'digest',
    'RCurl',
    'rlang',
    'parallel',
    'grid',
    'gridExtra',
    'knitr',
    'kableExtra',
    'ggplot2',
    'plotly',
    'ggfortify',
    'ggrepel',
    'gplots',
    'fastmatch',
    'pals',
    'Polychrome',
    'qdapTools',
    'beepr',
    'xtable',
    'pander',
    'grDevices',
    'graphics',
    'stats',
    'utils',
    'magick',
    'ade4',
    "MASS",
    "R.utils",
    "optparse",
    "targets",
    "usethis",
    "webshot",
    "reshape2",
    "data.table",
    "DT",
    "scales",
    "RColorBrewer",
    "readxl",
    "stringr",
    "tinytex"
  )

# GitHub Packages ----
gHubPkgs <-
  data.frame(
    docstring = 'dasonk/docstring',
    rstudioapi = 'rstudio/rstudioapi',
    easypackages = 'jakesherman/easypackages',
    redcapAPI = 'nutterb/redcapAPI',
    crayon = "r-lib/crayon",
    redcap = "epicentre-msf/redcap"
  )

# BioConductor Packages ----
biocPkgs <-
  c(
    'HDF5Array', 'rngtools', 'bumphunter','GEOquery', 'minfi', 'lumi', 'methylumi',
    'randomForest', 'glmnet','IlluminaHumanMethylation450kmanifest',
    'IlluminaHumanMethylation450kanno.ilmn12.hg19', 'Rtsne',
    'IlluminaHumanMethylationEPICmanifest', 'IlluminaHumanMethylationEPICanno.ilm10b2.hg19',
    'IlluminaHumanMethylationEPICanno.ilm10b4.hg19', 'MethylAid', 'conumee','BiocParallel',
    "Biobase","limma"
    )

# Extra Libraries ----
easyPkgs <- c('tidyverse','sjmisc','stringi','digest','RCurl','gridExtra','needs')

# Helper Functions ----
sw <- function(pkgOb){try(return(suppressMessages(suppressWarnings(pkgOb))),silent=T)}
rq <- function(pkg){return(sw(!require(pkg, character.only=T)))}
ld <- function(libName){sw(suppressPackageStartupMessages(library(libName,character.only=T)))}
up <- function(){update.packages(repos='http://cran.rstudio.com/',type = "source", ask=F, checkBuilt=T)}

#' pk.inst will load a package library or install it if it does not exist
#'
#' @param pkg a character string package name
#' @param ... additional parameters passed to install.packages()
#'
#' @return Message if loading the library was successful
pk.inst <- function(pkg){
    pk.opt <- list(pkg,dependencies=T,verbose=T,Ncpus = 6)
    message("Checking ", pkg, "...")
    tryCatch(
        expr={if(rq(pkg)){
            do.call(install.packages, c(pk.opt,list(type="both")))
        }else{message(pkg," ...load successful")}
        },
        warning=function(cond) {message("Warning caught:", cond, "\n-----------")
          do.call(install.packages, c(pk.opt,list(type="source")))
        },
        error=function(cond) {message("Error1 caught: ", cond, "\n-----------")
            do.call(install.packages, c(pk.opt,list(type="binary")))
        },
        custom_error=function(cond){message("Error2 caught: ", cond, "\n-----------")
          easypackages::packages(pkg,prompt=F, Ncpus = 6);message(cond)}
    )
}
# FUN: Installs package from github link
gh.inst <- function(pkNam, ghLnk,...){
    message("Checking ", pkNam, "...")
    gh.opt=list(ghLnk, dependencies=T, upgrade="never",...)
    tryCatch(
        expr={if(rq(pkNam)){
            params <- c(gh.opt,list(type="source"))
            do.call(devtools::install_github,params)
        }else{message(pkNam," ...load successful")}
        },
        error=function(cond){
            params  <- c(gh.opt,list(type="binary"))
            do.call(devtools::install_github, params)
        }
    )
}

# FUN: Installs package from Source link
srcInst <- function(fn,...){
  message("Checking ", basename(fn), "...")
  #downLoc = file.path("~/Desktop",basename(fn))
  #download.file(url=fn,destfile=downLoc)
    tryCatch(
        expr={install.packages(fn, repos=NULL, type="source", Ncpus = 6)},
        error=function(cond){
            install.packages(fn,repos=NULL, method="libcurl", type = "source")
          },
        warning=function(cond){
          install.packages(fn,repos=NULL, method="auto", type = "source")
          }
        )
}

# FUN: Installs package from Bioconductor
bioInst <- function(pknm){
    message("Checking ", pknm, "...")
    bio.opt <- list(pkgs=pknm, update=T, ask=F)
    tryCatch(
        expr={if(rq(pknm)){sw(do.call(BiocManager::install, c(bio.opt)))
        }else{message(pknm," ...load successful")}},
        error=function(cond){message("Package already loaded")}
    )
}

fixProf <- function(){
    txt1 <- "^[:blank:]*autoload\\(\"needs\", \"needs\"\\)"
    txt2 <- "\n\nautoload(\"needs\", \"needs\")\n\n"
    siteProf <- if(is.na(Sys.getenv("R_PROFILE", unset=NA))) {
        file.path(Sys.getenv("R_HOME"),"etc","Rprofile.site")
    } else {Sys.getenv("R_PROFILE")}
    if (!file.exists(siteProf)) {file.create(siteProf)}
    cxn <- file(siteProf); lines <- readLines(cxn)
    if (!any(grepl(txt1,lines))) {write(txt2, file=siteProf, append=T)}
    close(cxn)
}

fixNeeds <- function (){
    sysfile <- system.file("extdata", "promptUser", package="needs")
    write(0, file=sysfile); options(needs.promptUser=FALSE); invisible(needs:::autoload(TRUE))
}

checkNeeds <- function(){
    tryCatch(
        expr={
            if(!("needs" %in% rownames(installed.packages()))){
                install.packages("needs",dependencies=T,verbose=T, Ncpus = 6)
                fixNeeds();fixProf()
            }else{
                fixNeeds();fixProf()
            }
        },
        error=function(cond){
            devtools::install_github("joshkatz/needs", ref = "development",dependencies=T,verbose=T,upgrade="always")
            fixNeeds();fixProf()
        },
        warning=function(cond){
            devtools::install_github("joshkatz/needs", ref = "development",dependencies=T,verbose=T,upgrade="always")
            fixNeeds();fixProf()
        }
    )
}

setOptions <- function(){
    options(needs.promptUser=F); Sys.setenv(R_ENABLE_JIT=T)
    options("needs.promptUser"=F); options("promptUser"=F)
    options("device.ask.default" = FALSE)
    options("install.packages.compile.from.source"="Yes")
    options("install.packages.check.source"="yes")
    compiler::enableJIT(3); compiler::compilePKGS(enable=T);
    compiler::setCompilerOptions(suppressAll=T, optimize=3)
}

cpMakeV <- function(editFile=F){
    mkFi <- "/Volumes/CBioinformatics/jonathan/Rprojects/Methylation-scripts/Makevars"
    message("Copying Makevars to local dotfiles: ", mkFi)
    if(editFile){file.edit(mkFi)}
    if(file.exists(mkFi)){file.copy(mkFi,"~/.R/", overwrite=T )}
}

# FUN: Load Classifier Versions -
install.or.load <- function(pathtoFile=NULL, instNew=T, rmpkg=F) {
    cbio="/Volumes/CBioinformatics"
    if(is.null(pathtoFile)){
        mnp.pk.loc = file.path(cbio, "Methylation/in_house/mnp.v116/mnp.v11b6")
    }else{mnp.pk.loc<-pathtoFile}
    if(rmpkg){remove.packages("mnp.v11b6", lib=.libPaths()[[1]][1])}
    if(instNew){
        install.packages(mnp.pk.loc, repos=NULL, type="source", force=T, Ncpus = 6)
        install.or.load(instNew=F)
    }else{ld("mnp.v11b6")}
}

# Loads default packages or custom if input provided
#' @return package installs or loads
#' @export
loadPacks <- function(pkgs=cranPkgs, ezLibs=easyPkgs, ghPk=gHubPkgs, bioPks=biocPkgs) {
    library("base")
    rlis = getOption("repos")
    rlis["CRAN"] = "http://cran.us.r-project.org"
    options(repos = rlis)
    tryCatch(
        expr={
          setOptions(); fixProf()
            if (rq("devtools")) {install.packages("devtools",dependencies=T,verbose=T)}else{library("devtools")}
            checkNeeds()
            gh.inst("easypackages","jakesherman/easypackages")
            if (rq("tidyverse")) {pk.inst("tidyverse")}else{ld("tidyverse")}
            easypackages::packages("parallel","doSNOW","doParallel", "foreach","compiler", prompt=F)
            if(suppressWarnings(!require("BiocManager"))){
              install.packages("BiocManager", Ncpus = 6)
              BiocManager::install(version="3.10",update=T, ask=F, type="source")
            }
            #if(BiocManager::version()!='3.13'){BiocManager::install(version="3.13",update=T, ask=F, type="source")}
            if(rq("zip")){install.packages("zip", dependencies=T, type="binary")}
            invisible(lapply(1:length(pkgs),function(x){pk.inst(pkgs[x])}))
            invisible(lapply(names(ghPk),function(x){gh.inst(x,ghPk[1,x])}))
            invisible(lapply(1:length(bioPks),function(x){bioInst(bioPks[x])}))
            # Loads libraries
            ready <- ezLibs %in% sessioninfo::package_info()[["package"]]
            if(any(!ready)){lapply(ezLibs[!ready], require, character.only=TRUE)}
            fldx <- function(nmn = "foldAll") {return(rstudioapi::executeCommand(nmn))}
            assign(x = "fldx", value = fldx, envir = .GlobalEnv)
        },
        error= function(cond){
            message("\n~~~You encountered the following error during install:\n", cond,"\n")
            message("If there is a compile error, try running fixCompiles() and then try loadPacks() again")
        }
    )

}

fixCompiles <- function(){
    #system("/bin/bash -c $(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)")
    cpMakeV()
    system("brew install llvm gcc")
    system("defaults write org.R-project.R NSAppSleepDisabled -bool YES")
    system('echo "export PATH="/usr/local/sbin:$PATH"" >> ~/.zshrc >')
    system('echo "export PATH="/usr/local/sbin:$PATH"" >> ~/.bashrc >')
    system("export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES")
    system('echo "export OBJC_DISABLE_INITIALIZE_FORK_SAFETY="YES"" >> ~/.zshrc')
    system('echo "export OBJC_DISABLE_INITIALIZE_FORK_SAFETY="YES"" >> ~/.bashrc')
    system("xcode-select --install")
    system("xcodebuild -runFirstLaunch")
    system("brew install dpkg openssl libffi go python")
    system('export PATH="/usr/local/opt/python/libexec/bin:$PATH"')
    system('export GOROOT=/usr/local/opt/go/libexec')
    system("export PATH=$PATH:$GOROOT/bin")
    system("brew install cmake gsl gmp libopenmpt open-mpi")
    # export RSTUDIO_WHICH_R=/usr/local/bin/R
    system("brew install aspell gdal autoconf automake bash cairo cmake coreutils dos2unix exa fd fontconfig freetype fribidi fzf gawk gcc gdbm gdk-pixbuf gettext giflib git glib gmp gnu-sed gnu-tar gnutls graphite2 harfbuzz htop icu4c ilmbase imagemagick isl jpeg libcroco libde265 libevent libffi libgcrypt libgpg-error libheif libidn2 libmpc libomp libpng librsvg libtasn1 libtermkey libtiff libtool libunistring libuv libvterm libyaml little-cms2 llvm lua lua@5.1 luajit luarocks lzo mpfr msgpack multimarkdown ncurses neovim nettle nspr nss oniguruma openconnect openexr openjpeg openssl@1.1 p11-kit pandoc pandoc-citeproc pandoc-crossref pango pcre pcre2 pdfgrep perl pixman pkg-config poppler qt readline rename ripgrep ruby shared-mime-info shellcheck sqlite stoken swig texinfo tmux tree unbound unibilium vim webp wget x265")
    system('export CFLAGS+=-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
  export CCFLAGS+=-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
  export CXXFLAGS+=-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
  export CPPFLAGS+=-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk')
    system("brew install golang libstfl")
    if(dir.exists("/Library/Developer/CommandLineTools")){
        if(length(list.dirs("/Library/Developer/CommandLineTools", recursive = F))>0){
            system("ln -s /Applications/Xcode.app/Contents/Developer /Library/Developer/CommandLineTools")}}
    system("brew unlink gcc && brew link gcc --force")
    system('brew link --overwrite gcc')
    system('brew doctor')
}

startmsg <- function(){
  cbio = "/Volumes/CBioinformatics/"; zdriv = "/Volumes/molecular/Molecular"
  wmm = "You do not have this path mounted:\n"
  if (!dir.exists(cbio)) {warning(wmm, cbio)}
  if (!dir.exists(zdriv)) {warning(wmm, zdriv)}
  message("You have the following drives mounted:")
  system("ls /Volumes")
  stopifnot(dir.exists(cbio) | dir.exists(zdriv))
  #system("gcc --version")
}

colorMsg <- function(){
  mkred <- function(strMsg) {return(crayon::white$bgRed$bold(strMsg))}
  mkblu <- function(strMsg) {return(crayon::white$bgBlue$bold(strMsg))}
  ms1 <- paste0(
    crayon::white$bgGreen("Your mnp.v11b6 package is installed and loading"),"\n",
    crayon::white$bgGreen("To update in-house classifier to current version run:"),"\n",
    mkblu("install.or.load(pathtoFile=NULL, instNew=F, rmpkg=F)"),"\n"
  )
  ms2 <- paste0(
    mkred("mnp.v11b6 package is not installed executing the function to install:"),"\n",
    mkblu("install.or.load(instNew=T)"),"\n")
  return(c(ms1,ms2))
}

# Load all Functions ---------------------
startLoadingAll <- function() {
    sexEst = "https://github.com/jungch/sest/raw/master/sest.tar"
    mgmtLn = "https://git.io/JWKTo"
    cbioLn = "/Volumes/CBioinformatics/Methylation/"
    mnqDir = file.path(cbioLn,"in_house/mnp.v116/mnpqc_0.1.0.tar.gz")
    startmsg(); loadPacks()
    ms <- colorMsg()
    gh.inst("rmarkdown", 'rstudio/rmarkdown')
    if (rq("mgmtstp27")) {
      dLoc <- "~/Desktop/temp.tar.gz"
      download.file(url = mgmtLn, destfile = dLoc)
      sw(install.packages(dLoc, repos = NULL, type = "source"))
    }
    if (rq("sest")) {sw(srcInst(sexEst))}
    if (rq("mnpqc")) {install.packages(mnqDir, repos = NULL, type="source", verbose=T)}
    if (rq("mnp.v11b6")) {
      cat(ms[2]);install.or.load(instNew = T)}else{
      cat(ms[1]);install.or.load(instNew = F)
    }
}
startLoadingAll()
