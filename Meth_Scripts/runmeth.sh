#!/bin/bash
## ---------------------------
## Script name: runmeth.sh
## Purpose: cUrl source and execute Rscripts locally for clinical methylation report generation
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------

# RedCap API Token: Saved Locally (Hardcoded)
methAPI='XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' 

# INPUT ARGS: if args are empty assign NULL as default below ---------------------------------------------------
methRun=${1-NULL}  # methylation run id e.g. 22-MGDM17
PRIORITY=${2-NULL} # string of prioritized RD-numbers
runPath=${3-NULL}  # any custom directory to copy/run the idat files
redcapUp=${4-NULL} # to upload to redcap or not if server down single char i.e. "T" or "F"
runLocal=${5-NULL} # If the run directory should be executed without shared drives locally i.e. "T" or "F"

# HARDCODED GITHUB PATHS ---------------------------------------------------------------------------------------
GITHUBLINK="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/dockerScripts/"
GITHUBMAIN="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/"

# Hardcoded Terminal Color Variables ---------------------------------------------------------------------------
BG_BLU="$(tput setab 4)" # makes text background blue
BG_CYA="$(tput setab 6)" # makes text background cyan
BG_MAG="$(tput setab 5)" # makes text background magenta
BG_GRN="$(tput setab 2)" # makes text background green
BG_RED="$(tput setab 1)" # makes text background green
bold=$(tput bold)        # makes text bold
NORMAL=$(tput sgr0)      # resets default text

cd "$HOME"

# FUNCTION: EXIT when any command fails ------------------------------------------------------------------------
error_handled=false

# Function to exit when any command fails
trap_error() {
   local last_command=$1
   local exit_code=$2
   local line_no=$3

   if [[ $exit_code -ne 0 && $error_handled == false ]]; then
      printf "\n%sERROR on LAST COMMAND%s:\n%s\n" "${BG_RED}" "${NORMAL}" "${last_command}"
      printf "Command failed at line %s with exit code %s.\n" "${line_no}" "${exit_code}"
      error_handled=true
   fi
}

# Enable strict mode
set -Eeuo pipefail

trap 'last_command=$BASH_COMMAND; last_line=$LINENO' DEBUG
trap 'trap_error "${BASH_COMMAND}" $? "${LINENO}"' ERR EXIT

R_LOCAL="/usr/local/bin/R"
R_CURR="$(which R)"

if [[ "$R_CURR" != "$R_LOCAL" ]]; then
   printf "\n%sYour system is not using R located in:%s $R_LOCAL" "${BG_RED}" "${NORMAL}"
   printf "\nInstead R is located in: $R_CURR"
   exit 0
fi

function checkSystem {
   APPNAME="$1"
   if ! command -v "$APPNAME" &> /dev/null; then
      printf "\n${BG_RED}You did not install $APPNAME!${NORMAL} Stopping script\n"
      exit 0
   fi

}

checkSystem "latex"
checkSystem "pdflatex"
checkSystem "xquartz"
checkSystem "brew"
checkSystem "pandoc"
checkSystem "go"

# CHECK: Prints instruction if methRun parameter is empty ------------------------------------------------------
if [ -z "$methRun" ]; then
   echo "${BG_RED}You did not provide a Methylation Run ID name${NORMAL}. Ex. '21-MGDM30'"
   exit 1
fi

if [ -z "$redcapUp" ]; then
   redcapUp="T"
fi

# SANITY CHECK: Prints parameters ------------------------------------------------------------------------------
echo "~~~~~~~~~~~~~~~~~~~~~Parameters passed to RScript~~~~~~~~~~~~~~~~~~~~~"
echo "${BG_BLU}methAPI:${NORMAL} $methAPI"
echo "${BG_CYA}(arg1) methRun:${NORMAL} $methRun"
echo "${BG_MAG}(arg2) PRIORITY samples:${NORMAL} $PRIORITY"
echo "${BG_GRN}(arg3) Custom runPath:${NORMAL} $runPath"
echo "${BG_BLU}(arg4) To upload redcapUp:${NORMAL} $redcapUp"
printf "\nCurl files downloaded:\n"

# FUNCTIONS: Curl download and mkdir if not exists -------------------------------------------------------------
message_curl() {
   linkname=$1
   filename=$2
   printf "${BG_GRN}$HOME/$filename${NORMAL}\n"
   curl -k -# -L "$linkname$filename" >"$HOME/$filename"
}

check_directory() {
   CURRDIR=$1
   if [ ! -d "$CURRDIR" ]; then
      printf "\nCreating new Directory:\n$CURRDIR"
      mkdir -p "$CURRDIR"
   fi
}

message_curl ${GITHUBMAIN} "report.Rmd"
message_curl ${GITHUBLINK} "Methyl_QC.Rmd"
message_curl ${GITHUBMAIN} "methylExpress.R"

printf "\nExecuting Rscript with parameters input:\n"
color_params="${BG_MAG}$methAPI${NORMAL} ${BG_CYA}$methRun${NORMAL} ${BG_GRN}$PRIORITY${NORMAL} ${BG_BLU}$runPath${NORMAL} ${BG_MAG}$redcapUp${NORMAL} ${BG_GRN}$runLocal${NORMAL}\n"
printf "${BG_BLU}${bold}Rscript --verbose $HOME/methylExpress.R${NORMAL} $color_params"

# EXECUTE: Main command using user input  ----------------------------------------------------------------------
Rscript --verbose "$HOME/methylExpress.R" $methAPI $methRun $PRIORITY $runPath $redcapUp $runLocal

# ASSIGN: Default output paths ---------------------------------------------------------------------------------
if [ -n "$runPath" ]; then
   DEFAULTPATH="/Volumes/CBioinformatics/Methylation/Clinical_Runs/${methRun}/"
   runPath=$DEFAULTPATH
fi

LOCALDIR="$HOME/${methRun}"
RUNYEAR=${methRun%%-*}
MOLECDIR="/Volumes/molecular/MOLECULAR LAB ONLY/NYU-METHYLATION/Results/20${RUNYEAR}"
MOLECOUT="/Volumes/molecular/Molecular/MethylationClassifier/20${RUNYEAR}"

# SANITY CHECK: Prints output directories ----------------------------------------------------------------------
echo "${BG_GRN}runPath:${NORMAL} $runPath"
echo "${BG_GRN}MOLECDIR:${NORMAL} $MOLECDIR"
echo "${BG_GRN}LOCALDIR:${NORMAL} $LOCALDIR"

# EXECUTE: Rsync html to HOME directory then cp from HOME to Z-Drive -------------------------------------------
check_directory "$LOCALDIR"

printf "\n${BG_BLU}Rsyncing files:${NORMAL}\n"
rsync -vrthP --include '*.html' --exclude '*' "${runPath}" "$LOCALDIR"

printf "\nExecuting copy command:\n"
printf "${BG_GRN}cp -RvX ${LOCALDIR} ${MOLECDIR}${NORMAL}\n"

printf "${BG_BLU}Copying files to Z-drive from HOME:${NORMAL}\n"
cp -RvfX "$LOCALDIR" "$MOLECDIR"
cp -RvfX "$LOCALDIR" "$MOLECOUT"

echo "METHYLATION RUN COMPLETE"
