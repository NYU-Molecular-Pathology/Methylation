#!/bin/bash
## ---------------------------
## Script name: runMeth.sh
## Purpose: Download methylation pipeline scripts from Github with curl and execute Rscripts
## Date Created: August 5, 2021
## Date Last Modified: February 1, 2024
## Author: Jonathan Serrano
## Version: 1.0.0
## Copyright (c) NYULH Jonathan Serrano, 2024
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
   echo -e "\n%sYour system is not using R located in:%s $R_LOCAL" "${BG_RED}" "${NORMAL}"
   echo -e "\nInstead R is located in: $R_CURR"
   exit 0
fi

# Check if system requirements installed
function checkSystem {
   APPNAME="$1"
   if ! command -v "$APPNAME" &> /dev/null; then
      echo -e "\n${BG_RED}You did not install $APPNAME!${NORMAL} Stopping script\n"
      exit 0
   fi

}

unameOut="$(uname -s)"

if [[ "$unameOut" == "Darwin" ]]; then
   checkSystem "latex"
   checkSystem "pdflatex"
   checkSystem "xquartz"
   checkSystem "brew"
   checkSystem "pandoc"
   checkSystem "go"
fi

# CHECK: Prints instruction if methRun parameter is empty ------------------------------------------------------
[ -z "$methRun" ] && { echo "${BG_RED}Missing Input: Methylation Run ID ${NORMAL}. Ex. '21-MGDM30'"; exit 1; }
[ -z "$redcapUp" ] && { redcapUp="T"; }

# SANITY CHECK: Prints parameters ------------------------------------------------------------------------------
echo "~~~~~~~~~~~~~~~~~~~~~Parameters passed to RScript~~~~~~~~~~~~~~~~~~~~~"
echo "${BG_BLU}methAPI:${NORMAL} $methAPI"
echo "${BG_CYA}(arg1) methRun:${NORMAL} $methRun"
echo "${BG_MAG}(arg2) PRIORITY:${NORMAL} $PRIORITY"
echo "${BG_GRN}(arg3) runPath:${NORMAL} $runPath"
echo "${BG_BLU}(arg4) redcapUp:${NORMAL} $redcapUp"

# FUNCTIONS: Curl download and mkdir if not exists -------------------------------------------------------------
message_curl() {
   linkname=$1
   filename=$2
   echo -e "${BG_GRN}Downloading file from Github:${NORMAL}\n$HOME/$filename\n"
   curl -k -# -L "$linkname$filename" >"$HOME/$filename"
   chmod +rwx "$HOME/$filename"
}

check_directory() {
   CURRDIR=$1
   echo -e "${BG_GRN}Creating Directory:${NORMAL}\nmkdir -p $CURRDIR\n"
   mkdir -p "$CURRDIR"
}

message_red() {
   CMD_CODE=$1
   echo -e "${BG_RED}${CMD_CODE}${NORMAL}\n"
}

message_curl ${GITHUBMAIN} "report.Rmd"
message_curl ${GITHUBMAIN} "Methyl_QC.Rmd"
message_curl ${GITHUBMAIN} "methylExpress.R"
message_curl ${GITHUBMAIN} "PostRunSarc.R"

printf "\nExecuting Rscript with parameters input:\n"
color_params="${BG_MAG}$methAPI${NORMAL} ${BG_CYA}$methRun${NORMAL} ${BG_GRN}$PRIORITY${NORMAL} ${BG_BLU}$runPath${NORMAL} ${BG_MAG}$redcapUp${NORMAL} ${BG_GRN}$runLocal${NORMAL}\n"

echo -e "${BG_BLU}${bold}Rscript --verbose $HOME/methylExpress.R${NORMAL} $color_params"

# EXECUTE: Main command using user input  ----------------------------------------------------------------------
Rscript --verbose "$HOME/methylExpress.R" $methAPI $methRun $PRIORITY $runPath $redcapUp $runLocal
Rscript --verbose "$HOME/PostRunSarc.R" $methAPI $methRun

# ASSIGN: Default output paths ---------------------------------------------------------------------------------
[ -n "$runPath" ] && runPath="/Volumes/CBioinformatics/Methylation/Clinical_Runs/${methRun}/"

LOCALDIR="$HOME/Desktop/Html_${methRun}/${methRun}/"

RUNYEAR=${methRun%%-*}
MOLECDIR="/Volumes/molecular/MOLECULAR LAB ONLY/NYU-METHYLATION/Results/20${RUNYEAR}/"
MOLECOUT="/Volumes/molecular/Molecular/MethylationClassifier/20${RUNYEAR}/"

# SANITY CHECK: Prints output directories ----------------------------------------------------------------------
echo -e "\n----Output Directories----\n"
echo "${BG_GRN}runPath:${NORMAL}  $runPath"
echo "${BG_GRN}MOLECDIR:${NORMAL} $MOLECDIR"
echo "${BG_GRN}LOCALDIR:${NORMAL} $LOCALDIR"

copy_htmls() {
   COPY_DIR=$1
   OUTPUT_QC=$2
   rsync -avWP --progress --include '*.html' --include='*/' --exclude='*' "${LOCALDIR}" "$COPY_DIR"
   rsync -vWPh --progress "$OUTPUT_QC" "$COPY_DIR"
}

# COPY: Rsync html to HOME directory then cp from HOME to Z-Drive -------------------------------------------
check_directory "$LOCALDIR"

echo -e "\nRsyncing ${methRun} html files to Home directory:\n"
echo -e "${BG_GRN}rsync -vrthP --include '*.html' --exclude '*' '${runPath}' '$LOCALDIR'${NORMAL}\n"

rsync -av --include '*.html' --exclude '*' "${runPath}" "$LOCALDIR"
rsync -av "${runPath}${methRun}_QC_and_Classifier_Scores.csv" "$LOCALDIR"

OUTPUT_DIR_1="${MOLECDIR}${methRun}/"
OUTPUT_DIR_2="${MOLECOUT}${methRun}/"
OUTPUT_QC_FI="${LOCALDIR}${methRun}_QC_and_Classifier_Scores.csv"

check_directory "${OUTPUT_DIR_1}"
check_directory "${OUTPUT_DIR_2}"

HTMLCOPY="rsync -avWP --progress --include '*.html' --include='*/' --exclude='*'"

echo -e "\n${BG_BLU}Ensure output files are copied from your HOME to Z-drive:${NORMAL}\n"
message_red "${HTMLCOPY} \"${LOCALDIR}\" \"${OUTPUT_DIR_1}\""
message_red "rsync -vWPh --progress \"${OUTPUT_QC_FI}\" \"${OUTPUT_DIR_1}\""

message_red "${HTMLCOPY} \"${LOCALDIR}\" \"${OUTPUT_DIR_2}\""
message_red "rsync -vWPh --progress \"${OUTPUT_QC_FI}\" \"${OUTPUT_DIR_2}\""

copy_htmls "${OUTPUT_DIR_1}" "${OUTPUT_QC_FI}"
copy_htmls "${OUTPUT_DIR_2}" "${OUTPUT_QC_FI}"

echo -e "\n${BG_BLU}Verify files are copied here:${NORMAL}\n${MOLECDIR}\n"
open "${MOLECDIR}"
echo -e "\n${BG_GRN}METHYLATION RUN ${methRun} COMPLETE${NORMAL}\n"
