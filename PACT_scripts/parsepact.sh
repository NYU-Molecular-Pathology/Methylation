#!/bin/bash
## ---------------------------
## Script name: parsepact.sh
## Purpose: Initiate Rscript to parse PACT xlsm worksheet and rsync output to HPC
## Date Created: February 9, 2023
## Author: Jonathan Serrano
## Version: 1.2.0
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

# Hardcoded Variables --------------------------------------------------------------------------------
APITOKEN='8XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
methAPI='5XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
CSVOUTDIR="/gpfs/data/molecpathlab/production/samplesheets/LG-PACT/"
GITHUBLINK="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts/"

# Input Arguments ------------------------------------------------------------------------------------
PACTID=${1-NULL}           # if arg $1 is empty assign NULL as default PACT-22-12
RUN_ID=${2-NULL}           # if arg $2 is empty assign NULL as default 250817_NB501073_0999_ABCD2TBGXM
VAL_KWD=${3-$ILMN}         # if arg $3 is empty assign NULL as default as "-ILMNVAL"
KERBEROS=${4-$USER}        # if arg 4 is empty assign USER as default
MOLECAPI=${5-$API_DEFAULT} # molecular lab REDCap Database API Token default

kerbero="$KERBEROS"
CURR_USER="$USER"

# CHECK: Ensure correct R version path is used -------------------------------------------------------
RMAINDIR=$(which R)
[ "$RMAINDIR" != "/usr/local/bin/R" ] && {
   echo -e "R not loaded from /usr/local/bin/R\nCheck ~/.Rprofile"
   exit 1
}

# Text Color Vars ------------------------------------------------------------------------------------
bold=$(tput bold)        # makes console output text bold
BG_BLU="$(tput setab 4)" # makes console output text background blue
BG_GRN="$(tput setab 2)" # makes console output text background green
FG_YLW="$(tput setaf 3)" # makes console output text yellow color
NORM=$(tput sgr0)        # resets default console output text

set -Eeuo pipefail

# FUN: trap notify ERR -------------------------------------------------------------------------------
function notify {
   echo "Bash script exited!"
   echo "$(caller): ${BASH_COMMAND}"
}

trap notify ERR

# FUN: curl file from link & message -----------------------------------------------------------------
message_curl() {
   linkname=$1
   filename=$2
   printf "\n${FG_YLW}curl -# -L $linkname$filename >$HOME/$filename${NORM}\n"
   curl -# -L "$linkname$filename" >"$HOME/$filename"
   chmod ugo+rwx "$HOME/$filename"
   chmod +rwx "$HOME/$filename"
}

# FUN: Message which command is executing ------------------------------------------------------------
message_exe() {
   exeOrder=$1
   filename=$2
   extraMsg=$3
   printf "\n${BG_GRN}Executing Script #${exeOrder}:${NORM}"
   printf "Rscript $HOME/$filename ${BG_BLU}${extraMsg}${NORM}\n\n"
}

# FUN: Message a step in color blue ------------------------------------------------------------------
message_print() {
   stringMsg=$1
   extraMsg=$2
   printf "\n${bold}${BG_BLU}${stringMsg}:${NORM}\n${extraMsg}\n\n"
}

# FUN: Checks if BigPurple directory exists ----------------------------------------------------------
check_directory() {
   kerbero=$1
   PACTDIR=$2
   if ssh "$kerbero@bigpurple.nyumc.org" "[ ! -d $PACTDIR ]"; then
      printf "Creating new Directory:\n%s" "$PACTDIR"
      printf "\n%sssh \"%s@bigpurple.nyumc.org\" \"mkdir -p %s\"%s" "${BG_BLU}" "$kerbero" "$PACTDIR" "${NORM}"
      ssh "$kerbero@bigpurple.nyumc.org" "mkdir -p $PACTDIR"
   fi
}

# Check Args: Exit if input is missing ---------------------------------------------------------------
[ -z "$PACTID" ] && {
   echo "Missing argument #1: You did not provide a PACTID name (i.e. PACT-22-12)"
   exit 1
}
[ -z "$RUN_ID" ] && {
   echo "Missing argument #2: PACT Run name (i.e. 250817_NB501073_0999_ABCD2TBGXM)"
   exit 1
}

# Go to Home Folder and print Input Args -------------------------------------------------------------
cd "$HOME"
message_print "Input Kerberos ID" "$kerbero"
message_print "Input RUN_ID" "$RUN_ID"
message_print "Input PACTID" "$PACTID"

# CURL: Download latest R scripts --------------------------------------------------------------------
message_curl ${GITHUBURL} "pactParse.R"
message_curl ${GITHUBURL} "PactMethMatch.R"

# EXE 1: R Script for Samplesheet Generation ---------------------------------------------------------
message_exe 1 "pactParse.R" "$MOLECAPI $PACTID $RUN_ID"
Rscript pactParse.R $MOLECAPI "$PACTID" "$RUN_ID" "$VAL_KWD"

NEWCSV="$HOME/Desktop/${RUN_ID}-SampleSheet.csv"
PACTDIR="$CSVOUTDIR${RUN_ID}"
SHEETDIR="$CURR_USER@bigpurple.nyumc.org:$PACTDIR"

NEWCSV_VAL="$HOME/Desktop/${RUN_ID}-${VAL_KWD}-SampleSheet.csv"

if [ -f "${NEWCSV_VAL}" ]; then
   HAS_VALIDATION=true
else
   HAS_VALIDATION=false
fi

# Use the variable later in your script
if $HAS_VALIDATION; then
   echo "Validation File exists."
   PACTDIR_VAL="$CSVOUTDIR${RUN_ID}-${VAL_KWD}"
   SHEETDIR_VAL="$CURR_USER@bigpurple.nyumc.org:$PACTDIR_VAL"
fi

# Check if BigPurple Directory LG-PACT Directory Exists ----------------------------------------------
message_print "SampleSheet is on your Desktop" "$NEWCSV"
check_directory "$CURR_USER" "$PACTDIR"

if $HAS_VALIDATION; then
   check_directory "$CURR_USER" "$PACTDIR_VAL"

fi

# Rsync CSV from Desktop to BigPurple and chmod ------------------------------------------------------
message_print "Executing the following" "rsync -vrthP -e ssh $NEWCSV $SHEETDIR"
rsync -vrthP -e ssh "$NEWCSV" "$SHEETDIR"

if $HAS_VALIDATION; then
   rsync -vrthP -e ssh "$NEWCSV_VAL" "$SHEETDIR_VAL"
fi

message_print "Changing permissions" "ssh $CURR_USER@bigpurple.nyumc.org chmod -R g+rwx $PACTDIR"
ssh "$CURR_USER@bigpurple.nyumc.org" "chmod -R g+rwx $PACTDIR"

if $HAS_VALIDATION; then
   ssh "$CURR_USER@bigpurple.nyumc.org" "chmod -R g+rwx $PACTDIR_VAL"
fi

# CURL: Download latest shell scripts ----------------------------------------------------------------
message_curl ${GITHUBURL} "printPactCommands.sh"
message_curl ${GITHUBURL} "make_consensus.sh"

pactRunID=$(basename "${PACTID%.*}")

message_print "If needed, modify Desktop samplesheet and rsync again" "rsync -vrthP -e ssh $NEWCSV $SHEETDIR"

# EXE 2: HTML commands Generation --------------------------------------------------------------------
message_print "Saving Html File" "$HOME/printPactCommands.sh $RUN_ID ${pactRunID} NULL ${kerbero} >$HOME/${pactRunID}.html && open $HOME/${pactRunID}.html"

"$HOME/printPactCommands.sh" "$RUN_ID" "${pactRunID}" NULL "${kerbero}" >"$HOME/${pactRunID}.html" && open "$HOME/${pactRunID}.html"

if $HAS_VALIDATION; then
   "$HOME/printPactCommands.sh" "$RUN_ID-${VAL_KWD}" "${pactRunID}-${VAL_KWD}" NULL "${kerbero}" >"$HOME/${pactRunID}-${VAL_KWD}.html" && open "$HOME/${pactRunID}-${VAL_KWD}.html"
fi

open "$HOME/Desktop/${RUN_ID}-SampleSheet.csv"
