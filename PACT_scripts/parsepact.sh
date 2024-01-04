#!/bin/bash
## ---------------------------
## Script name: parsepact.sh
## Purpose: Initiate Rscript to parse PACT xlsm worksheet and rsync output to HPC
## Author: Jonathan Serrano
## Version: 1.0.0
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

# Hardcoded Variables --------------------------------------------------------------------------------
APITOKEN='8XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
methAPI='5XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
CSVOUTDIR="/gpfs/data/molecpathlab/production/samplesheets/LG-PACT/"
GITHUBLINK="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts/"

# Input Arguments ------------------------------------------------------------------------------------
pactID=${1-NULL}   # if arg $1 is empty assign NULL as default PACT-22-12
runID=${2-NULL}    # if arg $2 is empty assign NULL as default 250817_NB501773_0999_ABCD2TBGXM
kerbero=${3-$USER} # if arg $3 is empty assign USER as default
kerbero="$kerbero"

# Text Color Vars ------------------------------------------------------------------------------------
bold=$(tput bold)         # makes console output text bold
BG_BLUE="$(tput setab 4)" # makes console output text background blue
BG_GRN="$(tput setab 2)"  # makes console output text background green
FG_YLW="$(tput setaf 3)"  # makes console output text yellow color
normal=$(tput sgr0)       # resets default console output text

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
   printf "\n${FG_YLW}curl -# -L $linkname$filename >$HOME/$filename${normal}\n"
   printf "\n${BG_GRN}Executing:\n$HOME/$filename${normal}\n\n"
   curl -# -L $linkname$filename >$HOME/$filename
   chmod ugo+rwx "$HOME/$filename"
   chmod +rwx "$HOME/$filename"
}

# FUN: Message a step in color blue ------------------------------------------------------------------
message_print() {
   stringMsg=$1
   extraMsg=$2
   printf "\n${bold}${BG_BLUE}${stringMsg}:${normal}\n${extraMsg}\n\n"
}

# FUN: Checks if BigPurple directory exists ----------------------------------------------------------
check_directory() {
   kerbero=$1
   PACTDIR=$2
   if ssh "$kerbero@bigpurple.nyumc.org" "[ ! -d $PACTDIR ]"; then
      printf "Creating new Directory:\n$PACTDIR"
      echo "ssh \"$kerbero@bigpurple.nyumc.org\" \"mkdir -p $PACTDIR\""
      ssh "$kerbero@bigpurple.nyumc.org" "mkdir -p $PACTDIR"
   fi
}

# Print helpFunction in case parameter is empty ------------------------------------------------------
if [ -z "$pactID" ]; then
   echo "Missing argument #1: You did not provide a pactID name (i.e. PACT-22-12)"
   exit 1
fi
if [ -z "$runID" ]; then
   echo "Missing argument #2: PACT Run name (i.e. 250817_NB501073_0999_ABCD2TBGXM)"
   exit 1
fi

# Go to Home Folder and print out inputs, download R scripts -----------------------------------------
cd $HOME
message_print "Input Kerberos ID" "$kerbero"
message_print "Input runID" "$runID"
message_print "Input pactID" "$pactID"
message_curl ${GITHUBLINK} "pactParse.R"
message_curl ${GITHUBLINK} "PactMethMatch.R"

# EXE 1: R Script for Samplesheet Generation ---------------------------------------------------------
Rscript pactParse.R $APITOKEN "$pactID" "$runID"

NEWCSV=$(ls -t ~/Desktop/*-SampleSheet.csv | awk '{printf($0);exit}')
CSVBASE=$(basename ${NEWCSV%%-*})
PACTDIR="$CSVOUTDIR$CSVBASE"
SHEETDIR="$kerbero@bigpurple.nyumc.org:$PACTDIR"

# Rsync CSV from Desktop to BigPurple and chmod ------------------------------------------------------
message_print "SampleSheet is on your Desktop" "$NEWCSV"
message_print "Copying the file" "$NEWCSV\nTo: $SHEETDIR"
check_directory "$kerbero" "$PACTDIR"
message_print "Executing the following" "rsync -e ssh $NEWCSV $SHEETDIR"
rsync -e ssh "$NEWCSV" "$SHEETDIR"
message_print "Changing permissions" "ssh $kerbero@bigpurple.nyumc.org chmod -R g+rwx $PACTDIR"
ssh "$kerbero@bigpurple.nyumc.org" "chmod -R g+rwx $PACTDIR"

# EXE 2: R Script for MethylMatch Generation ---------------------------------------------------------
message_print "Executing Rscript" "Rscript PactMethMatch.R $methAPI $pactID"
Rscript PactMethMatch.R $methAPI "$pactID"

# EXE 3: HTML commands Generation --------------------------------------------------------------------
message_curl ${GITHUBLINK} "printPactCommands.sh"
message_curl ${GITHUBLINK} "make_consensus.sh"
message_print "Saving Html File" "$HOME/printPactCommands.sh $runID ${pactID} NULL ${kerbero} >$HOME/${pactID}.html && open $HOME/${pactID}.html"
$HOME/printPactCommands.sh $runID ${pactID} NULL ${kerbero} >$HOME/${pactID}.html && open $HOME/${pactID}.html
