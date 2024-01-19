#!/bin/bash
## ---------------------------
## Script name: BigPurpleCopy.sh
## Purpose: Copy a file from a local directory to a BigPurple HPC directory
## Date Created: October 23, 2023
## Date Last Modified: January 19, 2024
## Author: Jonathan Serrano
## Version: 1.0.0
## Copyright (c) NYULH Jonathan Serrano, 2024
## ---------------------------

INPUTFI_PATH=${1-NULL} # Full path to the file to be copied (i.e. /Users/myName/Desktop/myfile.csv)
OUTDIR=${2-NULL}       # Full path to the directory on HPC where the file should be copied
KERBEROS=${3-$USER}    # if arg $3 is empty assign USER as default

kerbero="$KERBEROS" # force console input as string

# Check if OUTDIR ends with a file separator slash ("/"). Append a slash if it doesn't.
if [[ "$OUTDIR" != */ ]]; then
   OUTDIR="${OUTDIR}/"
fi

# Input Arguments ------------------------------------------------------------------------------------
FIBASE=$(basename "${INPUTFI_PATH}")
HPC_FILE_DIR="${OUTDIR}${FIBASE}"
NEW_FILE_PATH="$kerbero@bigpurple.nyumc.org:$HPC_FILE_DIR"

# Text Color Vars ------------------------------------------------------------------------------------
bold=$(tput bold)         # makes console output text bold
BG_BLUE="$(tput setab 4)" # makes console output text background blue
normal=$(tput sgr0)       # resets default console output text

# FUN: Message a step in color blue ------------------------------------------------------------------
message_print() {
   stringMsg=$1
   extraMsg=$2
   printf "\n${bold}${BG_BLUE}${stringMsg}:${normal}\n${extraMsg}\n\n"
}

# FUN: Checks if BigPurple directory exists ----------------------------------------------------------
echo_hpc_command() {
    OUTDIR=$1
    kerbero=$2
    hpc_command='ssh "'${kerbero}@bigpurple.nyumc.org'" "mkdir -p '${OUTDIR}'"'
    echo -e "\nCreating new Directory:\n${OUTDIR}"
    echo -e "\n${BG_BLUE}${hpc_command}${normal}\n"
}

check_directory() {
    kerbero=$1
    OUTDIR=$2
    if ssh "$kerbero@bigpurple.nyumc.org" "[ ! -d $OUTDIR ]"; then
        echo_hpc_command "$OUTDIR" "$kerbero"
        ssh "$kerbero@bigpurple.nyumc.org" "mkdir -p $OUTDIR"
    fi
}

# Check if BigPurple Directory Directory Exists ----------------------------------------------
check_directory "$kerbero" "$OUTDIR"

# Rsync file from local directory to BigPurple ------------------------------------------------------
message_print "Executing the following" "rsync -vrthP -e ssh $INPUTFI_PATH $NEW_FILE_PATH"
rsync -vrthP -e ssh "$INPUTFI_PATH" "$NEW_FILE_PATH"

# Change permissions of the file on BigPurple ------------------------------------------------------
message_print "Changing permissions" "ssh $kerbero@bigpurple.nyumc.org chmod -R g+rwx $HPC_FILE_DIR"
ssh "$kerbero@bigpurple.nyumc.org" "chmod -R g+rwx $HPC_FILE_DIR"
