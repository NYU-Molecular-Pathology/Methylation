#!/bin/bash
## Script name: BigPurpleCopy.sh
## Purpose: Copy a file from a local directory to a BigPurple HPC directory
## Date Created: October 23, 2023
## Author: Jonathan Serrano
## Version: 1.0.0
## Copyright (c) NYULH Jonathan Serrano, 2024

INPUTFI_PATH=${1-NULL} # Full path to the file or directory to be copied
OUTDIR=${2-NULL}       # Full path to the destination on HPC

kerbero="$USER"

# Check if OUTDIR ends with a file separator ("/") & append if not
if [[ "$OUTDIR" != */ ]]; then
    OUTDIR="${OUTDIR}/"
fi

# Determine if INPUTFI_PATH is a file or directory
if [[ -f "$INPUTFI_PATH" ]]; then
    COPY_TYPE="file"
elif [[ -d "$INPUTFI_PATH" ]]; then
    COPY_TYPE="directory"
else
    echo "Error: The input '$INPUTFI_PATH' is not a file or directory."
    exit 1
fi

# Input Arguments -------------------------------------------------------------
FIBASE=$(basename "${INPUTFI_PATH}")
if [[ "$COPY_TYPE" == "directory" ]]; then
    HPC_FILE_DIR="${OUTDIR}${FIBASE}"  # Keep folder structure for directories
else
    HPC_FILE_DIR="${OUTDIR}"          # Copy file directly into the destination
fi
NEW_FILE_PATH="$kerbero@bigpurple.nyumc.org:$HPC_FILE_DIR"

# Text Color Vars -------------------------------------------------------------
bold=$(tput bold)         # makes console output text bold
BG_BLUE="$(tput setab 4)" # makes console output text background blue
normal=$(tput sgr0)       # resets default console output text

# FUN: Message a step in color blue -------------------------------------------
message_print() {
    stringMsg=$1
    extraMsg=$2
    echo -e "\n${bold}${BG_BLUE}${stringMsg}:${normal}\n${extraMsg}\n\n"
}

# FUN: Checks if BigPurple directory exists -----------------------------------
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

# Check if BigPurple Directory Exists -----------------------------------------
check_directory "$kerbero" "$OUTDIR"

# Rsync file or directory from local to BigPurple -----------------------------
message_print "Executing the following" "rsync -vrthP -e ssh $INPUTFI_PATH $NEW_FILE_PATH"
rsync -vrthP -e ssh "$INPUTFI_PATH" "$NEW_FILE_PATH"

# Change permissions of the file or directory on BigPurple --------------------
message_print "Changing permissions" "ssh $kerbero@bigpurple.nyumc.org chmod -R g+rwx $HPC_FILE_DIR"
ssh "$kerbero@bigpurple.nyumc.org" "chmod -R g+rwx $HPC_FILE_DIR"
