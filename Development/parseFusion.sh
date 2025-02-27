#!/bin/bash

# Input Args --------------------------------------------------------------------------------------------------------------
SHEETPATH=${1-$NULL} # if arg $1 is empty will be NULL i.e. /Users/Jonathan/Desktop/22-MGFS25.xlsx
kerbero=${2-$USER}   # OPTIONAL if arg $2 is empty assign $USER as default else i.e. whoami kerberosid

# Print helpFunction in case parameter is empty -----------------------------------------------------
if [ -z "$SHEETPATH" ]; then
   echo "Missing argument #1: You did not provide a SHEETPATH (i.e. path to the Fusion Excel Workbook /Users/kerberos/Downloads/23-MGFS5.xlsx)"
   exit 1
fi

bold=$(tput bold)         # makes console output text bold
BG_BLUE="$(tput setab 4)" # makes console output text background blue
BG_GRN="$(tput setab 2)"  # makes console output text background green
FG_YLW="$(tput setaf 3)"  # makes console output text yellow color
BG_MAG="$(tput setab 5)"  # makes console output text background magenta
FG_CYA="$(tput setaf 6)"  # makes console output text color cyan
normal=$(tput sgr0)       # resets default console output text

samSheetDest="/gpfs/data/molecpathlab/production/samplesheets/archer/"
GITHUBMAIN="https://github.com/NYU-Molecular-Pathology/Methylation/raw/main/Development/"
FSID=$(basename ${SHEETPATH%%.*}) # from 23-MGFS##.xlsx

# Download Python Script and Archer Sequences --------------------------------------------------------------------------------------------------------------
message_curl() {
   linkname=$1
   filename=$2
   printf "${bold}${BG_GRN}$HOME/$filename${normal}\n"
   curl -# -L $linkname$filename >$HOME/$filename
   chmod +rwx $HOME/$filename
}

# FUN: Message a step in color blue -----------------------------------------------------------------
message_print() {
   stringMsg=$1
   extraMsg=$2
   printf "\n${bold}${BG_BLUE}${stringMsg}:${normal}\n${extraMsg}\n\n"
}

cd $HOME

message_print "Input Kerberos ID" "$kerbero"
message_print "Input SHEETPATH" "$SHEETPATH"
message_print "Input FSID" "$FSID"

message_curl ${GITHUBMAIN} "Archer_Index2_Sequences.xlsx"
message_curl ${GITHUBMAIN} "gen_sample_sheet_fixed.py"

# Execute Python Script to Generate Samplesheet ---------------------------------------------------------------------------------------------------------
message_print "Executing" "python3 gen_sample_sheet_fixed.py -t \"${SHEETPATH}\" -o \"${HOME}/Desktop/\" --I5_index \"${HOME}/Archer_Index2_Sequences.xlsx\""
python3 gen_sample_sheet_fixed.py -t "$SHEETPATH" -o "$HOME/Desktop/" --I5_index "$HOME/Archer_Index2_Sequences.xlsx"

NEWCSV=$(ls -t ~/Desktop/*-SampleSheet.csv | awk '{printf($0);exit}')
message_print "SampleSheet output is on your Desktop" "$NEWCSV"

# Save Run String Variable Information -----------------------------------------------------------------------------------------------------------------
FUSIONRUNID=$(basename ${NEWCSV%%-*})
DEMUXDIR="/gpfs/data/molecpathlab/production/Demultiplexing/${FG_YLW}$FUSIONRUNID${normal}"
userID="$kerbero"
samSheetPath="$userID@bigpurple.nyumc.org:$samSheetDest$FUSIONRUNID"

# Create Directory Samplesheet Run Driectory in BigPurple -----------------------------------------------------------------------------------------------
if ssh "$userID@bigpurple.nyumc.org" "[ ! -d $samSheetDest$FUSIONRUNID ]"; then
   message_print "Creating new Directory" "\n$samSheetDest$FUSIONRUNID"
   message_print "Executing" "ssh ${FG_CYA}$userID${normal}@bigpurple.nyumc.org mkdir -p $samSheetDest$FUSIONRUNID"
   ssh "$userID@bigpurple.nyumc.org" "mkdir -p $samSheetDest$FUSIONRUNID"
fi

# Copy Output File with Rsync to BigPurple --------------------------------------------------------------------------------------------------------------
message_print "Copying the file" "\n${NEWCSV}\nto\n${samSheetPath}"
message_print "Executing the following command" "rsync -e ssh $NEWCSV ${samSheetPath}"
rsync -e ssh "$NEWCSV" "$samSheetPath"

# Edit new Directory Permissions on BigPurple --------------------------------------------------------------------------------------------------------------
message_print "Changing directory permissions" "ssh \"$kerbero@bigpurple.nyumc.org\" chmod -R g+rwx $samSheetDest$FUSIONRUNID"
ssh "$kerbero@bigpurple.nyumc.org" "chmod -R g+rwx $samSheetDest$FUSIONRUNID"

# Print Out Commands to complete the run ----------------------------------------------------------------------------------------------------------------
message_curl ${GITHUBMAIN} "printFusionCommands.sh"
chmod ugo+rwx "$HOME/printFusionCommands.sh"
chmod +rwx "$HOME/printFusionCommands.sh"

message_print "Saving Html File" "$HOME/printFusionCommands.sh $FUSIONRUNID ${FSID} ${kerbero} >${HOME}/${FSID}.html && open $HOME/${FSID}.html"

$HOME/printFusionCommands.sh $FUSIONRUNID ${FSID} ${kerbero} >${HOME}/${FSID}.html && open $HOME/${FSID}.html
