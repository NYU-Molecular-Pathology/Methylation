#!/bin/bash
runID=${1-NULL} # if arg $1 is empty assign NULL as default else i.e. 123456_NB501073_0212_AHT3V7BGXK
pactRun=${2-NULL} # if arg $2 is empty assign NULL as default else i.e. PACT-21-28
kerbero=${3-$USER} # if arg $3 is empty assign $USER as default else i.e. whoami kerberosid

FILE="$HOME/demuxQC.sh"

bold=$(tput bold) # makes console output text bold
BG_GRN="$(tput setab 2)" # makes console output text background green
normal=$(tput sgr0) # resets default console output text format

if [ -f "$FILE" ]; then
   sh "$HOME/demuxQC.sh" $runID $pactRun > "$HOME/$pactRun_stages.txt"
else
   curl -o "$HOME/demuxQC.sh" -L https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts/demuxQC.sh
   chmod gu+rwx "$HOME/demuxQC.sh"
   sh "$HOME/demuxQC.sh" $runID $pactRun > "$HOME/$pactRun_stages.txt"
fi

shopt -s nocasematch
echo "${BG_GRN}Enter an NGS stage to view steps and press ${bold}[return]:${normal}"
echo "1: Create SampleSheet.csv to BigPurple & Demultiplexing"
echo "2: In-House and Philips Uploads"
echo "3: QC Generation"
echo "4: Copy the QC and Output to Zdrive"
echo "5: Consensus Files"
echo ""
read NGSTAGE
echo -n "The steps for stage $NGSTAGE are: "
echo ""
case $NGSTAGE in
1)
sed -n '1,27p' "$HOME/$pactRun_stages.txt"
;;
2)
sed -n '28,54p' "$HOME/$pactRun_stages.txt"
;;
3)
sed -n '55,67p' "$HOME/$pactRun_stages.txt"
;;
4)
sed -n '68,113p' "$HOME/$pactRun_stages.txt"
;;
5)
sed -n '114,150p' "$HOME/$pactRun_stages.txt"
;;
*)
echo -n "Invalid Stage Number"
esac
echo ""
