#!/bin/bash

methAPI='XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' # RedCap API Token

# INPUT ARGS: if args are empty assign NULL as default below
methRun=${1-NULL}  # methylation run id e.g. 22-MGDM17
PRIORITY=${2-NULL} # string of prioritized RD-numbers
runPath=${3-NULL}  # any custom directory to copy/run the idat files
redcapUp=${4-NULL} # to upload to redcap or not if server down single char i.e. "T" or "F"
runLocal=${5-NULL} # If the run directory should be executed without shared drives locally i.e. "T" or "F"

# HARDCODED GITHUB PATHS ---------------------------------------------------------------------------------------
GITHUBLINK="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/dockerScripts/"
GITHUBMAIN="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/"

# Hardcoded Terminal Color Variables ---------------------------------------------------------------------------
# tput changes the console text output display
BG_BLUE="$(tput setab 4)" # makes text background blue
BG_CYAN="$(tput setab 6)" # makes text background cyan
BG_MAG="$(tput setab 5)"  # makes text background magenta
BG_GRN="$(tput setab 2)"  # makes text background green
BG_RED="$(tput setab 1)"  # makes text background green
bold=$(tput bold)         # makes text bold
normal=$(tput sgr0)       # resets default text

cd $HOME

# FUNCTION: EXIT when any command fails ------------------------------------------------------------------------
trap_error() {
   LSTCMD=$1
   EXCODE=$2
   if [ $EXCODE -ne 0 ]; then
      printf "\n${BG_RED}ERROR on LAST COMMAND${normal}:\n${LSTCMD}\n"
      echo "Command filed with exit code $EXCODE."
   fi
}

set -e
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'trap_error "${last_command}" $?' EXIT
# CHECK: Prints instruction if methRun parameter is empty
if [ -z "$methRun" ]; then
   echo "You did not provide a Methylation Run ID name. Ex. '21-MGDM30'"
   exit 1
fi

if [ -z "$redcapUp" ]; then
   redcapUp="T"
fi

echo "~~~~~~~~~~~~~~~~~~~~~Parameters passed to RScript~~~~~~~~~~~~~~~~~~~~~"
echo "${BG_BLUE}methAPI:${normal} $methAPI"
echo "${BG_CYAN}(arg1) methRun:${normal} $methRun"
echo "${BG_MAG}(arg2) PRIORITY samples:${normal} $PRIORITY"
echo "${BG_GRN}(arg3) Custom runPath:${normal} $runPath"
echo "${BG_BLUE}(arg4) To upload redcapUp:${normal} $redcapUp"
printf "\nCurl files downloaded:\n"

message_curl() {
   linkname=$1
   filename=$2
   printf "${BG_GRN}$HOME/$filename${normal}\n"
   curl -k -# -L $linkname$filename >$HOME/$filename
}

check_directory() {
   CURRDIR=$1
   if [ ! -d "$CURRDIR" ]; then
      printf "Creating new Directory:\n$CURRDIR"
      mkdir -p "$CURRDIR"
   fi
}

message_curl ${GITHUBMAIN} "report.Rmd"
message_curl ${GITHUBLINK} "Methyl_QC.Rmd"
message_curl ${GITHUBMAIN} "methylExpress.R"

printf "\nExecuting Rscript with parameters input:\n"
color_params="${BG_MAG}$methAPI${normal} ${BG_CYAN}$methRun${normal} ${BG_GRN}$PRIORITY${normal} ${BG_BLUE}$runPath${normal} ${BG_MAG}$redcapUp${normal} ${BG_GRN}$runLocal${normal}\n"
printf "${BG_BLUE}${bold}Rscript --verbose $HOME/methylExpress.R${normal} $color_params"
Rscript --verbose $HOME/methylExpress.R $methAPI $methRun $PRIORITY $runPath $redcapUp $runLocal

if [ -n "$runPath" ]; then
   DEFAULTPATH="/Volumes/CBioinformatics/Methylation/Clinical_Runs/${methRun}/"
   runPath=$DEFAULTPATH
fi

LOCALDIR="$HOME/${methRun}"
RUNYEAR=${methRun%%-*}
MOLECDIR="/Volumes/molecular/MOLECULAR LAB ONLY/NYU-METHYLATION/Results/20${RUNYEAR}"

echo "${BG_GRN}runPath:${normal} $runPath"
echo "${BG_GRN}MOLECDIR:${normal} $MOLECDIR"
echo "${BG_GRN}LOCALDIR:${normal} $LOCALDIR"

check_directory "$LOCALDIR"
printf "\nRsyncing files:\n"
rsync -vrthP --include '*.html' --exclude '*' "${runPath}" "$LOCALDIR"
#check_directory "$MOLECDIR"
#chmod +rwx "$LOCALDIR"
#cd "$MOLECDIR"
printf "${BG_GRN}cp -RvX ${LOCALDIR} ${MOLECDIR}${normal}"
#rsync -vrthP --include '*.html' "${LOCALDIR}/" "$PWD"
cp -RvX "$LOCALDIR" "$MOLECDIR"
echo "METHYLATION RUN ENDED"
