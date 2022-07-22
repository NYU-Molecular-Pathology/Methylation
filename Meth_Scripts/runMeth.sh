#!/bin/bash

methAPI='XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' # RedCap API Token

# if args are empty assign NULL as default below
methRun=${1-NULL}  # methylation run id e.g. 22-MGDM17
PRIORITY=${2-NULL} # string of prioritized RD-numbers
runPath=${3-NULL}  # any custom directory to copy/run the idat files
redcapUp=${4-NULL} # to upload to redcap or not if server down single char i.e. "T" or "F"

cd $HOME

GITHUBLINK="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/dockerScripts/"
GITHUBMAIN="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/"

# Print instruction if methRun parameter is empty
if [ -z "$methRun" ]; then
   echo "You did not provide a Methylation Run ID name. Ex. '21-MGDM30'"
   exit 1
fi

if [ -z "$redcapUp" ]; then
   redcapUp="T"
fi

# tput changes the console text output display
BG_BLUE="$(tput setab 4)" # makes text background blue
BG_CYAN="$(tput setab 6)" # makes text background cyan
BG_MAG="$(tput setab 5)"  # makes text background magenta
BG_GRN="$(tput setab 2)"  # makes text background green
bold=$(tput bold)         # makes text bold
normal=$(tput sgr0)       # resets default text

echo "~~~~~~~~~~~~~~~~~~~~~Parameters passed to RScript~~~~~~~~~~~~~~~~~~~~~"
echo "${BG_BLUE}methAPI:${normal} $methAPI"
echo "${BG_CYAN}(arg1) methRun:${normal} $methRun"
echo "${BG_MAG}(arg2) PRIORITY samples:${normal} $PRIORITY"
echo "${BG_GRN}(arg3) Custom runPath:${normal} $runPath"
echo "${BG_BLUE}(arg4) To upload redcapUp:${normal} $redcapUp"
printf "\nCurl files downloaded:\n"

[ ! -d "$HOME/logos" ] && mkdir $HOME/logos

message_curl() {
   linkname=$1
   filename=$2
   printf "${BG_GRN}$HOME/$filename${normal}\n"
   curl -# -L $linkname$filename >$HOME/$filename
}

message_curl ${GITHUBLINK} "report.Rmd"
message_curl ${GITHUBLINK} "logos.html"
message_curl ${GITHUBLINK} "logos/logos.css"
message_curl ${GITHUBLINK} "logos/logos.png"
message_curl ${GITHUBLINK} "logos/logos.jpeg"

message_curl ${GITHUBLINK} "Methyl_QC.Rmd"
message_curl ${GITHUBMAIN} "methylExpress.R"
printf "\nExecuting Rscript with parameters input:\n"
color_params="${BG_MAG}$methAPI${normal} ${BG_CYAN}$methRun${normal} ${BG_GRN}$PRIORITY${normal} ${BG_BLUE}$runPath${normal} ${BG_MAG}$redcapUp${normal}"
printf "${BG_BLUE}${bold}Rscript --verbose $HOME/methylExpress.R${normal} $color_params"
Rscript --verbose $HOME/methylExpress.R $methAPI $methRun $PRIORITY $runPath $redcapUp
