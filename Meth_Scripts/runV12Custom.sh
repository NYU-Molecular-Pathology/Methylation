#!/bin/bash

METHAPI='XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' # RedCap API Token

# INPUT ARGS: if args are empty assign NULL as default below
CSVPATH=${1-NULL}
MAKESARC=${2-NULL}

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
# CHECK: Prints instruction if CSVPATH parameter is empty
if [ -z "$CSVPATH" ]; then
    echo "You did not provide a csv File Path"
    exit 1
fi

echo "~~~~~~~~~~~~~~~~~~~~~Parameters passed to RScript~~~~~~~~~~~~~~~~~~~~~"
echo "${BG_BLUE}(inFile) METHAPI:${normal} $METHAPI"
echo "${BG_CYAN}(arg1) CSVPATH:${normal} $CSVPATH"
echo "${BG_BLUE}(arg2) MAKESARC:${normal} $MAKESARC"

# Curl function with message -----------------------------------------------------------------------------------
printf "\n${BG_GRN}Curl files downloaded:${normal}\n"
message_curl() {
    linkname=$1
    filename=$2
    printf "${BG_GRN}$HOME/$filename${normal}\n"
    curl -# -L $linkname$filename >$HOME/$filename
}

# Download Required files with curl to $HOME -------------------------------------------------------------------
message_curl ${GITHUBMAIN} "report.Rmd"
message_curl ${GITHUBLINK} "Methyl_QC.Rmd"
message_curl ${GITHUBMAIN} "ReRunV12_express.R"
message_curl ${GITHUBMAIN} "SarcReport.Rmd"

printf "\nExecuting Rscript with parameters input:\n"
color_params="${BG_MAG}$METHAPI${normal} ${BG_CYAN}$CSVPATH${normal} ${BG_GRN}$MAKESARC${normal}\n"
printf "${BG_BLUE}${bold}Rscript --verbose $HOME/ReRunV12_express.R${normal} $color_params"

Rscript --verbose $HOME/ReRunV12_express.R $METHAPI $CSVPATH $MAKESARC

echo "METHYLATION RUN ENDED"
