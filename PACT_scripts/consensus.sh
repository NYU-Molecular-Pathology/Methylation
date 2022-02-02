#!/bin/bash
pactRun=${1-NULL} # if arg $1 is empty assign NULL as default else i.e. PACT-21-28
pactUser=${2-$USER} # if arg $2 is empty assign $USER as default else i.e. whoami kerberosid
workDir=${3-NULL} # if arg $3 is empty assign NULL as default  

bold=$(tput bold) # makes console output text bold
undLN=$(tput smul) # makes console output text underline
BG_BLUE="$(tput setab 4)" # makes console output text background blue
BG_CYAN="$(tput setab 6)" # makes console output text background cyan
BG_MAG="$(tput setab 5)" # makes console output text background magenta
BG_GRN="$(tput setab 2)" # makes console output text background green
FG_YLW="$(tput setaf 3)" # makes console output text yellow color
FG_RED="$(tput setaf 1)" # makes console output text red color
FG_GRN="$(tput setaf 2)" # makes console output text color green
FG_BLU="$(tput setaf 4)" # makes console output text color blue
FG_CYA="$(tput setaf 6)" # makes console output text color cyan
normal=$(tput sgr0) # resets default console output text

case "$workDir" in */)
    echo "${FG_CYA}workDir:${normal}"
    echo "${FG_CYA}$workDir${normal}";;
*)
    echo "${FG_CYA}workDir:${normal}"
    workDir="$workDir/"
    echo "${FG_CYA}$workDir${normal}";;
esac

concenDir="${workDir}${pactRun}_consensus/"

csvFile="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts/PACT_desc.csv"
rmdFile="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts/PACT_consensus.Rmd"
echo " "
echo "${bold}${BG_BLUE}Checking the directory for your PACT consensus:${normal}"
echo "[ -d \"$workDir${FG_RED}$pactRun${normal}_consensus\" ] ||"
echo "mkdir \"$workDir${FG_RED}$pactRun${normal}_consensus\""

checkDir=`[ -d "${workDir}${pactRun}_consensus" ] || mkdir "${workDir}${pactRun}_consensus"`
echo "$checkDir"

echo "cd \"$workDir${FG_RED}$pactRun${normal}_consensus\""
echo `cd $concenDir`

echo "${bold}${BG_BLUE}To this directory, copy all the .cnv.plot.pdf facets files from:${normal}"
echo "/gpfs/data/molecpathlab/production/NGS607/${FG_YLW}runID${normal}/output/cnv/FACETS/*.pdf"
echo " "
echo "${bold}${BG_BLUE}Copy the following files to the concenus Rmd directory:${normal}"
echo "${FG_RED}$pactRun${normal}_MethylMatch.xlsx"
echo "${FG_YLW}$runID${normal}-SampleSheet.csv"
echo " "

echo "${bold}${BG_BLUE}Copy the QC file below to the Rmd directory${normal}"
echo "/gpfs/data/molecpathlab/production/NGS607/${FG_YLW}runID${normal}/${FG_RED}$pactRun${normal}-QC.tsv"
echo " "
echo "${bold}${BG_BLUE}Downloading files to working directory:${normal}"
descFi="$concenDir${pactRun}_desc.csv"
rmdFi="$concenDir${pactRun}_consensus.Rmd"

if [ -f "$descFi" ]; then
    echo "${FG_RED}$pactRun${normal}_desc.csv exists."
else 
    echo "Copying ${FG_RED}$pactRun${normal}_desc.csv from GitHub..."
    echo "curl -o ${FG_RED}$pactRun${normal}_desc.csv -L $csvFile -s"
    echo `curl -o ${descFi} -L $csvFile -s`
fi

if [ -f "$rmdFi" ]; then
    echo "${FG_RED}$pactRun${normal}_consensus.Rmd exists."
else 
    echo "Copying ${FG_RED}$pactRun${normal}_consensus.Rmd from GitHub..."
    echo "curl -o ${FG_RED}$pactRun${normal}_consensus.Rmd -L $rmdFile -s"
    echo `curl -o $rmdFi -L $rmdFile -s`
fi
echo " "
echo "${bold}${BG_BLUE}Knitting Rmd file with description csv:${normal}"
echo "Rscript -e \"rmarkdown::render('${FG_RED}$pactRun${normal}_consensus.Rmd', params = list(pactName=${FG_RED}$pactRun${normal}, userName=${FG_GRN}$pactUser${normal}, workDir=${FG_CYA}$workDir${normal}${normal}))\""
echo `Rscript -e "rmarkdown::render('$rmdFi', params = list(pactName='${pactRun}', userName='$pactUser', workDir='$workDir'))"`
