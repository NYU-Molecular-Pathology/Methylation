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

consenDir="${workDir}${pactRun}_consensus/"

csvFile="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts/PACT_desc.csv"
rmdFile="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts/PACT_consensus.Rmd"
echo " "
echo "${bold}${BG_BLUE}Checking the directory for your PACT consensus:${normal}"
echo "[ -d \"$workDir${FG_RED}$pactRun${normal}_consensus\" ] ||"
echo "mkdir \"$workDir${FG_RED}$pactRun${normal}_consensus\""

checkDir=`[ -d "${workDir}${pactRun}_consensus" ] || mkdir "${workDir}${pactRun}_consensus"`
echo "$checkDir"

echo "cd \"$workDir${FG_RED}$pactRun${normal}_consensus\""
changeDir=`cd $consenDir`
#echo $changeDir

echo "${bold}${BG_BLUE}To this directory, copy all the .cnv.plot.pdf facets files from:${normal}"
echo "/gpfs/data/molecpathlab/production/NGS607/${FG_YLW}runID${normal}/output/cnv/FACETS/*.pdf"
echo " "
echo "${bold}${BG_BLUE}Copy the following files to the concenus Rmd directory:${normal}"
echo "${FG_RED}$pactRun${normal}_MethylMatch.xlsx"
echo "${FG_YLW}$runID${normal}-SampleSheet.csv"
echo " "

echo "${bold}${BG_BLUE}Copy the QC file below to the Rmd directory${normal}"

echo " "
echo "${bold}${BG_BLUE}Checking if files are in the directory directory:${normal}"

descFi="$consenDir${pactRun}_desc.csv"
rmdFi="$consenDir${pactRun}_consensus.Rmd"


checkSam=`ls ${consenDir}*-SampleSheet.csv`

if [ -f "$descFi" ]; then
    echo "${FG_GRN}${pactRun}_desc.csv${normal} exists."
else 
    echo "Copying ${FG_RED}$pactRun${normal}_desc.csv from GitHub..."
    echo "curl -o ${FG_RED}$pactRun${normal}_desc.csv -L $csvFile -s"
    echo `curl -o ${descFi} -L $csvFile -s`
fi

if [ -f "$rmdFi" ]; then
    echo "${FG_GRN}${pactRun}_consensus.Rmd${normal} exists."
else 
    echo "~~~~~"
    echo "Copying ${FG_RED}$pactRun${normal}_consensus.Rmd from GitHub..."
    echo "curl -o ${FG_RED}$pactRun${normal}_consensus.Rmd -L $rmdFile -s"
    echo `curl -o $rmdFi -L $rmdFile -s`
fi

if [ -f "$checkSam" ]; then
    echo "${FG_GRN}-SampleSheet.csv${normal} exists."
else
   echo "~~~~~"
   echo "${FG_RED}-SampleSheet.csv${normal} is missing from the knit directory."
   echo "Generating Samplesheet with makePactSheet.sh:"
   echo `/Volumes/CBioinformatics/PACT/makePactSheet.sh ${pactRun}`
   echo "Copy -SampleSheet.csv file to the knit directory and run again."
   echo "Exiting"
   exit 0
fi

if [ -f "$consenDir${pactRun}_MethylMatch.xlsx" ]; then
    echo "${FG_GRN}${pactRun}_MethylMatch.xlsx${normal} exists."
else 
    echo "~~~~~"
    echo "${FG_RED}${pactRun}_MethylMatch.xlsx${normal} is missing from the knit directory."
    echo "Writing xlsx using /Volumes/CBioinformatics/PACT/methylMatch.sh:"
    echo `/Volumes/CBioinformatics/PACT/methylMatch.sh ${pactRun}`
    echo "Exiting"
    exit 0
fi

if [ -f "$consenDir${pactRun}-QC.tsv" ]; then
    echo "${FG_GRN}${pactRun}-QC.tsv${normal} exists."
else 
    echo "~~~~~"
    echo "${FG_RED}${pactRun}-QC.tsv${normal} is missing from the knit directory."
    echo "~~~~~"
    echo "Copy the tsv file below:"
    echo "/gpfs/data/molecpathlab/production/NGS607/${FG_YLW}runID${normal}/${FG_RED}$pactRun${normal}-QC.tsv"
    echo "To your directory:"
    echo "$consenDir"
    echo "Exiting"
    exit 0
fi

echo " "
echo "${bold}${BG_BLUE}Knitting Rmd file with description csv:${normal}"
echo "Rscript -e \"rmarkdown::render('${FG_RED}$pactRun${normal}_consensus.Rmd', params = list(pactName=${FG_RED}$pactRun${normal}, userName=${FG_GRN}$pactUser${normal}, workDir=${FG_CYA}$workDir${normal}${normal}))\""
echo `Rscript -e "rmarkdown::render('$rmdFi', params = list(pactName='${pactRun}', userName='$pactUser', workDir='$workDir'))"`
