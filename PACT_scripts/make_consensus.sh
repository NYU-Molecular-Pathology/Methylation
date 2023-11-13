#!/bin/bash
## ---------------------------
## Script name: make_consensus.sh
## Purpose: Generate PACT Consensus HTML report with given input
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------

# Input args
runID=${1-NULL}   # if arg $1 is empty assign NULL as default else i.e. 123456_NB501073_0212_AHT3V7BGXK
pactRun=${2-NULL} # if arg $2 is empty assign NULL as default else i.e. PACT-21-28
DEFAULTVALUE="/Volumes/CBioinformatics/jonathan/pact/consensus/"
consensusDir="${3:-$DEFAULTVALUE}"
kerbero=${4-$USER} # if arg $3 is empty assign $USER as default else i.e. whoami kerberosid

# Current Year Calculation
currYear=$(date +'%Y')

# Hardcoded paths
pactGithub="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts"
outputDir="/molecular/MOLECULAR LAB ONLY/NYU PACT Patient Data/Results/Bioinformatics/"
#clinicalOuput="/Volumes/molecular/Molecular/NGS607/${currYear}/${runID}/output/clinical/"

mkdir -p "${consensusDir}${pactRun}_consensus" 

cd "${consensusDir}${pactRun}_consensus" || exit

curl -# -L ${pactGithub}/PACT_consensus.Rmd >"${consensusDir}${pactRun}_consensus/${pactRun}_consensus.Rmd"

/Volumes/CBioinformatics/PACT/getMethylMatch.sh "${pactRun}" "${runID}"

cd "${consensusDir}${pactRun}_consensus" || exit

volMolec="/Volumes/molecular/Molecular/"
deskDir="$HOME/Desktop/${pactRun}/"

mkdir -p "${deskDir}"

echo "Copying files: rsync -vrthP --include='*.png' '${volMolec}REDCap/cnv_facets/${pactRun}/' '${deskDir}'"

rsync -vrthP --include="*.png" "${volMolec}REDCap/cnv_facets/${pactRun}/" "${deskDir}"
#rsync -vrthP "${clinicalOuput}${pactRun}-Somatic_Variants.html" "${deskDir}"
#rsync -vrthP "${clinicalOuput}${pactRun}.html" "${deskDir}"
rsync -vrthP "${volMolec}REDCap/cnv_facets/${pactRun}/${pactRun}-QC.tsv" "${deskDir}"
rsync -vrthP "${volMolec}NGS607/${currYear}/${runID}/${pactRun}_Hotspots.tsv" "${deskDir}"
mkdir -p "${consensusDir}${pactRun}_consensus/cnvpng/"
rsync -vrthP --include="*.png" "${deskDir}" "${consensusDir}${pactRun}_consensus/cnvpng/"
rsync -vrthP --include="*.tsv" "${deskDir}" "${consensusDir}${pactRun}_consensus"
rsync -vrthP "$HOME/Desktop/${runID}-SampleSheet.csv" "${consensusDir}${pactRun}_consensus"
rsync -vrthP "$HOME/Desktop/${pactRun}_MethylMatch.xlsx" "${consensusDir}${pactRun}_consensus"
#rsync -vrthP "${deskDir}${pactRun}-Somatic_Variants.html" "${consensusDir}${pactRun}_consensus"
#rsync -vrthP "${deskDir}${pactRun}.html" "${consensusDir}${pactRun}_consensus"

cd "${HOME}" && curl -# -L ${pactGithub}/MakeIndelList.R >"${HOME}/MakeIndelList.R"

chmod +rwx "${HOME}/MakeIndelList.R"

RScript --verbose "${HOME}/MakeIndelList.R" "${pactRun}"

cd "${consensusDir}${pactRun}_consensus/" || exit

Rscript --verbose -e "rmarkdown::render('${consensusDir}${pactRun}_consensus/${pactRun}_consensus.Rmd', params=list(pactName='${pactRun}', userName='${kerbero}'))" && open "${consensusDir}${pactRun}_consensus/${pactRun}_consensus.html"

echo "Trying to copy files:"
echo "rsync -vrthP ${consensusDir}${pactRun}_consensus/${pactRun}_consensus.html \"/Volumes${outputDir}${currYear}/${pactRun}/\""

rsync -vrthP "${consensusDir}${pactRun}_consensus/${pactRun}_consensus.html" "/Volumes${outputDir}${currYear}/${pactRun}/"
