#!/bin/bash
## ---------------------------
## Script name: make_consensus.sh
## Purpose: Generate PACT Consensus HTML report with given input
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------

runID=${1-NULL}   # if arg $1 is empty assign NULL as default else i.e. 123456_NB501073_0212_AHT3V7BGXK
pactRun=${2-NULL} # if arg $2 is empty assign NULL as default else i.e. PACT-21-28
DEFAULTVALUE="/Volumes/CBioinformatics/jonathan/pact/consensus/"
consensusDir="${3:-$DEFAULTVALUE}"
kerbero=${4-$USER} # if arg $3 is empty assign $USER as default else i.e. whoami kerberosid

currYear=$(date +'%Y')
pactGithub="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts"

mkdir -p "${consensusDir}${pactRun}_consensus" 

cd "${consensusDir}${pactRun}_consensus" || exit

curl -# -L ${pactGithub}/PACT_consensus.Rmd >"${consensusDir}${pactRun}_consensus/${pactRun}_consensus.Rmd"

/Volumes/CBioinformatics/PACT/getMethylMatch.sh "${pactRun}" "${runID}"

cd "${consensusDir}${pactRun}_consensus" || exit

echo "Copying files: cp -fvX '/Volumes/molecular/Molecular/REDCap/cnv_facets/${pactRun}/*.pdf' '$HOME/Desktop/${runID}-SampleSheet.csv' '$HOME/Desktop/${pactRun}_MethylMatch.xlsx /Volumes/molecular/Molecular/NGS607/${currYear}/${runID}/output/clinical/${pactRun}-Somatic_Variants.html" "/Volumes/molecular/Molecular/NGS607/${currYear}/${runID}/output/clinical/${pactRun}.html' '/Volumes/molecular/Molecular/REDCap/cnv_facets/${pactRun}/${pactRun}-QC.tsv' '/Volumes/molecular/Molecular/NGS607/${currYear}/${runID}/output/${pactRun}_Hotspots.tsv' ./
"
cp -fvX "/Volumes/molecular/Molecular/REDCap/cnv_facets/${pactRun}/*.pdf" "$HOME/Desktop/${runID}-SampleSheet.csv" "$HOME/Desktop/${pactRun}_MethylMatch.xlsx /Volumes/molecular/Molecular/NGS607/${currYear}/${runID}/output/clinical/${pactRun}-Somatic_Variants.html" "/Volumes/molecular/Molecular/NGS607/${currYear}/${runID}/output/clinical/${pactRun}.html" "/Volumes/molecular/Molecular/REDCap/cnv_facets/${pactRun}/${pactRun}-QC.tsv" "/Volumes/molecular/Molecular/NGS607/${currYear}/${runID}/output/${pactRun}_Hotspots.tsv" ./

cd "${HOME}" && curl -# -L ${pactGithub}/MakeIndelList.R >"${HOME}/MakeIndelList.R"

chmod +rwx "${HOME}/MakeIndelList.R"

RScript --verbose "${HOME}/MakeIndelList.R" "${pactRun}"

cd "${consensusDir}${pactRun}_consensus/" || exit

Rscript --verbose -e "rmarkdown::render('${consensusDir}${pactRun}_consensus/${pactRun}_consensus.Rmd', params=list(pactName='${pactRun}', userName='${kerbero}'))"
