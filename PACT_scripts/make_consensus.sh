#!/bin/bash
## ---------------------------
## Script name: make_consensus.sh
## Purpose: Generate PACT Consensus HTML report with given input
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------

# Input args -------------------------------------------------------------------------------------------
runID=${1-NULL}   # if arg $1 is empty assign NULL as default else i.e. 123456_NB501073_0212_AHT3V7BGXK
pactRun=${2-NULL} # if arg $2 is empty assign NULL as default else i.e. PACT-21-28
DEFAULTVALUE="/Volumes/CBioinformatics/jonathan/pact/consensus/"
consensusDir="${3:-$DEFAULTVALUE}"
kerbero=${4-$USER} # if arg $3 is empty assign $USER as default else i.e. whoami kerberosid

# Hardcoded paths
pactGithub="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts"
outputDir="/molecular/MOLECULAR LAB ONLY/NYU PACT Patient Data/Results/Bioinformatics/"

# Color Variables --------------------------------------------------------------------------------------
BG_GRN="$(tput setab 2)" # makes text background green
BG_RED="$(tput setab 1)" # makes text background green
NORMAL=$(tput sgr0)      # resets default text

# Function to exit when any command fails
trap_error() {
   local last_command=$1
   local exit_code=$2
   local line_no=$3

   if [[ $exit_code -ne 0 && $error_handled == false ]]; then
      printf "\n%sERROR on LAST COMMAND%s:\n%s\n" "${BG_RED}" "${NORMAL}" "${last_command}"
      printf "Command failed at line %s with exit code %s.\n" "${line_no}" "${exit_code}"
      error_handled=true
   fi
}

# Enable strict mode
set -Eeuo pipefail
trap 'last_command=$BASH_COMMAND; last_line=$LINENO' DEBUG
trap 'trap_error "${BASH_COMMAND}" $? "${LINENO}"' ERR EXIT

# Current Year Calculation -----------------------------------------------------------------------------
currYear=$(date +'%Y')

# Function to display and execute code
msg_code() {
    printf "\n${BG_GRN}Executing:${NORMAL}\n"
    printf "$*\n"
    eval "$*"
}

# Function to create directories
create_dir() {
    mkdir -p "$1"
}

# Main execution
create_dir "${consensusDir}${pactRun}_consensus"
cd "${consensusDir}${pactRun}_consensus" || exit
curl -# -L ${pactGithub}/PACT_consensus.Rmd >"${consensusDir}${pactRun}_consensus/${pactRun}_consensus.Rmd"

/Volumes/CBioinformatics/PACT/getMethylMatch.sh "${pactRun}" "${runID}"

cd "${consensusDir}${pactRun}_consensus" || exit

# Consolidated rsync commands
MOLEC_VOL="/Volumes/molecular/Molecular/"
DESK_DIR="$HOME/Desktop/${pactRun}/"
FACET_DIR="${consensusDir}${pactRun}_consensus/cnvpng/"

create_dir "${DESK_DIR}"
create_dir "${FACET_DIR}"

msg_code rsync -vrthP --include=\"*.png\" \"${MOLEC_VOL}REDCap/cnv_facets/${pactRun}/\" \"${DESK_DIR}\"
msg_code rsync -vrthP \"${MOLEC_VOL}REDCap/cnv_facets/${pactRun}/${pactRun}-QC.tsv\" \"${DESK_DIR}\"
msg_code rsync -vrthP \"${MOLEC_VOL}NGS607/${currYear}/${runID}/${pactRun}_Hotspots.tsv\" \"${DESK_DIR}\"

printf "\n${BG_GRN}Executing:${NORMAL}\n"
printf "cp ${DESK_DIR}/*.png ${FACET_DIR}/"
cp "${DESK_DIR}"/*.png "${FACET_DIR}/"

msg_code rsync -vrthP --include=\"*.tsv\" \"${DESK_DIR}\" \"${consensusDir}${pactRun}_consensus\"
msg_code rsync -vrthP \"$HOME/Desktop/${runID}-SampleSheet.csv\" \"${consensusDir}${pactRun}_consensus\"
msg_code rsync -vrthP \"$HOME/Desktop/${pactRun}_MethylMatch.xlsx\" \"${consensusDir}${pactRun}_consensus\"

# Final operations for generating html
cd "${HOME}" && curl -# -L ${pactGithub}/MakeIndelList.R >"${HOME}/MakeIndelList.R"
chmod +rwx "${HOME}/MakeIndelList.R"

RScript --verbose "${HOME}/MakeIndelList.R" "${pactRun}"
cd "${consensusDir}${pactRun}_consensus/" || exit

echo "Running Rscript:"

Rscript --verbose -e "rmarkdown::render('${consensusDir}${pactRun}_consensus/${pactRun}_consensus.Rmd', params=list(pactName='${pactRun}', userName='${kerbero}'))" && open "${consensusDir}${pactRun}_consensus/${pactRun}_consensus.html"

echo "Trying to copy html file to ouput folder:"
echo "/Volumes${outputDir}${currYear}/${pactRun}/"

#rsync -vrthP "${consensusDir}${pactRun}_consensus/${pactRun}_consensus.html" "/Volumes${outputDir}${currYear}/${pactRun}/"

CONSENSUS_FILE="${consensusDir}${pactRun}_consensus/${pactRun}_consensus.html"
FINAL_DEST="/Volumes${outputDir}${currYear}/${pactRun}/"

# Copy HTML file from consensus directory to DESK_DIR
cp "${CONSENSUS_FILE}" "${DESK_DIR}"

# Copy HTML file from DESK_DIR to final destination
cp -X "${DESK_DIR}${pactRun}_consensus.html" "${FINAL_DEST}"
