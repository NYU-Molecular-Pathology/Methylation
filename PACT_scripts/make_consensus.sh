#!/bin/bash
## ---------------------------
## Script name: make_consensus.sh
## Purpose: Generate PACT Consensus HTML report with given input
## Author: Jonathan Serrano
## Date Created: August 21, 2023
## Author: Jonathan Serrano
## Version: 1.0.0
## Copyright (c) NYULH Jonathan Serrano, 2024
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
error_handled=true

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
yearPart=${pactRun:5:2}
# Concatenate "20" with the extracted part
currYear="20${yearPart}"
#currYear=$(date +'%Y')

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
WORK_DIR="${consensusDir}${pactRun}_consensus"

create_dir "${WORK_DIR}"
cd "${WORK_DIR}" || exit
curl -# -L ${pactGithub}/PACT_consensus.Rmd >"${WORK_DIR}/${pactRun}_consensus.Rmd"

/Volumes/CBioinformatics/PACT/getMethylMatch.sh "${pactRun}" "${runID}"

cd "${WORK_DIR}" || exit

# Consolidated rsync paths --------------------------------------------------------------------------------
MOLEC_VOL="/Volumes/molecular/Molecular/"
DESK_DIR="$HOME/Desktop/${pactRun}/"
PNG_OUT_DIR="${WORK_DIR}/cnvpng/"
VAF_DIR="${MOLEC_VOL}REDCap/cnv_facets/${pactRun}/VAF_plots"
DESK_VAF="${DESK_DIR}VAF_plots"

create_dir "${DESK_DIR}"
create_dir "${PNG_OUT_DIR}"

# COPY FROM: Z-drive TO: DESKTOP --------------------------------------------------------------------------
msg_code rsync -vrthP --include=\"*.png\" \"${MOLEC_VOL}REDCap/cnv_facets/${pactRun}/\" \"${DESK_DIR}\"
msg_code rsync -vrthP \"${MOLEC_VOL}REDCap/cnv_facets/${pactRun}/${pactRun}-QC.tsv\" \"${DESK_DIR}\"
msg_code rsync -vrthP \"${MOLEC_VOL}NGS607/${currYear}/${runID}/${pactRun}_Hotspots.tsv\" \"${DESK_DIR}\"

# Copy demux-samplesheet.csv to desktop
msg_code rsync -vrthP \"${MOLEC_VOL}REDCap/cnv_facets/${pactRun}/demux-samplesheet.csv\" \"${DESK_DIR}\"

# Copy VAF QC output files if availible
if [ -d "${VAF_DIR}" ]; then
    msg_code rsync -vrthP \"${VAF_DIR}\" \"${DESK_DIR}\"
fi

printf "\n${BG_GRN}Executing:${NORMAL}\n"
printf "cp ${DESK_DIR}*.png ${PNG_OUT_DIR}"

# COPY FROM: Desktop TO: Consensus Directory --------------------------------------------------------------
cp "${DESK_DIR}"*.png "${PNG_OUT_DIR}"
msg_code rsync -vrthP --include=\"*.tsv\" \"${DESK_DIR}\" \"${WORK_DIR}\"
msg_code rsync -vrthP \"$HOME/Desktop/${runID}-SampleSheet.csv\" \"${WORK_DIR}\"

# Check if _MethylMatch.xlsx file exists ------------------------------------------------------------------
METH_MATCH="$HOME/Desktop/${pactRun}_MethylMatch.xlsx"
if [ -f "${METH_MATCH}" ]; then
    msg_code rsync -vrthP \"${METH_MATCH}\" \"${WORK_DIR}\"
else
    echo -e "${BG_RED}NO METHYLMATCH FILE COPIED${NORMAL}"
fi

# Copy demux-samplesheet.csv to consensus
msg_code rsync -vrthP \"$HOME/Desktop/${pactRun}/demux-samplesheet.csv\" \"${WORK_DIR}\"

# Copy VAF QC output files if availible
if [ -d "${VAF_DIR}" ]; then
    msg_code rsync -vrP \"${DESK_VAF}\" \"${WORK_DIR}/\"
fi

# EXECUTE: Rscripts for generating HTML Report ------------------------------------------------------------
cd "${HOME}" && curl -# -L ${pactGithub}/MakeIndelList.R >"${HOME}/MakeIndelList.R"
chmod +rwx "${HOME}/MakeIndelList.R"

RScript --verbose "${HOME}/MakeIndelList.R" "${pactRun}"
cd "${WORK_DIR}/" || exit

echo "Running Rscript:"
Rscript --verbose -e "rmarkdown::render('${WORK_DIR}/${pactRun}_consensus.Rmd', params=list(pactName='${pactRun}', userName='${kerbero}'))" && open "${WORK_DIR}/${pactRun}_consensus.html"

# COPY FINAL Output Report -----------------------------------------------------------------------------
echo "Trying to copy html file to ouput folder:"
echo "/Volumes${outputDir}${currYear}/${pactRun}/"

CONSENSUS_FILE="${WORK_DIR}/${pactRun}_consensus.html"
FINAL_DEST="/Volumes${outputDir}${currYear}/${pactRun}/"

# Copy HTML file from consensus directory to DESK_DIR
cp "${CONSENSUS_FILE}" "${DESK_DIR}"
chmod +rwx "${DESK_DIR}${pactRun}_consensus.html"

# Copy HTML file from DESK_DIR to final destination
cd "${FINAL_DEST}" || exit
rsync -vrthP "${DESK_DIR}${pactRun}_consensus.html" .

echo -e "Opening folder:\n${FINAL_DEST}"
open "${FINAL_DEST}"
echo -e "Ensure ${pactRun}_consensus.html is viewable in the output directory"
