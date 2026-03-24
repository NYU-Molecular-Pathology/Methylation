#!/bin/bash
## Script name: make_consensus.sh
## Purpose: Generate PACT Consensus HTML report with given input
## Author: Jonathan Serrano
## Date Created: August 21, 2023
## Date Modified: March 23, 2026
## Author: Jonathan Serrano
## Version: 1.0.1
## Copyright (c) NYULH Jonathan Serrano, 2024

# Hardcoded paths
GITHUB_URL="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts"
RESULTS_DIR="/Volumes/molecular/MOLECULAR LAB ONLY/NYU PACT Patient Data/Results/Bioinformatics"
SCRIPTS_DIR="/Volumes/CBioinformatics/Bash_Scripts"
DEFAULT_DIR="/Volumes/CBioinformatics/jonathan/pact/consensus"

# Input args ------------------------------------------------------------------
RUN_ID=${1-NULL}    # i.e. 123456_NB501073_0212_AHT3V7BGXK
PACT_ID=${2-NULL}   # i.e. PACT-YY-28
CONSENSUS_DIR="${3:-$DEFAULT_DIR}"
kerbero=${4-$USER}

WORK_DIR="${CONSENSUS_DIR}/${PACT_ID}_consensus"
CURRENT_YR="20${RUN_ID:0:2}"
METH_MATCH="$HOME/Desktop/${PACT_ID}_MethylMatch.xlsx"

# Consolidated rsync paths ----------------------------------------------------
MOLEC_VOL="/Volumes/molecular/Molecular"
TMB_MSI_DIR="${RESULTS_DIR}/${CURRENT_YR}/${PACT_ID}/TMB_MSI"
MSI_TSV="msi_validation.tsv"
TMB_TSV="annotations.paired.tmb.validation.2callers.tsv"

DESK_DIR="$HOME/Desktop/${PACT_ID}"
PNG_OUT_DIR="${WORK_DIR}/cnvpng/"
FACETS_DIR="${MOLEC_VOL}/REDCap/cnv_facets/${PACT_ID}"
VAF_DIR="${FACETS_DIR}/VAF_plots"
DESK_VAF="${DESK_DIR}/VAF_plots"
TMB_MSI_OUT="${WORK_DIR}/TMB_MSI/"

NGS607_DIR="${MOLEC_VOL}/NGS607/${CURRENT_YR}/${RUN_ID}"

# Color Variables -------------------------------------------------------------
BG_GRN="$(tput setab 2)" # makes text background green
BG_RED="$(tput setab 1)" # makes text background red
BG_BLU="$(tput setab 4)" # makes text background blue
NORMAL=$(tput sgr0)      # resets default text

# Function to exit when any command fails
set -Eeuo pipefail; last_cmd="" err_line=0 err_fi=""
msg_red() { echo "" && echo -e "${BG_RED}$1${NORMAL}\n"; }
msg_err() { local ec=$?; msg_red "ERROR: ${last_cmd}\nFile: ${err_fi}:${err_line} (Exit: ${ec})"; }
trap 'last_cmd=$BASH_COMMAND; err_line=$LINENO; err_fi=${BASH_SOURCE[0]}' DEBUG
trap 'msg_err' ERR

# Function to display and execute code
msg_rsync() {
    local flags=${3:-}
    echo -e "\n${BG_BLU}rsync:${NORMAL}\nrsync -vrhP ${flags} \"$1\" \"$2\"\n"
    rsync -vrhP "${flags}" "$1" "$2"
}
msg_curl() {
    out_fi=$1; out_dir=${2:-$HOME}; outpath="${out_dir}/${out_fi}"
    cd "$out_dir" || exit
    echo -e "${BG_GRN}Downloading file from Github:${NORMAL}\n$outpath\n"
    curl -k -# -L "$GITHUB_URL/$out_fi" >"$outpath" && chmod +rwx "$outpath"
    msg_red "--------------END--------------"
}

# Main execution --------------------------------------------------------------
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}" || exit

curl -# -L ${GITHUB_URL}/PACT_consensus.Rmd >"${WORK_DIR}/${PACT_ID}_consensus.Rmd"
msg_curl "MakeIndelList.R"
msg_curl "hs_metric_consensus.R"

/Volumes/CBioinformatics/PACT/getMethylMatch.sh "${PACT_ID}" "${RUN_ID}" & wait

cd "${WORK_DIR}" || exit

mkdir -p "${DESK_DIR}"
mkdir -p "${DESK_DIR}/TMB_MSI"
mkdir -p "${PNG_OUT_DIR}"
mkdir -p "${TMB_MSI_OUT}"

cd "$HOME" || exit
# Z-drive TO Desktop
msg_rsync "${FACETS_DIR}/" "${DESK_DIR}" "--include=*.png"
msg_rsync "${FACETS_DIR}/${PACT_ID}-QC.tsv" "${DESK_DIR}"
msg_rsync "${NGS607_DIR}/${PACT_ID}_Hotspots.tsv" "${DESK_DIR}"
msg_rsync "${TMB_MSI_DIR}/${MSI_TSV}" "${DESK_DIR}/TMB_MSI/${MSI_TSV}"
msg_rsync "${TMB_MSI_DIR}/${TMB_TSV}" "${DESK_DIR}/TMB_MSI/${TMB_TSV}"

# Copy VAF QC output files if availible
[ -d "${VAF_DIR}" ] && msg_rsync "${VAF_DIR}" "${DESK_DIR}"

# Desktop TO Consensus
cp "${DESK_DIR}/"*.png "${PNG_OUT_DIR}"
msg_rsync "${DESK_DIR}/" "${WORK_DIR}" "--include=*.tsv"
msg_rsync "${DESK_DIR}/TMB_MSI/" "${TMB_MSI_OUT}" "--include=*.tsv"

# Check if _MethylMatch.xlsx file exists --------------------------------------
if [ -f "${METH_MATCH}" ]; then
    msg_rsync "${METH_MATCH}" "${WORK_DIR}"
else
    msg_red "NO METHYLMATCH FILE COPIED!"
fi

# Copy demux-samplesheet.csv to consensus
msg_rsync "$HOME/Desktop/${PACT_ID}/demux-samplesheet.csv" "${WORK_DIR}"

# Copy VAF QC output files if availible
[ -d "${VAF_DIR}" ] && rsync -vrP "${DESK_VAF}" "${WORK_DIR}/"

# EXECUTE: Rscripts for generating HTML Report --------------------------------
msg_red "Checking if GOS idat files need to be copied..."
source "${SCRIPTS_DIR}/copy_marcin_idats.sh" "${PACT_ID}" || true

RScript --verbose "${HOME}/MakeIndelList.R" "${PACT_ID}"
RScript --verbose "${HOME}/hs_metric_consensus.R" "${RUN_ID}" "${PACT_ID}"

cd "${WORK_DIR}/" || exit

echo "Running rmarkdown:"
Rscript --verbose -e "rmarkdown::render('${WORK_DIR}/${PACT_ID}_consensus.Rmd', params=list(pactName='${PACT_ID}', userName='${kerbero}'))" && open "${WORK_DIR}/${PACT_ID}_consensus.html"

# COPY FINAL Output Report ----------------------------------------------------
CONSENSUS_OUT="${WORK_DIR}/${PACT_ID}_consensus.html"
FINAL_DEST="${RESULTS_DIR}/${CURRENT_YR}/${PACT_ID}/"
echo "Copying to: ${FINAL_DEST}"
[ ! -f "$CONSENSUS_OUT" ] && { echo "Consensus failed to generate: ensure all input files are copied to ${WORK_DIR}"; exit 1; }

# Copy HTML file from consensus directory to DESK_DIR
cp "${CONSENSUS_OUT}" "${DESK_DIR}/"
chmod +rwx "${DESK_DIR}/${PACT_ID}_consensus.html"

# Copy HTML file from DESK_DIR to final destination
cd "${FINAL_DEST}" || exit
rsync -vP "${DESK_DIR}/${PACT_ID}_consensus.html" .

echo -e "Opening folder:\n${FINAL_DEST}"
open "${FINAL_DEST}"
echo -e "Ensure ${PACT_ID}_consensus.html is viewable in the output directory"

# Check if PHILIP_LOG exists and open it if it does
PHILIP_LOG="$HOME/Desktop/${PACT_ID}_philips_missing.txt"

if [ -f "$PHILIP_LOG" ]; then
    echo "File $PHILIP_LOG exists. Some samples missing from data dump..."
    open "$PHILIP_LOG"
else
    echo "File $PHILIP_LOG does not exist. All data dump files downloaded"
fi
