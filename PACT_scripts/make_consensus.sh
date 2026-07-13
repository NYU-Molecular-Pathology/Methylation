#!/bin/bash
## Script name: make_consensus.sh
## Purpose: Generate PACT Consensus HTML report with given input
## Author: Jonathan Serrano
## Date Created: August 21, 2023
## Author: Jonathan Serrano
## Version: 1.0.1
## Copyright (c) NYULH Jonathan Serrano, 2026

# shellcheck disable=SC1091

# Hardcoded paths
GIT_URL="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts"
RESULTS_DIR="/Volumes/molecular/MOLECULAR LAB ONLY/NYU PACT Patient Data/Results/Bioinformatics"
SCRIPTS_DIR="/Volumes/CBioinformatics/Bash_Scripts"
DEFAULT_DIR="/Volumes/CBioinformatics/jonathan/pact/consensus"

# Input args ------------------------------------------------------------------
RUN_ID=${1-NULL}
PACT_ID=${2-NULL}
CONSENSUS_DIR="${3:-$DEFAULT_DIR}"
kerbero=${4-$USER}

[[ "$RUN_ID" == "NULL" || "$PACT_ID" == "NULL" ]] && 
	{ echo "<RUN_ID> or <PACT_ID> cannot be NULL" && exit 1; }

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
DEMUX_CSV="${RESULTS_DIR}/${CURRENT_YR}/${PACT_ID}/demux-samplesheet.csv"

# Color Variables -------------------------------------------------------------
BG_RED="$(tput setab 1)" # makes text background red
BG_BLU="$(tput setab 4)" # makes text background blue
BG_YLW="$(tput setab 3)$(tput setaf 0)" # yellow background, black text
NORMAL=$(tput sgr0)      # resets default text

# trap and notify for ERRORS --------------------------------------------------
set -Eeuo pipefail
msg_red() { echo -e "\n${BG_RED}$1${NORMAL}\n"; }
trap 'msg_red "ERROR:$BASH_COMMAND\nFILE: ${BASH_SOURCE[0]}\nLINE:$LINENO"' ERR

# Function to display and execute code
msg_code() { echo -e "Executing:\n${BG_BLU}$*${NORMAL}\n"; "$@"; }
msg_rsync() { msg_code rsync -vrhP ${3:+$3 }"$1" "$2"; }
msg_curl() {
    out_fi=$1; out_dir=${2:-$HOME}; dest="${out_dir}/${out_fi}"
    cd "$out_dir" || exit
    msg_code curl -k -# -L "$GIT_URL/$out_fi" -o"$dest" && chmod +rwx "$dest"
}

# Main execution --------------------------------------------------------------
mkdir -p "${WORK_DIR}" "${DESK_DIR}/TMB_MSI"
cd "${WORK_DIR}" || exit
mkdir -p "${PNG_OUT_DIR}" "${TMB_MSI_OUT}"

# Download Rmd and scripts from Github
RMD_FILE="${WORK_DIR}/${PACT_ID}_consensus.Rmd"

msg_code curl -# -L ${GIT_URL}/PACT_consensus.Rmd -o "${RMD_FILE}"
msg_curl "MakeIndelList.R"
msg_curl "hs_metric_consensus.R"

msg_code /Volumes/CBioinformatics/PACT/getMethylMatch.sh "${PACT_ID}" "${RUN_ID}"

# Z-drive TO Desktop
cd "$HOME" || exit
msg_rsync "${FACETS_DIR}/" "${DESK_DIR}" "--include=*.png --exclude=*"
msg_rsync "${FACETS_DIR}/${PACT_ID}-QC.tsv" "${DESK_DIR}"
msg_rsync "${NGS607_DIR}/${PACT_ID}_Hotspots.tsv" "${DESK_DIR}"
msg_rsync "${TMB_MSI_DIR}/${MSI_TSV}" "${DESK_DIR}/TMB_MSI/${MSI_TSV}"
msg_rsync "${TMB_MSI_DIR}/${TMB_TSV}" "${DESK_DIR}/TMB_MSI/${TMB_TSV}"

# Copy VAF QC output files if availible
[ -d "${VAF_DIR}" ] && msg_rsync "${VAF_DIR}" "${DESK_DIR}"

# Desktop TO Consensus
if [ ! -w "${PNG_OUT_DIR}" ]; then
    msg_red "Permission denied! Remove or rename the previous consenus folder:\n${WORK_DIR}" && exit 1
fi
cp "${DESK_DIR}/"*.png "${PNG_OUT_DIR}"
msg_rsync "${DESK_DIR}/" "${WORK_DIR}" "--include=*.tsv"
msg_rsync "${DESK_DIR}/TMB_MSI/" "${TMB_MSI_OUT}" "--include=*.tsv"

# Check if _MethylMatch.xlsx file exists --------------------------------------
[ -f "${METH_MATCH}" ] && msg_rsync "${METH_MATCH}" "${WORK_DIR}"

# Copy demux-samplesheet.csv to consensus
msg_rsync "${DEMUX_CSV}" "$HOME/Desktop/${PACT_ID}/"
msg_rsync "$HOME/Desktop/${PACT_ID}/demux-samplesheet.csv" "${WORK_DIR}"

# Copy VAF QC output files if availible
[ -d "${VAF_DIR}" ] && rsync -vrP "${DESK_VAF}" "${WORK_DIR}/"

# EXECUTE: Rscripts for generating HTML Report --------------------------------
echo -e "${BG_YLW}Checking if GOS idat files need to be copied...${NORMAL}"
( source "${SCRIPTS_DIR}"/copy_marcin_idats.sh "${RUN_ID}" "${PACT_ID}" ) || true

msg_code RScript --verbose "${HOME}/MakeIndelList.R" "${PACT_ID}"
msg_code RScript --verbose "${HOME}/hs_metric_consensus.R" "${RUN_ID}" "${PACT_ID}"

msg_code cd "${WORK_DIR}/" || exit 1
msg_code Rscript --verbose -e "rmarkdown::render('${RMD_FILE}', params=list(pactName='${PACT_ID}', userName='${kerbero}'))" 

# COPY FINAL Output Report ----------------------------------------------------
CONSENSUS_OUT="${WORK_DIR}/${PACT_ID}_consensus.html"
FINAL_DEST="${RESULTS_DIR}/${CURRENT_YR}/${PACT_ID}/"

[ ! -f "$CONSENSUS_OUT" ] && { echo "Consensus failed to generate: ensure all input files are copied to ${WORK_DIR}"; exit 1; }

open "${CONSENSUS_OUT}"

# Copy HTML file from consensus directory to DESK_DIR
cp "${CONSENSUS_OUT}" "${DESK_DIR}/"
chmod +rwx "${DESK_DIR}/${PACT_ID}_consensus.html"

# Copy HTML file from DESK_DIR to final destination
cd "${FINAL_DEST}" || exit
msg_code rsync -vP "${DESK_DIR}/${PACT_ID}_consensus.html" .

echo -e "Ensure ${PACT_ID}_consensus.html is viewable in the output directory"
open "${FINAL_DEST}"
