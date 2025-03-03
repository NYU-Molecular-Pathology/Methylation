#!/bin/bash
## Script name: zdrive_copier.sh
## Purpose: Copies PACT output files to Molecular Z-drive from DMN
## Date Created: August 28, 2023
## Author: Jonathan Serrano
## Version: 1.0.0
## Copyright (c) NYULH Jonathan Serrano, 2025

# MAIN ARGS: if args are empty assign NULL as default below -------------------
RUN_ID=${1-NULL}   # PACT RUN_ID i.e. 123456_NB501073_0999_AHT3V7BGXK
PACT_ID=${2-NULL}  # PACT Batch ID i.e. PACT-YY-##
KERBERO=${3-$USER} # assign $USER as default

# Console Output Color Variables ----------------------------------------------
BG_GRN="$(tput setab 2)" # makes text background green
BG_RED="$(tput setab 1)" # makes text background red
NORM=$(tput sgr0)        # resets default text color and background

# Enable strict mode, trap errors, and validate inputs ------------------------
set -Eeuo pipefail
trap_error() {
   local exit_code=$?
   local last_li=$1
   if [[ $exit_code -ne 0 ]]; then
      echo -e "\n${BG_RED}ERROR on LAST COMMAND${NORM}:\n$BASH_COMMAND\n"
      echo -e "Failed at line ${last_li} with exit code ${exit_code}\n"
   fi
}

trap 'trap_error $LINENO' ERR

if [[ "$RUN_ID" == "NULL" || "$PACT_ID" == "NULL" ]]; then
   echo "Error: both RUN_ID and PACT_ID must be provided."
   echo "Usage: $0 <RUN_ID> <PACT_ID> [KERBERO]"
   exit 1
fi

# PATHS and Directories -------------------------------------------------------
YEAR_DIR="20${PACT_ID:5:2}"
ZDRIVE_MNT="/mnt/${KERBERO}/molecular"
HPC_RUN_DIR="/gpfs/data/molecpathlab/production/NGS607/${RUN_ID}"
HPC_OUT_DIR="${HPC_RUN_DIR}/output"
HPC_CNV_DIR="${HPC_OUT_DIR}/cnv/FACETS/"
HPC_VAF_DIR="${HPC_OUT_DIR}/VAF_plots"
RESULTS_DIR="MOLECULAR LAB ONLY/NYU PACT Patient Data/Results/Bioinformatics"

ZDRIVE_OUT="${ZDRIVE_MNT}/Molecular/NGS607/${YEAR_DIR}/${RUN_ID}/output/"
ZDRIVE_RES="${ZDRIVE_MNT}/${RESULTS_DIR}/${YEAR_DIR}/${PACT_ID}"
FACETS_DIR="${ZDRIVE_MNT}/Molecular/REDCap/cnv_facets/${PACT_ID}/"
ZDRIVE_RUN="${ZDRIVE_MNT}/Molecular/NGS607/${YEAR_DIR}/${RUN_ID}/"

TMB_OUTPUT="${HPC_OUT_DIR}/annotations.paired.tmb.validation.2callers.tsv"
MSI_OUTPUT="${HPC_OUT_DIR}/clinical/msi_validation.tsv"

# Helper Functions: Message line being executed and colorized -----------------
msg_code() {
   echo -e "\nExecuting:\n${BG_GRN}$1${NORM}\n"
}

msg_dir(){
   NEW_DIR=$1
   msg_code "mkdir -p \"${NEW_DIR}\""
   mkdir -p "${NEW_DIR}"
}

# FUNCTION: Message code and rsync directories --------------------------------
print_rsync() {
   INPUT_DIR=$1
   OUT_PATH=$2
   msg_code "rsync -vrthP ${INPUT_DIR} ${OUT_PATH}"
   rsync -vrthP "${INPUT_DIR}" "${OUT_PATH}"
}

# Once mounted, create the output directories in /MOLECULAR/NGS607/ -----------
msg_dir "${ZDRIVE_OUT}alignments"
msg_dir "${ZDRIVE_RES}"
msg_dir "${ZDRIVE_RES}/TMB_MSI"

# RSYNC: files from BigPurple to /MOLECULAR/NGS607/ ---------------------------
print_rsync "${HPC_OUT_DIR}/alignments/recalibrated" "${ZDRIVE_OUT}alignments/"
print_rsync "${HPC_OUT_DIR}/alignments/raw" "${ZDRIVE_OUT}alignments/"
print_rsync "${HPC_OUT_DIR}/annotations" "${ZDRIVE_OUT}"
print_rsync "${HPC_OUT_DIR}/clinical" "${ZDRIVE_OUT}"
print_rsync "${HPC_OUT_DIR}/CollectHsMetrics" "${ZDRIVE_OUT}"

# COPY MAD values from output if directory exists and ReConCNV 
ERR_MSG="\n${BG_RED}${HPC_OUT_DIR}/MAD_CNV does not exist! Skipping...${NORM}"
print_rsync "${HPC_OUT_DIR}/MAD_CNV" "${ZDRIVE_RES}" || echo -e "${ERR_MSG}"
print_rsync "${HPC_OUT_DIR}/${PACT_ID}_reconCNV" "${ZDRIVE_RES}"
print_rsync "${TMB_OUTPUT}" "${ZDRIVE_RES}/TMB_MSI"
print_rsync "${MSI_OUTPUT}" "${ZDRIVE_RES}/TMB_MSI" || echo -e "${BG_RED}MSI tsv output is missing!${NORM}"

# COPY TMB VARIANTS FOLDER FROM OUTPUT TO CLINICAL ----------------------------
print_rsync "${HPC_OUT_DIR}/TMB_Variants" "${ZDRIVE_RES}"

# COPY TSV files from HPC output to Molecular Run directory -------------------
rsync -vrthP --include="*.tsv" --exclude="*" "${HPC_OUT_DIR}/" "${ZDRIVE_RUN}"

# COPY remaining files in clinical, .PNG, -QC.tsv, and samplesheet ------------
cd "${ZDRIVE_RES}/" || exit 1
print_rsync "${HPC_OUT_DIR}/clinical/" "./"
rsync -vrthP --include="*.png" --exclude="*" "${HPC_CNV_DIR}" "${FACETS_DIR}"
print_rsync "${HPC_RUN_DIR}/${PACT_ID}-QC.tsv" "${FACETS_DIR}"
print_rsync "${HPC_RUN_DIR}/demux-samplesheet.csv" "${FACETS_DIR}"

# Change permissions and remove any duplicate msi_validation.tsv --------------
if [ -d "${HPC_VAF_DIR}" ]; then
   print_rsync "${HPC_VAF_DIR}" "${FACETS_DIR}"
   chmod -R g+rwx "${FACETS_DIR}"
fi

chmod -R g+rwx "${ZDRIVE_RUN}"
chmod -R g+rwx "${ZDRIVE_RES}"
rm -rf "${ZDRIVE_RES}"/msi_validation.tsv
