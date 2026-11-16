#!/bin/bash
## Script name: runMeth.sh
## Purpose: Download methylation pipeline scripts from Github and execute run
## Date Created: August 5, 2021
## Date Modified: March 20, 2026
## Author: Jonathan Serrano
## Version: 1.0.2
## Copyright (c) NYULH Jonathan Serrano, 2026

# RedCap API Token ------------------------------------------------------------
IFS= read -r REDCAP_API_TOKEN < "/Volumes/CBioinformatics/scripts/METH_DB_API.txt"

# INPUT ARGS: if args are empty assign NULL as default below ------------------
METH_ID=${1-NULL} # methylation run id e.g. 22-MGDM17

if [[ -z "$METH_ID" || "$METH_ID" == "NULL" ]]; then
    msg_red "Missing input: Run ID!" && exit 1
fi

# HARDCODED GITHUB PATHS ------------------------------------------------------
GH_MAIN="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R"

# Hardcoded Terminal Color Variables ------------------------------------------
BG_BLU="$(tput setab 4)"; BG_CYA="$(tput setab 6)"
BG_GRN="$(tput setab 2)"; BG_RED="$(tput setab 1)"
NORMAL="$(tput sgr0)"

# Helper functions for printing colored messages and error handling -----------
msg_red() { echo "" && echo -e "${BG_RED}$1${NORMAL}\n"; }
msg_grn() { echo "" && echo -e "${BG_GRN}$1${NORMAL}\n"; }
msg_code() { echo -e "\n${2:-Executing:\n}${BG_GRN}${1}${NORMAL}\n"; }
fail_msg() { if [ -z "${2+x}" ]; then msg_red "$1 Continuing..."; else msg_red "$1  $2"; fi; }

version_gte() { [ "$(printf '%s\n' "$1" "$2" | sort -V | head -n 1)" = "$2" ]; }
check_sys() { if ! command -v "$1"; then msg_red "$1 is not installed!" && exit 0; fi; }
check_dir() { if [ ! -d "$1" ]; then msg_grn "Creating Dir: $1" && mkdir -p "$1"; fi; }

# FUNCTIONS: Curl download and mkdir if not exists ----------------------------
msg_curl() {
	local git_url="$1" out_fi="$2"
	echo -e "${BG_GRN}Downloading file from Github:${NORMAL}\n$HOME/$out_fi\n"
	curl -k -# -L "$git_url/$out_fi" >"$HOME/$out_fi" && chmod +rwx "$HOME/$out_fi"
}

set -Eeuo pipefail
last_cmd="" err_line=0 err_fi=""
trap_error() { local ec=$?; fail_msg "ERR: ${last_cmd}\nFile: ${err_fi}:${err_line} ${ec}" ""; }

trap 'last_cmd=$BASH_COMMAND; err_line=$LINENO; err_fi=${BASH_SOURCE[0]}' DEBUG
trap 'trap_error' ERR

[ -z "$METH_ID" ] && msg_red "Missing Input: Run ID! Ex. '21-MGDM30'" && exit 1

RUN_DIR="/Volumes/CBioinformatics/Methylation/Clinical_Runs/${METH_ID}/"

cd "$HOME" || { msg_red "Failed to cd to: $HOME" && exit 1; }

DESK_DIR="$HOME/Desktop/Html_${METH_ID}/${METH_ID}/"
R_LOCAL="/usr/local/bin/R"
R_CURR="$(which R)"

if [[ "$R_CURR" != "$R_LOCAL" ]]; then
	msg_red "Your system is not using R located in: $R_LOCAL"
	msg_red "\nInstead R is located in: $R_CURR" && exit 0
fi

# Check if it is a Clinical or Research Run
IS_RESEARCH="FALSE"
RESEARCH_OUT="/Volumes/snudem01labspace/FINAL_PDF_Reports_Brain"

[[ "$METH_ID" == MR* ]] && IS_RESEARCH="TRUE"

RUNYEAR=${METH_ID%%-*}
[ "$IS_RESEARCH" = "TRUE" ] && RUNYEAR=${METH_ID:2:2}

CLINICAL_OUT="/Volumes/molecular/Molecular/MethylationClassifier/20${RUNYEAR}"
MOL_RESULTS_DIR="/Volumes/molecular/MOLECULAR LAB ONLY/NYU-METHYLATION/Results/20${RUNYEAR}"

# Get the R version to compare versions, returns 0 if $1 >= $2
version=$(R --version | grep -o 'R version [0-9.]*' | grep -o '[0-9.]*')

# Compare the current R version with 4.3.2
if version_gte "$version" "4.3.2"; then
	echo "R version is $version, which is greater than or equal to 4.3.2"
else
	echo "R version is $version, which is less than 4.3.2"
	echo "Update R version: https://cran.r-project.org/bin/macosx/" && exit 1
fi

# Check if system requirements installed
unameOut="$(uname -s)"

if [[ "$unameOut" == "Darwin" ]]; then
	check_sys "latex"
	check_sys "pdflatex"
	check_sys "xquartz"
	check_sys "brew"
	check_sys "pandoc"
	check_sys "go"
fi

msg_col() { echo -e "$1$2:${NORMAL} $3"; }

# SANITY CHECK: Prints parameters ---------------------------------------------
echo "~~~~~~~~~~~~~~~~~~~~~Parameters passed to RScript~~~~~~~~~~~~~~~~~~~~~"
msg_col "$BG_BLU" "REDCAP_API_TOKEN" "$REDCAP_API_TOKEN"
msg_col "$BG_CYA" "(arg1) METH_ID"   "$METH_ID"

msg_curl ${GH_MAIN} "Methyl_QC.Rmd"
msg_curl ${GH_MAIN} "methylExpress.R"

echo -e "\nExecuting Rscript with the following args:\n"
r_params="$REDCAP_API_TOKEN $METH_ID\n"
msg_col "$BG_BLU" "Rscript --verbose $HOME/methylExpress.R" "$r_params"

#########################################
# EXECUTE: Main command using user input --------------------------------------
#########################################
Rscript --verbose "$HOME/methylExpress.R" "$REDCAP_API_TOKEN" "$METH_ID" 

# SANITY CHECK: Prints output directories -------------------------------------
echo -e "\n----Output Directories----\n"
echo "${BG_GRN}RUN_DIR:${NORMAL}  $RUN_DIR"
echo "${BG_GRN}MOL_RESULTS_DIR:${NORMAL} $MOL_RESULTS_DIR"
echo "${BG_GRN}DESK_DIR:${NORMAL} $DESK_DIR"

# --include '*.html' --exclude='.*' --exclude='*'
copy_htmls() {
	DEST_DIR=$(dirname "$2")
	echo -e "COPY_DIR:\n$1\nOUTPUT_DIR:\n$2"
	cd "$DEST_DIR" || { msg_red "Failed to cd to ${DEST_DIR}"; exit 1; }
	mkdir -p "${METH_ID}" && cd "$2"
	rsync -vrh --include '*html' --include '*csv' --exclude '.*' --exclude '*' "$1" "."

	QC_FILE1="${2%/}/${METH_ID}_qc_data.csv"
	QC_FILE2="${2%/}/${METH_ID}_QC_and_Classifier_Scores.csv"

	if [[ ! -f "${QC_FILE1}" || ! -f "${QC_FILE2}" ]]; then
        cp -v "$1"/*.{html,csv} "$2"
    fi

	if [[ ! -f "${QC_FILE1}" || ! -f "${QC_FILE2}" ]]; then
        msg_red "Missing QC files: ${QC_FILE1} or ${QC_FILE2}"
        exit 1
    fi
}

rename_htmls() {
    local src="${1%/}/" dst="${2%/}/"
    local line sample mp file base name
    local -a fields
	sam_sheet="${src}samplesheet.csv"
    if [[ ! -f "$sam_sheet" ]]; then msg_red "Missing CSV: $sam_sheet" && exit 1; fi

    echo -e "${BG_BLU}Renaming HTML files to Desktop: ${dst}:${NORMAL}"
    mkdir -p "$dst"

    while IFS= read -r line; do
        IFS=, read -r -a fields <<< "$line"
        sample="${fields[0]//$'\r'/}"
        [[ "$sample" == RD-* ]] || continue

        mp="${fields[7]//$'\r'/}"
        mp="${mp//\//-}"
        mp="${mp//[[:space:]]/}"

        file="${src}${sample}.html"
        [[ -f "$file" ]] || file="${src}${sample}_QC_FAILED.html"
        [[ -f "$file" ]] || continue

        base="${file##*/}"
        name="${sample}_${mp}.html"

        if [[ -e "${dst}${name}" ]]; then
            echo "${name} already exists, skipping ${base}"
        else
            cp "$file" "${dst}${name}" && echo "Copied ${base} -> ${name}"
        fi
    done < <(tail -n +2 "$sam_sheet")

	qc_file="${src}${METH_ID}_QC.html"

    if [[ ! -f "$qc_file" ]]; then
		msg_red "Missing QC file: $qc_file"
	else
		cp "$qc_file" "$dst" && echo "Copied $qc_file -> $dst"
	fi
}

require_smb_mount() {
	local mountpoint="$1"
	if [[ ! -d "$mountpoint" ]]; then msg_red "SMB mount missing: $mountpoint"; exit 1; fi
	if [[ ! -w "$mountpoint" ]]; then msg_red "Mount is not writable: $mountpoint"; exit 1; fi
}

mountpoint="/Volumes/molecular/Molecular"
[[ "$IS_RESEARCH" == "TRUE" ]] && mountpoint="/Volumes/snudem01labspace/idats"
require_smb_mount "$mountpoint"

# COPY: Rsync html to DESK_DIR directory then cp from DESK_DIR to Z-Drive -----
check_dir "$DESK_DIR"

OUTPUT_DIR_1="${MOL_RESULTS_DIR}/${METH_ID}/"
OUTPUT_DIR_2="${CLINICAL_OUT}/${METH_ID}/"

htmls_out=""
if [[ -d "$OUTPUT_DIR_1" ]]; then
    htmls_out=$(find "$OUTPUT_DIR_1" -maxdepth 1 -type f -name '*.html' -print -quit)
fi

if [[ -n "$htmls_out" ]]; then
    echo -e "Output already exists in:\n${OUTPUT_DIR_1}\nExiting..." && exit 0
fi

echo -e "\nRsyncing ${METH_ID} html files to Desktop:\n$DESK_DIR"
rsync -rv --include='*_QC_and_Classifier_Scores.csv' --exclude='*' "${RUN_DIR}" "$DESK_DIR"

# If a run is clinical same results copy to two output directories
if [ "$IS_RESEARCH" != "TRUE" ]; then
	rename_htmls "$RUN_DIR" "$DESK_DIR"
	rsync -rv  --include='*control*.html' --exclude='*' "${RUN_DIR}" "$DESK_DIR"
	rsync -rv --include='*QC*.csv' --include='*qc*.csv' --exclude='*' "${RUN_DIR}" "$DESK_DIR"
	
	copy_htmls "${DESK_DIR}" "${OUTPUT_DIR_1}"
	copy_htmls "${DESK_DIR}" "${OUTPUT_DIR_2}"

else
	rsync -rv --include '*.html' --exclude='*' "${RUN_DIR}" "$DESK_DIR"
	copy_htmls "${DESK_DIR}" "${RESEARCH_OUT}/${METH_ID}/"
fi

echo -e "${BG_BLU}Verify files are copied here:${NORMAL}\n${MOL_RESULTS_DIR}\n"
echo -e "\n${BG_GRN}METHYLATION RUN ${METH_ID} COMPLETE${NORMAL}\n"

open "${MOL_RESULTS_DIR}"

redcap_logs="${RUN_DIR}${METH_ID}_import_log.tsv"
upload_logs="${RUN_DIR}${METH_ID}_upload_log.tsv"

if [[ -f "$redcap_logs" ]]; then open "$redcap_logs"; else echo "No REDCap import warnings!"; fi
if [[ -f "$upload_logs" ]]; then open "$upload_logs"; else echo "No upload warnings!"; fi
