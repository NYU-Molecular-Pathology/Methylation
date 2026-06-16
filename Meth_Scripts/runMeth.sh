#!/bin/bash
## ---------------------------
## Script name: runMeth.sh
## Purpose: Executes methylation pipeline scripts and copies reports to output directories
## Date Created: August 5, 2021
## Date Last Modified: June 15, 2026
## Author: Jonathan Serrano
## Version: 1.1.0
## Copyright (c) NYULH Jonathan Serrano, 2026
## ---------------------------

# RedCap API Token: Saved Locally (Hardcoded)
methAPI='XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' 

# INPUT ARGS: if args are empty assign NULL as default below ------------------
RUN_ID=${1-NULL}   # methylation run id e.g. 22-MGDM17
PRIORITY=${2-NULL} # string of prioritized RD-numbers
PROJ_DIR=${3-NULL} # any custom directory to copy/run the idat files
redcapUp=${4-NULL} # to upload to redcap or not single char i.e. "T" or "F"
runLocal=${5-NULL} # If the run is local or on shared drive i.e. "T" or "F"

# HARDCODED GITHUB PATHS ------------------------------------------------------
GITHUBMAIN="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R"

# Hardcoded Terminal Color Variables ------------------------------------------
BG_BLU="$(tput setab 4)" # makes text background blue
BG_CYA="$(tput setab 6)" # makes text background cyan
BG_MAG="$(tput setab 5)" # makes text background magenta
BG_GRN="$(tput setab 2)" # makes text background green
BG_RED="$(tput setab 1)" # makes text background red
NORMAL=$(tput sgr0)      # resets default text output

# Helper functions for printing colored messages and error handling -----------
msg_red() { echo "" && echo -e "${BG_RED}$1${NORMAL}\n"; }
msg_grn() { echo "" && echo -e "${BG_GRN}$1${NORMAL}\n"; }
msg_code() { echo -e "\n${2:-Executing:\n}${BG_GRN}${1}${NORMAL}\n"; }
fail_msg() { if [ -z "${2+x}" ]; then msg_red "$1 Continuing execution..."; else msg_red "$1  $2"; fi; }
version_gte() { [ "$(printf '%s\n' "$1" "$2" | sort -V | head -n 1)" = "$2" ]; }
checkSystem() { if ! command -v "$1"; then
	msg_red "$1 is not installed or not in PATH."
	exit 0
fi; }
check_directory() { if [ ! -d "$1" ]; then msg_grn "Creating Dir: $1" && mkdir -p "$1"; fi; }

# FUNCTIONS: Curl download and mkdir if not exists ----------------------------
msg_curl() {
	git_url=$1
	out_fi=$2
	echo -e "${BG_GRN}Downloading file from Github:${NORMAL}\n$HOME/$out_fi\n"
	curl -k -# -L "$git_url/$out_fi" >"$HOME/$out_fi" && chmod +rwx "$HOME/$out_fi"
}

set -Eeuo pipefail
last_cmd="" err_line=0 err_fi=""
trap_error() {
	local ec=$?
	fail_msg "ERROR: ${last_cmd}\nFile: ${err_fi}:${err_line} (Exit: ${ec})" ""
}
trap 'last_cmd=$BASH_COMMAND; err_line=$LINENO; err_fi=${BASH_SOURCE[0]}' DEBUG
trap 'trap_error' ERR

[ -z "$RUN_ID" ] && msg_red "Missing Input: Run ID! Ex. '21-MGDM30'" && exit 1
[ -z "$redcapUp" ] && { redcapUp="F"; }

RUN_DIR="/Volumes/CBioinformatics/Methylation/Clinical_Runs/${RUN_ID}/"
[ ! -n "$PROJ_DIR" ] && RUN_DIR="${PROJ_DIR}/${RUN_ID}/"

cd "$HOME" || {
	msg_red "Failed to cd to: $HOME"
	exit 1
}

LOCALDIR="$HOME/Desktop/Html_${RUN_ID}/${RUN_ID}/"
R_LOCAL="/usr/local/bin/R"
R_CURR="$(which R)"

if [[ "$R_CURR" != "$R_LOCAL" ]]; then
	msg_red "Your system is not using R located in: $R_LOCAL"
	msg_red "\nInstead R is located in: $R_CURR" && exit 0
fi

# Check if it is a Clinical or Research Run
IS_RESEARCH="FALSE"
[[ "$RUN_ID" == MR* ]] && IS_RESEARCH="TRUE"
RUNYEAR=${RUN_ID%%-*}
[ "$IS_RESEARCH" = "TRUE" ] && RUNYEAR=${RUN_ID:2:2}
MOLECOUT="/Volumes/molecular/Molecular/MethylationClassifier/20${RUNYEAR}/"
OUT_DIR="/Volumes/molecular/MOLECULAR LAB ONLY/NYU-METHYLATION/Results/20${RUNYEAR}/"
[ "$IS_RESEARCH" = "TRUE" ] && OUT_DIR="/Volumes/snudem01labspace/FINAL_PDF_Reports_Brain/"

# Get the R version to compare versions, returns 0 if $1 >= $2
version=$(R --version | grep -o 'R version [0-9.]*' | grep -o '[0-9.]*')

# Compare the current R version with 4.3.2
if version_gte "$version" "4.3.2"; then
	echo "R version is $version, which is greater than or equal to 4.3.2"
else
	echo "R version is $version, which is less than 4.3.2"
	echo "Update R version to continue at: https://cran.r-project.org/bin/macosx/" && exit 1
fi

# Check if system requirements installed
unameOut="$(uname -s)"

if [[ "$unameOut" == "Darwin" ]]; then
	checkSystem "latex"
	checkSystem "pdflatex"
	checkSystem "xquartz"
	checkSystem "brew"
	checkSystem "pandoc"
	checkSystem "go"
fi

msg_col() { echo -e "$1$2:${NORMAL} $3"; }

# SANITY CHECK: Prints parameters ---------------------------------------------
echo "~~~~~~~~~~~~~~~~~~~~~Parameters passed to RScript~~~~~~~~~~~~~~~~~~~~~"
msg_col "$BG_BLU" "methAPI"          "$methAPI"
msg_col "$BG_CYA" "(arg1) RUN_ID"    "$RUN_ID"
msg_col "$BG_MAG" "(arg2) PRIORITY"  "$PRIORITY"
msg_col "$BG_GRN" "(arg3) PROJ_DIR"  "$PROJ_DIR"
msg_col "$BG_BLU" "(arg4) redcapUp"  "$redcapUp"

msg_curl ${GITHUBMAIN} "Methyl_QC.Rmd"
msg_curl ${GITHUBMAIN} "methylExpress.R"

echo -e "\nExecuting Rscript with parameters input:\n"
r_params="$methAPI $RUN_ID $PRIORITY $PROJ_DIR $redcapUp $runLocal\n"
msg_col "$BG_BLU" "Rscript --verbose $HOME/methylExpress.R" "$r_params"

#########################################
# EXECUTE: Main command using user input --------------------------------------
#########################################
Rscript --verbose "$HOME/methylExpress.R" $methAPI "$RUN_ID" "$PRIORITY" "$PROJ_DIR" "$redcapUp" "$runLocal"

# SANITY CHECK: Prints output directories -------------------------------------
echo -e "\n----Output Directories----\n"
echo "${BG_GRN}RUN_DIR:${NORMAL}  $RUN_DIR"
echo "${BG_GRN}OUT_DIR:${NORMAL} $OUT_DIR"
echo "${BG_GRN}LOCALDIR:${NORMAL} $LOCALDIR"

# --include '*.html' --exclude='.*' --exclude='*'
copy_htmls() {
	DEST_DIR=$(dirname "$2")
	DIR_NAME=$(basename "$2")
	echo -e "COPY_DIR:\n$1\nOUTPUT_DIR:\n$2"
	cd "$DEST_DIR" || { msg_red "Failed to cd to ${DEST_DIR}"; exit 1; }
	mkdir -p "${DIR_NAME}" && cd "${DIR_NAME}"
	rsync -vrhP --include '*.html' --include '*.csv' --exclude '.*' --exclude '*' "$1" "."
}

# Copy HTML files to results directory, renaming RD-* matches with MP number.
rename_htmls() {
    local src_dir="${1%/}/"
    local dst_dir="${2%/}/"
    local csv_file="${src_dir}samplesheet.csv"

    if [[ ! -f "$csv_file" ]]; then msg_red "CSV not found: ${csv_file}";
        return 1
    fi

    echo -e "${BG_BLU}Copying HTML files to ${dst_dir}:${NORMAL}"
    mkdir -p "$dst_dir"

    # Indexed arrays:
    #   samples[i] -> mps[i]
    #   used[i]    -> whether that samplesheet row has already been matched
    local -a samples mps used fields
    local line sample mp

    # CSV columns: 1 = sample, 8 = MP number
    while IFS= read -r line; do
        IFS=, read -r -a fields <<< "$line"

        sample="${fields[0]//$'\r'/}"

        if [[ "$sample" != RD-* ]]; then continue; fi

        mp="${fields[7]//$'\r'/}"
        mp="${mp//\//-}"
        mp="${mp//[[:space:]]/}"

        samples+=("$sample")
        mps+=("$mp")
        used+=(0)
    done < <(tail -n +2 "$csv_file")

    local html_file base stem suffix new_name
    local i best_i best_mp best_sample best_len

    for html_file in "${src_dir}"*.html; do
        if [[ ! -f "$html_file" ]]; then continue; fi

        base="$(basename "$html_file")"
        stem="${base%.html}"

        best_i=-1
        best_mp=""
        best_sample=""
        best_len=-1

        for ((i = 0; i < ${#samples[@]}; i++)); do
            if [[ "${used[$i]}" == 1 ]]; then continue; fi

            # Match only if:
            #   1. stem is exactly the sample name, or
            #   2. stem starts with the sample name followed by ., _, or -
            #
            # Prevents RD-26-MGDM24-1 from incorrectly matching
            # RD-26-MGDM24-10 or RD-26-MGDM24-12.
            if [[ "$stem" == "${samples[$i]}" || "$stem" == "${samples[$i]}"[._-]* ]]; then
                if (( ${#samples[$i]} > best_len )); then
                    best_i=$i
                    best_len=${#samples[$i]}
                    best_sample="${samples[$i]}"
                    best_mp="${mps[$i]}"
                fi
            fi
        done

        if [[ "$best_i" -ge 0 && -n "$best_mp" ]]; then
            suffix="${stem#"$best_sample"}"
            new_name="${best_sample}_${best_mp}${suffix}.html"
            used[best_i]=1
        else
            new_name="$base"
        fi

        if [[ -e "${dst_dir}${new_name}" ]]; then
            echo "WARN: ${new_name} already exists, skipping ${base}"
        else
            cp "$html_file" "${dst_dir}${new_name}" &&
                echo "Copied ${base} -> ${new_name}"
        fi
    done
}

# COPY: Rsync html to LOCALDIR directory then cp from LOCALDIR to Z-Drive -----
check_directory "$LOCALDIR"

OUTPUT_DIR_1="${OUT_DIR}${RUN_ID}/"
OUTPUT_DIR_2="${MOLECOUT}${RUN_ID}/"

# Check if OUTPUT_DIR_1 exists
if [ -d "$OUTPUT_DIR_1" ]; then
	echo -e "Output files already exist in directory:\n${OUTPUT_DIR_1}"
	echo "Exiting script..." && exit 0
fi

echo -e "\nRsyncing ${RUN_ID} html files to Desktop:\n$LOCALDIR"
rsync -rv --include='*_QC_and_Classifier_Scores.csv' --exclude='*' "${RUN_DIR}" "$LOCALDIR"

if [ "$IS_RESEARCH" != "TRUE" ]; then
	rename_htmls "$RUN_DIR" "$LOCALDIR"
	copy_htmls "${LOCALDIR}" "${OUTPUT_DIR_1}"
	copy_htmls "${LOCALDIR}" "${OUTPUT_DIR_2}"
else
	rsync -rv --include '*.html' --exclude='*' "${RUN_DIR}" "$LOCALDIR"
	copy_htmls "${LOCALDIR}" "${OUTPUT_DIR_1}"
fi

echo -e "${BG_BLU}Verify files are copied here:${NORMAL}\n${OUT_DIR}\n"
echo -e "\n${BG_GRN}METHYLATION RUN ${RUN_ID} COMPLETE${NORMAL}\n"

open "${OUT_DIR}"

redcap_logs="${RUN_DIR}${RUN_ID}_import_log.tsv"
upload_logs="${RUN_DIR}${RUN_ID}_upload_log.tsv"

[ -f "$redcap_logs" ] && open "${redcap_logs}" || echo "No REDCap warnings!"
[ -f "$upload_logs" ] && open "$upload_logs" || echo "No upload warnings!"
