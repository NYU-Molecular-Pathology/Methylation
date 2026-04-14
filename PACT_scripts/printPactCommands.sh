#!/bin/bash
## ---------------------------
## Script name: printPactCommands.sh
## Purpose: Print out all the copiable commands used for executing the PACT pipeline
## Author: Jonathan Serrano
## Version: 1.0.1
## Date Updated: March 25, 2026
## Copyright (c) NYULH Jonathan Serrano, 2026
## ---------------------------

DEFAULT_DIR="/Volumes/CBioinformatics/jonathan/pact/consensus/"

RUN_ID=${1-NULL}  		# i.e. 123456_NB501073_0212_AHT3V7BGXK
PACT_ID=${2-NULL}		# i.e. 26-PACT33
CONSENSUS_DIR=${3-NULL} # path to where consensus generates
kerbero=${4-$USER} 		# i.e. whoami kerberosid

[[ -n "${3}" ]] && CONSENSUS_DIR="${DEFAULT_DIR}"
IS_SOPHIA=$([[ "${PACT_ID:0:2}" =~ ^[0-9]{2}$ ]] && echo true || echo false)
BASH_HELPERS="/gpfs/data/molecpathlab/scripts/bash_helpers"

if [[ "$IS_SOPHIA" == "true" ]]; then
	year_part=${PACT_ID:0:2}
	BASH_HELPERS="/gpfs/data/molecpathlab/scripts/bash_helpers_SG"
else
	year_part=${PACT_ID:5:2}
fi

VERS="1.0.1"
FG_YLW='<span style="color:#cf6a00;">' # makes text darkorange color
FG_RED='<span style="color:#cc0000">'  # makes text red color
FG_GRN='<span style="color:#6aa84f">'  # makes text color green
FG_BLU='<span style="color:dodgerblue">'     # makes text color blue
FG_CYA='<span style="color:cyan;">'    # makes text color #baffc9
FG_MAG='<span style="color:#ff00ff;">' # makes text color magenta
WHT_BG="<span style='background-color:white;margin-left:30px;padding:3px;margin-top:5px;margin-bottom:5px;line-height:1.2!important;'>"
NORMAL="</span>" # resets default text
BOX1=' <div class="boxed"> '
BOX2=' </div> '

if [[ "$IS_SOPHIA" == "true" ]]; then
	year_part=${PACT_ID:0:2}
else
	year_part=${PACT_ID:5:2}
fi

TD_DATE="$(date +"%B %d, %Y %-I:%M%P %Z")"
PACT_ID="${FG_RED}${PACT_ID}${NORMAL}"
kerbero="${FG_MAG}${kerbero}${NORMAL}"
RUN_ID="${FG_YLW}${RUN_ID}${NORMAL}"

YEAR_DIR="${FG_BLU}20${year_part}${NORMAL}"
PROD_DIR="/gpfs/data/molecpathlab/production"
outputDir="/molecular/MOLECULAR LAB ONLY/NYU PACT Patient Data/Results/Bioinformatics/"
DEMUXDIR="${PROD_DIR}/Demultiplexing/${RUN_ID}"

echo "
<style>
:root {
    --bg: #f5f7fb; --surface: #ffffff; --surface-2: #f8fafc;
    --text: #18212f; --muted: #5b6472; --border: #d9e1ea;
    --accent: #2563eb; --accent-hover: #1d4ed8; --success: #15803d;
    --shadow: 0 8px 24px rgba(15, 23, 42, 0.08);
    --radius: 14px; --radius-sm: 10px;
    --code-bg: #111827; --code-text: #f3f4f6;
}

*, *:before, *:after {box-sizing: border-box;}

html {font-size: 16px;}

html body {
    width: auto; color: var(--text);
    font-family: 'Helvetica Neue', Helvetica, 'Segoe UI', Arial, freesans, sans-serif;
}

body {
    margin: 0; padding: 24px; width: auto;
    background: #626c79 !important; color: var(--text) !important;
    display: block; white-space: normal;
}

main {display: grid; max-width: 1100px; margin: 20px auto;}

.markdown-preview {
    margin: 0 auto !important; max-width: 1100px !important;
    width: 100% !important; height: 100% !important;
    box-sizing: border-box !important; padding-top: 4px !important;
}

h1, h2, h3 {
    width: auto; color: var(--text);
    line-height: 1.2 !important; padding-bottom: 0 !important;
    margin-top: 0 !important; margin-bottom: 0.6rem !important;
}

h1 {font-size: 1.3rem;}

h5 {
    padding: 2px !important; margin-top: 4px !important;
    margin-bottom: 4px !important; color: var(--text);
}

ol {line-height: 1.4 !important; margin-bottom: 0 !important;}

hr {
    border: 0; border-top: 1px solid #e5e7eb;
    margin: 0.8rem 0 1rem 0;
}

.boxed {
    background: rgb(90, 90, 90) !important; border: 3px solid var(--border);
    padding: 18px 20px; border-radius: var(--radius);
    width: 100%; display: block; white-space: normal;
    font-size: 12px !important; box-shadow: var(--shadow);
    margin: 0.9rem 0 1.25rem 0;
}

.tocbox {
    background: #eef3f7 !important; border: 1px solid #b8c4cf;
    padding: 14px 18px; border-radius: 12px;
    width: auto; max-width: 720px; display: inline-block;
    white-space: normal; font-size: 15px !important; line-height: 1.4;
    box-shadow: 0 4px 14px rgba(15, 23, 42, 0.08);
    margin: 0.9rem auto 1.25rem auto;
}

.toc-title {
    margin: 0 0 8px 0; font-size: 24px; font-weight: 800;
    letter-spacing: 0.02em; color: #24384a;
}

.tocbox ol, .tocbox ul {margin: 0; padding-left: 1.2rem;}
.tocbox li {
    margin: 0.28rem 0; line-height: 1.4; color: #1f2937;
    word-break: break-word; overflow-wrap: anywhere;
}
.tocbox a {color: #1e3a5f; text-decoration: none; font-weight: 600;}
.tocbox a:hover {color: #2563eb; text-decoration: underline;}

.input-table-wrap {
    display: inline-block; margin-top: 8px; margin-bottom: 12px;
    background: #2b2521; border: 1px solid #5a4c42;
    border-radius: 10px; box-shadow: 0 4px 14px rgba(0, 0, 0, 0.22);
    overflow: hidden;
}

.input-table {
    border-collapse: collapse; width: auto; min-width: 480px;
    font-size: 13px; line-height: 1.25;
}

.input-table th, .input-table td {
    padding: 6px 10px; text-align: left; vertical-align: top;
    border-bottom: 1px solid #4c4037;
}

.input-table tr:last-child th, .input-table tr:last-child td {border-bottom: none;}

.input-table th {
    background: #3a322d; color: #eadfce; font-weight: 700;
    white-space: nowrap; width: 160px; border-right: 1px solid #5a4c42;
}

.input-table td {
    background: #2b2521; color: #f5eee6;
    word-break: break-word; overflow-wrap: anywhere; min-width: 260px;
}

.input-table-title {
    background: #3a322d !important; color: #f2e7d5 !important; text-align: center !important;
    font-weight: 800; font-size: 16px; letter-spacing: 0.02em;
    padding: 10px 12px !important; border-bottom: 1px solid #5a4c42 !important;
}

code {
    background-color: transparent !important; color: var(--code-text) !important;
    margin-left: 0 !important; padding: 0 !important; line-height: 1.6 !important;
    font-family: Menlo, Monaco, Consolas, 'Courier New', monospace !important;
    font-size: 13px !important;
}

pre, pre[class*=\"language-\"] {
    position: relative; overflow: auto; margin: 8px 0 14px 0;
    padding: 1rem 1rem 1rem 7.5rem;
    border-radius: 12px; background: var(--code-bg) !important;
    border: 1px solid #1f2937;
    box-shadow: inset 0 1px 0 rgba(255,255,255,0.03);
}

button {
    align-items: center; appearance: none; box-sizing: border-box;
    background: rgba(255,255,255,0.08); border: 1px solid rgba(255,255,255,0.14);
    border-radius: 10px; box-shadow: none; color: #ffffff;
    cursor: pointer; display: inline-flex;
    font-family: 'Helvetica Neue', Helvetica, 'Segoe UI', Arial, sans-serif;
    font-size: 12px; font-weight: 600; height: auto;
    justify-content: center; line-height: 1.2; padding: 6px 10px;
    position: absolute; top: 10px; left: 10px;
    text-align: center; text-decoration: none;
    transition: background-color 0.18s ease, border-color 0.18s ease, transform 0.18s ease;
    user-select: none; -webkit-user-select: none; touch-action: manipulation;
    vertical-align: top; white-space: nowrap;
}

button:hover {background: rgba(255,255,255,0.16); border-color: rgba(255,255,255,0.22);}
button:active, button:focus {outline: none;}
button:focus:not(:active) {box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.28);}

.pressed {
    background-color: #06a96e;
    background-image: linear-gradient(1deg, #00aa6c, #14C667 99%);
    border-color: rgba(34, 197, 94, 0.35) !important;
}

.text-block {
    position: relative; margin: 8px 0 14px 0;
    padding: 14px 16px 14px 16px; border-radius: 10px;
    background: #ffffff !important; border: 1px solid #d9e1ea;
    color: #000000 !important;
    font-family: Calibri, Aptos, 'Segoe UI', Arial, sans-serif !important;
    font-size: 14px; line-height: 1.45; white-space: pre-wrap;
    word-break: break-word; overflow-wrap: anywhere;
    box-shadow: 0 2px 8px rgba(15, 23, 42, 0.04);
}

.text-copy-btn {
    position: absolute; top: 10px; left: 10px;
    align-items: center; appearance: none; display: inline-flex;
    background: #e2e8f0; border: 1px solid #cbd5e1; border-radius: 8px;
    color: #0f172a; cursor: pointer;
    font-family: 'Helvetica Neue', Helvetica, 'Segoe UI', Arial, sans-serif;
    font-size: 12px; font-weight: 600; line-height: 1.2;
    padding: 5px 10px; white-space: nowrap;
}

.text-copy-btn:hover {background: #cbd5e1;}
.text-copy-btn:focus, .text-copy-btn:active {outline: none;}

.text-copy-btn.pressed {
    background: #dcfce7 !important; border-color: #86efac !important;
    color: #166534 !important;
}

.text-copy-pad {padding-top: 38px;}

.page-title {
    padding-top: 8px !important; margin: 0 0 6px 0 !important; text-align: center;
    font-size: 42px !important; font-weight: 700 !important; letter-spacing: 0.02em;
    color: #1e3a5f !important;
}

.page-ti-ru {
    width: 180px; height: 2px; background: #cbd5e1;
    margin: 0 auto 12px auto; border-radius: 999px;
}

.stagehead {
    display: inline-block; margin-top: 15px !important; margin-bottom: 8px !important;
    padding: 0.4rem 1rem 0.4rem 1.1rem; font-size: 32px; font-weight: 700 !important;
    letter-spacing: 0.14em; text-transform: uppercase;
    color: var(--accent) !important; background: #eff6ff;
    border: 4px solid #20262b; border-radius: 4px;
    font-family: 'Helvetica Neue', Helvetica, 'Segoe UI', Arial, sans-serif;
    -webkit-text-stroke-width: 0;
}

@media (max-width: 700px) {
    body {padding: 14px;}
    .boxed, .tocbox {padding: 16px;}
}
</style>
"

msg_note() {
    wordString=$1
    xtraString=$2
    echo "<h4 style='margin:0;'>${FG_YLW}${wordString}${NORMAL} <u style='color:white;font-style:italic;'>${xtraString}</u></h4>"
    echo " "
}

# Stages -----------------------
msg_stage() {
    echo -e " \n<a class='stagehead' id='stage-${1}'>STAGE-${1}</a>\n"
    echo -e " <h3 style='margin:5px 0 8px 0!important; font-style:italic; font-weight:800; font-size:24px; letter-spacing:0.01em; font-family:Aptos, Calibri, \"Segoe UI\", Arial, sans-serif; color:#ffffff; text-shadow:0 1px 1px rgba(0,0,0,0.18);'>${2}</h3>\n${BOX1}\n "
}

msg_code() {
    echo -e "<pre>\n<code class=\"language-bash\" id=\"copy\">${1}</code>\n</pre>"
}

msg_txt() {
    echo -e "<div class='text-block text-copy-pad'><button type='button' class='text-copy-btn'>Copy Text</button>"
    echo -e "<div class='text-content'>${1}</div></div>"
}

msg_step() {
    stepNumb=$1
    stepColor=$2
    stepString=$3
    bg_color="<span style=\"background-color:$stepColor;border-radius:5px!important;padding:5px 8px!important;line-height:1.5;\"> "
    echo -e "<h5 style='font-size:14;margin-top:10px!important;'>${bg_color}${stepNumb}. ${stepString}${NORMAL}</h5>\n<hr>"
}

white_bg() {
    echo "${WHT_BG}<u style='padding:2px;margin-top:5px;margin-bottom:5px;line-height:2;font-size:14pt;'>${1}</u></span></br>"
}

msg_bold() { echo "<u><b>${1}</b></u></br>"; }

msg_white() {
    echo "<p style=\"color:white!important;line-height:1.0;padding:3px;margin:0;\">${1}</p>"
}

print_toc() {
    liststart="<li><a href=\"#stage-${1}\">"
    echo "${liststart}${2}</a></li>"
}

msg_input() { echo "<tr><th>${1}</th><td>${2}</td></tr>";}
msg_meta() { echo "<span style='font-weight:bold'>${1}</span>: ${2}</br>"; }
msg_title() { echo "<h1 class='page-title'>${1}</h1><div class='page-ti-ru'></div>";}
msg_meta "Author" "Jonathan Serrano"; msg_meta "Script Version" "${VERS}"; msg_meta "Current Date" "${TD_DATE}"; echo "</br>"

msg_title "LG-PACT Commands"

echo "<div class='input-table-wrap'><table class='input-table'>"
echo "<tr><th class='input-table-title' colspan='2'>Script Inputs</th></tr>"
msg_input "RUN_ID" "${RUN_ID}"
msg_input "PACT_ID" "${PACT_ID}"
msg_input "Consensus Directory" "${CONSENSUS_DIR}"
msg_input "Kerberos ID" "${kerbero}"
echo -e "</table></div><br></br>" && echo ' <div class="tocbox"> '

# Table of Contents -----------------------
echo -e "<div class='toc-title'>Table of Contents</div>\n"
echo '<ol start="1" style="font-size:14;">'
print_toc 1 "Demultiplexing"
print_toc 2 "Execute In-House Pipeline and Philips Uploads"
print_toc 3 "Generate In-House QC and Copy the Output Files and QC to Zdrive"
print_toc 4 "Generate Consensus"
echo -e "</ol>\n$BOX2\n<br></br>"

# Stage 1 -----------------------
msg_stage 1 "Begin Demultiplexing"
msg_step 1 "#acacff" "Once sequencing is complete, start demux by logging into HPC"
msg_code "ssh -Y ${kerbero}@bigpurple.nyumc.org"
msg_step 2 "#acacff" "Go into the demux-nf2 directory and execute the deploy command"
msg_code "${BASH_HELPERS}/start_demux.sh ${RUN_ID}"
echo "$BOX2"

# Stage 2 -----------------------
if [[ "$IS_SOPHIA" == "true" ]]; then
	msg_stage 2 "Deploy and execute the NGS607 Pipeline"
	msg_step 1 "#ffffba" "After Demux finishes, check the QC by pasting the link below in a web browser to open in a SFTP client like CyberDuck:"
	msg_code "sftp://bigpurple.nyumc.org${DEMUXDIR}/output/${RUN_ID}.report.html"
	msg_step "2" "#ffffba" "Begin the SG pipeline by executing the script below:"
	msg_code "${BASH_HELPERS}/start_pipeline.sh ${RUN_ID}"
	echo "$BOX2"
else
	msg_stage 2 "Execute In-House Pipeline and Philips Uploads"
	msg_step 1 "#ffffba" "After Demux finishes, check the QC by pasting the link below in a web browser to open in a SFTP client like CyberDuck:"
	msg_code "sftp://bigpurple.nyumc.org${DEMUXDIR}/output/${RUN_ID}.report.html"
	msg_step "2" "#ffffba" "Return to the Demux directory, execute the pass and upload *make* commands to upload to Philips, then deploy the in-house pipeline with EITHER 2a or 2b"
	msg_code "ssh ${kerbero}@bigpurple.nyumc.org"
	msg_step "2a" "#acacff" "If all samples in the run are production samples then execute the following command 2a"
	msg_code "${BASH_HELPERS}/start_pipeline.sh ${PACT_ID} ${RUN_ID}"
	msg_step "2b" "#acacff" "If the run has any filler or previous validation cases that need to be changed in the sample sheet, execute the steps below separately after editing the run"
	msg_code "${BASH_HELPERS}/make_passed_uploads.sh ${PACT_ID} ${RUN_ID}"
	msg_code "${BASH_HELPERS}/deploy_607_submit.sh ${PACT_ID} ${RUN_ID}"
	msg_step 3 "#ffffba" "ssh pgm@pgmapllcdcpvm01.nyumc.org to the isg-uploads folder to ensure the files are accessible and check the pgm log to verify uploading to Philips every 30 min"
	msg_code "cat pgm/log/uploads.log"
	echo "$BOX2"
fi

# Stage 3 -----------------------
msg_stage 3 "Copy QC and output data to the Molecular Z-drive"
msg_step 1 "#baffc9" "Copy the calibrated BAMs and run the post-scripts, then ssh to the data mover node by executing below:"
msg_code "${BASH_HELPERS}/bam_copier.sh ${RUN_ID} ${PACT_ID}"
msg_step 2 "#baffc9" "From the data mover node, mount /mnt/${kerbero}/molecular, and execute zdrive_copier.sh to copy the data from HPC to the Molecular share drive."
msg_code "/mnt/${kerbero}/molecular/Molecular/Validation/Scripts/zdrive_copier.sh ${RUN_ID} ${PACT_ID}"
msg_step 3 "#baffc9" "Once the files are copied, email the PACT team to notify them by replying to the PACTers email thread with the following message template:"
if [[ "$IS_SOPHIA" == "true" ]]; then
	msg_txt "Hi all,\nThe in-house pipeline completed for ${PACT_ID}. The data for this week's PACT run is copied here:
	smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular/Molecular/NGS607/${YEAR_DIR}/${RUN_ID}/\n
	The QC and output is copied here:
	smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular/MOLECULAR LAB ONLY/NYU PACT Patient Data/Results/Bioinformatics/${YEAR_DIR}/${PACT_ID}/${PACT_ID}.html
	QC for upload to Sophia can be reviewed by logging into: https://pactqc.nyumc.org/
	"
else
	msg_txt "Hi all,\nThe in-house pipeline completed for ${PACT_ID}. The data for this week's PACT run is copied here:
	smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular/Molecular/NGS607/${YEAR_DIR}/${RUN_ID}/\n
	The QC and output is copied here:
	smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular/MOLECULAR LAB ONLY/NYU PACT Patient Data/Results/Bioinformatics/${YEAR_DIR}/${PACT_ID}/${PACT_ID}.html"
fi
echo "$BOX2"

# Stage 4 -----------------------
if [[ "$IS_SOPHIA" == "true" ]]; then
	msg_stage 4 "Generate Methylation CNV Consensus Report"
	msg_step 1 "#bae1ff" "Once the SG results are ready, open your LOCAL terminal and execute the following script that was downloaded to your home folder by parseSophia.sh"
	msg_code "${HOME}/make_consensus.sh ${RUN_ID} ${PACT_ID}"
	msg_step 2 "#bae1ff" "Send an email to notify the html report has been generated and copied using the following as a template:"
	msg_txt "Hi all,\nThe methylation CNV consensus is copied here:\nsmb://shares-cifs.nyumc.org/apps/acc_pathology${outputDir}${YEAR_DIR}/${PACT_ID}/${PACT_ID}_consensus.html"
	echo "$BOX2"
else
	msg_stage 4 "Methylation and Philips CNV Consensus Report"
	msg_step 1 "#bae1ff" "In your LOCAL terminal execute the make_consensus.sh script Monday after the Philips data dumps into /molecular/Molecular/Philips_SFTP"
	msg_code "${HOME}/make_consensus.sh ${RUN_ID} ${PACT_ID}"
	msg_step 2 "#bae1ff" "Send an email to notify the file is ready"
	msg_txt "Hi all,\nThe methylation CNV consensus is copied here:\nsmb://shares-cifs.nyumc.org/apps/acc_pathology${outputDir}${YEAR_DIR}/${PACT_ID}/${PACT_ID}_consensus.html"
	echo "$BOX2"
fi

echo "
<script>
document.addEventListener('DOMContentLoaded', () => {
    const copyToClipboard = async (text) => {
        if (navigator.clipboard && window.isSecureContext) {
            await navigator.clipboard.writeText(text);
            return true;
        }

        const ta = document.createElement('textarea');
        ta.value = text;
        ta.setAttribute('readonly', '');
        ta.style.position = 'absolute';
        ta.style.left = '-9999px';
        document.body.appendChild(ta);
        ta.select();

        let ok = false;
        try {
            ok = document.execCommand('copy');
        } catch (err) {
            ok = false;
        }

        document.body.removeChild(ta);
        return ok;
    };

    document.querySelectorAll('.text-copy-btn').forEach((button) => {
        button.addEventListener('click', async () => {
            const block = button.closest('.text-block');
            if (!block) {return;}

            const content = block.querySelector('.text-content');
            if (!content) {return;}

            const ok = await copyToClipboard(content.textContent);

            if (ok) {
                button.innerText = 'Text Copied';
                button.classList.add('pressed');
            } else {
                button.innerText = 'Copy Failed';
            }
        });
    });

    document.querySelectorAll('pre').forEach((block) => {
        const button = document.createElement('button');
        button.type = 'button';
        button.innerText = 'Copy Code';

        button.addEventListener('click', async () => {
            const code = block.querySelector('code');
            if (!code) {return;}

            const ok = await copyToClipboard(code.textContent);

            if (ok) {
                button.innerText = 'Code Copied';
                button.classList.add('pressed');
            } else {
                button.innerText = 'Copy Failed';
            }
        });

        block.prepend(button);
    });
});
</script>
"
