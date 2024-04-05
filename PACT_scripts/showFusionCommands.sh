#!/bin/bash
## ---------------------------
## Script name: printFusionCommands.sh
## Purpose: Print out all the copiable commands used for executing the PACT pipeline
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------

# INPUT ARGS -----------------------------------------------------------------------------------------
RUN_ID=${1-NULL}     # if arg $1 is empty assign NULL as default i.e. 220819_NB551999_0999_AHVTN5XYZ3
SHEETPATH=${2-$NULL} # if arg $2 is empty will be NULL i.e. "/Users/$USER/Downloads/22-MGFS25.xlsx"
FUSION_ID=${3-$NULL} # if arg $3 is empty, value is set to basename of SHEETPATH i.e. "22-MGFS25"
KERBERO=${4-$USER}   # if arg $4 is empty, value is set to $USER as default i.e. kerberos ID

FG_YLW='<span style="color:#cf6a00;">' # makes text darkorange color
FG_GRN='<span style="color:#6aa84f">'  # makes text color green
FG_CYA='<span style="color:cyan;">'    # makes text color #baffc9
WHT_BG="<span style='background-color:white;margin-left:30px;padding:3px;margin-top:5px;margin-bottom:5px;line-height:1.2!important;'>"
normal="</span>" # resets default text
BOX1=' <div class="boxed"> '
BOX2=' </div> '

# If FUSION_ID is NULL use `basename` to get FUSION_ID by removing the extension -----------------------------
if [ "$FUSION_ID" == "NULL" ]; then
    FUSION_ID=$(basename -- "$SHEETPATH")
    FUSION_ID="${FUSION_ID%.*}"
fi

KERBERO="${FG_CYA}${KERBERO}${normal}"
RUN_ID="${FG_YLW}${RUN_ID}${normal}"
SHEETDIR="/gpfs/data/molecpathlab/production/samplesheets/archer"
XLSXPATH="${SHEETDIR}/worksheets/${FUSION_ID}.xlsx"

echo "
<link href='https://fonts.googleapis.com/css?family=Allerta Stencil' rel='stylesheet'>


<style>
h1, h2, h3 {
    line-height: 40% !important;
    padding-bottom: 0px !important;
    margin-top: 0px !important;
    width: auto;
}

html body {
    color: darkblue;
    width: auto;
}

html body {
    font-family: 'Helvetica Neue',Helvetica,'Segoe UI',Arial,freesans,sans-serif;
}


.pressed {
    background-color: #06a96e;
    background-image: linear-gradient(1deg, #00aa6c, #14C667 99%);
}


.boxed {
    background: rgb(90, 90, 90) !important;
    border: 3px solid black;
    padding: 10px;
    border-radius: 10px;
    width: auto;
    display: inline-block;
    white-space: nowrap;
    font-size:12 !important;
}

.tocbox {
    background: rgb(255, 255, 255) !important;
    border: 3px solid black;
    padding: 20px;
    border-radius: 10px;
    width: auto;
    display: inline-block;
    white-space: nowrap;
    font-size:16 !important;
}

h5 {
    padding:2px !important;
    margin-bottom:4px !important;
    margin-top:4px !important;
}

body {
    background-color: grey !important;
    width: 100%;
    margin-left: 10px;
    width: auto;
    display: inline-block;
    white-space: nowrap;
}

.markdown-preview {
    margin-right: 10px !important;
    width: 100% !important;
    height: 100% !important;
    box-sizing: content-box !important;
    margin-left: 20px !important;
    padding-top: 4px !important;
    margin-top: 4px !important;
}

ol {
    line-height: 1.2 !important;
    margin-bottom:0!important;
}

code {
    background-color: black !important;
    color: white !important;
    margin-left: 20px !important;
    padding: 3px !important;
    line-height: 1.75 !important;
    font-family: Menlo,Monaco,Consolas,'Courier New',monospace!important;
    font-size:13 !important;
}

*, *:before, *:after {
    box-sizing: border-box;
}

pre[class*=\"language-\"] {
    position:relative;
    overflow: auto;
    margin:5px 0;
    padding:1.75rem 0 1.75rem 1rem;
    border-radius:10px;
}

button {
    align-items: center;
    appearance: none;
    background-color: #3EB2FD;
    background-image: linear-gradient(1deg, #4F58FD, #149BF3 99%);
    background-size: calc(100% + 20px) calc(100% + 20px);
    border-radius: 100px;
    border-width: 0;
    box-shadow: none;
    box-sizing: border-box;
    color: #FFFFFF;
    cursor: pointer;
    display: inline-flex;
    font-family: CircularStd,sans-serif;
    font-size: 14;
    height: auto;
    justify-content: center;
    line-height: 1.0;
    padding: 3px 10px;
    position: relative;
    text-align: center;
    text-decoration: none;
    transition: background-color .2s,background-position .2s;
    user-select: none;
    -webkit-user-select: none;
    touch-action: manipulation;
    vertical-align: top;
    white-space: nowrap;
}

button:active,
button:focus {
    outline: none;
}

button:hover {
    background-position: -20px -20px;
}

button:focus:not(:active) {
    box-shadow: rgba(40, 170, 255, 0.25) 0 0 0 .125em;
}

main {
    display: grid;
    max-width: 600px;
    margin: 20px auto;
}


h1{
    font-size:1.3rem;
}

.stagehead {
    font-size: 40px;
    font-weight: 600!important;
    background-image: linear-gradient(to bottom, #CA4246, #E16541, #F18F43, #8B9862)!important;
    color: transparent!important;
    display: block;
    background-clip: text!important;
    -webkit-background-clip: text!important;
    font-family: 'Allerta Stencil';
    margin-bottom: 0px !important;
    margin-top: 15px !important;
    -webkit-text-stroke-width: 0.05px;
    -webkit-text-stroke-color: black;
}

</style>
"

msg_note() {
    wordString=$1
    xtraString=$2
    echo "<h4 style='margin:0;'>${FG_YLW}${wordString}${normal} <u style='color:white;font-style:italic;'>${xtraString}</u></h4>"
    echo " "
}

# Stages -----------------------
msg_stage() {
    stageString=$1
    stageTitle=$2
    echo " "
    echo "<a class='stagehead' id='stage-${stageString}'>STAGE-${stageString}</a>"
    echo " "
    echo "<h3 style='font-style:italic;color:#004c00'>${stageTitle}</h3>"
    echo "$BOX1"
    echo " "
}

msg_code() {
    codeString=$1
    echo "<pre>"
    echo "<code class=\"language-bash\" id=\"copy\">${codeString}</code>"
    echo "</pre>"
}

msg_step() {
    stepNumb=$1
    stepColor=$2
    stepString=$3
    bg_color="<span style=\"background-color:$stepColor;border-radius:5px!important;padding:3px!important;line-height:1.1;\"> "
    echo "<h5 style='font-size:13;margin-top:10px!important;'>${bg_color}${stepNumb}. ${stepString}${normal}</h5>"
    echo "<hr>"
}

white_bg() {
    msgString=$1
    echo "${WHT_BG}<u style='padding:2px;margin-top:5px;margin-bottom:5px;line-height:2;font-size:14pt;'>${msgString}</u></span></br>"
}

msg_bold() {
    msgString=$1
    echo "<u><b>${msgString}</b></u></br>"
}

msg_white() {
    msgString=$1
    echo "<p style=\"color:white!important;line-height:1.0;padding:3px;margin:0;\">${msgString}</p>"
}

print_toc() {
    stepNumb=$1
    stepString=$2
    liststart="<li><a href=\"#stage-${stepNumb}\">"
    echo "${liststart}${stepString}</a></li>"
}

echo "<span style='font-weight: bold'>Author</span>: Jonathan Serrano</br>"
echo "<span style='font-weight: bold'>Current Date</span>: $(date)</br>"
echo "<span style='font-weight: bold'>Created with</span>: <pre>$HOME/showFusionCommands.sh $RUN_ID ${SHEETPATH} ${FUSION_ID} >$HOME/${FUSION_ID}.html</pre></br>"

echo "</br>"
echo "<h2 style='padding-top: 5px !important; -webkit-text-stroke-width: 1px; -webkit-text-stroke-color: black;'>${FG_GRN}LG-PACT Commands${normal}</h2>"
msg_step 1 "white" "FUSION RUN: ${RUN_ID}</br>"
msg_step 2 "white" "FUSION ID: ${FUSION_ID}</br>"
msg_step 3 "white" "KERBEROs ID: ${KERBERO}</br>"
echo "</br>"
echo ' <div class="tocbox"> '

# Table of Contents -----------------------
echo "<h2 style='margin-top: 0px;font-size:20;'> Table of Contents </h2>"
echo " "
echo '<ol start="0" style="font-size:14;">'
print_toc 1 "Demultiplexing"
print_toc 2 "Copy Output"
echo "</ol>"
echo "$BOX2"

# Stage 1 -----------------------
msg_stage 1 "Demux Automation Steps"
msg_step 1 "#ffdfba" "Before 8pm on Friday, the Fusion worksheet has been copied to BigPurple with the script below"
msg_code "/Volumes/CBioinformatics/FUSION/copyFusionSheet.sh ${SHEETPATH}"
msg_step 2 "#ffdfba" "If automation fails, the commands to create the SampleSheet and manually deploy are below"
msg_code "ssh ${KERBERO}@bigpurple.nyumc.org"
msg_code "module load python/cpu/3.6.5"
msg_code "mkdir -p ${SHEETDIR}/${RUN_ID} && chmod -R g+rwx ${SHEETDIR}/${RUN_ID}"
msg_code "cd ${SHEETDIR}/${RUN_ID}"
msg_code "python /gpfs/data/molecpathlab/bin/ArcherDX/gen_sample_sheet_fixed.py -t ${XLSXPATH} -o ."
msg_code "cd /gpfs/data/molecpathlab/pipelines/demux-nf2/"
msg_code "make deploy RUNID=${RUN_ID} SAMPLESHEET=${SHEETDIR}/${RUN_ID}/${RUN_ID}-SampleSheet.csv SEQTYPE=Archer"
msg_code "cd /gpfs/data/molecpathlab/production/Demultiplexing/${RUN_ID}"
msg_code "make update && make submit"
msg_code "cat ./logs/slurm*.out"
echo "$BOX2"

# Stage 2 -----------------------
msg_stage 2 "Copying Output Files to Z-drive"
msg_step 1 "#ffbfda" "After Demux finishes, ssh to data mover node"
msg_code "ssh ${KERBERO}@bigpurple.nyumc.org"
msg_code "ssh ${KERBERO}@dmn-0002"
msg_code "mount /mnt/${KERBERO}/molecular"
msg_step 2 "#ffbfda" "This script will copy the files to the Z-drive and email if rsync successful:"
msg_code "/mnt/${KERBERO}/molecular/Molecular/Validation/Scripts/fusion_run_copy.sh ${RUN_ID} ${FUSION_ID}"
msg_step 3 "#ffbfda" "Exit the data mover node and start V7 automation with the command"
msg_code "/gpfs/data/molecpathlab/production/Demultiplexing/240401_NB501073_0321_AH57KHAFX7/bin/UploadArcherFastqs.py -j ${FUSION_ID} -d /gpfs/data/molecpathlab/production/Demultiplexing/${RUN_ID}/output/reads/ArcherDx_Run"
echo "$BOX2"

echo "
<script>
    function copyEvent(id) {
        var str = document.getElementById(id);
        window.getSelection().selectAllChildren(str);
        document.execCommand(\"Copy\")
    }

    const copyButtonLabel = \"Copy Code\";
    let blocks = document.querySelectorAll(\"pre\");

    blocks.forEach((block) => {
        if (navigator.clipboard) {
            let button = document.createElement(\"button\");
            button.innerText = copyButtonLabel;
            button.addEventListener(\"click\", copyCode);
            button.onclick = copyEvent('copy');
            block.prepend(button);
        }
    });

    async function copyCode(event) {
        const button = event.srcElement;
        const pre = button.parentElement;
        let code = pre.querySelector(\"code\");
        let text = code.innerText;
        await navigator.clipboard.writeText(text);
        button.innerText = \"Code Copied\";
        button.className = \"pressed\";
    }

    function clearSelection(){
    if (window.getSelection) {window.getSelection().removeAllRanges();}
    else if (document.selection) {document.selection.empty();}
    }

    window.onload = clearSelection;

</script>
"
