#!/bin/bash
FUSIONRUNID=${1-NULL} # if arg $1 is empty assign NULL as default else i.e. 123456_NB501073_0212_AHT3V7BGXK
FSID=${2-NULL}        # if arg $2 is empty assign NULL as default else i.e. PACT-21-28
kerbero=${3-$USER}    # if arg $3 is empty assign $USER as default else i.e. whoami kerberosid

python_version=$(python3 --version)
python_path=$(which python3)

FG_YLW='<span style="color:#cf6a00;">' # makes text darkorange color
FG_RED='<span style="color:#cc0000">'  # makes text red color
FG_GRN='<span style="color:#6aa84f">'  # makes text color green
FG_BLU='<span style="color:blue">'     # makes text color blue
FG_CYA='<span style="color:cyan;">'    # makes text color #baffc9
WHT_BG="<span style='background-color:white;margin-left:30px;padding:3px;margin-top:5px;margin-bottom:5px;line-height:1.2!important;'>"
normal="</span>" # resets default text
BOX1=' <div class="boxed"> '
BOX2=' </div> '
currYear=$(date +'%Y')
#currYear=$(date +"%Y") #date2022
productionDir="/gpfs/data/molecpathlab/production"
outputDir="/molecular/MOLECULAR LAB ONLY/NYU PACT Patient Data/Results/Bioinformatics/"
reprtDir="/gpfs/data/molecpathlab/bin/QC_reprot/"
pactGithub="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts"

SHEETDIR="${productionDir}/samplesheets/archer/${FUSIONRUNID}"
DEMUXDIR="${productionDir}/Demultiplexing/${FUSIONRUNID}"
molecDir="${productionDir}/NGS607/"

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

button:active, button:focus {
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
  font-size:40px !important;
  font-weight: 600!important;
  background-image: linear-gradient(to bottom, #2c3e50, orange, purple )!important;
  color: transparent!important;
  display: block;
  background-clip: text!important;
  -webkit-background-clip: text!important;
  font-family: 'Allerta Stencil';
  margin-bottom: 0px !important;
  -webkit-text-stroke-width: 1px;
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
    #echo "<code>${codeString}</code></br>"
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
lastTwo=${FSID: -2}
runMid=${FSID:0:2}
FSID="${FG_RED}${FSID}${normal}"
kerbero="${FG_CYA}${kerbero}${normal}"
FUSIONRUNID="${FG_YLW}${FUSIONRUNID}${normal}"
currYear="${FG_BLU}20$runMid${normal}"
rsyncDir="/gpfs/home/${kerbero}/molecpathlab/production/NGS607/${FUSIONRUNID}/output"
zdrive="/mnt/${kerbero}/molecular/Molecular"

echo "<span style='font-weight: bold'>Author</span>: Jonathan Serrano</br>"
echo "<span style='font-weight: bold'>Current Date</span>: $(date)</br>"
echo "</br>"
echo "<h2 style='padding-top: 10px !important; -webkit-text-stroke-width: 1px; -webkit-text-stroke-color: black;'>${FG_GRN}FUSION Seq Commands${normal}</h2>"
echo "<h2 style='padding-top: 10px !important; -webkit-text-stroke-width: 1px; -webkit-text-stroke-color: black;'>${FG_BLU}Your Input Args${normal}</h2>"
msg_step 1 "white" "Sequencer Run ID FUSIONRUNID: ${FUSIONRUNID}</br>"
msg_step 2 "white" "FUSION Run Name FSID: ${FSID}</br>"
msg_step 3 "white" "Kerberos ID: ${kerbero}</br>"
msg_step 4 "white" "Python3 path: ${python_version}</br>"
echo "</br>"
echo ' <div class="tocbox"> '
#echo " "
# Table of Contents -----------------------
echo "<h2 style='margin-top: 0px;font-size:20;'> Table of Contents </h2>"
echo " "
echo '<ol start="0" style="font-size:14;">'
print_toc 0 "Create the SampleSheet.csv & Copy to BigPurple"
print_toc 1 "Demultiplexing Steps"
print_toc 2 "Transfer demux output data to Z-Drive"
print_toc 3 "Notify the Lab and Generate HeatMap and QC"
echo "</ol>"
echo "$BOX2"
# Stage 0 -----------------------
msg_stage 0 "Create the SampleSheet.csv & Copy to BigPurple"
msg_step 1 "#ffb3ba" "First mount the Molecular drive and execute the script. For example"
msg_code "diskutil mountdisk /Volumes/CBioinformatics/"
msg_code "/Volumes/CBioinformatics/FUSION/parseFusion.sh /Users/${kerbero}/Downloads/${FSID}.xlsx"
msg_step 2 "#ffb3ba" "Review the generated SampleSheet.csv ouput below and make sure it contains no errors"
msg_code "${HOME}/Desktop/${FUSIONRUNID}-SampleSheet.csv"
msg_step 3 "#ffb3ba" "Check the script copied with group read/write permissions to the folder"
msg_code "${productionDir}/samplesheets/archer/${FUSIONRUNID}/${FUSIONRUNID}-SampleSheet.csv"
echo "$BOX2"

# Stage 1 -----------------------
msg_stage 1 "Demultiplexing Steps"
msg_step 1 "#ffdfba" "Once sequencing is finished, start demultiplexing by logging into BigPurple"
msg_code "ssh -Y ${kerbero}@bigpurple.nyumc.org"
msg_step 2 "#ffdfba" "Go into demux-nf2"
msg_code "cd /gpfs/data/molecpathlab/pipelines/demux-nf2/"
msg_step 3 "#ffdfba" "Execute the deploy command"
msg_code "make deploy RUNID=${FUSIONRUNID} SAMPLESHEET=${productionDir}/samplesheets/archer/${FUSIONRUNID}/${FUSIONRUNID}-SampleSheet.csv SEQTYPE=Archer"
msg_step 4 "#ffdfba" "Go into the Demultiplexing run folder:"
msg_code "cd $DEMUXDIR"
msg_step 5 "#ffdfba" "Update/submit and cat logs"
msg_code "make update && make submit"
echo "$BOX2"

# Stage 2 -----------------------
msg_stage 2 "Transfer demux output data to Z-Drive"
msg_step 1 "#ffffba" "SSH into BigPurple and change permissions for the output run directory"
msg_code "ssh ${kerbero}@bigpurple.nyumc.org"
msg_code "chmod -R g+rwx ${DEMUXDIR}/output"
msg_step 2 "#ffffba" "ssh to data mover node and transfer file to Z drive"
msg_code "ssh ${kerbero}@dmn-0002"
msg_code "mount /mnt/${kerbero}/molecular"
msg_note "NOTE:" "Some times the 0001 will be down and the 0002 will take over.  We usually just use 0002"
msg_step 3 "ffffba" "Go to the output Directory and Rsync the data"
msg_code "mkdir /mnt/${kerbero}/molecular/Molecular/Demultiplexing/${FUSIONRUNID}"
msg_code "rsync -vrthP /gpfs/home/${kerbero}/molecpathlab/production/Demultiplexing/${FUSIONRUNID}/output /mnt/$kerbero${normal}/molecular/Molecular/Demultiplexing/$FUSIONRUNID"
msg_step 4 "ffffba" "Email the FUSION team once the QC files are generated and copied to the Z-drive:"
msg_code "Demultiplexing completed for ${FSID}.
The data for this week’s Fusion run is copied here:
smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular/Molecular/Demultiplexing/${FUSIONRUNID}/"
echo "$BOX2"

# Stage 3 -----------------------
msg_stage 3 "Generate HeatMap and QC"
msg_step 1 "#baffc9" "After Xiaojun says the data is downloaded to the Z-drive, download the github repo for the QC"
msg_code "gh repo clone NYU-Molecular-Pathology/FusionSeq_QC"
msg_step 2 "#baffc9" "Execute the Python script in that directory using your sample inputs"
msg_code "cd FusionSeq_QC && ${python_path} control_QC_both.py -d '/Volumes/molecular/MOLECULAR LAB ONLY/NYU FUSION SEQer/${currYear} reports' -r ${FSID} -o $HOME/FusionSeq_QC"
msg_step 3 "#baffc9" "After data is output, generate the QC and email"
msg_code "Heatmaps are generated and are copied here:
smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular/MOLECULAR LAB ONLY/NYU FUSION SEQer/${currYear} reports/${FSID}/"
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
