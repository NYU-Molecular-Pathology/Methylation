#!/bin/bash
## ---------------------------
## Script name: printPactCommands.sh
## Purpose: Print out all the copiable commands used for executing the PACT pipeline
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2023
## ---------------------------

runID=${1-NULL}   # if arg $1 is empty assign NULL as default else i.e. 123456_NB501073_0212_AHT3V7BGXK
pactRun=${2-NULL} # if arg $2 is empty assign NULL as default else i.e. PACT-21-28
DEFAULTVALUE="/Volumes/CBioinformatics/jonathan/pact/consensus/"
consensusDir="${3:-$DEFAULTVALUE}"
kerbero=${4-$USER} # if arg $3 is empty assign $USER as default else i.e. whoami kerberosid

FG_YLW='<span style="color:#cf6a00;">' # makes text darkorange color
FG_RED='<span style="color:#cc0000">'  # makes text red color
FG_GRN='<span style="color:#6aa84f">'  # makes text color green
FG_BLU='<span style="color:blue">'     # makes text color blue
FG_CYA='<span style="color:cyan;">'    # makes text color #baffc9
WHT_BG="<span style='background-color:white;margin-left:30px;padding:3px;margin-top:5px;margin-bottom:5px;line-height:1.2!important;'>"
normal="</span>" # resets default text
BOX1=' <div class="boxed"> '
BOX2=' </div> '
year_part=${pactRun:5:2}
currYear="20${year_part}"
#currYear=$(date +"%Y") #date2022
#evernoteLink='https://www.evernote.com/shard/s331/sh/5416e425-83c7-5aeb-0683-6667fb3d6f8e/07ef3e8f603ecdff3afe5da18f0204f2'
consensusDir='/Volumes/CBioinformatics/jonathan/pact/consensus/'
productionDir="/gpfs/data/molecpathlab/production"
outputDir="/molecular/MOLECULAR LAB ONLY/NYU PACT Patient Data/Results/Bioinformatics/"
#reprtDir="/gpfs/data/molecpathlab/bin/QC_reprot/"
pactGithub="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts"
SHEETDIR="${productionDir}/samplesheets/LG-PACT/${runID}"
DEMUXDIR="${productionDir}/Demultiplexing/${runID}"
#molecDir="${productionDir}/NGS607/"

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
  background-image: linear-gradient(
        to bottom,
        #CA4246, 
        #E16541, 
        #F18F43, 
        #8B9862)!important;
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
#lastTwo=${pactRun: -2}
runMid=${pactRun:5:2}
pactRun="${FG_RED}${pactRun}${normal}"
kerbero="${FG_CYA}${kerbero}${normal}"
runID="${FG_YLW}${runID}${normal}"
currYear="${FG_BLU}20$runMid${normal}"
gpfsHome="/gpfs/home/${kerbero}/"
rsyncDir="${gpfsHome}molecpathlab/production/NGS607/${runID}/output"
zdrive="/mnt/${kerbero}/molecular/Molecular"

echo "<span style='font-weight: bold'>Author</span>: Jonathan Serrano</br>"
echo "<span style='font-weight: bold'>Current Date</span>: $(date)</br>"
echo "</br>"
echo "<h2 style='padding-top: 5px !important; -webkit-text-stroke-width: 1px; -webkit-text-stroke-color: black;'>${FG_GRN}LG-PACT Commands${normal}</h2>"
msg_step 1 "white" "PACT RUNID: ${runID}</br>"
msg_step 2 "white" "PACT Run Name: ${pactRun}</br>"
msg_step 3 "white" "Consensus Directory: ${consensusDir}</br>"
msg_step 4 "white" "Kerberos ID: ${kerbero}</br>"
echo "</br>"
echo ' <div class="tocbox"> '
#echo " "
# Table of Contents -----------------------
echo "<h2 style='margin-top: 0px;font-size:20;'> Table of Contents </h2>"
echo " "
echo '<ol start="0" style="font-size:14;">'
print_toc 1 "Demultiplexing"
print_toc 2 "Execute In-House Pipeline and Philips Uploads"
print_toc 3 "Generate In-House QC and Copy the Output Files and QC to Zdrive"
print_toc 4 "Generate Consensus"
echo "</ol>"
echo "$BOX2"
# Stage 0 -----------------------
# msg_stage 0 "Create SampleSheet.csv & Copy to BigPurple"
# msg_step 1 "#ffb3ba" "First mount the Molecular drive and execute the shell script"
# msg_code "diskutil mountdisk /Volumes/CBioinformatics/"
# msg_code "/Volumes/CBioinformatics/PACT/parsepact.sh ${pactRun} ${runID}"
# msg_step 2 "#ffb3ba" "Review the generated SampleSheet.csv ouput below and notify the lab if it contains any errors"
# msg_code "${HOME}/Desktop/${runID}-SampleSheet.csv"
# msg_step 3 "#ffb3ba" "Verify the script copied the samplesheet to BigPurple and added group read/write permissions to the directory:"
# msg_code "${productionDir}/samplesheets/LG-PACT/${runID}/${runID}-SampleSheet.csv"
# echo "$BOX2"

# Stage 1 -----------------------
msg_stage 1 "Demux Steps"
msg_step 1 "#ffdfba" "Once sequencing is finished, start demultiplexing by logging into BigPurple"
msg_code "ssh -Y ${kerbero}@bigpurple.nyumc.org"
msg_step 2 "#ffdfba" "Go into the demux-nf2 directory and execute the deploy command"
msg_code "/gpfs/data/molecpathlab/scripts/bash_helpers/start_demux.sh ${pactRun} ${runID}"
#msg_code "cd /gpfs/data/molecpathlab/pipelines/demux-nf2/ && make deploy RUNID=${runID} SAMPLESHEET=$SHEETDIR/${runID}-SampleSheet.csv SEQTYPE=NGS607"
#msg_step 3 "#ffdfba" "Execute the deploy command for your run"
#msg_code "make deploy RUNID=${runID} SAMPLESHEET=$SHEETDIR/${runID}-SampleSheet.csv SEQTYPE=NGS607"
#msg_step 3 "#ffdfba" "Go into the newly deployed Demultiplexing run directory then make update and submit"
#msg_code "cd $DEMUXDIR && make update && make submit"
#msg_code "make update && make submit"
echo "$BOX2"

# Stage 2 -----------------------
msg_stage 2 "Execute In-House Pipeline and Philips Uploads"
msg_step 1 "#ffffba" "After Demux finishes, check the QC by pasting the link below in a web browser to open in a SFTP client like CyberDuck:"
msg_code "sftp://bigpurple.nyumc.org${DEMUXDIR}/output/${runID}.report.html"
#msg_code "chmod -R g+rwx ${DEMUXDIR}/output"
msg_step "2" "#ffffba" "Return to the Demux directory, execute the pass and upload *make* commands to upload to Philips, then deploy the in-house pipeline with EITHER 2a or 2b"
msg_code "ssh ${kerbero}@bigpurple.nyumc.org"
msg_step "2a" "#acacff" "If all samples in the run are production samples then execute the following command 2a"
msg_code "/gpfs/data/molecpathlab/scripts/bash_helpers/start_pipeline.sh ${pactRun} ${runID}"
msg_step "2b" "#acacff" "If the run has any filler or previous validation cases that need to be changed in the sample sheet, execute the steps below separately after editing the run"
msg_code "/gpfs/data/molecpathlab/scripts/bash_helpers/make_passed_uploads.sh ${pactRun} ${runID}"
msg_code "/gpfs/data/molecpathlab/scripts/bash_helpers/deploy_607_submit.sh ${pactRun} ${runID}"
#msg_code "make passed && make uploads && make deploy-NGS607 && cd ${productionDir}/NGS607/${runID} && chmod -R g+rwx ${productionDir}/isg-uploads/${runID}"
#msg_code "cd ${productionDir}/NGS607/${runID}"
#msg_step 3 "#ffffba" "Update and then submit the slurm job and check your logs using squeue -u ${kerbero} && lt logs/*"
#msg_code "make update && make submit"
#msg_code "squeue -u ${kerbero} && lt logs/*"
msg_step 3 "#ffffba" "ssh pgm@pgmapllcdcpvm01.nyumc.org to the isg-uploads folder to ensure the files are accessible and check the pgm log to verify uploading to Philips every 30 min"
#msg_code "chmod -R g+rwx ${productionDir}/isg-uploads/${runID}"
#msg_code "ssh pgm@pgmlcdcpvm01.nyumc.org"
#msg_note "NOTE:" "If you forget the pgm password, we have it saved in Evernote"
msg_code "cat pgm/log/uploads.log"
msg_note "Reminder:" "If you forget the pgm password, check Evernote. Ensure Philips ISPM uploads are not pending or missing. Nextflow will email you any errors and when the pipeline completes"
echo "$BOX2"

# Stage 3 -----------------------
# msg_stage 3 "Generate QC"
# msg_step 1 "#d9d2e9" "After NextFlow emails you the pipeline completed successfully, return to BigPurple and generate the QC"
# msg_code "ssh ${kerbero}@bigpurple.nyumc.org"
# msg_step 2 "#d9d2e9" "Go to the run production directory, edit permissions, and then execute the Python Automations"
# msg_code "module load python/cpu/3.8.11 && cd ${molecDir}${runID} && mkdir -p ${productionDir}/NGS607/${runID}/output/clinical/ && chmod -R g+rwx ${productionDir}/NGS607/${runID}/output/"
# msg_code "python ${reprtDir}snp_overlap.py -o ${molecDir}${runID}/output/ && python ${reprtDir}generate_html_report.py -o ${molecDir}${runID} -p ${pactRun} -r ${runID} && python3 ${reprtDir}variants_qc.py -rid ${runID} -rdir ${molecDir}${runID}/output/ -pactid ${pactRun}"
# msg_code "python /gpfs/data/molecpathlab/bin/QC_reprot/detect_hotspots.py -rid ${runID} -pactid ${pactRun}"
# msg_step 3 "#d9d2e9" "Change permissions for group access to the new QC files ouput"
# msg_code "chmod -R g+rwx ${productionDir}/NGS607/${runID}/output/"
# msg_note "If this step breaks, try re-installing miniconda and the requirements.txt with the outlined commands in step 4"
# msg_step 4 "#a980ff" "If you have MiniConda3 skip to step 5, otherwise install it and type \"no\" when prompted if you wish the installer to initialize Miniconda3 by running conda init"
# msg_code "cd /gpfs/home/${kerbero} && wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh && bash Miniconda3-py38_4.10.3-Linux-x86_64.sh"
# msg_step 4a "#a980ff" "To complete installation, cd into your bin directory, then CLOSE the terminal window (exit BigPurple) to restart it"
# msg_code "cd ${gpfsHome}miniconda3/bin/ && ./conda init"
# msg_step 4b "#a980ff" "After closing the console, ssh back into BigPurple and use conda to install pip"
# msg_code "ssh ${kerbero}@bigpurple.nyumc.org"
# msg_code "conda install pip"
# msg_code "pip install -r ${gpfsHome}molecpathlab/development/NGS_QC_xf/requirements.txt"
# msg_step 5 "#d9d2e9" "Go to the QC python directory and generate the QC file"
# msg_code "conda activate && python ${gpfsHome}molecpathlab/development/NGS_QC_xf/xf_pactqc.py -rdir ${molecDir}${runID}/output -pactid ${pactRun} && conda deactivate"
# msg_code "chmod -R g+rwx ${productionDir}/NGS607/${runID}/output/"
# echo "$BOX2"

BAMSDIR="/gpfs/data/molecpathlab/production/NGS607/${runID}/output/alignments/recalibrated"

# Stage 3 -----------------------
msg_stage 3 "Copy the QC files and Output data to the Molecular Z-drive"
msg_step 1 "#baffc9" "Create the output BAM directory and copy calibrated and then ssh to the data mover node by executing below:"
msg_code "/gpfs/data/molecpathlab/scripts/bash_helpers/bam_copier.sh ${runID} ${pactRun}"
#msg_code "mkdir -p \"/gpfs/data/clinpathlab/external/${pactRun}\" && rsync -vrthP ${BAMSDIR}/*.dd.ra.rc.bam \"/gpfs/data/clinpathlab/external/${pactRun}/\" && rsync -vrthP ${BAMSDIR}/*.dd.ra.rc.bam.bai \"/gpfs/data/clinpathlab/external/${pactRun}/\" && chmod -R ag+rwx \"/gpfs/data/clinpathlab/external/${pactRun}\""
msg_step 2 "#baffc9" "From the data mover node, mount /mnt/${kerbero}/molecular, and execute zdrive_copier.sh below"
#msg_code "ssh ${kerbero}@dmn-0002"
#msg_code "mount /mnt/${kerbero}/molecular"
#msg_note "NOTE:" "Occasionally, the 0001 will be down and the 0002 will take over.  We usually just use 0002"
#msg_step 2 "#baffc9" "Once mounted, create the output directories in /MOLECULAR/NGS607/"
msg_code "/mnt/${kerbero}/molecular/Molecular/Validation/Scripts/zdrive_copier.sh ${runID} ${pactRun}"
# msg_code "mkdir -p \"${zdrive}/NGS607/${currYear}/${runID}/output/alignments\" \"/mnt/${kerbero}${outputDir}${currYear}/${pactRun}\""
# msg_code "chmod -R g+rwx \"/mnt/${kerbero}${outputDir}${currYear}/${pactRun}\" && chmod -R g+rwx ${gpfsHome}molecpathlab/production/NGS607/${runID}/output"
# msg_step 3 "#baffc9" "Once created, rsync the files from BigPurple to /MOLECULAR/NGS607/"
# msg_white "<b>Copy Alignments and Annotations then theClinical Folder and .tsv files:</b>"
# msg_code "rsync -vrthP ${rsyncDir}/alignments/${FG_GRN}recalibrated${normal} ${zdrive}/NGS607/${currYear}/${runID}/output/${FG_GRN}alignments/${normal} && rsync -vrthP ${rsyncDir}/${FG_GRN}annotations${normal} ${zdrive}/NGS607/${currYear}/${runID}/output/${FG_GRN}annotations/${normal} && rsync -vrthP ${rsyncDir}/${FG_GRN}clinical${normal} ${zdrive}/NGS607/${currYear}/${runID}/output/ && rsync -vrthP ${rsyncDir}/${FG_GRN}*.tsv${normal} ${zdrive}/NGS607/${currYear}/${runID}/"
# msg_step 4 "#baffc9" "Next, go to /MOLECULAR LAB ONLY/NYU PACT Patient Data/ on the Z-drive and copy the remaining clinical output files, PDF facets, and QC tsv "
# msg_code "cd \"/mnt/${kerbero}${outputDir}${currYear}/${pactRun}/\" && rsync -vrthP ${rsyncDir}/${FG_GRN}clinical/${normal} ./"
# msg_code "rsync -vrthP ${productionDir}/NGS607/${runID}/output/cnv/FACETS/*.pdf ${productionDir}/NGS607/${runID}/${pactRun}-QC.tsv ${zdrive}/REDCap/cnv_facets/${pactRun}/"
msg_step 3 "#baffc9" "Email the PACT team once the QC files are copied to notify them the following"
msg_code "The in-house pipeline completed for ${pactRun}. The data for this week’s PACT run is copied here:
smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular/Molecular/NGS607/${currYear}/${runID}/
The QC and output is copied here:
smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular/MOLECULAR LAB ONLY/NYU PACT Patient Data/Results/Bioinformatics/${currYear}/${pactRun}/${pactRun}.html"
echo "$BOX2"

# Stage 4 -----------------------
msg_stage 4 "Methylation and Philips CNV Consensus Report"
msg_step 1 "#bae1ff" "In your LOCAL terminal execute the make_consensus.sh script Monday after the Philips data dumps into /molecular/Molecular/Philips_SFTP" #, create a new directory for your PACT consensus in CBioinformatics drive and curl the template RMD file from GitHub"
msg_code "${HOME}/make_consensus.sh ${runID} ${pactRun}"
# msg_code "mkdir -p \"${consensusDir}${pactRun}_consensus\" && cd \"${consensusDir}${pactRun}_consensus\""
# msg_code "curl -# -L ${pactGithub}/PACT_consensus.Rmd >${consensusDir}${pactRun}_consensus/${pactRun}_consensus.Rmd"
# msg_step 1 "#bae1ff" "Re-generate the latest MethylMatch data for concensus if there was a new run last Friday"
# msg_code "/Volumes/CBioinformatics/PACT/getMethylMatch.sh ${pactRun} ${runID}"
# msg_step 2 "#bae1ff" "From the concensus directory, copy the .cnv.plot.pdf facets and QC from the Z-drive, and MethylMatch.xlsx from the Desktop"
# msg_code "cd \"${consensusDir}${pactRun}_consensus\" && cp /Volumes/molecular/Molecular/REDCap/cnv_facets/${pactRun}/*.pdf ~/Desktop/${runID}-SampleSheet.csv ~/Desktop/${pactRun}_MethylMatch.xlsx /Volumes/molecular/Molecular/NGS607/${currYear}/${runID}/output/clinical/${pactRun}-Somatic_Variants.html /Volumes/molecular/Molecular/NGS607/${currYear}/${runID}/output/clinical/${pactRun}.html /Volumes/molecular/Molecular/REDCap/cnv_facets/${pactRun}/${pactRun}-QC.tsv /Volumes/molecular/Molecular/NGS607/${currYear}/${runID}/output/${pactRun}_Hotspots.tsv ./"
# # msg_step 3 "#bae1ff" "Use the Evernote Guide to save the Somatic Variants as a .csv file to the downloads folder: ${HOME}/Downloads/export_mytable_MM_DD_${currYear}.csv"
# # msg_code "open ${consensusDir}${pactRun}_consensus/${pactRun}-Somatic_Variants.html && open ${evernoteLink}"
# msg_step 3 "#bae1ff" "Download and run the R script to save the ${pactRun}_desc.csv file"
# msg_code "cd ${HOME} && curl -# -L ${pactGithub}/MakeIndelList.R >${HOME}/MakeIndelList.R && chmod +rwx ${HOME}/MakeIndelList.R"
# msg_code "RScript --verbose ${HOME}/MakeIndelList.R ${pactRun}"
# msg_step 4 "#bae1ff" "Knit the rMarkdown file in your consensus directory"
# #msg_code "cp ${HOME}/Desktop/${pactRun}_desc.csv ${consensusDir}${pactRun}_consensus/"
# msg_code "cd ${consensusDir}${pactRun}_consensus/ && Rscript --verbose -e \"rmarkdown::render('${consensusDir}${pactRun}_consensus/${pactRun}_consensus.Rmd', params=list(pactName='${pactRun}', userName='${kerbero}'))\""
# msg_step 2 "#bae1ff" "Once the CNV concensus html is created, copy it to \"/Volumes${outputDir}${currYear}/${pactRun}/\""
# msg_code "open ${consensusDir}${pactRun}_consensus/${pactRun}_consensus.html"
# msg_code "cp ${consensusDir}${pactRun}_consensus/${pactRun}_consensus.html \"/Volumes${outputDir}${currYear}/${pactRun}/\""
msg_step 2 "#bae1ff" "Send an email to notify the file is ready"
msg_code "Hi all,
The methylation CNV consensus is copied here:
smb://shares-cifs.nyumc.org/apps/acc_pathology${outputDir}${currYear}/${pactRun}/${pactRun}_consensus.html"
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
