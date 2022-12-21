#!/bin/bash
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
currYear=$(date +'%Y')
#currYear=$(date +"%Y") #date2022
evernoteLink='https://www.evernote.com/shard/s331/sh/5416e425-83c7-5aeb-0683-6667fb3d6f8e/07ef3e8f603ecdff3afe5da18f0204f2'
consensusDir='/Volumes/CBioinformatics/jonathan/pact/consensus/'
productionDir="/gpfs/data/molecpathlab/production"
outputDir="/molecular/MOLECULAR LAB ONLY/NYU PACT Patient Data/Results/Bioinformatics/"
reprtDir="/gpfs/data/molecpathlab/bin/QC_reprot/"
pactGithub="https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts"

SHEETDIR="${productionDir}/samplesheets/LG-PACT/${runID}"
DEMUXDIR="${productionDir}/Demultiplexing/${runID}"
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
  font-size:40px !important;
  font-weight: 600!important;
  background-image: linear-gradient(to bottom, #2c3e50, blue, purple )!important;
  color: transparent!important;
    display: block;
  background-clip: text!important;
  -webkit-background-clip: text!important;
  font-family: 'Allerta Stencil';
  margin-bottom: 0px !important;
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
lastTwo=${pactRun: -2}
runMid=${pactRun:5:2}
pactRun="${FG_RED}${pactRun}${normal}"
kerbero="${FG_CYA}${kerbero}${normal}"
runID="${FG_YLW}${runID}${normal}"
currYear="${FG_BLU}20$runMid${normal}"
rsyncDir="/gpfs/home/${kerbero}/molecpathlab/production/NGS607/${runID}/output"
zdrive="/mnt/${kerbero}/molecular/Molecular"

echo "<h2 style='padding-top: 10px !important;'>${FG_BLU}Your Input${normal}</h2>"
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
print_toc 0 "Setup SampleSheet.csv & Copy to BigPurple"
print_toc 1 "Demultiplexing"
print_toc 2 "Execute In-House Pipeline and Philips Uploads"
print_toc 3 "Generate In-House QC"
print_toc 4 "Copy the Output Files and QC to Zdrive"
print_toc 5 "Generate Consensus"
echo "</ol>"
echo "$BOX2"
# Stage 0 -----------------------
msg_stage 0 "Create SampleSheet.csv & Copy to BigPurple"
msg_step 1 "#ffb3ba" "First mount the Molecular drive and execute the script"
msg_code "diskutil mountdisk /Volumes/CBioinformatics/"
msg_code "/Volumes/CBioinformatics/PACT/parsepact.sh ${pactRun} ${kerbero}"
msg_step 2 "#ffb3ba" "Review the generated SampleSheet.csv ouput below and make sure it contains no errors"
msg_code "~/Desktop/${runID}-SampleSheet.csv"
msg_step 3 "#ffb3ba" "Check the script copied with group read/write permissions to the folder"
msg_code "${productionDir}/samplesheets/LG-PACT/${runID}/${runID}-SampleSheet.csv"
echo "$BOX2"

# Stage 1 -----------------------
msg_stage 1 "Demux Steps"
msg_step 1 "#ffdfba" "Once sequencing is finished, start demultiplexing by logging into BigPurple"
msg_code "ssh -Y ${kerbero}@bigpurple.nyumc.org"
msg_step 2 "#ffdfba" "Go into demux-nf2"
msg_code "cd ~/molecpathlab/pipelines/demux-nf2/"
msg_step 3 "#ffdfba" "Execute the deploy command"
msg_code "make deploy RUNID=${runID} SAMPLESHEET=$SHEETDIR/${runID}-SampleSheet.csv SEQTYPE=NGS607"
msg_step 4 "#ffdfba" "Go into the Demultiplexing run folder"
msg_code "cd $DEMUXDIR"
msg_step 5 "#ffdfba" "Update/submit and cat logs"
msg_code "make update && make submit"
echo "$BOX2"

# Stage 2 -----------------------
msg_stage 2 "Execute In-House Pipeline and Philips Uploads"

msg_step 1 "#ffffba" "After Demux finishes, check the output folder for QC file.  An SFTP client like CyberDuck can be used to view:"
#echo "Use a SFTP client like CyberDuck to view ~/molecpathlab/production/quicksilver/${runID}"
white_bg "sftp://bigpurple.nyumc.org${DEMUXDIR}/output/${runID}.report.html"
msg_step 2 "#ffffba" "Open directory permissions for others"
msg_code "chmod -R g+rwx ${DEMUXDIR}/output"
msg_step 3 "#ffffba" "Now, within the same directory, run the pass and upload *make* commands"
msg_code "cd $DEMUXDIR"
msg_code "make passed && make uploads"
msg_step 4 "#ffffba" "Next, deploy the main PACT pipeline and go to that directory"
msg_code "make deploy-NGS607"
msg_code "cd ${productionDir}/NGS607/${runID}"
msg_step 5 "#ffffba" "Update and then submit the slurm job"
msg_code "make update && make submit"
msg_step 6 "#ffffba" "Cat slurm logs to check the run started"
msg_code "squeue -u ${kerbero} && lt logs/*"
msg_step 7 "#ffffba" "Ensure the isg-uploads folder is accessible:"
msg_code "chmod -R g+rwx ${productionDir}/isg-uploads/${runID}"
msg_step 8 "#ffffba" "Check the pgm log to make sure the files are being uploaded to Philips every 30 min"

msg_code "ssh pgm@pgmlcdcpvm01.nyumc.org"
msg_code "(Enter pgm password from Evernote)"
msg_code "cat pgm/log/uploads.log"
msg_step 9 "#ffffba" "Check Philips IntelliSpace to make sure file uploads are not pending review or missing"
echo "$BOX2"

# Stage 3 -----------------------
msg_stage 3 "Generate QC"
msg_step 1 "#d9d2e9" "Execute the following commands in BigPurple after pipeline completes"
msg_code "module load python/cpu/3.8.11"
msg_step 2 "#d9d2e9" "Go to the run production directory:"
msg_code "cd ${molecDir}${runID}"
msg_step 3 "#d9d2e9" "Execute Python Automations"
msg_code "python ${reprtDir}snp_overlap.py -o ${molecDir}${runID}/output && python ${reprtDir}generate_html_report.py -o ${molecDir}${runID} -p ${pactRun} -r ${runID}"
msg_code "python3 ${reprtDir}variants_qc.py -rid ${runID} -rdir ${molecDir}${runID}/output/ -pactid ${pactRun}"
msg_step 4 "#d9d2e9" "Change permissions for group access"
msg_code "chmod -R g+rwx ${productionDir}/NGS607/${runID}/output"
echo "$BOX2"

# Stage 4 -----------------------
msg_stage 4 "Copy the QC and Output to Zdrive"
msg_step 1 "#baffc9" "In BigPurple, mount the molecular drive to your data mover node"
msg_code "ssh ${kerbero}@dmn-0002"
msg_code "mount /mnt/${kerbero}/molecular"
msg_note "NOTE:" "Some times the 0001 will be down and the 0002 will take over."
msg_step 2 "#baffc9" "Once mounted, create the directories in /MOLECULAR/NGS607/"
msg_code "mkdir -p \"${zdrive}/NGS607/${currYear}/${runID}/output/alignments\" \"/mnt/${kerbero}${outputDir}${currYear}/${pactRun}\""
msg_code "chmod -R g+rwx \"/mnt/${kerbero}${outputDir}${currYear}/${pactRun}\" && chmod -R g+rwx /gpfs/home/${kerbero}/molecpathlab/production/NGS607/${runID}/output"
msg_step 3 "#baffc9" "Once created, copy files to /MOLECULAR/NGS607/"
msg_white "<b>Copy Alignments and Annotations:</b>"
msg_code "rsync -vrthP ${rsyncDir}/alignments/${FG_GRN}deduplicated${normal} ${zdrive}/NGS607/${currYear}/${runID}/output/${FG_GRN}alignments/${normal} && rsync -vrthP ${rsyncDir}/${FG_GRN}annotations${normal} ${zdrive}/NGS607/${currYear}/${runID}/output/${FG_GRN}annotations/${normal}"
msg_white "<b>Copy Clinical Folder and .tsv files:</b>"
msg_code "rsync -vrthP ${rsyncDir}/${FG_GRN}clinical${normal} ${zdrive}/NGS607/${currYear}/${runID}/output/ && rsync -vrthP ${rsyncDir}/${FG_GRN}*.tsv${normal} ${zdrive}/NGS607/${currYear}/${runID}/"
msg_step 4 "#baffc9" "Next, copy files to /MOLECULAR LAB ONLY/NYU PACT Patient Data/"
msg_code "cd \"/mnt/${kerbero}${outputDir}${currYear}/${pactRun}/\" && rsync -vrthP ${rsyncDir}/${FG_GRN}clinical/${normal} ./"
msg_step 5 "#baffc9" "<b>Optional--Rsync the PDF facet files and  QC tsv for concensus</b>"
msg_code "rsync -vrthP ${productionDir}/NGS607/${runID}/output/cnv/FACETS/*.pdf ${productionDir}/NGS607/${runID}/${pactRun}-QC.tsv ${zdrive}/REDCap/cnv_facets/${pactRun}/"
msg_step 6 "#baffc9" "Email notify the PACT team once the QC files are copied"
white_bg "From ${productionDir}/NGS607/${runID}/output/clinical/:"
white_bg "${pactRun}-Somatic_Variants.html <b>AND</b> ${pactRun}-Germline_Variants.html <b>AND</b> ${pactRun}.html"
msg_step 7 "#baffc9" "Check the following directories copied into Z-drive ${runID} folder"
white_bg "From ${productionDir}/NGS607/${runID}/output"
white_bg "/clinical"
white_bg "/annotations/"
white_bg "/alignments/deduplicated"
echo "<h4 style='font-size:14;margin:2px;padding:2px;'>${FG_YLW}Only transfer all .tsv, alignments, annotations, clinical${normal}</h4>"
msg_note "NOTE:" "The command runs as \"assume we are in the destination directory\""
msg_white "-v=verbose"
msg_white "-r=recrusive"
msg_white "-t=preserve modification times"
msg_white "-h=human readable outputs"
msg_white "-P=preserve partially transferred files (--partial) + print information showing the progress of the transfer (--progress)"
echo "$BOX2"
# Stage 5 -----------------------
msg_stage 5 "Consensus Files"
msg_step 1 "#bae1ff" "Create a new directory for your PACT consensus"
msg_code "mkdir \"${consensusDir}${pactRun}_consensus\""
msg_code "cd \"${consensusDir}${pactRun}_consensus\""
msg_step 2 "#bae1ff" "Download the template RMD file and csv description file from GitHub to that folder"
msg_code "curl -o ${pactRun}_desc.csv -L ${pactGithub}/PACT_desc.csv -s && curl -o ${pactRun}_consensus.Rmd -L ${pactGithub}/PACT_consensus.Rmd -s"
msg_step 3 "#bae1ff" "To this directory, copy all the .cnv.plot.pdf facets files and Samplesheets"
msg_code "cp /Volumes/molecular/Molecular/REDCap/cnv_facets/${pactRun}/*.pdf ~/Desktop/${runID}-SampleSheet.csv ~/Desktop/${pactRun}_MethylMatch.xlsx ./"
msg_step 4 "#bae1ff" "Next copy the QC file below to the consensus Rmd directory${normal}"
msg_code "cp /Volumes/molecular/Molecular/NGS607/${currYear}/${runID}/output/clinical/${pactRun}-Somatic_Variants.html /Volumes/molecular/Molecular/NGS607/${currYear}/${runID}/output/clinical/${pactRun}.html /Volumes/molecular/Molecular/REDCap/cnv_facets/${pactRun}/${pactRun}-QC.tsv ./"
msg_step 5 "#bae1ff" "Use the Evernote Guide to download the indels from our in-house QC html file as a .csv file and then run the Rscript:"
white_bg "<b>${evernoteLink}</b>"
msg_code "RScript devtools::source_url('https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts/MakeIndelList.R')"
msg_step 6 "#bae1ff" "Open Philips Intellispace and the ${pactRun}_desc.csv file and annotate the variant calls"
white_bg "<b>${consensusDir}${pactRun}_consensus/${pactRun}_desc.csv</b>"
msg_note "NOTE -" "Make sure you Rsync the facet files for the run from the data mover node to the Z-drive:"
msg_code "rsync -vrthP ${productionDir}/NGS607/${runID}/output/cnv/FACETS/*.pdf ${zdrive}/REDCap/cnv_facets/${pactRun}/"
msg_step 7 "#bae1ff" "Once you have completely filled the consensus description file, knit the Rmd in Rstudio and view it"
msg_code "Rscript --verbose -e \"rmarkdown::render('${pactRun}_consensus.Rmd', params=list(pactName='${pactRun}',userName='Jonathan'))\""
msg_step 8 "#bae1ff" "After the concensus html file is created, copy to the output folder"
msg_code "cp ${consensusDir}${pactRun}_consensus/${pactRun}_consensus.html \"/Volumes${outputDir}${currYear}/${pactRun}/\""
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

        setTimeout(() => {
            button.innerText = copyButtonLabel;
        }, 1000)

    }
</script>
"
