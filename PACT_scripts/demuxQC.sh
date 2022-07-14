#!/bin/bash
runID=${1-NULL} # if arg $1 is empty assign NULL as default else i.e. 123456_NB501073_0212_AHT3V7BGXK
pactRun=${2-NULL} # if arg $2 is empty assign NULL as default else i.e. PACT-21-28
kerbero=${3-$USER} # if arg $3 is empty assign $USER as default else i.e. whoami kerberosid


bold=$(tput bold) # makes console output text bold
undLN=$(tput smul) # makes console output text underline
BG_BLUE="$(tput setab 4)" # makes console output text background blue
BG_CYAN="$(tput setab 6)" # makes console output text background cyan
BG_MAG="$(tput setab 5)" # makes console output text background magenta
BG_GRN="$(tput setab 2)" # makes console output text background green
FG_YLW="$(tput setaf 3)" # makes console output text yellow color
FG_RED="$(tput setaf 1)" # makes console output text red color
FG_GRN="$(tput setaf 2)" # makes console output text color green
FG_BLU="$(tput setaf 4)" # makes console output text color blue
FG_CYA="$(tput setaf 6)" # makes console output text color cyan
normal=$(tput sgr0) # resets default console output text


SHEETDIR="/gpfs/data/molecpathlab/production/samplesheets/LG-PACT/${FG_YLW}$runID${normal}"
DEMUXDIR="/gpfs/data/molecpathlab/production/Demultiplexing/${FG_YLW}$runID${normal}"
molecDir="/gpfs/data/molecpathlab/production/NGS607/"
reprtDir="/gpfs/data/molecpathlab/bin/QC_reprot/"
currYear=`date +"%Y"` #2022

lastTwo=${pactRun: -2}
runMid=${pactRun:5:2}

echo " "
echo "${bold}${FG_YLW}------------Your Input------------${normal}"
echo "${bold}1. ${BG_MAG}PACT RUNID:${normal} $runID"
echo "${bold}2. ${BG_GRN}PACT Run Name:${normal} $pactRun"
echo "${bold}3. ${BG_BLUE}Kerberos ID:${normal} $kerbero"
echo " "
echo " "

echo "${bold}${FG_YLW}------------STAGE 0: Create SampleSheet.csv & Copy to BigPurple------------${normal}"

echo "${bold}${BG_CYAN}1. First mount the Molecular drive and execute the script:${normal}"
echo "diskutil mountdisk /Volumes/CBioinformatics/"
echo "/Volumes/CBioinformatics/PACT/parsepact.sh ${FG_RED}$pactRun${normal} ${FG_CYA}$kerbero${normal}"
echo "${bold}${BG_CYAN}2. Review the generated SampleSheet.csv ouput below and email to wetlab:${normal}"
echo "~/Desktop/${FG_YLW}$runID${normal}-SampleSheet.csv"
echo "${bold}${BG_CYAN}3. Check the script copied with group read/write permissions to the folder:${normal}"
echo "/gpfs/data/molecpathlab/production/samplesheets/LG-PACT/${FG_YLW}$runID${normal}/$runID-SampleSheet.csv"

echo "${bold}${FG_YLW}------------STAGE 1: Demux Steps------------${normal}"

echo "${bold}${BG_BLUE}1. Once sequencing is finished, start demultiplexing by logging into BigPurple:${normal}"
echo "ssh -Y ${FG_CYA}$kerbero${normal}@bigpurple.nyumc.org"
echo "${bold}${BG_BLUE}2. Go into demux-nf2:${normal}"
echo "cd ~/molecpathlab/pipelines/demux-nf2/"
echo "${bold}${BG_BLUE}3. Execute the deploy command:${normal}"
echo "make deploy RUNID=${FG_YLW}$runID${normal} SAMPLESHEET=$SHEETDIR/${FG_YLW}$runID${normal}-SampleSheet.csv SEQTYPE=NGS607"
echo "${bold}${BG_BLUE}4. Go into the Demultiplexing run folder:${normal}"
echo "cd $DEMUXDIR"
echo "${bold}${BG_BLUE}5. Update/submit and cat logs:${normal}"
echo "make update"
echo "make submit"

echo "${bold}${FG_YLW}------------STAGE 2: In-House and Philips Uploads------------${normal}"

echo "${bold}${BG_MAG}1. Once Demux finishes, check output dir for QC file:${normal}"
#echo "Use a SFTP client like CyberDuck to view ~/molecpathlab/production/quicksilver/${FG_YLW}$runID${normal}"
echo "sftp://bigpurple.nyumc.org/gpfs/data/molecpathlab/production/Demultiplexing/${FG_GRN}$runID${normal}/output/"
echo "${bold}${BG_MAG}Open directory permissions for others:${normal}"
echo "chmod -R g+rwx /gpfs/data/molecpathlab/production/Demultiplexing/${FG_GRN}$runID${normal}/"
echo "${FG_GRN}$runID${normal}.report.html"
echo "${bold}${BG_MAG}2. Now, from the same demultiplexing folder, run a pass command and upload:${normal}"
echo "cd $DEMUXDIR"
echo "make passed"
echo "make uploads"
echo "${bold}${BG_MAG}3. Next, deploy the main PACT pipeline and go to that directory:${normal}"
echo "make deploy-NGS607"
echo "cd /gpfs/data/molecpathlab/production/NGS607/${FG_YLW}$runID${normal}"
echo "${bold}${BG_MAG}4. Update and then submit the slurm job:${normal}"
echo "make update"
echo "make submit"
echo "${bold}${BG_MAG}5. Cat slurm logs to check the run started:${normal}"
echo "squeue -u ${FG_CYA}$kerbero${normal}"
echo "lt logs/*"
echo "${bold}${BG_MAG}6. Make Sure the isg-uploads folder is accessible:${normal}"
echo "chmod -R g+rwx /gpfs/data/molecpathlab/production/isg-uploads/${FG_YLW}$runID${normal}"
echo "${bold}${BG_MAG}7. Check the pgm log to make sure the files are being uploaded to Philips every 30 min:${normal}"
echo "ssh pgm@pgmlcdcpvm01.nyumc.org"
echo "enter password from Evernote"
echo "cat pgm/log/uploads.log"
echo "${bold}${BG_MAG}8. Check Philips IntelliSpace to make sure file uploads are not pending review or missing${normal}"
echo " "

echo "${bold}${FG_YLW}------------STAGE 3: QC Generation------------${normal}"

echo "${bold}${BG_GRN}1. Execute the following commands in BigPurple after pipeline completes:${normal}"
echo "module load python/cpu/3.8.11"
echo "${bold}${BG_GRN}2. Change permissions to the output directory and then cd into it:${normal}"
echo "chmod -R g+rwx /gpfs/data/molecpathlab/production/NGS607/${FG_YLW}$runID${normal}/"
printf "cd $molecDir${FG_YLW}$runID${normal}\n"
#echo "${bold}${BG_GRN}2. Deploy the QC make:${normal}"
#echo "make update"
#printf "make deploy PRODDIR=$molecDir${FG_YLW}$runID${normal} RUNID=${FG_YLW}$runID${normal} FASTQDIR=$molecDir${FG_YLW}$runID${normal}/output/reads/trimmed\n"
echo "${bold}${BG_GRN}3. Execute Python Automations:${normal}"
printf "python ${reprtDir}snp_overlap.py -o $molecDir${FG_YLW}$runID${normal}/output\n\n"
printf "python ${reprtDir}generate_html_report.py -o $molecDir${FG_YLW}$runID${normal} -p ${FG_RED}$pactRun${normal} -r ${FG_YLW}$runID${normal}\n\n"
printf "python3 ${reprtDir}variants_qc.py -rid ${FG_YLW}$runID${normal} -rdir $molecDir${FG_YLW}$runID${normal}/output -pactid ${FG_RED}$pactRun${normal}\n\n"

echo "${bold}${FG_YLW}------------STAGE 4: Copy the QC and Output to Zdrive------------${normal}"

echo "${bold}${BG_CYAN}  1. Create the following folder in the Molecular Drive:${normal}"
echo "mkdir \"/Volumes/molecular/MOLECULAR LAB ONLY/NYU PACT Patient Data/Results/Bioinformatics/${FG_BLU}20$runMid${normal}/${FG_RED}$pactRun${normal}\""
echo "${bold}${BG_CYAN}  2. Create the following folder in the Zdrive:${normal}"
echo "mkdir -p \"/Volumes/molecular/Molecular/NGS607/${FG_BLU}20$runMid${normal}/${FG_YLW}$runID${normal}/output/alignments\""
echo "${bold}${BG_CYAN}  3. In BigPurple, mount the molecular drive to your data mover node:${normal}"
echo "ssh ${FG_CYA}$kerbero${normal}@dmn-0002"
echo "mount /mnt/${FG_CYA}$kerbero${normal}/molecular"
echo "${bold}${FG_YLW}NOTE:${normal} ${undLN}Some times the 0001 will be down and the 0002 will take over.${normal}"
echo " "
echo "${bold}${BG_CYAN}  4. Once mounted, copy files to ${normal}${BG_CYAN}/MOLECULAR/NGS607/:${normal}"
echo "${undLN}Copy Alignments${normal}"
echo "rsync -vrthP /gpfs/home/${FG_CYA}$kerbero${normal}/molecpathlab/production/NGS607/${FG_YLW}$runID${normal}/output/alignments/${FG_GRN}deduplicated${normal} /mnt/${FG_CYA}$kerbero${normal}/molecular/MOLECULAR/NGS607/${FG_BLU}20$runMid${normal}/${FG_YLW}$runID${normal}/output/${FG_GRN}alignments/${normal}"
echo " "
echo "${undLN}Copy Annotations${normal}"
echo "rsync -vrthP /gpfs/home/${FG_CYA}$kerbero${normal}/molecpathlab/production/NGS607/${FG_YLW}$runID${normal}/output/${FG_GRN}annotations${normal} /mnt/${FG_CYA}$kerbero${normal}/molecular/MOLECULAR/NGS607/${FG_BLU}20$runMid${normal}/${FG_YLW}$runID${normal}/output/${FG_GRN}annotations/${normal}"
echo " "
echo "${undLN}Copy Clinical${normal}"
echo "rsync -vrthP /gpfs/home/${FG_CYA}$kerbero${normal}/molecpathlab/production/NGS607/${FG_YLW}$runID${normal}/output/${FG_GRN}clinical${normal} /mnt/${FG_CYA}$kerbero${normal}/molecular/MOLECULAR/NGS607/${FG_BLU}20$runMid${normal}/${FG_YLW}$runID${normal}/output/"
echo " "
echo "${undLN}Copy .tsv${normal}"
echo "rsync -vrthP /gpfs/home/${FG_CYA}$kerbero${normal}/molecpathlab/production/NGS607/${FG_YLW}$runID${normal}/output/${FG_GRN}*.tsv${normal} /mnt/${FG_CYA}$kerbero${normal}/molecular/MOLECULAR/NGS607/${FG_BLU}20$runMid${normal}/${FG_YLW}$runID${normal}/"
echo " "
echo "${bold}${BG_CYAN}  5. Next, copy files to ${normal}${BG_CYAN}/MOLECULAR LAB ONLY/NYU PACT Patient Data/:${normal}"
echo "cd \"/mnt/${FG_CYA}$kerbero${normal}/molecular/MOLECULAR LAB ONLY/NYU PACT Patient Data/Results/Bioinformatics/${FG_BLU}20$runMid${normal}/${FG_RED}$pactRun${normal}/\""
echo "rsync -vrthP /gpfs/home/${FG_CYA}$kerbero${normal}/molecpathlab/production/NGS607/${FG_YLW}$runID${normal}/output/${FG_GRN}clinical/${normal} ./"
echo ""
echo "Optional: Rsync the facets files for concensus"
echo "rsync -vrthP /gpfs/data/molecpathlab/production/NGS607/${FG_YLW}$runID${normal}/output/cnv/FACETS/*.pdf /mnt/serraj10/molecular/Molecular/REDCap/cnv_facets/${FG_RED}$pactRun${normal}/"
echo "/gpfs/data/molecpathlab/production/NGS607/${FG_YLW}$runID${normal}/${FG_RED}$pactRun${normal}-QC.tsv /mnt/serraj10/molecular/Molecular/REDCap/cnv_facets/${FG_RED}$pactRun${normal}/"
echo "${bold}${BG_CYAN}  6. Email notify the PACT team once the QC files are copied:${normal}"
echo "From /gpfs/data/molecpathlab/production/NGS607/${FG_YLW}$runID${normal}/output/clinical/"
echo "${FG_RED}$pactRun${normal}-Somatic_Variants.html and ${FG_RED}$pactRun${normal}-Germline_Variants.html"
echo "${FG_RED}$pactRun${normal}.html"
echo "${bold}${BG_CYAN}  7. Check the following directories copied into Z-drive ${FG_YLW}$runID${normal}${bold}${BG_CYAN} folder:${normal}"
echo "From /gpfs/data/molecpathlab/production/NGS607/${FG_YLW}$runID${normal}/output"
echo "/clinical"
echo "/annotations"
echo "/alignments/deduplicated"
echo "${bold}${FG_YLW}Only transfer all .tsv, alignments, annotations, clinical${normal}"
echo "${bold}${FG_YLW}NOTE:${normal} ${undLN}The command runs as  \"assume we are in the destination directory\"${normal}
-v verbose 
-r recrusive 
-t preserve modification times 
-h human readable outputs
-P preserve partially transferred files (--partial) + print information showing the progress of the transfer (--progress)"

echo "${bold}${FG_YLW}------------STAGE 5: Consensus Files------------${normal}"

echo "${bold}${BG_BLUE}1. Create a new directory for your PACT consensus:${normal}"
echo "mkdir \"/Volumes/CBioinformatics/jonathan/pact/consensus/${FG_RED}$pactRun${normal}_consensus\""
echo "cd \"/Volumes/CBioinformatics/jonathan/pact/consensus/${FG_RED}$pactRun${normal}_consensus\""
echo "${bold}${BG_BLUE}2. Download the template RMD file and csv description file from GitHub to that folder:${normal}"
echo "curl -o ${FG_RED}$pactRun${normal}_desc.csv -L https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts/PACT_desc.csv -s"
echo "curl -o ${FG_RED}$pactRun${normal}_consensus.Rmd -L https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts/PACT_consensus.Rmd -s"
echo "${bold}${BG_BLUE}3. To this directory, copy all the .cnv.plot.pdf facets files from:${normal}"

echo "cp /Volumes/molecular/Molecular/REDCap/cnv_facets/${FG_RED}$pactRun${normal}/*.pdf ./"
echo "${bold}${BG_BLUE}4. Copy the following files to this concenus Rmd directory from the Desktop:${normal}"
echo "cp ~/Desktop/${FG_RED}$pactRun${normal}_MethylMatch.xlsx ./"
echo "cp ~/Desktop/${FG_YLW}$runID${normal}-SampleSheet.csv ./"
echo "${bold}${BG_BLUE}5. Next copy the QC file below to the consensus Rmd directory${normal}"
echo "cp /Volumes/molecular/Molecular/REDCap/cnv_facets/${FG_RED}$pactRun${normal}/*.tsv ./"
echo "cp /Volumes/molecular/Molecular/NGS607/$currYear/${FG_YLW}$runID${normal}/output/clinical/${FG_RED}$pactRun${normal}.html ./"
echo "cp /Volumes/molecular/Molecular/NGS607/$currYear/${FG_YLW}$runID${normal}/output/clinical/${FG_RED}$pactRun${normal}-Somatic_Variants.html ./"
echo "${bold}${BG_BLUE}6. Open Philips Intellispace and the in-house ${FG_RED}$pactRun${normal}${bold}${BG_BLUE}.html file and annotate the abberations description file with your samples in:${normal}"
echo "/Volumes/CBioinformatics/jonathan/pact/consensus/${FG_RED}$pactRun${normal}_consensus/${FG_RED}$pactRun${normal}_desc.csv"
echo "${bold}${BG_BLUE}7. Once you have completely filled the consensus description file, knit the Rmd and view it:${normal}"
echo "Rscript --verbose -e \"rmarkdown::render('${FG_RED}$pactRun${normal}_consensus.Rmd', params=list(pactName = '${FG_RED}$pactRun${normal}', userName = 'Jonathan Serrano'))\""
echo "${bold}${BG_BLUE}8. After the concensus html file is created, copy to the output folder:${normal}"
echo "cp /Volumes/CBioinformatics/jonathan/pact/consensus/${FG_RED}$pactRun${normal}_consensus/${FG_RED}$pactRun${normal}_consensus.html \"/Volumes/molecular/MOLECULAR LAB ONLY/NYU PACT Patient Data/Results/Bioinformatics/${FG_BLU}20$runMid${normal}/${FG_RED}$pactRun${normal}/\""
echo "Test with:"
echo "Rscript --verbose -e \"rmarkdown::render('/Volumes/CBioinformatics/jonathan/pact/consensus/${FG_RED}$pactRun${normal}_consensus/${FG_RED}$pactRun${normal}_consensus.Rmd', params=list(pactName = '${FG_RED}$pactRun${normal}'))\""
