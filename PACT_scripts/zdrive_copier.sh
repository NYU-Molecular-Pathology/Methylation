#!/bin/bash

# MAIN ARGS
runID=${1-NULL}   # if arg $1 is empty assign NULL as default else i.e. 123456_NB501073_0212_AHT3V7BGXK
pactRun=${2-NULL} # if arg $2 is empty assign NULL as default else i.e. PACT-21-28
kerberos=${4-$USER} # if arg $3 is empty assign $USER as default else i.e. whoami kerberosid

# DERIVED ARGS
kerbero="${kerberos}"
gpfsHome="/gpfs/home/${kerbero}/"
productionDir="/gpfs/data/molecpathlab/production"
outputDir="/molecular/MOLECULAR LAB ONLY/NYU PACT Patient Data/Results/Bioinformatics/"
rsyncDir="${gpfsHome}molecpathlab/production/NGS607/${runID}/output"
zdrive="/mnt/${kerbero}/molecular/Molecular"
runMid=${pactRun:5:2}
currYear="20$runMid"

#Once mounted, create the output directories in /MOLECULAR/NGS607/
mkdir -p "${zdrive}/NGS607/${currYear}/${runID}/output/alignments" "/mnt/${kerbero}${outputDir}${currYear}/${pactRun}"
chmod -R g+rwx "/mnt/${kerbero}${outputDir}${currYear}/${pactRun}"
chmod -R g+rwx "${gpfsHome}molecpathlab/production/NGS607/${runID}/output"

# Once created, rsync the files from BigPurple to /MOLECULAR/NGS607/
rsync -vrthP "${rsyncDir}/alignments/recalibrated" "${zdrive}/NGS607/${currYear}/${runID}/output/alignments/"
rsync -vrthP "${rsyncDir}/annotations" "${zdrive}/NGS607/${currYear}/${runID}/output/annotations/"
rsync -vrthP "${rsyncDir}/clinical" "${zdrive}/NGS607/${currYear}/${runID}/output/"
rsync -vrthP "${rsyncDir}/*.tsv" "${zdrive}/NGS607/${currYear}/${runID}/"

#Go to /MOLECULAR LAB ONLY/NYU PACT Patient Data/ on the Z-drive and copy the remaining clinical output files, PDF facets, and QC tsv "
cd "/mnt/${kerbero}${outputDir}${currYear}/${pactRun}/"
rsync -vrthP "${rsyncDir}/clinical/" ./
rsync -vrthP "${productionDir}/NGS607/${runID}/output/cnv/FACETS/*.pdf" "${productionDir}/NGS607/${runID}/${pactRun}-QC.tsv" "${zdrive}/REDCap/cnv_facets/${pactRun}/"
