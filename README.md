# How To Setup Clinical Methylation Classifier

## Methylation Pipline Overview
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/1bcd4fcdb6fb8c1908cb2d38fcfc7cd2ffffe8a2/screenshots/meth_pipeline_uml.png" alt="drawing" width="100%"/><br/>
<span style="color:blue">
## Essential Downloads
The classifier runs on R version 3.6.3 and up.  It is not compatible with 3.3.3<br />
1. **R 4.1** from CRAN: [https://cran.r-project.org/bin/macosx/](https://cran.r-project.org/bin/macosx/)<br />
2. **RStudio 1.4**: [https://download1.rstudio.org/desktop/macos/RStudio-1.4.1717.dmg](https://download1.rstudio.org/desktop/macos/RStudio-1.4.1717.dmg)<br />
3. **XQuartz**: [https://github.com/XQuartz/XQuartz/releases/download/XQuartz-2.8.1/XQuartz-2.8.1.dmg](https://github.com/XQuartz/XQuartz/releases/download/XQuartz-2.8.1/XQuartz-2.8.1.dmg)<br />
4. **LaTeX** for Mac: [https://www.tug.org/mactex/mactex-download.html](https://www.tug.org/mactex/mactex-download.html)<br />
5. **Homebrew**: https://brew.sh/ you can install using the following line in terminal:<br />
`/bin/bash -c $(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh`<br />
6. **Pandoc**: https://pandoc.org/installing.html<br />
7. **XCode 12.0** or higher for Mac OS: [https://developer.apple.com/download/all/?q=xcode](https://developer.apple.com/download/all/?q=xcode)<br>
8. **Java 8 JDK**:
https://www.oracle.com/java/technologies/downloads/#java8-mac
9. **[Rswitch](https://rud.is/rswitch/)**<br />
</span>
### NOTE
- R v4.1 includes compile and Tckl dependencies. brew can install libomp and cairo if needed.
- After downloading R & RStudio **do not install** until you have unlocked the
[System Preferences Privacy & Security Panel](https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/Notes/SystemPermissions.md).
___
## Network Drive Mount Paths
### To install & run the pipeline, it is critical to mount the following network smb shared drives:

Open Finder and press **⌘(CMD) + K** then paste each of these directories, login name and password is your NYUMC\KerberosID
<br>
`smb://research-cifs.nyumc.org/Research/CBioinformatics/`<br />
`smb://research-cifs.nyumc.org/Research/snudem01lab/snudem01labspace`<br />
`smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular`<br />
___
<span style="color:green">
# Install pipeline and start a Run in Terminal
### To run the pipeline from your terminal, simply execute the following command:<br />
`/Volumes/CBioinformatics/Methylation/runMeth.sh 21-MGDM_TEST`<br />

If you have issues with the automation, you can open methylExpress.R which downloads to your home directory in RStudio<br />
</span>
___
## Input Paths
### Files are copied to the work directory by their RUNID name and YEAR, including the worksheet and idats for example:<br />
#### Worksheet
/Volumes/molecular/MOLECULAR LAB ONLY/NYU-METHYLATION/WORKSHEETS/2022/22-MGDM17.xlsm
#### .idat files input directory:
/Volumes/molecular/Molecular/iScan/
## Output Paths
### Files are saved to MethylationClassifier/YEAR/RUNID, for example:<br/>
/Volumes/molecular/Molecular/MethylationClassifier/2022/22-MGDM17<br>
/Volumes/molecular/Molecular/MethylationClassifier/2020/20-MGDM5
___
### There are two system Rscript to run methylExpress.R with the arguments in order:<br />
`arg[1]` is the **token** for the API call ('#######################')<br />
`arg[2]` is the **RunID** which if NULL runs the latest Clinical Worksheet ('MR21-099')<br />
`arg[3]` is the **selectRds** parameter which is to prioritize samples being run (NULL)<br />
`arg[4]` is the **baseFolder** parameter which is optional if you want to run/save output to a different directory (NULL)<br />

Alternatively, you can source then run the github script locally using [methylExpress.R](https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/R/methylExpress.R)

## Run the Test Case after installation with the following command:<br />
`/Volumes/CBioinformatics/Methylation/runMeth.sh 21-MGDM_TEST`<br />

# Email List of Methylation samples which are PACT run

`curl -o PactMethMatch.R -L https://git.io/J41Wp; Rscript --verbose PactMethMatch.R 'TOKENAPI12345667891011' '/Users/PATH/TO/CSV/Desktop/210715_NB501073_9999_ABCDEFGHIJRLK-SampleSheet.csv'`

# Email PACT csv file:
`curl -o pactParse.R -L https://git.io/J0kfR; Rscript --verbose pactParse.R 'TOKENAPI12345667891011' 'PACT-21-##'`

## **RUNNING METHYLATION CLI**
To run the Clinical or Research Methylation pipeline, simply use the locally stored Shell Script in:<br>
`/Volumes/CBioinformatics/Methylation/runMeth.sh`
You can copy runMeth.sh and create an alias or symlink to execute more easily.  For example:<br>
`alias runmeth='bash ~/script/runmeth.sh'`
or
`echo "alias runmeth='bash ~/script/runmeth.sh'" >> ~/.bashrc`

The shell script takes the following argument parameters:<br>
`#!/bin/bash`<br>
`methAPI='XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' # RedCap API Token`<br>
`methRun=${1-NULL} # methylation run id e.g. MDGM22-3`<br>
`PRIORITY=${2-NULL} # string of prioritized RD-numbers`<br>
`runPath=${3-NULL} # any custom directory to copy/run the idat files`<br>
`redcapUp=${4-NULL} # to upload to redcap or not if server down single char i.e. "T" or "F"`<br>

`curl -o methylExpress.R -L https://git.io/JWujj; Rscript --verbose methylExpress.R $methAPI $methRun $PRIORITY $runPath`<br>

You can locally copy or symlink the runMeth.sh file to execute more easily<br>

## Troubleshooting

### How to upload manually to REDCap
1. Login with your kerberos ID to https://redcap.nyumc.org/
2. On the left-hand sidebar scroll all the way down the Reports Bookmarks until you see the folder: >>>>CURRENT Runs~~~~~ and 3) >>>>>CLINICAL Current Run
3. Here, you can click on the RD-number of choice and then select "Upload html file" under the methylation menue
4. Optionally, you can also select "  Add / Edit Records" menu in the left sidebar and find your RD-number in the "Search query" field

### How to check worksheet
1. In your current working directory or in the path /Volumes/CBioinformatics/Methylation/Clinical_Runs/22-MGDM##, open the RUNID.xlsm file
2. On the Review ribbon, click unprotect worksheet and unprotect tab
3. Right-click the "worksheet" tab at the bottom and unhide... raw_labels tab
4. If any "#ref" errors either drag the formula down to correct or type "=" and select the cell in the first tab "worksheet" and press return.  For example "=worksheet!B25" references cell B25 in the tab named "worksheet"
