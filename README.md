# Clinical Methylation Classifier Setup

## Methylation Pipline Overview
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/1bcd4fcdb6fb8c1908cb2d38fcfc7cd2ffffe8a2/screenshots/meth_pipeline_uml.png" alt="drawing" width="100%"/><br/>

## Essential Downloads
Download and install the following packages:<br>
1. **R 4.1** from CRAN: https://cran.r-project.org/bin/macosx/<br />
2. **RStudio 1.4 or later**: https://www.rstudio.com/products/rstudio/download/#download<br />
3. **XQuartz**: https://www.xquartz.org/<br />
4. **LaTeX** for Mac: https://www.tug.org/mactex/mactex-download.html<br />
5. **Pandoc**: https://pandoc.org/installing.html<br />
6. **XCode 12.0** or higher for Mac OS: https://developer.apple.com/download/all/?q=xcode<br>
7. **Java 8 JDK**: https://www.oracle.com/java/technologies/downloads/#java8-mac
8. **Rswitch**: https://rud.is/rswitch/<br />
9. **Homebrew**: https://brew.sh/ you can install using the following line in terminal:<br />
`/bin/bash -c $(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh`<br />

**NOTE**<br>
- R v4.1 includes compile and Tckl dependencies. brew can install libomp and cairo if needed.
- After downloading R & RStudio unlocked the
[System Preferences Privacy & Security Panel](https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/Notes/SystemPermissions.md) before installing packages.
___
## Network Drive Mount Paths
### To install & run the pipeline, it is critical to mount the following network smb shared drives:

Open Finder and press **⌘(CMD) + K** then paste each of these directories, login name and password is your NYUMC\KerberosID
<br>
`smb://research-cifs.nyumc.org/Research/CBioinformatics/`<br />
`smb://research-cifs.nyumc.org/Research/snudem01lab/snudem01labspace`<br />
`smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular`<br />

# Install pipeline and start a Run in Terminal
To run the pipeline from your terminal, simply execute the following command:<br />
`/Volumes/CBioinformatics/Methylation/runMeth.sh 21-MGDM_TEST`<br />
___
### Input Paths
*Files are copied to the work directory by their RUNID name and YEAR, including the worksheet and idats for example:*
- **Worksheets** `/Volumes/molecular/MOLECULAR LAB ONLY/NYU-METHYLATION/WORKSHEETS/2022/22-MGDM17.xlsm`<br />
- **.idat files**`/Volumes/molecular/Molecular/iScan/`
### Default Working Directory
- **Input files are copied and report files are generated on the Cbioinformatics drive:**<br />
`/Volumes/CBioinformatics/Methylation/Clinical_Runs/22-MGDM17`<br />
### Output Paths
- **Html report files saved to the working directory are copied to the Z-drive**<br/>
For example, run 22-MGDM17 report files would be output in the following directories:<br/>
`/Volumes/molecular/Molecular/MethylationClassifier/2022/22-MGDM17`<br>
`/Volumes/molecular/MOLECULAR LAB ONLY/NYU-METHYLATION/Results/2022/22-MGDM17`
___
## **Executing Methylation CLI**
To run the Clinical or Research Methylation pipeline, simply use the locally stored Shell Script in:<br>
`/Volumes/CBioinformatics/Methylation/runMeth.sh`<br/>
This shell script uses Curl to download the files from this repo and executes methylExpress.R in the terminal.
<br>
You can copy runMeth.sh and create an alias or symlink to execute more easily.  For example:<br>
`alias runmeth='bash ~/script/runmeth.sh'`
or
`echo "alias runmeth='bash ~/script/runmeth.sh'" >> ~/.bashrc`

### The shell script takes the following argument parameters:<br>
**`methAPI='XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'`** # RedCap API Token <br>
**`methRun=${1-NULL}`** # methylation run id e.g. MDGM22-3<br>
`PRIORITY=${2-NULL}` # string of prioritized RD-numbers<br>
`runPath=${3-NULL}` # any custom directory to copy/run the idat files<br>
`redcapUp=${4-NULL}` # to upload to redcap or not if server down single char i.e. "T" or "F"<br>
<br>
**runmeth.sh downloads methylExpress.R and other files using curl:**<br>
`curl -o methylExpress.R -L https://git.io/JWujj`<br>
**Finally runmeth.sh passes your input parameters as args to methylExpress.R:**<br>
`Rscript --verbose methylExpress.R` **`$methAPI` `$methRun` `$PRIORITY` `$runPath`**<br>

### The four positional arguments from *runmeth.sh* are passed to the Rscript *methylExpress.R*:<br />
`arg[1]` is the **token** for the API call ('#######################')<br />
`arg[2]` is the **RunID** which if NULL runs the latest Clinical Worksheet ('MR21-099')<br />
`arg[3]` is the **selectRds** parameter which is to prioritize samples being run (NULL)<br />
`arg[4]` is the **baseFolder** parameter which is optional if you want to run/save output to a different directory (NULL)<br />

Alternatively, instead of passing the RunID to runmeth.sh, you can source and download this repository and then locally edit args in [methylExpress.R](https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/R/methylExpress.R) to run manually.

## Run the Test Case after installation with the following command:<br />
`/Volumes/CBioinformatics/Methylation/runMeth.sh 21-MGDM_TEST`<br />

## Troubleshooting

### Pipeline Installation Issues
If you have issues with package installation or dependencies:<br />
Make sure compilers are installed by opening Xcode.app or executing `sudo xcode-select --install`<br>
Then, execute the all_installer.R script by copy and pasting the raw contents of the script below into Rstudio before running runmeth.sh again<br />
https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/Research/all_installer.R<br />
To resolve any problems during automation, you can open methylExpress.R in RStudio which is downladed by runmeth.sh to your home directory.<br />

### REDCap errors
Once your run completes check in your run directory if there is any *upload_log.tsv* file or *redcaperrors.txt*.  If these files exist, they may note any files or data which would have been over-written in the database.  Check with the wet lab if any RD-numbers were duplicated or previously used for the samples listed in the upload_log.tsv file. <br>
**NOTE** When running the test case run *21-MGDM_TEST* you may notice an error with the upload log as these reports would already exist in REDCap.  It is normal for the test case to fail uploading since the REDCap database already contains the data and files for the test run.

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
