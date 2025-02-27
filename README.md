<h1>
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/fd470d3a25816e77789a0741dfa62d584bd93b7d/screenshots/clinical_pipeline.png" width="100%"/><br/> </h1>

### Table of Contents 
- [📖 Methylation Pipline Overview](#-methylation-pipline-overview)
  * [💻 Essential Downloads](#-essential-downloads)
  * [🌐 Network Drive Mount Paths](#-network-drive-mount-paths)
- [⚡️ Quickstart](#%EF%B8%8F-quickstart)</br>
     [<sup>-  Input Paths</sup>](#input-paths)</br>
     [<sup>-  Default Working Directory</sup>](#default-working-directory)</br>
     [<sup>-  Output Paths</sup>](#output-paths)</br>
- [⚙️ Executing Methylation CLI](#%EF%B8%8F-executing-methylation-cli)
    + [🤖 runmeth.sh parameters](#-runmethsh-parameters)
    + [🧮 Passing Arguments to R](#-passing-arguments-to-r)
- [🧪 Run the Test Case](#-run-the-test-case)
- [⚠️ Troubleshooting](#%EF%B8%8F-troubleshooting)


## 📖 Methylation Pipline Overview
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/1bcd4fcdb6fb8c1908cb2d38fcfc7cd2ffffe8a2/screenshots/meth_pipeline_uml.png" alt="drawing" width="100%"/><br/>

## 💻 Essential Downloads
You can automatically install these requirements using:
```
/Volumes/CBioinformatics/Methylation/install_requirements.sh
```
Note: you may need to enter your password to enable sudo permissions.
<br>
Use **ARM** (-arm64.pkg) package downloads for *M1/M2 Macs* & **Intel** (-x86_64.pkg) for older non-Apple Silicon Chip Based Macs)

- [X] **R 4.3 or higher**: https://cran.r-project.org/bin/macosx/ <br />
- [X] **Gfortran** from GitHub use dmg installer: https://github.com/fxcoudert/gfortran-for-macOS/releases/ <br />
- [X] **RStudio 2023.06.0+421 or later**: https://www.rstudio.com/products/rstudio/download/#download<br />
- [X] **XQuartz**: https://www.xquartz.org/<br />
- [X] **LaTeX** for Mac: https://www.tug.org/mactex/mactex-download.html [Direct DL](https://mirror.ctan.org/systems/mac/mactex/MacTeX.pkg) or `brew install --cask basictex`<br />
- [X] **Pandoc**: https://pandoc.org/installing.html<br />
- [X] **XCode command line tools** for Mac OS: in iTerm or Terminal enter `xcode-select --install`<br>
- [X] **Homebrew**: https://brew.sh/ <br />
- [X] **Library Magic, Sqlite and Proj**: `brew install libmagic sqlite proj tcl-tk`<br />
- [X] **Compilers+**: `brew install llvm aspell gdal autoconf automake gcc libgit2 openssl@3 zlib go pandoc git libffi`<br />

<summary>NOTE</summary>

You may need to unlock permissions before installing packages in the Mac's **System Preferences Privacy & Security Panel**:<br>
https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/Notes/SystemPermissions.md 

___

# ❗First Time Running Classifier Pre-install packages
- You can install all the dependencies above by executing the script on the CBioinformatics shared drive: <code style="color : green">text</code>/Volumes/CBioinformatics/Methylation/install_requirements.sh</code>
- After you have installed all the required system dependencies above in Essential Downloads above, you must install all the R packages needed to install and run the classifiers.
- Before running the classifier for the first time run the Rscript below, `all_installer.R`, to install any R-package dependencies.  The script only needs to be run the first time installing the classifier on new systems.<br />
- For better debugging, paste the raw code from the URL into RStudio:<br />
https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/all_installer.R <br />
___

## 🌐 Network Drive Mount Paths
- To install & run the pipeline, it is critical to mount the following network smb shared drives:
- Open Finder and press **⌘(CMD) + K** then paste each of the directories below, using NYUMC\KerberosID as the login name and password is your kerberos password. <br>
```
smb://research-cifs.nyumc.org/Research/CBioinformatics/
```
```
smb://research-cifs.nyumc.org/Research/snudem01lab/snudem01labspace
```
```
smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular
```
# ⚡️ Quickstart
#### 1. Download the shell script to your home folder or another directory:
  + You can download [runMeth.sh](https://github.com/NYU-Molecular-Pathology/Methylation/blob/32c3b043bd2fd27de4106bc56b8d4f13ac42d48d/Meth_Scripts/runMeth.sh) in this repo under Methylation/Meth_Scripts/ or use curl/wget:
 ```
 curl -# -L https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/Meth_Scripts/runMeth.sh >$HOME/runMeth.sh
 ```
#### 2. Open the shell script, paste your REDCap API token in the **methAPI** field on line 3, and save it.</br>  
 You can use `nano $HOME/runMeth.sh` </br>
 ```bash
 methAPI="XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"  #Paste your API Token here
 ```
  + Note: Your API Token can be found in "All Samples DataBase" on the left-side panel called "API" in REDCap.  Explained [here](https://redcap.nyumc.org/apps/redcap/api/help/?content=tokens)
 
#### 3. Add permissions to the script to be executable:
 ```ruby
 chmod +rwx $HOME/runMeth.sh
 ```

 + If install of any packages fail, be sure to check the troubleshooting section at the bottom of this page
___
### Input Paths
*Files are copied to the work directory by their RUNID name and YEAR, including the worksheet and idats for example:*
- **Worksheets** `/Volumes/molecular/MOLECULAR LAB ONLY/NYU-METHYLATION/WORKSHEETS/2022/22-MGDM17.xlsm`<br />
- **.idat files** `/Volumes/molecular/Molecular/iScan/`
### Default Working Directory
- **Input files are copied and report files are generated on the Cbioinformatics drive:**<br />
`/Volumes/CBioinformatics/Methylation/Clinical_Runs/22-MGDM17`<br />
### Output Paths
- **Html report files saved to the working directory are copied to the Z-drive**<br/>
For example, run 22-MGDM17 report files would be output in the following directories:<br/>
`/Volumes/molecular/Molecular/MethylationClassifier/2022/22-MGDM17`<br>
`/Volumes/molecular/MOLECULAR LAB ONLY/NYU-METHYLATION/Results/2022/22-MGDM17`
___

# ⚙️ Executing Methylation CLI
To run the Clinical or Research Methylation pipeline, simply use the locally stored Shell Script in:<br>
`/Volumes/CBioinformatics/Methylation/runMeth.sh`<br/>
- This shell script uses Curl to download the files from this repo and takes four positional argument inputs to execute methylExpress.R in the terminal.
- The bash script stores your REDCap API token locally and only requires the methylation run ID to be entered.
- You can copy runMeth.sh and create an alias or symlink to execute more easily.  For example:<br>
`alias runmeth='bash $HOME/runMeth.sh'` or `echo "alias runmeth='bash $HOME/runMeth.sh'" >> ~/.bashrc`
---
### 🤖 runMeth.sh parameters
The shell script takes the following positional arguments:<br>
```R
methAPI='XXXXXXXX' # (hardcoded) Your REDCap API Token
methRun=${1-NULL}  # methylation run id e.g. 22-MGDM17
PRIORITY=${2-NULL} # string of prioritized RD-numbers
runPath=${3-NULL}  # any custom directory to copy/run the idat files
redcapUp=${4-NULL} # to upload to redcap or not if server down single char i.e. "T" or "F"
runLocal=${5-NULL} # If the run directory should be executed without shared drives locally i.e. "T" or "F"
```
___

### 🧮 Passing Arguments to R
The four positional arguments from *runmeth.sh* are passed to the Rscript *methylExpress.R*:<br />
`arg[1]` is the **`token`** for the API call ('#######################')<br />
`arg[2]` is the **`RunID`** which if NULL runs the latest Clinical Worksheet *22-MGDM17*<br />
`arg[3]` is the **`selectRds`** parameter which is to prioritize samples being run (NULL)<br />
`arg[4]` is the **`baseFolder`** parameter which is optional if you want to run/save output to a different directory (NULL)<br />

Alternatively, instead of passing the RunID to runmeth.sh, you can source and download this repository and then locally edit args in [methylExpress.R](https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/R/methylExpress.R) to run manually.

# 🧪 Run the Test Case
After installation, test the pipeline from your terminal, by executing the test case:<br />
```ruby
$HOME/runMeth.sh 21-MGDM_TEST
```
or if you have not saved the runMeth.sh script locally:
```ruby
/Volumes/CBioinformatics/Methylation/runMeth.sh 21-MGDM_TEST
```
You can then check the output to confirm each html report was generated in the output directory: <br/>
`/Volumes/CBioinformatics/Methylation/Clinical_Runs/21-MGDM_TEST`
`ls -lha "$HOME/Desktop/html_21-MGDM_TEST/21-MGDM_TEST/"`

**NOTE**: When running the test case (*21-MGDM_TEST*) you may notice an **error** with the upload log as these reports would already exist in REDCap. It is normal for the test case html files to **fail uploading** since the REDCap database already contains the data and files for the test run, 21-MGDM_TEST.<br>
<br>

## To run the Sarcoma Classifier or re-Run Individual Samples
 - For Individual Cases: Execute the script directly with RD-numbers, for example:<br>
 `Rscript --verbose /Volumes/CBioinformatics/Methylation/Clinical_Runs/Sarcoma_runs/methylExpress_sarcoma.R RD-15-123 RD-16-1234 RD-17-321`<br>
 - For Several/Bulk Cases: Execute the script by passing the path to a csv file containing a list of RD-numbers in the first column, for example:<br>
`Rscript --verbose /Volumes/CBioinformatics/Methylation/Clinical_Runs/Sarcoma_runs/methylExpress_sarcoma.R /Path/To/Desktop/MyListRDs.csv`<br>
In the event the shared drive is not accessible, the script without the API token is availible [here](https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/methylExpress_sarcoma.R)

# ⚠️ Troubleshooting
<details>
<summary>Pipeline Installation Issues</summary>
 1. If you have issues with package installation or dependencies:<br />
Make sure compilers are installed by opening Xcode.app or executing `sudo xcode-select --install`<br>
 <br>
 2. Then, execute the all_installer.R script by copy and pasting the raw contents of the script below into Rstudio before running runmeth.sh again: https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/Research/all_installer.R<br />
<br>
 3. To resolve any problems during automation, you can open methylExpress.R in RStudio which is downladed by runmeth.sh to your home directory.<br />
<br>
 4. Try to run `sudo xcode-select -s /Library/Developer/CommandLineTools` and `brew install gdal proj` then install the package **rgdal** in Rstudio.<br />
 <br>
 5. Download the libraries below from their sources:<br />
 (a) sqlite-autoconf-3330000.tar.gz from "https://www.sqlite.org/download.html".<br />
 (b) tiff-4.1.0.tar.gz from "https://download.osgeo.org/libtiff/"<br />
 (c) proj-7.2.0.tar.gz from "https://proj.org/download.html#current-release"<br />
 (d) libgeotiff-1.6.0.tar.gz from "https://download.osgeo.org/geotiff/libgeotiff/"<br />
 (e) geos-3.8.1.tar.bz2 from "https://trac.osgeo.org/geos"<br />
 (f) gdal-3.2.0.tar.gz from "https://gdal.org/download.html"<br />

</details>

<details>
<summary>REDCap errors</summary>
 1. Once your run completes check in your run directory if there is any *upload_log.tsv* file or *redcaperrors.txt*.  If these files exist, they may note any files or data which would have been over-written in the database.
 2. Check with the wet lab if any RD-numbers were duplicated or previously used for the samples listed in the upload_log.tsv file.
 3. Make sure your API token is not NULL and that REDCap is not down for maintenence here: https://redcap.nyumc.org/apps/redcap/<br>
 4. Check if any of the urls in the notification or API calls have been broken by a new version of REDCap. For Example, the link: https://redcap.nyumc.org/apps/redcap/redcap_v13.1.35/API/project_api.php?pid=24752 if broken, modify the URL to match REDCap Version i.e. /redcap/**redcap_v13.2.57**/) </br>
 Additional resources are here: https://redcap.nyumc.org/apps/redcap/index.php?action=help&newwin=1
</details>

<details>
<summary>REDCap Email Notification issues</summary>
The automatic email notifications are located on the left-side panel called "Alerts & Notifications".  If you need to change an output path in the email or change the year in the email, click on edit for Alert #1:Research Run Complete or Alert #2:Clinical Run Complete.<br>  
View the "Applications Overview" video here: https://redcap.nyumc.org/apps/redcap/index.php?action=training <br>
A detailed guide for Alerts is availible here: https://www.ctsi.ufl.edu/wordpress/files/2019/06/REDCap-Alerts-Notifications-User-Guide.pdf <br>
Additional resources are here: https://redcap.nyumc.org/apps/redcap/index.php?action=help&newwin=1 <br>
</details>

<details>
<summary>How to upload manually to REDCap</summary>
 1. Login with your kerberos ID to https://redcap.nyumc.org/<br />
 2. On the left-hand sidebar scroll all the way down the Reports Bookmarks until you see the folders:<br />
 `>>>>CURRENT Runs~~~~~ and 3) >>>>>CLINICAL Current Run`<br />
 3. Here, you can click on the RD-number of choice and then select "Upload html file" under the methylation menue<br />
 4. Optionally, you can also select "Add / Edit Records" menu in the left sidebar and find your RD-number in the "Search query" field<br />
 5. To upload the sample classifier details, such as the values and scores, a csv file named &lt;run_id&gt;_Redcap.csv is saved on the Desktop in a folder a run folder created named with the &lt;run_id&gt;. This file can be uploaded in the import tab of REDCap under *Data Import Tool* in the sidebar. The folder will also contain a &lt;run_id&gt;_samplesheet.csv file used in the run derived from the .xlsm file.<br /> 
 Additionally, this file is copied to: `/Volumes/CBioinformatics/Methylation/Clinical_Runs/csvRedcap/&lt;run_id&gt;/&lt;run_id&gt;_Redcap.csv`<br />
</details>

<details>
<summary>Issues with installing or running packages</summary>
 1. If you are getting compiler errors or all_installer.R fails, try installing additional system dependencies with brew and restart your R session: https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/Development/brewFix.sh<br />
 <br>
 2. If you still have errors with compiling or installing a package, try removing you MakeVars directory in: <br/>
 `rm -rf $HOME/.R/Makevars`<br/>
</details>

<details>
<summary>Fix wet lab worksheet</summary>
 1. In your run path /Volumes/CBioinformatics/Methylation/Clinical_Runs/22-MGDM##, open the RUNID.xlsm file<br />
 2. On the Review ribbon, click unprotect worksheet and unprotect tab<br />
 3. Right-click the "worksheet" tab at the bottom and unhide... raw_labels tab<br />
 4. If any "#ref" errors either drag the formula down to correct or type "=" and select the cell in the first tab "worksheet" and press return.<br />
 For example "=worksheet!B25" references cell B25 in the tab named "worksheet"
</details>
