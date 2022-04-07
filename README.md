# How To Setup Clinical Methylation Classifier 
## Essential Downloads
The classifier runs on R version 3.6.3 and up.  It is not compatible with 3.3.3<br />
Download **R 4.1** from CRAN: [https://cran.r-project.org/bin/macosx/](https://cran.r-project.org/bin/macosx/)<br />
Download **RStudio 1.4**: [https://download1.rstudio.org/desktop/macos/RStudio-1.4.1717.dmg](https://download1.rstudio.org/desktop/macos/RStudio-1.4.1717.dmg)<br />
Download **XQuartz**: [https://github.com/XQuartz/XQuartz/releases/download/XQuartz-2.8.1/XQuartz-2.8.1.dmg](https://github.com/XQuartz/XQuartz/releases/download/XQuartz-2.8.1/XQuartz-2.8.1.dmg)<br />
Download **LaTeX** for Mac: [https://www.tug.org/mactex/mactex-download.html](https://www.tug.org/mactex/mactex-download.html)<br />
Download **Homebrew**: https://brew.sh/ you can install using the following line in terminal:<br />
`/bin/bash -c $(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh`<br />
Download and Install XCode 12.0 or higher for Mac OS: [https://developer.apple.com/download/all/?q=xcode](https://developer.apple.com/download/all/?q=xcode)<br>
Install Java 8 JDK:
https://www.oracle.com/java/technologies/downloads/#java8-mac

- R v4.1 includes compile and Tckl dependencies. brew can install libomp and cairo if needed.
- After downloading R & RStudio **do not install** until you have unlocked the
[System Preferences Privacy & Security Panel](https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/Notes/SystemPermissions.md).
- Download the [Rswitch](https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/Rswitch/) app, which allows you to switch between R version libraries, for example:
`/Library/Frameworks/R.framework/Versions/3.6/Resources/library/minfi`
`/Library/Frameworks/R.framework/Versions/4.1/Resources/library/minfi`
___
## Network Drive Mount Paths
### To install & run the pipeline, it is critical to mount the following network smb shared drives:

Open Finder and press **⌘(CMD) + K** then paste each of these directories, login name and password is your NYUMC\KerberosID
<br>
`smb://research-cifs.nyumc.org/Research/CBioinformatics/`<br />
`smb://research-cifs.nyumc.org/Research/snudem01lab/snudem01labspace`<br />
`smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular`<br />
___
## Install pipeline and start a Run in Terminal
### To run the pipeline from your terminal, simply execute the following command:<br />
`wget -L https://git.io/JWujj -O methylExpress.R; Rscript --verbose methylExpress.R '#######################' '21-MGDM_TEST' NULL`<br />

If you have issues with wget, you can also use curl to run the script:

`curl -o methylExpress.R -L https://git.io/JWujj; Rscript --verbose methylExpress.R '#######################' '21-MGDM_TEST' NULL`<br />

### There are two system Rscript to run methylExpress.R with the arguments in order:<br />
`arg[1]` is the **token** for the API call ('#######################')<br />
`arg[2]` is the **RunID** which if NULL runs the latest Clinical Worksheet ('MR21-099')<br />
`arg[3]` is the **selectRds** parameter which is to prioritize samples being run (NULL)<br />
`arg[4]` is the **baseFolder** parameter which is optional if you want to run/save output to a different directory (NULL)<br />

Alternatively, you can source then run the github script locally using [methylExpress.R](https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/methylExpress.R)

# Email List of Methylation samples which are PACT run

`curl -o PactMethMatch.R -L https://git.io/J41Wp; Rscript --verbose PactMethMatch.R 'TOKENAPI12345667891011' '/Users/PATH/TO/CSV/Desktop/210715_NB501073_9999_ABCDEFGHIJRLK-SampleSheet.csv'`

# Email PACT csv file:
`curl -o pactParse.R -L https://git.io/J0kfR; Rscript --verbose pactParse.R 'TOKENAPI12345667891011' 'PACT-21-##'`

## Pipeline Functions Overview
### First function is LoadInstall_new.R which lists all the dependencies and required packages 
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/421ecaf2793e8260d83fba35ce6d130e58fc3f0c/screenshots/loadinstall.png" alt="drawing" width="450"/><br />

### Next, the functions inside SetRunParams.R load the default file paths and parameter values
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/421ecaf2793e8260d83fba35ce6d130e58fc3f0c/screenshots/setparam.png" alt="drawing" width="850"/><br />

### The functions inside CopyInputs.R help to setup the methylation run by copying idats and parsing the worksheet to the working directory
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/421ecaf2793e8260d83fba35ce6d130e58fc3f0c/screenshots/copyin.png" alt="drawing" width="850"/><br />

### The functions inside CopyOutput.R copy the html reports to the shared drives and copy them to REDCap also
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/be9bea1e0ee763d964aa18caab8f444369dc1030/screenshots/copyout.png" alt="drawing" width="850"/><br />

### The functions inside pipelineHelper.R generate the reports and save the html to the current directory and csv files on the desktop
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/a880a5bcba2f9c6242286116715f6a2e10c0440a/screenshots/pipeline.png" alt="drawing" width="850"/><br />

### The main script that is run from the command line is methylExpress.R which sources each of these scripts into the global environment and takes in the user parameters
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/5e54893b0cc7d54212508f78d3239e787cf09c26/screenshots/methEx.png" alt="drawing" width="850"/><br />

## **RUNNING METHYLATION CLI**
To run the Clinical or Research Methylation pipeline, simply use the locally stored Shell Script in:<br>
`/Volumes/CBioinformatics/Methylation/runMeth.sh`
The shell script takes the following argument parameters:<br>
`#!/bin/bash`

`methAPI='XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' # RedCap API Token`
`methRun=${1-NULL} # if arg $1 is empty assign NULL as default`
`PRIORITY=${2-NULL} # if arg $2 is empty assign NULL as default`
`runPath=${3-NULL} # if arg $3 is empty assign NULL as default`
`BIuser="$(whoami)"`

`# Print helpFunction in case methRun parameter is empty`
`if [ -z "$methRun" ]`
`then`
`   echo "You did not provide a Methylation Run ID name. Ex. '21-MGDM30'";`
`   exit 1`
`fi`

`# Begin script in case all parameters are entered`
`curl -o methylExpress.R -L https://git.io/JWujj; Rscript --verbose methylExpress.R $methAPI $methRun $PRIORITY $runPath`<br>

You can locally copy or symlink the runMeth.sh file to execute more easily
