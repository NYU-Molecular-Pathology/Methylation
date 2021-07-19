# How To Setup Clinical Methylation Classifier 
## Essential Downloads
The classifier runs on R version 3.6.3 and up.  It is not compatible with 3.3.3<br />
Download **R 4.1** from CRAN: [https://cran.r-project.org/bin/macosx/](https://cran.r-project.org/bin/macosx/)<br />
Download **RStudio 1.4**: [https://download1.rstudio.org/desktop/macos/RStudio-1.4.1717.dmg](https://download1.rstudio.org/desktop/macos/RStudio-1.4.1717.dmg)<br />
Download **XQuartz**: [https://github.com/XQuartz/XQuartz/releases/download/XQuartz-2.8.1/XQuartz-2.8.1.dmg](https://github.com/XQuartz/XQuartz/releases/download/XQuartz-2.8.1/XQuartz-2.8.1.dmg)<br />
Download **LaTeX** for Mac: [https://www.tug.org/mactex/mactex-download.html](https://www.tug.org/mactex/mactex-download.html)<br />
Download **Homebrew**: https://brew.sh/ you can install using the following line in terminal:<br />
`/bin/bash -c $(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh`

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
### There are two system Rscript to run methylExpress.R with the arguments in order:<br />
`arg[1]` is the **token** for the API call ('#######################')<br />
`arg[2]` is the **RunID** which if NULL runs the latest Clinical Worksheet ('MR21-099')<br />
`arg[3]` is the **selectRds** parameter which is to prioritize samples being run (NULL)<br />

Alternatively, you can source then run the github script locally using [methylExpress.R](https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/methylExpress.R)

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
