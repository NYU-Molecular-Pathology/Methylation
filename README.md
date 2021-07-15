# How To Setup Clinical Methylation Classifier ---

### Essential Downloads
The classifier runs on R version 3.6.3 and up.  It is not compatible with 3.3.3

- Download **R 4.1** from CRAN: [https://cran.r-project.org/bin/macosx/](https://cran.r-project.org/bin/macosx/)
- Download **RStudio 1.4**: [https://download1.rstudio.org/desktop/macos/RStudio-1.4.1717.dmg](https://download1.rstudio.org/desktop/macos/RStudio-1.4.1717.dmg)
- Download **XQuartz**: [https://github.com/XQuartz/XQuartz/releases/download/XQuartz-2.8.1/XQuartz-2.8.1.dmg](https://github.com/XQuartz/XQuartz/releases/download/XQuartz-2.8.1/XQuartz-2.8.1.dmg)
- Download **LaTeX** for Mac: [https://www.tug.org/mactex/mactex-download.html](https://www.tug.org/mactex/mactex-download.html)
- Download **Homebrew**, you can install using the following line in terminal:
`/bin/bash -c $(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh`
- R v4.1 includes compile and Tckl dependencies. brew can install libomp and cairo if needed.
- After downloading R & RStudio **do not install** until you have unlocked the
[System Preferences Privacy & Security Panel](https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/Notes/SystemPermissions.md).
- Download the [Rswitch](https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/Rswitch/) app, which allows you to switch between R version libraries, for example:
`/Library/Frameworks/R.framework/Versions/3.6/Resources/library/minfi`
`/Library/Frameworks/R.framework/Versions/4.1/Resources/library/minfi`

## Network Drive Mount Paths

---

### To install & run the pipeline, it is critical to mount the following network smb shared drives:

Open Finder and press **⌘(CMD) + K** then paste each of these directories, login name and password is your NYUMC\KerberosID
<br>
`smb://research-cifs.nyumc.org/Research/CBioinformatics/`<br />
`smb://research-cifs.nyumc.org/Research/snudem01lab/snudem01labspace`<br />
`smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular`<br />

## Install and run the pipeline

---

Run [LoadInstall_new.R](https://github.com/NYU-Molecular-Pathology/Methylation/blob/8b32e9a3c90ebf6b568a6c6582a020f6729867ef/LoadInstall_new.R) to install all package dependencies and the pipeline on a new system.  This script is also called when running the pipeline.

## Start a Run in Terminal

### To run the pipeline from your terminal, simply execute the following command:<br />
`[ ! -f methylExpress.R ] && wget -L https://git.io/JWujj; Rscript --verbose methylExpress.R '12456789abcdefghijklmnop''MR21-099' NULL`<br />

### Args
There are two system Rscript to run methylExpress.R with the arguments in order:<br />

`arg[1]` is the token for the API call ('12456789abcdefghijklmnop')<br />
`arg[2]` is the RunID which if NULL runs the latest Clinical Worksheet ('MR21-099')<br />
`arg[3]` is the selectRds parameter which is to prioritize samples being run (NULL)<br />

Alternatively, you can source then run the github script:
