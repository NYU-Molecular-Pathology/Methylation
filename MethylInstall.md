# How To Setup Clinical Methylation Classifier

First, download RStudio 4.1 from:
[https://download1.rstudio.org/desktop/macos/RStudio-1.4.1717.dmg](https://download1.rstudio.org/desktop/macos/RStudio-1.4.1717.dmg)

- The classifier runs on R version 3.6.3 and up.  It is not compatible with 3.3.3  
- Classifier Package install times are below:
<br> R 3.6 approximately 36.4 mintues
<br> R 4.0 approximately 15.2 mintues
<br> R 4.1 approximately 6.33 mintues

Install R from Cran: [https://cran.r-project.org/bin/macosx/](https://cran.r-project.org/bin/macosx/)

Download XQuartz with any new install of R or Rstudio you have:
[https://github.com/XQuartz/XQuartz/releases/download/XQuartz-2.8.1/XQuartz-2.8.1.dmg](https://github.com/XQuartz/XQuartz/releases/download/XQuartz-2.8.1/XQuartz-2.8.1.dmg)

Next, go to Mac OS System Preferences by pressing **⌘(CMD) + Space** and typing "System Pref" and click on Security & Privacy:<br>
<img src="/Users/serraj10/Desktop/syspref1.png" alt="drawing" width="650"/>

<br>
1. Click on the "General" Tab and click the lock to make changes
<br>
<img src="/Users/serraj10/Desktop/syspref2.png" alt="drawing" width="650"/>

2. Keep this window open as you install RStudio and R as Mac OS will block any app install.
<br>
3. Once you install R and Rstudio, click the "Privacy" tab in system preferences
<br>
4. Select "Full Disk Access", and add R and Rstudio to prevent any slow file reads or permission issues.
<br>
<img src="/Users/serraj10/Desktop/syspref3.png" alt="drawing" width="650"/>


If you don't have homebrew, install brew by running the following command in R:

`system("/bin/bash -c $(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)`

The latest version of brew installs any casks safely in $USER/local and symlinks.  It is reversable as it does not overwrite an existing components nor does it overwrite any Mac OS System Components such as Clang which are part of Xcode.

To avoid any install or compile errors, you must have Xcode installed and accept the licence.  You can run the following commands in terminal or in R:

`system("xcode-select --install")`

`system("xcodebuild -runFirstLaunch")`

## Getting the Source Code

To install and run the pipeline, you must mount two drives:
<br><br>
Open Finder and press **⌘(CMD) + K** then add three directories, login name is NYUMC\KerberosID:
<br><br>
`smb://research-cifs.nyumc.org/Research/CBioinformatics/`

`smb://research-cifs.nyumc.org/Research/snudem01lab/snudem01labspace`

`smb://shares-cifs.nyumc.org/apps/acc_pathology/molecular`

Next press **⌘(CMD) + Shift + G** or click GO > Go to Folder... in the Finder menubar and paste: 
<br>`/Volumes/CBioinformatics/jonathan/Rprojects/Clinical_Methylation/Clinical_Methylation_Run.Rmd`

`/Volumes/CBioinformatics/Methylation/load_install_test/LoadInstall.R`

