# How To Setup Clinical Methylation Classifier

First, download RStudio 4.1 from:
[https://download1.rstudio.org/desktop/macos/RStudio-1.4.1717.dmg](https://download1.rstudio.org/desktop/macos/RStudio-1.4.1717.dmg)

- The classifier runs on R version 3.6.3 and up.  It is not compatible with 3.3.3  
- Classifier Package install times are below:
<br> R 3.6 approximately 36.4 mintues
<br> R 4.0 approximately 15.2 mintues
<br> R 4.1 approximately 6.33 mintues

Download the following Mac Packages, but **do not install** until you have opened the system preferences Security Pane as below.

I recommend to download 4.1 and run on a clean install to avoid changing any other packages.  4.1 includes compile and Tckl dependencies.  Attached in this directory are menubar apps for Rswitch which allows you to seamlessly switch between package libraries:
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/25020451c5d320ff6aa70e65421f7a828f4a6905/screenshots/rswitch.png" alt="drawing" width="250"/>

Download R from Cran: [https://cran.r-project.org/bin/macosx/](https://cran.r-project.org/bin/macosx/)

Download XQuartz with any new install of R or Rstudio you have:
[https://github.com/XQuartz/XQuartz/releases/download/XQuartz-2.8.1/XQuartz-2.8.1.dmg](https://github.com/XQuartz/XQuartz/releases/download/XQuartz-2.8.1/XQuartz-2.8.1.dmg)

Download LaTeX for Mac in order to be able to properly render the HTML report outputs
[https://www.tug.org/mactex/mactex-download.html](https://www.tug.org/mactex/mactex-download.html)


Next, go to Mac OS System Preferences by pressing **⌘(CMD) + Space** and typing "System Pref" and click on Security & Privacy:<br>
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/0e124ca5b0b06278ccfd03569e0b8cd769e9fd2b/screenshots/syspref1.png" alt="drawing" width="650"/>

<br>
### 1. Click on the "General" Tab and click the lock to make changes
<br>
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/0e124ca5b0b06278ccfd03569e0b8cd769e9fd2b/screenshots/syspref2.png" alt="drawing" width="650"/>

### 2. Keep this window open as you install RStudio and R as Mac OS will block any app install.
<br>
### 3. Once you install R and Rstudio, click the "Privacy" tab in system preferences
<br>
### 4. Select "Full Disk Access", and add R and Rstudio to prevent any slow file reads or permission issues.
<br>
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/0e124ca5b0b06278ccfd03569e0b8cd769e9fd2b/screenshots/syspref3.png" alt="drawing" width="650"/>

### To avoid any install or compile errors, you must have Xcode installed and accept the licence.  You can run the following commands in R or paste into terminal the content in quotes:

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

Next press **⌘(CMD) + Shift + G** or click the Finder Menubar menu GO > Go to Folder... in the Finder menubar and paste: 

<br>`/Volumes/CBioinformatics/jonathan/Rprojects/Clinical_Methylation/Clinical_Methylation_Run.Rmd`

`/Volumes/CBioinformatics/Methylation/load_install_test/LoadInstall.R`

If you don't have homebrew, install brew by running the following command in R:

`system("/bin/bash -c $(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)`

Next run the LoadInstall.R to install all the package dependencies.

The latest version of brew installs any casks safely in $USER/local and symlinks.  It is reversable as it does not overwrite an existing components nor does it overwrite any Mac OS System Components such as Clang which are part of Xcode. To stop a brew cask from loading simply execute:
<br>
`brew unlink [packageName]`

