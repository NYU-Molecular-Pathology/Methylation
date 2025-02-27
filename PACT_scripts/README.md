# Shell Scripts
### 📃 makePactSheet.sh (generates -samplesheet.csv only)
 - Executing the shell script command below generates -samplesheet.csv for the PACT run:
   + `/Volumes/CBioinformatics/PACT/makePactSheet.sh PACT-YY-##`
   + Here *YY* stands for the year, for example, 2022 will be 22, and *##* is an integer for the run from 01 to 99+
 - The csv file output is saved on the local Desktop and then copied to the BigPurple run folder
 - Once the script has created the sample sheet csv file, it is emailed using the REDCap API. 
 - Alternatively, you can also pass the path a .xlsm or .xlsx worksheet emailed from the wetlab. For example:
   + `makePactSheet.sh /Users/Jonathan/Desktop/PACT-22-15.xlsm`

### 🖼️ methylMatch.sh (Generates matched .xlsx of RD-numbers matching NGS case)
 - To generate a list of Methylation cases matched to the current pact run execute the following shell commmand:<br>
   + `/Volumes/CBioinformatics/PACT/methylMatch.sh PACT-YY-##`
 - methylMatch.sh will save the pact sheet PACT-YY-##_MethylMatch.xlsx to the local Desktop, then email it using the REDCap API
 - If any RD-numbers are found matching the MRNs on the PACT worksheet, the script will then generate the CNV png files for those RD-numbers, save them to the local Desktop and copy them to the Z: drive folder here:
```ruby
/molecular/Molecular/MethylationClassifier/CNV_PNG/
```
 - The Shell scripts should be saved locally or symlinked to make the script easier to execute without Network file path dependencies.

### ⚙️ parsePact.sh (Generates both csv and methylmatch .xslx)
 - **parsePact.sh** executes both makePactSheet.sh and methylMatch.sh
 - ***makePactSheet.sh*** and ***methylMatch.sh*** have one arg passed: The experiment name without quotes, for example `PACT-21-32`
 - Both curl calls can be run consecutively using parsePact.sh with the same argument parameter as makePactSheet.sh:
```ruby
/Volumes/CBioinformatics/PACT/parsepact.sh PACT-YY-##
```
The API tokens are saved within the shell files where $pactID is the experiment name arg1 input by the user.

The provided script `parsepact.sh` is designed to automate the process of parsing a PACT XLSM worksheet using R, followed by synchronizing the output to a high-performance computing (HPC) environment. Here is a detailed breakdown of the script:

## Script Metadata
- **Name:** `parsepact.sh`
- **Purpose:** Initiate an R script to parse PACT XLSM worksheet and sync the output to an HPC.
- **Date Created:** February 9, 2023
- **Author:** Jonathan Serrano
- **Version:** 1.2.0
- **Copyright:** NYULH Jonathan Serrano, 2024

## Hardcoded Variables
- **APITOKEN:** `'8XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'`
- **methAPI:** `'5XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'`
- **CSVOUTDIR:** `"/gpfs/data/molecpathlab/production/samplesheets/LG-PACT/"`
- **GITHUBLINK:** `"https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts/"`

## Input Arguments
- **PACTID:** `($1)` - Default: `NULL`
- **RUN_ID:** `($2)` - Default: `NULL`
- **VAL_KWD:** `($3)` - Default: `"-ILMNVAL"`
- **KERBEROS:** `($4)` - Default: `$USER`
- **MOLECAPI:** `($5)` - Default: `$API_DEFAULT`

## Variable Assignments
- **kerbero:** `$KERBEROS`
- **CURR_USER:** `$USER`

## R Version Check
Ensures the correct R version is used:
```bash
RMAINDIR=$(which R)
[ "$RMAINDIR" != "/usr/local/bin/R" ] && {
   echo -e "R not loaded from /usr/local/bin/R
Check ~/.Rprofile"
   exit 1
}
```

## Text Color Variables
Defines variables for colored and bold text output for better readability:
```bash
bold=$(tput bold)
BG_BLU="$(tput setab 4)"
BG_GRN="$(tput setab 2)"
FG_YLW="$(tput setaf 3)"
NORM=$(tput sgr0)
```

## Error Handling
Sets up a trap to notify if the script exits due to an error:
```bash
set -Eeuo pipefail

function notify {
   echo "Bash script exited!"
   echo "$(caller): ${BASH_COMMAND}"
}

trap notify ERR
```

## Utility Functions
- **`message_curl`**: Downloads a file using `curl` and changes its permissions.
- **`message_exe`**: Prints which command is executing.
- **`message_print`**: Prints a message in blue.
- **`check_directory`**: Checks if a directory exists on the remote HPC and creates it if it does not.

## Argument Validation
Checks if required arguments are provided, exits if missing:
```bash
[ -z "$PACTID" ] && {
   echo "Missing argument #1: You did not provide a PACTID name (i.e. PACT-22-12)"
   exit 1
}
[ -z "$RUN_ID" ] && {
   echo "Missing argument #2: PACT Run name (i.e. 250817_NB501073_0999_ABCD2TBGXM)"
   exit 1
}
```

## Main Script Execution

1. **Navigate to Home Folder and Print Input Arguments**
   ```bash
   cd "$HOME"
   message_print "Input Kerberos ID" "$kerbero"
   message_print "Input RUN_ID" "$RUN_ID"
   message_print "Input PACTID" "$PACTID"
   ```

2. **Download Latest R Scripts**
   ```bash
   message_curl ${GITHUBURL} "pactParse.R"
   message_curl ${GITHUBURL} "PactMethMatch.R"
   ```

3. **Execute R Script for Samplesheet Generation**
   ```bash
   message_exe 1 "pactParse.R" "$MOLECAPI $PACTID $RUN_ID"
   Rscript pactParse.R $MOLECAPI "$PACTID" "$RUN_ID" "$VAL_KWD"
   ```

4. **Define and Check CSV Files**
   ```bash
   NEWCSV="$HOME/Desktop/${RUN_ID}-SampleSheet.csv"
   PACTDIR="$CSVOUTDIR${RUN_ID}"
   SHEETDIR="$CURR_USER@bigpurple.nyumc.org:$PACTDIR"

   NEWCSV_VAL="$HOME/Desktop/${RUN_ID}-${VAL_KWD}-SampleSheet.csv"

   if [ -f "${NEWCSV_VAL}" ]; then
      HAS_VALIDATION=true
   else
      HAS_VALIDATION=false
   fi
   ```

5. **Check and Create Directories on HPC**
   ```bash
   check_directory "$CURR_USER" "$PACTDIR"
   if $HAS_VALIDATION; then
      PACTDIR_VAL="$CSVOUTDIR${RUN_ID}-${VAL_KWD}"
      SHEETDIR_VAL="$CURR_USER@bigpurple.nyumc.org:$PACTDIR_VAL"
      check_directory "$CURR_USER" "$PACTDIR_VAL"
   fi
   ```

6. **Synchronize CSV Files to HPC**
   ```bash
   message_print "Executing the following" "rsync -vrthP -e ssh $NEWCSV $SHEETDIR"
   rsync -vrthP -e ssh "$NEWCSV" "$SHEETDIR"
   if $HAS_VALIDATION; then
      rsync -vrthP -e ssh "$NEWCSV_VAL" "$SHEETDIR_VAL"
   fi
   ```

7. **Change Permissions on HPC**
   ```bash
   message_print "Changing permissions" "ssh $CURR_USER@bigpurple.nyumc.org chmod -R g+rwx $PACTDIR"
   ssh "$CURR_USER@bigpurple.nyumc.org" "chmod -R g+rwx $PACTDIR"
   if $HAS_VALIDATION; then
      ssh "$CURR_USER@bigpurple.nyumc.org" "chmod -R g+rwx $PACTDIR_VAL"
   fi
   ```

8. **Download and Execute Additional Scripts**
   ```bash
   message_curl ${GITHUBURL} "printPactCommands.sh"
   message_curl ${GITHUBURL} "make_consensus.sh"

   pactRunID=$(basename "${PACTID%.*}")

   message_print "If needed, modify Desktop samplesheet and rsync again" "rsync -vrthP -e ssh $NEWCSV $SHEETDIR"
   ```

9. **Generate HTML Commands**
   ```bash
   message_print "Saving Html File" "$HOME/printPactCommands.sh $RUN_ID ${pactRunID} NULL ${kerbero} >$HOME/${pactRunID}.html && open $HOME/${pactRunID}.html"

   "$HOME/printPactCommands.sh" "$RUN_ID" "${pactRunID}" NULL "${kerbero}" >"$HOME/${pactRunID}.html" && open "$HOME/${pactRunID}.html"
   if $HAS_VALIDATION; then
      "$HOME/printPactCommands.sh" "$RUN_ID-${VAL_KWD}" "${pactRunID}-${VAL_KWD}" NULL "${kerbero}" >"$HOME/${pactRunID}-${VAL_KWD}.html" && open "$HOME/${pactRunID}-${VAL_KWD}.html"
   fi

   open "$HOME/Desktop/${RUN_ID}-SampleSheet.csv"
   ```

This script is structured to ensure the correct execution of R scripts, synchronization of results to an HPC, and appropriate notification and handling of any errors that may occur.


# 📑 Printing PACT Demultiplexing instructions

1. Save the shell script to a local directory with execute permissions, for example:
 ```ruby
curl -L https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/PACT_scripts/PrintNGS.sh -o $HOME/PrintNGS.sh | chmod +rwx $HOME/PrintNGS.sh
 ```
2. Next, add an alias to your bash or zsh to execute the script from anywhere:
  ```ruby
  alias printngs=$HOME/PrintNGS.sh
  echo -e "\nalias printngs='$HOME/PrintNGS.sh'" >> .zshrc
  ```
  - PrintNGS.sh will take two parameters: the PACT run name and the RUNID. For example:
  ```ruby
  printngs 220907_NB501073_012345678_ABCDEFG1234 PACT-22-99
  ```
  - The script will download **demuxQC.sh** to your $HOME directory and will generate a text file named **"PACT-YY-##_stages.txt"** in your $HOME folder from the output of the parameters passed.  You can then cat PACT-YY-##_stages.txt to see all commands for each stage in the NGS pipeline.
## ⚡ Printing NGS Stage commands
 - Execute the following shell script passing the PACT run ID and run name:
 ```ruby
./PrintNGS.sh 220101_NB501073_0123_ABCDEFG1234 PACT-22-XX
 ```
 - The Shell script is also availible on BigPurple here:
 ```ruby
 "$HOME/molecpathlab/development/bash_scripts/PrintNGS.sh"
 ```
# TroubleShooting Tips

<details>
  <summary>PrintNGS.sh not found</summary>
Ensure you have permissions to read, write, and execute the directory where you are using curl to download the script as well as add execute permissions to the aliased bash script.
</details>

<details>
  <summary>Unknown Font Error</summary>
If you have an unknown system font error, make sure you have XQuartz installed for Mac and MacTex (LaTeX) installed<br>
You may need cairo if not installed which can be installed using `brew install cairo`
</details>

<details>
  <summary>Aborted/Fatal error during CNV PNG creation</summary>
If you encounter a system abort or memory allocation error, it may just be a timeout issue with generating the segments for a complex CNV plot<br>
The cnv script works by generating the CNV ggplot from the mnpv11b6 package and then saving it as a ggplot widget and rendering the html as a png file.<br>
</details>


<details>
  <summary>System Memory or output error</summary>
Check if you have an existing temp file in ~/Desktop/temp.html and delete it.  The temp file is the ggplotly object saved as html which is then opened in a headless chromium browser and then webshot into a PNG file.  It may be because the temp.html file could not be overwritten <br>


For example, where the following variables:<br>
`
sampleName = "RD-21-322"
tempPathFi = "~/Desktop/temp.html"
fn = "~/Desktop/RD-21-322_cnv.png"
sex = "male"`

Are being passed into the function gen.cnv.png from<br>
https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/PACT_scripts/generateCNV.R<br>
After lines 96:

```R
  xx <- mnp.v11b6::MNPcnv(Mset,sex = sex, main = sampleID) # Get cnv object
  thePlot <- mnp.v11b6::MNPcnvggplotly(xx, getTables = F) # saves cnv ggplot
  p <- plotly::ggplotly(thePlot)                          # converts plot to ggplotly
  htmlwidgets::saveWidget(                                # saves as widget in temp.html
       widget=plotly::as.widget(p), 
       file=tempPathFi)
  # Below saves a screenshot of the temp.html file as a PNG image
  webshot2::webshot(url=tempPathFi, file = fn, cliprect = "viewport", delay = 2.5, vwidth = 2340, vheight = 1344)
```
</details>

<details>
  <summary>File Copy or Mount Path issue</summary>
 - The mounted drive paths are checked to idat folders in Molecular and Snuderlabspace shares so that they are accessible.
 - You may see in red if the path to idat folder in /Volumes/snudem01labspace is not mounted; however, if your idats are all in the molecular drive, the script will continue if not issues occured.
 - If any of the cnv png output files on the desktop are not copied to the output folder, check that the files do not already exist as copying will be skipped or if the network molecular drive is not accessible, the png files will remain on the desktop.
</details>
