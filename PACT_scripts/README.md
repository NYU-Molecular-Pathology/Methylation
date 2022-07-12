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

# 📑 Printing PACT Demultiplexing instructions

1. Save the shell script to a local directory with execute permissions, for example:
 ```ruby
curl -o ./PrintNGS.sh -L https://github.com/NYU-Molecular-Pathology/Methylation/PACT_scripts/PrintNGS.sh
chmod gu+rwx ./PrintNGS.sh
 ```
2. Next, execute the shell script to print the steps and commands.
  - PrintNGS.sh will take two parameters: the PACT run name and the RUNID. For example:
  ```ruby
  PrintNGS.sh 220907_NB501073_012345678_ABCDEFG1234 PACT-22-99
  ```
  - The script will download **demuxQC.sh** to your $HOME directory and will generate a text file named **"PACT-YY-##_stages.txt"** in your $HOME folder from the output of the parameters passed.  You can then cat PACT-YY-##_stages.txt to see all commands for each stage in the NGS pipeline.
## ⚡ Printing NGS Stage commands
 - Execute the following shell script passing the PACT run ID and run name:
 ```ruby
PrintNGS.sh 220101_NB501073_0123_ABCDEFG1234 PACT-22-XX
 ```
 - The Shell script is also availible on BigPurple here:
 ```ruby
 /gpfs/home/serraj10/molecpathlab/development/bash_scripts/PrintNGS.sh
 ```
# TroubleShooting Tips

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
