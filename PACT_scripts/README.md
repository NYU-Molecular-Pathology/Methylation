# Shell Scripts
### makePactSheet.sh (generates -samplesheet.csv only)
 - Executing the shell script command below generates -samplesheet.csv for the PACT run:
`/Volumes/CBioinformatics/PACT/makePactSheet.sh PACT-YY-##`<br>
 - Here *YY* stands for the year, for example, 2022 will be 22, and *##* is an integer for the run from 01 to 99+
 - The csv file output is saved on the local Desktop and then copied to the BigPurple run folder
 - Once the script has created the sample sheet csv file, it is emailed using the REDCap API. 

Alternatively, you can also pass the path to the .xlsm worksheet emailed from the wetlab.

### methylMatch.sh (Generates matched .xlsx of RD-numbers matching NGS case)
 - To generate a list of Methylation cases matched to the current pact run execute the following shell commmand:<br>
`/Volumes/CBioinformatics/PACT/methylMatch.sh PACT-YY-##`
 - methylMatch.sh will save the pact run PACT-YY-##_MethylMatch.xlsx email and save CNV png files execute:

makePactSheet.sh and methylMatch.sh can be saved locally or symlinked to make the script easier to execute.

### parsePact.sh (Generates both csv and methylmatch .xslx)
 - **parsePact.sh** executes both makePactSheet.sh and methylMatch.sh
 - ***makePactSheet.sh*** and ***methylMatch.sh*** have one arg passed: The experiment name without quotes, for example `PACT-21-32`
 - Both curl calls can be run consecutively using parsePact.sh with the same argument parameter as makePactSheet.sh:
```ruby
/Volumes/CBioinformatics/PACT/parsepact.sh PACT-YY-##
```
The API tokens are saved within the shell files where $pactID is the experiment name arg1 input by the user.

# Printing PACT Demultiplexing instructions

1. First copy the shell script to a local directory, for example:
 ```ruby
 cp /Volumes/CBioinformatics/PACT/demuxQC.sh ~/
 ```
2. Next, execute the shell script and save it as a text file.
  - demuxQC.sh will take two parameters: the PACT run name and the RUNID
  - For example:
  ```ruby
  demuxQC.sh 220907_NB501073_012345678_ABCDEFG1234 PACT-22-99 > PACT-22-99.txt
  ```
3. Now that you have the output saved locally as a text file you can reference it to copy and paste the commands for each stage.
  - For example, to view the first two stages, print the first 27 (n) lines:
  ```ruby
  head -n 27 PACT-22-99.txt
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
