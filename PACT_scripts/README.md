# Shell Scripts
## makePactSheet.sh
The shell script command bellow generate and email the samplesheet.csv and copy to BigPurple:<br>
`/Volumes/CBioinformatics/PACT/makePactSheet.sh PACT-YY-##`<br>

## methylMatch.sh
To generate the PACT-21-##_MethylMatch.xlsx email and save CNV png files execute:<br>
`/Volumes/CBioinformatics/PACT/methylMatch.sh PACT-22-##`

makePactSheet.sh and methylMatch.sh can be saved locally or symlinked to make the script easier to execute.

**methylMatch.sh has one arg**
1. The Experiment Name without quotes, for example PACT-21-32

**makePactSheet.sh has two args**
1. The Experiment Name without quotes, for example PACT-21-99
2. Your kerberos id is the second arg

## parsePact.sh (both makePactSheet.sh and methylMatch.sh)
Both curl calls can be run consecutively using parsePact.sh with the same two arguments as makePactSheet.sh:

`curl -o pactParse.R -L https://git.io/J0kfR -s; Rscript pactParse.R $APITOKEN "$pactID"`<br>
`curl -o PactMethMatch.R -L https://git.io/JEp7l -s; Rscript PactMethMatch.R $methAPI "$pactID"`<br>

The API tokens are saved within the shell files and $pactID is the experiment name arg1 input by the user.

# TroubleShooting Tips

## Unknown Font Error
If you have an unknown system font error, make sure you have XQuartz installed for Mac and MacTex (LaTeX) installed<br>
You may need cairo if not installed which can be installed using `brew install cairo`

## Aborted/fatal error during cnv png creation
If you encounter a system abort or memory allocation error, it may just be a timeout issue with generating the segments for a complex CNV plot<br>
The cnv script works by generating the CNV ggplot from the mnpv11b6 package and then saving it as a ggplot widget and rendering the html as a png file.<br>

## System Memory or output error
Check if you have an existing temp file in ~/Desktop/temp.html and delete it.  The temp file is the ggplotly object saved as html which is then opened in a headless chromium browser and then webshot into a PNG file.  It may be because the temp.html file could not be overwritten
<br>

For example, where the following variables:<br>
`
sampleName = "RD-21-322"
tempPathFi = "~/Desktop/temp.html"
fn = "~/Desktop/RD-21-322_cnv.png"
sex = "male"`

Are being passed into the function gen.cnv.png from<br>
https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/PACT_scripts/generateCNV.R<br>
After lines 96:
<br>
`message("~~~~~~~~~~~~~~~Generating ", sampleName, " cnv plot...")`<br>
`xx <- mnp.v11b6::MNPcnv(Mset,sex = sex,main = sampleID)` # Get cnv object<br>
`thePlot<-supM(mnp.v11b6::MNPcnvggplotly(xx, getTables = F))` # saves cnv ggplot<br>
`p<-plotly::ggplotly(thePlot)` # converts plot to ggplotly<br>
`htmlwidgets::saveWidget(widget=plotly::as.widget(p), file=tempPathFi)` # saves as widget in temp.html<br>

`# Below saves a screenshot of the temp.html file as a PNG image`<br>

`webshot2::webshot(url=tempPathFi, file = fn, cliprect = "viewport", delay = 2.5, vwidth = 2340, vheight = 1344)`<br>              

## File Copy or Mount Path issue
The mounted drive paths are checked to idat folders in Molecular and Snuderlabspace shares so that they are accessible.  You may see in red if the path to idat folder in /Volumes/snudem01labspace is not mounted; however, if your idats are all in the molecular drive, the script will continue if not issues occured. <br>
If any of the cnv png output files on the desktop are not copied to the output folder, check that the files do not already exist as copying will be skipped or if the network molecular drive is not accessible, the png files will remain on the desktop.
