# Shell Scripts
## makePactSheet.sh
The shell script command bellow generate and email the samplesheet.csv and copy to BigPurple:<br>
`/Volumes/CBioinformatics/PACT/makePactSheet.sh PACT-21-## kerberosID`<br>

## methylMatch.sh
To generate the PACT-21-##_MethylMatch.xlsx email and save CNV png files execute:<br>
`/Volumes/CBioinformatics/PACT/methylMatch.sh PACT-21-##`

makePactSheet.sh and methylMatch.sh can be saved locally or symlinked to make the script easier to execute.

**methylMatch.sh has one arg**
1. The Experiment Name without quotes, for example PACT-21-32

**makePactSheet.sh has two args**
1. The Experiment Name without quotes, for example PACT-21-99
2. Your kerberos id is the second arg

Both curl calls can be run consecutively using parsePact.sh with the same two arguments as makePactSheet.sh:

`curl -o pactParse.R -L https://git.io/J0kfR -s; Rscript pactParse.R $APITOKEN "$pactID"`<br>
`curl -o PactMethMatch.R -L https://git.io/JEp7l -s; Rscript PactMethMatch.R $methAPI "$pactID"`

The API tokens are saved within the shell files and $pactID is the experiment name arg1 input by the user.

# TroubleShooting Tips

## Unknown Font Error
If you have an unknown system font error, make sure you have XQuartz installed for Mac and MacTex (LaTeX) installed<br>
You may need cairo if not installed which can be installed using `brew install cairo`

## Aborted/fatal error during cnv png creation
If you encounter a system abort or memory allocation error, it may just be a timeout issue with generating the segments for a complex CNV plot<br>
The cnv script works by generating the CNV ggplot from the mnpv11b6 package and then saving it as a ggplot widget and rendering the html as a png file.<br>



