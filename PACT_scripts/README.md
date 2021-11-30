The shell script command bellow generate and email the samplesheet.csv and copy to BigPurple:<br>
`/Volumes/CBioinformatics/PACT/makePactSheet.sh PACT-21-## kerberosID`<br>

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

The API tokens are saved within the shell files
