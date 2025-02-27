# How to Selectively run Samples to Generate Reports

## Generate Classifier Reports for a list of RD-numbers or Single Sample
The [methylExpress_custom.R](https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/methylExpress_custom.R) script allows 
you to input RD-number(s) and output classifier report(s) for only those cases regardless of which run they are on. This script will:<br>
1. Parse RD-numbers as input<br>
2. Generate a "samplesheet.csv" file in the working directory using a REDCap API token<br>
3. Copy the .idat files ma to the working directory<br>
4. Generate either the Sarcoma or Brain Classifier reports for the RD-numbers input.<br>
5. Optionally, the output can be uploaded to REDCap.<br>

## Input Parameters
methylExpress_custom.R has two parameters that can be used to input RD-numbers: `rd_csv` or `rd_numbers`<br>
<br>
`rd_numbers`: *a character list of one more more RD-numbers.*<br>
Examples: `rd_numbers = c("RD-22-123", "RD-22-456", "RD-22-789")` OR `rd_numbers = c("RD-22-123")`<br>
<br>
`rd_csv`: *a file path to a comma separated values file (.csv) containing **two** or more RD-numbers.*<br>
Example: `rd_csv = "/Users/MyName/Desktop/my_rdnumbers.csv"`<br>

### NOTE
 - The CSV file should have **no** Header and the rd-numbers should be listed in the first column only.<br>
 - Do not use the csv file to input a single RD-number. You must have at least two rows of data in the .csv file column for it to parse as a dataframe. Use the *rd_numbers* variable to input a single RD-number instead.
 - You can use either `rd_numbers` or `rd_csv` to input RD-numbers. One of the variables must be NULL.<br>
 For example, if you use `rd_csv = "/Users/MyName/Desktop/my_rdnumbers.csv"`, then `rd_numbers = NULL`<br>

## Run Flags
methylExpress_custom.R has three boolean flags that can be changed depending on the desired output:<br>
<br>
`runSarcoma`: If TRUE, the html output will be from the Sarcoma Classifier and not Brain Classifier<br>
<br>
`makeNewSheet`: If TRUE, a new minfi samplesheet is created using the input parameters.<br>
**NOTE**: *Any existing "samplesheet.csv" in the current directory will be deleted and replaced if TRUE.*<br>
<br>
`forcedUpload`: If TRUE, **all** .html files in the current directory are imported into REDCap, replacing the previous report(s).<br>

## Running methylExpress_custom.R
You can download the script by executing the command below:
```
curl -k -# -L "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/methylExpress_custom.R" >"$HOME/methylExpress_custom.R"
```
After downloading, be sure to replace the placeholder ("XXXX") with your REDCap API token.<br>
<br>
To run methylExpress_custom.R with the API token and runPath pre-configured, you can use the file in the Shared Drive:<br>
`/Volumes/CBioinformatics/Methylation/Clinical_Runs/Custom_ReRun/methylExpress_custom.R`

## Demo
You can test the script by:
1. Be sure you have all the classifier [requirements](https://github.com/NYU-Molecular-Pathology/Methylation/blob/5696d877690d3165fcb489f53fd53fa023214058/README.md) and dependencies installed first with [all_installer.R](https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/all_installer.R)
2. Download [methylExpress_custom.R](https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/R/methylExpress_custom.R), and open it in Rstudio, setting the [REDCap API token](https://redcap.nyumc.org/apps/redcap/api/help/?content=tokens).
3. Set `rd_numbers = NULL`
4. Download the demo csv file here: https://github.com/NYU-Molecular-Pathology/Methylation/blob/main/R/CustomRuns/Demo/demo_input.csv
5. Assign `rd_csv` to the path where you saved demo_input.csv (i.e., `rd_csv = "/path/to/downloads/folder/demo_input.csv"`)
6. Modify `runFolder` if you want to copy idats and run in a different directory, then execute the script.

