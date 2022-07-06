## Pipeline Functions Overview
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/01d70d5ba1c9dd0e090cf02c225a5a909c78317e/screenshots/meth_pipeline_uml.png" width="100%"/><br/>

### First function is LoadInstall_new.R which lists all the dependencies and required packages 
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/421ecaf2793e8260d83fba35ce6d130e58fc3f0c/screenshots/loadinstall.png" alt="drawing" width="450"/><br />

### Next, the functions inside SetRunParams.R load the default file paths and parameter values
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/421ecaf2793e8260d83fba35ce6d130e58fc3f0c/screenshots/setparam.png" alt="drawing" width="850"/><br />

### The functions inside CopyInputs.R help to setup the methylation run by copying idats and parsing the worksheet to the working directory
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/421ecaf2793e8260d83fba35ce6d130e58fc3f0c/screenshots/copyin.png" alt="drawing" width="850"/><br />

### The functions inside CopyOutput.R copy the html reports to the shared drives and copy them to REDCap also
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/be9bea1e0ee763d964aa18caab8f444369dc1030/screenshots/copyout.png" alt="drawing" width="850"/><br />

### The functions inside pipelineHelper.R generate the reports and save the html to the current directory and csv files on the desktop
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/a880a5bcba2f9c6242286116715f6a2e10c0440a/screenshots/pipeline.png" alt="drawing" width="850"/><br />

### The main script that is run from the command line is methylExpress.R which sources each of these scripts into the global environment and takes in the user parameters
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/5e54893b0cc7d54212508f78d3239e787cf09c26/screenshots/methEx.png" alt="drawing" width="850"/><br />
