#!/usr/bin/env python
# coding=utf-8
from __future__ import unicode_literals

import json
import logging
import os
import subprocess

# Log file path and output Config
logging.basicConfig(filename='./bsapp.log', filemode='w', format='%(name)s - %(levelname)s - %(message)s', level=logging.INFO)

# Default Parameters
outPath = "/mnt/instruments/illumina/"
keyName = "PACT"  # or MGFS
fiFilter = "bin,xml,csv"
fieldFtr = "ExperimentName"
days = 7  # cutoff of days to download
bsPath = "./bs"  # ./bs default /usr/local/bin/bs local

# Commands called: load config, list runs as JSON format, download files
cmd1 = f'eval $({bsPath} load config default.cfg)' 
cmd2 = f'{bsPath} list runs --newer-than={days}d --filter-field={fieldFtr} --filter-term={keyName} -f json'
cmd3 = f"{bsPath} download run --extension={fiFilter} -i"  # --overwrite

# Call Commands
os.system(cmd1)
allExp = json.loads(subprocess.check_output(cmd2, shell=True))

def download_bs_info(allExp):
    """_summary_
    For each PACT experiment name, if run status is Complete AND folder does not exist in outDir, mkdir and download files

    Args:
        allExp (list): json string list of BaseSpace CLI filtered shell return
    """
    currCmd = ""
    isDone = [ex['Status'] for ex in allExp] # 'Complete' or 'Running' or 'Failed'
    runIds = [ex['Id'] for ex in allExp] # 236787594
    expNam = [ex["ExperimentName"] for ex in allExp] # PACT-22-99
    runNam = [ex['Name'] for ex in allExp] # 220616_NB551888_0199_AHMTVVBGXL 
    for index, runID in enumerate(runIds):
        outDir = outPath + expNam[index] # Change output name/dir here
        if isDone[index] == "Complete": # only downloads run status Complete
            if os.path.isdir(outDir) == False:
                try:
                    os.mkdir(outDir)
                    subprocess.call(['chmod', '-R', '+w', outDir])
                    currCmd = f'{cmd3} {runID} -o {outDir}'
                    os.system(currCmd)
                except Exception:
                    logging.error(f'An error occured with mkdir and chmod for {outDir}')
                    pass
        else:
            logging.info(f'Run {runNam[index]} staus is {isDone[index]}')

if not allExp:
    logging.info(f"No {keyName} runs to import from last {days} days")
else:
    download_bs_info(allExp)
