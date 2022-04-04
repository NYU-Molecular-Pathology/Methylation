#!/usr/bin/env python3
__author__ = "Jonathan Serrano"
__version__ = "1.1"

import argparse
import os
import shutil
import subprocess
import zipfile
from cmath import nan
from functools import reduce
from pathlib import Path

import numpy as np
import pandas as pd

philipsFtp = "/Volumes/molecular/Molecular/Philips_SFTP"  # path to Philips data dump

parser = argparse.ArgumentParser(description="Provide path to run directory's output dir")
parser.add_argument('-z', '--zipdir', help='Your directory with zip files', required=True)
parser.add_argument('-o', '--outdir', help='Your Genexus file output directory', required=False)
args = parser.parse_args()

#pd.set_option('display.max_columns', 100)
#pd.set_option('display.max_rows', 100)

def getFiles(searchDir, str_pattern):
    fiLi = list(Path(searchDir).rglob(str_pattern))
    return(fiLi)


def un_zipFiles(zipdir, outDir):
    zFiles = list(set(getFiles(zipdir, "*.zip")) -
                  set(getFiles(zipdir, "*copy.zip")))
    zipObj = list(map(zipfile.ZipFile, zFiles))  # create zipfile object
    def joinDir(zp): return os.path.join(
        outDir, os.path.splitext(os.path.basename(zp))[0])
    unzipDir = list(map(lambda zp: joinDir(zp), zFiles))
    totalZip = list(range(0, len(zipObj)))

    unzipped = np.asarray(list(map(lambda pth: os.path.exists(pth), unzipDir)))
    marr = np.ma.MaskedArray(totalZip, mask=unzipped)
    list(map(lambda zp: zipObj[zp].extractall(unzipDir[zp]),  totalZip))
    list(map(lambda zp: zp.close(),  zipObj))


def maskCol(df, kwd, colNam='Sample Name'):
    ngsMsk = df.loc[:, colNam].str.contains(kwd, na=False)
    return(ngsMsk)


def grabFirst(df, fltr, colName, delim='-NGS'):
    if colName == 'Sample Name':
        delim = '-CF2'
    first_str = df.loc[fltr, 'Sample Name'].str.split(delim, 1).str[0]
    df.loc[fltr, colName] = first_str
    return(df)


def grabEnd(df, fltr, colName, col2='Sample Name'):
    last_str = df.loc[fltr, col2].str.split('-', 2).str[-1]
    df.loc[fltr, colName] = last_str
    return(df)


def formatDash(df):
    """Masks Columns which contain NGS, but are not formatted as NGS-
    Args:
        df (dataframe): dataframe from formatNGS function with column named 'Sample Name'
    """
    msk = maskCol(df, 'NGS') & ~maskCol(df, 'NGS-')
    if any(msk):
        ngsNum = df.loc[msk, 'Sample Name'].str
        newNgs = ngsNum[0:3] + "-" + ngsNum[3:5] + "-0" + ngsNum[6:10]
        df.loc[msk, 'Sample Name'] = newNgs

    return(df)


def formatNGS(df):
    """Reformats NGS and CF numbers into separate columns
    Changes NGSYY-### to NGS-YY-0### and separates CF number

    Args:
        df (dataframe): a dataframe from Genexus csv output zip
    """
    df = df.copy()
    # Mask columns split into NGS and CF numbers
    ngsMsk = maskCol(df, '-NGS')
    cf2Msk = maskCol(df,  '-CF2')

    # Take names starting with CF then ending with -CF
    df = grabFirst(df, ngsMsk, 'Library Name')
    df = grabEnd(df, cf2Msk, 'Library Name')

    # Take names starting with NGS then ending with -NGS
    df = grabEnd(df, ngsMsk, 'Sample Name')
    df = grabFirst(df, cf2Msk, 'Sample Name')

    # Formatting NGS from NGSYY-### to NGS-YY-0###
    df = formatDash(df)
    return(df)


def mergeWrite(dfList, outFi, outDir):
    """Concats pandas dataframe list and saves to output directory as name input

    Args:
        dfList (list): list of all read dataframes from listed files
        outFi (string): string name of file to be saved i.e. 'myfileName.csv'
        outDir (string): full output directory path to save outFi
    """
    saveFi = os.path.join(outDir, outFi)
    df_merged = pd.concat(dfList)
    print("Saving file... " + str(os.path.basename(saveFi)))
    pd.DataFrame.to_csv(df_merged, saveFi, index=False)
    return(df_merged)


def chopSams(df):
    chop = df['Sample'].str.split('_').str
    df = df.assign(Sample=chop[0], RunID=chop[1], Result=chop[3])
    first_column = df.pop('Sample')
    df.insert(0, 'Sample Name', first_column)
    df['Library Name'] = df['Sample Name']
    return(df)


def parse_tsvFile(tf, custom=True):
    headKey = None
    if custom:
        headKey = list(range(0, 4))
    df = pd.read_csv(tf, sep="\t", dtype=str, names=headKey,
                     skip_blank_lines=True).dropna(thresh=1)
    df = df.assign(Sample=os.path.basename(Path(tf).parents[0]))
    df = df.copy()
    # CF22-0##-NGS22-### _ 22-GX-004 _ Result _ 270
    df = chopSams(df)
    df = formatNGS(df)
    return(df)


def parse_indels(tf):
    df = parse_tsvFile(tf, False)
    fltr = maskCol(df, '', 'Oncomine Variant Class')
    df2 = df.loc[fltr, :]

    if(df2.empty):
        df2 = df.iloc[0].to_frame().T
        df2.iloc[0, 1:21] = "None"
    return(df2)


def getUnique(df):
    swap = {".0": "", ".1": " RNA", ".2": " NTC_QC", ".3": " NTC_QC_RNA"}
    df = df.iloc[0:1, ]
    s = df.columns.to_series()
    df.columns = s.astype(str) + '.' + s.groupby(s).cumcount().astype(str)
    df.columns = reduce(lambda col, r: col.str.replace(
        r[0], r[1], regex=True), swap.items(),  df.columns)
    return(df)


def parse_QCFile(tf):
    headKey = list(range(0, 6))
    # Columns to pull from csv file containing varying duplicates
    colsToPull = [
        'Sample Name', 'Key Signal', 'Percent Loading', 'Raw Read Accuracy',
        'MAPD', 'Mapped Reads', 'Mean AQ20 Read Length (bp)', 'Mean Read Length (bp)', 'Mean Read Cov', 'Median Mol Cov', 'Uniformity Of Base Coverage',
        'RNA Expression Ctrls Detected', 'Average Reads Per Lane', 'Base Call Accuracy', 'Start Date', 'Library Name'
    ]

    # Final columns for re-index, inserting NAs for any samples missing column names
    finalCols = [
        'Sample Name', 'Key Signal', 'Percent Loading', 'Raw Read Accuracy',
        'MAPD', 'Mapped Reads', 'Mean AQ20 Read Length (bp)', 'Mean Read Length (bp)', 'Mean Read Cov', 'Median Mol Cov', 'Uniformity Of Base Coverage',
        'Mapped Reads RNA', 'Mean Read Length (bp) RNA', 'Mean AQ20 Read Length (bp) RNA', 'RNA Expression Ctrls Detected',
        'Average Reads Per Lane', 'Base Call Accuracy', 'Mapped Reads NTC_QC', 'Mean Read Length (bp) NTC_QC',
        'Mean Read Length (bp) NTC_QC_RNA', 'Start Date', 'Library Name'
    ]
    df = pd.read_csv(tf, sep=",", dtype=str, index_col=0,
                     names=headKey, lineterminator='\n', na_values="",
                     keep_default_na=False, squeeze=True, encoding='utf-8',
                     skip_blank_lines=True).dropna(thresh=1)

    try:
        df = df.T.dropna(thresh=1)[colsToPull]
        print()
    except:
        colsToPull[8] = 'Mean Red Cov'
        df = df.T.dropna(thresh=1)[colsToPull]

    df = getUnique(df.copy())
    df = df.reindex(columns=finalCols)
    df = formatNGS(df)
    return(df)


def readInput(outDir, spr):
    catFi = None
    if(spr == "\t"):
        fiLi = getFiles(outDir, 'Summary.tsv')
        dfLists = list(map(lambda tf: parse_tsvFile(tf), fiLi))
        if dfLists:
            catFi = mergeWrite(dfLists, "callerSum.csv", outDir)
        else:
            print("callerSum.csv is empty")
            
        fiLi = getFiles(outDir, 'Snvindel.tsv')
        dfLists = list(map(lambda tf: parse_indels(tf), fiLi))
        if dfLists:
            catFi = mergeWrite(dfLists, "snvIndels.csv", outDir)
        else:
            print("snvIndels.csv is empty")

    else:
        fiLi = getFiles(outDir, 'Info.csv')
        dfLists = list(map(lambda tf: parse_QCFile(tf), fiLi))
        if dfLists:
            catFi = mergeWrite(dfLists, "qcInfo.csv", outDir)
        else:
            print("Info.csv is empty")
    return(catFi)

def cnvPhilips(cnvDf, sam):
    cnvDf.assign(Test_Case=sam, Variant='CNV', Mutation_Type=cnvDf['AberrationType'],
                 Other=cnvDf['Chrom'], Comments=cnvDf['CopyNumber'], IGV='Copy Number')
    cnvDf['In.NYU'] = cnvDf['In.Philips'] = nan
    dropMut = cnvDf['Mutation.Type'] != "Hemizygous Loss"
    return(cnvDf[dropMut, ])


def getPercent(snvOut, colName):
    snvOut[colName] = snvOut[colName].str.strip("%")
    snvOut[colName] = snvOut[colName].astype('float')
    return(snvOut)


def snvPhilips(snvOut, sam):
    snvOut['Sample'] = sam
    colOrder = ['Sample', "HGNC_gene", "AberrationType", "Coordinate", "tumor_freq",
                "normal_freq", "tumor_dp", "normal_dp", "THERAPY_AVAILABILITY", "HGVSp_Short", "SomaticStatus"]
    # Canonical_HGVS_Protein
    varFilter = "splice|intron|UTR"  # splice synonymous #intron UTR variant
    synFilter = "synonymous variant"
    #colFilter = ['tumor_freq', 'normal_freq', 'tumor_dp']

    snvOut = snvOut[colOrder]  # reorder columns
    snvOut = snvOut.copy()

    snvOut = getPercent(snvOut, "tumor_freq")
    snvOut = getPercent(snvOut, "normal_freq")
    snvOut = getPercent(snvOut, "tumor_dp")

    rule1 = snvOut['tumor_freq'] >= 5.0
    rule2 = snvOut['normal_freq'] < 2.0
    rule3 = snvOut['tumor_dp'] > 200

    snvOut = snvOut.loc[rule1, :]
    snvOut = snvOut.loc[rule2, :]
    snvOut = snvOut.loc[rule3, :]
    
    ngsMsk = snvOut['AberrationType'].str.contains(varFilter, na=False)
    snvOut = snvOut[:][~ngsMsk]
    ngsMsk = snvOut['AberrationType'].str.startswith(synFilter)
    snvOut = snvOut[:][~ngsMsk]
    return(snvOut)


def dataDumpCsv(outDir, sam, csvFi):
    snvFi = os.path.join(outDir, "philipsNGS", sam, csvFi)
    if (os.path.isfile(snvFi)):
        snvInfo = pd.read_csv(snvFi, dtype=str)
        if (len(snvInfo.index) > 0):
            snvOut = snvPhilips(snvInfo, sam)
            return(snvOut)


def checkDataDump(sam, outDir):
    outPath = os.path.join(outDir, "zipfiles")  # path to copy zip files
    # path to search for NGS.zip
    dumpDir = os.path.join(philipsFtp, sam + ".zip")
    # destination to copy NGS.zip
    zipOutDir = os.path.join(outPath, sam + ".zip")
    # folder where to unzip files
    unzipDir = os.path.join(outDir, "philipsNGS")

    if (os.path.isfile(dumpDir)):
        print(dumpDir + " File Exists!")
        if (os.path.isdir(outPath) == False):
            os.mkdir(outPath)
        if (os.path.isfile(zipOutDir) == False):
            try:
                print("Copying file...")
                shutil.copy(dumpDir, zipOutDir)
            except:
                print(sam + " zip file failed to copy to: " + zipOutDir)
        if (os.path.isdir(os.path.join(unzipDir, sam)) == False):
            try:
                print("Unzipping file...")
                zipObj = zipfile.ZipFile(zipOutDir)  # create zipfile object
                zipObj.extractall(unzipDir)
                zipObj.close()
            except:
                print(sam + " zip file failed to extract")
        return(True)
    else:
        print(dumpDir + " does not exists in Philips data dump")
        return(False)


def parseFiles(outDir):
    infoReads = readInput(outDir, ',')
    summaryReads = readInput(outDir, '\t')
    samList = list(infoReads['Sample Name'])
    list(map(lambda sam: checkDataDump(sam, outDir), samList))
    ispmData = list(map(lambda sam: dataDumpCsv(
        outDir, sam, "aberration_snv.csv"), samList))
    if ispmData[0] is not None:
        mergeWrite(ispmData, "snvPhilips.csv", outDir)
    #nexusData = pd.read_csv(os.path.join(outDir, "snvIndel.csv"))
    #philipsData=pd.read_csv(os.path.join(outDir, "snvPhilips.csv"))


def downloadFi():
    # Download the RMD file and Knit using the output csv files
    url = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/Methylation/main/Research/genexus.Rmd"
    command = ["curl", "-o", "genexus.Rmd", "-L", url]
    subprocess.call(command)

def knitReport(outDir):
    runName = os.path.splitext(os.path.basename(outDir))[0] + ".html"
    render = "rmarkdown::render(\'genexus.Rmd\'" + " , " + "output_file="+"\'"+runName + "\'" + \
        " , " + "output_dir ="+"\'"+outDir+"\'" + " , " + \
        "params ="+"list(outDir=\'"+outDir+"\'"+"))"
    command2 = ["Rscript", "--verbose", "-e", render]
    subprocess.call(command2)


def main(zipdir, outDir=None):
    if "/" not in zipdir:
        runYr = "20" + zipdir.split("-")[0] 
        if not runYr.isnumeric():
            runYr = "20" + zipdir.split("-")[1]    
        outDir = "/Volumes/molecular/Molecular/Validation/Genexus/Results/" + runYr+ "/" + zipdir
        zipdir = "/Volumes/molecular/Molecular/Validation/Genexus/downloaded_data/" + zipdir + "/"
        print("outDir set to:")
        print(outDir)
        print("Zip directory is :")
        print(zipdir)
    if not os.path.isdir(zipdir):
        print("The directory is not found:")
        print(zipdir)
    else:
        if os.path.isdir(outDir):
            parseFiles(outDir)
        else:
            un_zipFiles(zipdir, outDir)
            parseFiles(outDir)
        if os.path.isfile("genexus.Rmd")==False:
            downloadFi()
        knitReport(outDir)

#"/Volumes/molecular/Molecular/Validation/Genexus/downloaded_data/22-GX-009/"
#"/Volumes/molecular/Molecular/Validation/Genexus/Results/2022/22-GX-009"
main(args.zipdir, args.outdir)
