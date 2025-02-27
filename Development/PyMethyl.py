import subprocess
import sys
import os
import numpy as np
import pandas as pd


def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

def check_and_install(package):
    try:
        __import__(package)
    except ImportError:
        install(package)


check_and_install("methylprep")
check_and_install("methylcheck")

def calculate_lrr_baf_tumor_purity(idat_path):
    # Process the .idat files using methylprep
    process.run_pipeline(data_dir=idat_path, save_uncorrected=True)

    # Load the methylprep-generated CSV files
    green_filepath = os.path.join(idat_path, 'export', 'uncorrected', 'green_channel.csv')
    red_filepath = os.path.join(idat_path, 'export', 'uncorrected', 'red_channel.csv')

    green = pd.read_csv(green_filepath, index_col=0)
    red = pd.read_csv(red_filepath, index_col=0)

    # Apply quantile normalization to green and red channel data
    green_norm = methylcheck.quantile_normalize(green)
    red_norm = methylcheck.quantile_normalize(red)

    # Calculate LRR
    lrr = (green_norm + red_norm).applymap(lambda x: 0 if x == 0 else np.log2(x))

    # Calculate BAF
    baf = green_norm / (green_norm + red_norm)

    # Calculate tumor purity (using Bcell as normal cell proportion)
    # This is a simplified estimation of tumor purity
    tumor_purity = 1 - baf.mean(axis=1)

    return lrr, baf, tumor_purity

# Set the path to the folder containing the .idat files
idat_path = "path/to/your/idat/files"

# Calculate LRR, BAF, and tumor purity
lrr, baf, tumor_purity = calculate_lrr_baf_tumor_purity(idat_path)

print("LRR:")
print(lrr)
print("BAF:")
print(baf)
print("Tumor Purity:")
print(tumor_purity)
#pip install methylprep methylcheck
