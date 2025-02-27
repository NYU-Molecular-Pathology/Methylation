__version__ = "1.2.1"
import pandas as pd
import os
import datetime
import argparse
import warnings

warnings.filterwarnings('ignore', category=UserWarning, module='openpyxl')


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--template", type=str, required=True,
                        help="Template Sheet")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output directory")
    parser.add_argument("--I5_index", type=str, required=False,
                        help="I5 Index sheet",
                        default='/gpfs/data/molecpathlab/data/archer/Archer_Index2_Sequences.xlsx')
    parser.add_argument("--debug", action="store_true",
                        help="debug mode print commands")
    return parser.parse_args()


def get_i5_index(row):
    try:
        idx = I5_Index_df.loc[I5_Index_df['P5 Adapter']
                              == row['I5_Index_ID']].index
        return I5_Index_df.iloc[idx]['Reverse Complement of P5 sequence for NextSeq sample sheets'].values[0]
    except Exception as e:
        print(e)


def get_i7_index(row):
    try:
        idx = I7_Index_df.loc[I7_Index_df['I7_Index_1_ID']
                              == row['I7_Index_ID']].index
        return I7_Index_df.iloc[idx]['index'].values[0]
    except Exception as e:
        print(e)


def format_template(template_file):
    ts = pd.read_excel(template_file, 'Worksheet', skiprows=5,
                       index_col=None, na_values=['NA'], engine='openpyxl')
    lastNrow = ts['Sample #'].index[ts['Sample #'].apply(pd.isnull)][0] - 1
    print("Last Row -------------------------", lastNrow)
    template_sheet = ts.loc[:lastNrow, :].copy()
    template_sheet['index'] = template_sheet.apply(
        lambda x: get_i7_index(x), axis=1)
    template_sheet['index2'] = template_sheet.apply(
        lambda x: get_i5_index(x), axis=1)
    template_sheet['Sample_ID'] = template_sheet['Accession#'].astype(
        str) + "-" + template_sheet['RNA #']
    template_sheet['Sample_Name'] = template_sheet['Sample_ID']
    template_sheet['Sample_Plate'] = "Plate"
    template_sheet['Sample_Well'] = "Well"
    template_sheet['Sample_Project'] = 'ArcherDx_Run'
    template_sheet['Description'] = 'Description'
    template_sheet['GenomeFolder'] = 'PhiX\Illumina\RTA\Sequence\WholeGenomeFASTA'
    sample_sheet = template_sheet[['Sample_ID', 'Sample_Name', 'Sample_Plate', 'Sample_Well',
                                   'I7_Index_ID', 'index', 'I5_Index_ID', 'index2', 'Sample_Project', 'Description', 'GenomeFolder']]
    return sample_sheet


def get_sheet_header(proj_name, date_str):
    header = """[Header]
IEMFileVersion	4
Investigator Name
Project Name	%s
Experiment Name	%s
Date	%s
Workflow	GenerateFASTQ
Application	FASTQ Only
Assay	Nextera
Description	Description
Chemistry	Amplicon
[Reads]
151
151
[Settings]
Adapter	CTGTCTCTTATACACATCT
,,,,,,,,,,
[Data]
""" % (proj_name, proj_name, date_str)
    return header


def print_debug_msg(I7_Index_df, I5_Index_df):
    print("======I7_Index======")
    print(I7_Index_df.columns.values)
    print(I7_Index_df.head(3).to_string())
    print("======I5_Index======")
    print(I5_Index_df.columns.values)
    print(I5_Index_df.head(3).to_string())


def main():
    global I7_Index_df
    global I5_Index_df
    args = get_options()
    today = datetime.date.today()
    date_str = "%s/%s/%s" % (today.month, today.day, today.year)
    I5_xls = pd.ExcelFile(args.I5_index, engine='openpyxl')
    I7_Index_df = pd.read_excel(I5_xls, 'Sheet1', engine='openpyxl')
    A_df = pd.read_excel(I5_xls, 'Set A', engine='openpyxl')
    B_df = pd.read_excel(I5_xls, 'Set B', engine='openpyxl')
    C_df = pd.read_excel(I5_xls, 'Set C', engine='openpyxl')
    I5_Index_df = pd.concat([A_df, B_df, C_df]).reset_index(drop=True)
    template_file = pd.ExcelFile(args.template, engine='openpyxl')
    run_info = pd.read_excel(template_file, 'Worksheet', skiprows=2, nrows=1,
                             index_col=1, na_values=['NA'], engine='openpyxl')['Sequencing Run ID']
    run_name = run_info[0]
    run_proj = run_info.index.values[0]
    print("RUN INFO -------------------------",
          "Project: ", run_proj, "Run ID: ", run_name)

    if args.debug:
        print_debug_msg(I7_Index_df, I5_Index_df)

    header = get_sheet_header(run_proj, date_str)
    sheet_out = os.path.join(args.output, '%s-SampleSheet.csv' % run_name)
    with open(sheet_out, 'w') as of:
        of.write(header.replace("\t", ","))

    sample_sheet = format_template(template_file)
    print(sample_sheet.to_string())
    sample_sheet.to_csv(sheet_out, sep=",", index=False, mode='a')


if __name__ == "__main__":
    main()
