#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse


"""
Syntax:
./create_wg0_samplesheet.py -h 
"""

parser = argparse.ArgumentParser(
    description="")
parser.add_argument("--work_dir", type=str, required=True, help="")
parser.add_argument("--method", type=str, required=True, default=["CellBender", "CellRanger"], help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import pandas as pd
import glob
import os

def load_fpaths(method):
    pool_index = None
    count_file = None
    barcode_file = None
    bam_file = os.path.join("outs", "possorted_genome_bam.bam")
    if method == "CellRanger":
        count_file = os.path.join("outs", "filtered_feature_bc_matrix.h5")
        barcode_file = os.path.join("outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
    elif method == "CellBender":
        count_file = os.path.join("cellbender_feature_bc_matrix_filtered.h5")
        barcode_file = os.path.join("cellbender_feature_bc_matrix_cell_barcodes.csv")
    else:
        print("Error, unexpected method '{}'".format(method))
        exit()

    data = []
    for folder, filename, colname in [(method, count_file, "Counts"),
                                      (method, barcode_file, "Barcodes"),
                                      ("CellRanger", bam_file, "Bam")]:
        pool_index = len(filename.split(os.sep)) + 1
        fpaths = glob.glob(os.path.join(args.work_dir, folder, "*", filename))

        column = {}
        for fpath in fpaths:
            pool = fpath.split(os.sep)[-pool_index].replace("Run1", "")
            column[pool] = fpath
        data.append(pd.Series(column, name=colname))
    df = pd.concat(data, axis=1)
    df.index.name = "Pool"
    df.reset_index(drop=False, inplace=True)
    df.dropna(axis=1, how='all', inplace=True)

    return df

print("Creating WGO samplesheet file")
df = load_fpaths(method=args.method)

df.to_csv(os.path.join(args.work_dir, "Combine_Results", "wg0_file_directories_{}.tsv".format(args.method)), sep="\t", header=True, index=False)
print("Saved file of shape: {}".format(df.shape))