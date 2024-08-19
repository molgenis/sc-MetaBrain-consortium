#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse
import os

"""
Syntax: 
./pivot_expression.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--data", type=str, required=True, help="")
parser.add_argument("--filter", type=str, required=False, help="")
parser.add_argument("--values", type=str, required=True, help="")
parser.add_argument("--index", type=str, nargs="+", required=True, help="")
parser.add_argument("--out_index", type=str, required=False, help="")
parser.add_argument("--columns", type=str, nargs="+", required=True, help="")
parser.add_argument("--remove_dupl", action='store_true', help="")
parser.add_argument("--out", type=str, required=True, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

if args.out_index is None:
    if len(args.index) == 1:
        args.out_index = args.index
    else:
        print("Error, please set --out_index")
        exit()

import numpy as np
import pandas as pd

print("Loading data")
df = pd.read_csv(args.data, sep="\t", header=0, index_col=None)
if args.filter == None:
    if "filter" in df.columns:
        print("Error")
        exit()
    df["filter"] = "none"
    args.filter = "filter"
    print("\tLoaded dataframe: {} with shape: {}".format(os.path.basename(args.data), df.shape))
print(df)

print("\nPivotting data.")
for filter_value in df[args.filter].unique():
    print("Filtering: {}".format(filter_value))

    pivot_df = pd.pivot_table(df.loc[df[args.filter] == filter_value, :], values=args.values, index=args.index, columns=args.columns)
    for key in args.index:
        if key != args.out_index:
            pivot_df = pivot_df.droplevel(key)
    print("\tData frame has shape: {}".format(pivot_df.shape))

    if len(set(pivot_df.columns)) != len(pivot_df.columns):
        print("\tError, columns are not unique.")
        continue

    # Check for duplicates in the genes or barcodes.
    n_genes = pivot_df.shape[0]
    index_mask = np.ones(n_genes, dtype=bool)
    if len(set(pivot_df.index)) != len(pivot_df.index):
        print("\tWarning, genes are not unique.")
        if args.remove_dupl:
            seen = set()
            duplicates = []
            for index in pivot_df.index:
                if index in seen:
                    duplicates.append(index)
                seen.add(index)

            index_mask = np.array([gene not in duplicates for gene in pivot_df.index], dtype=bool)
            print("\tRemoving '{:,}' duplicate genes having {:,} rows in total.".format(len(duplicates), np.size(index_mask) - np.sum(index_mask)))

    # Save.
    outpath = args.out + filter_value.replace(" / ", "").replace(" ", ".") + ".txt.gz"
    pivot_df.loc[index_mask, :].to_csv(outpath, sep="\t", header=True, index=True, compression="gzip")
    print("\tSaved dataframe: {} with shape: {}".format(os.path.basename(outpath), pivot_df.shape))

print("Done")