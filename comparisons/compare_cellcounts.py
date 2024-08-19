#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse

"""
Syntax: 
./compare_cellcounts.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--indir", type=str, required=False, default=None, help="")
parser.add_argument("--data1", type=str, required=False, default=None, help="")
parser.add_argument("--label1", type=str, required=False, default=None, help="")
parser.add_argument("--data2", type=str, required=False, default=None, help="")
parser.add_argument("--label2", type=str, required=False, default=None, help="")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import warnings
warnings.simplefilter("ignore", UserWarning)


import pandas as pd
import scanpy
import glob
import os

def get_count_matrices(self, type):
    pool_paths = glob.glob(os.path.join(self.input_dir, type, "*"))

    data = []
    for pool_path in pool_paths:
        pool = os.path.basename(pool_path)

        if type == "CellRanger":
            inpath = os.path.join(self.input_dir, type, pool, "outs", "filtered_feature_bc_matrix.h5")
            if os.path.exists(inpath):
                data.append([pool, inpath])
        elif type == "CellBender":
            inpath = os.path.join(self.input_dir, type, pool, "cellbender_remove_background_output_filtered.h5")
            if os.path.exists(inpath):
                data.append([pool, inpath])
        else:
            print("Error, unknown type")
            exit()

    return pd.DataFrame(data, columns=["Pool", type])

def load_coo_df(fpath):
    data = scanpy.read_10x_h5(fpath)
    barcodes = {index: barcode for index, barcode in enumerate(data.obs_names)}

    var_names_df = data.var['gene_ids'].to_frame().reset_index(names="gene_symbol")
    var_names_df["gene"] = var_names_df["gene_ids"] + "_" + var_names_df["gene_symbol"]
    if (data.var_names.to_numpy() != var_names_df["gene_symbol"]).all():
        print("Error in merging genes metadata.")
        exit()
    genes = {index: barcode for index, barcode in enumerate(var_names_df["gene"])}
    del var_names_df

    coo = data.X.tocoo(copy=False)
    df = pd.DataFrame({'index': coo.row, 'col': coo.col, 'data': coo.data})
    df["barcode"] = df["index"].map(barcodes)
    df["gene"] = df["col"].map(genes)
    return df[["barcode", "gene", "data"]]


df = pd.DataFrame({"Pool": ["Manual"], args.label1: [args.data1], args.label2: [args.data2]})
if args.data1 is None and args.data2 is None:
    df1 = args.get_count_matrices(type="CellRanger")
    df2 = args.get_count_matrices(type="CellBender")
    df = df1.merge(df2, on="Pool")
    del df1, df2
    args.label1 = "CellRanger"
    args.label2 = "CellBender"
print(df)

for _, row in df.iterrows():
    print("Processing Pool '{}'".format(row["Pool"]))
    df1 = load_coo_df(fpath=row[args.label1])
    print("\tLoaded {} h5 file with {:,} barcodes and {:,} features.".format(args.label1, len(df1["barcode"].unique()), len(df1["gene"].unique())))

    print(df1)
    df2 = load_coo_df(fpath=row[args.label2])
    print("\tLoaded {} h5 file with {:,} barcodes and {:,} features.".format(args.label2, len(df2["barcode"].unique()), len(df2["gene"].unique())))
    print(df2)

    # Merge the sparse matrix coordinates.
    df = df1.merge(df2, how="left", on=["barcode", "gene"], suffixes=["_" + args.label1, "_" + args.label2]).fillna(0)
    df["delta"] = df["data_" + args.label1] - df["data_" + args.label2]
    value_cols = [column for column in df.columns if column not in ["barcode", "gene"]]
    df.sort_values(by="delta", ascending=False, inplace=True)
    del df1, df2
    print("Merged coordinates data:")
    print(df)
    print("\tMerged h5 files with {:,} barcodes and {:,} features.".format(len(df["barcode"].unique()), len(df["gene"].unique())))

    print("\nGene level summary:")
    gene_df = df.groupby("gene").sum(value_cols).sort_values(by="delta", ascending=False)
    gene_df["delta%"] = gene_df["delta"] / gene_df[["data_" + args.label1, "data_" + args.label2]].max(axis=1)
    # gene_df.sort_values(by="delta%", ascending=False, inplace=True)
    print(gene_df)
    print(gene_df.iloc[:50, :])

    print("\nBarcode level summary:")
    barcode_df = df.groupby("barcode").sum(value_cols).sort_values(by="delta", ascending=False)
    barcode_df["delta%"] = barcode_df["delta"] / barcode_df[["data_" + args.label1, "data_" + args.label2]].max(axis=1)
    # barcode_df.sort_values(by="delta%", ascending=False, inplace=True)
    print(barcode_df)
    print(barcode_df.iloc[:50, :])

print("\nDone")