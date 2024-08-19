#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse
import os

"""
Syntax: 
./transform_cell_metadata.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--poolsheet", type=str, required=True, help="")
parser.add_argument("--wg1", type=str, required=True, help="")
parser.add_argument("--wg2", type=str, required=True, help="")
parser.add_argument("--cell_annot", type=str, required=True, help="")
parser.add_argument("--mode", type=str, required=True, choices=["cellbender_to_cellranger", "cellranger_to_cellbender"], help="")
args = parser.parse_args()

if not args.poolsheet.endswith(".tsv"):
    print("Error, expected 'tsv' suffix for --poolsheet.")
    exit()
if not args.wg1.endswith(".tsv.gz"):
    print("Error, expected 'tsv.gz' suffix for --wg1.")
    exit()
if not args.wg2.endswith(".tsv.gz"):
    print("Error, expected 'tsv.gz' suffix for --wg2.")
    exit()
if not args.cell_annot.endswith(".tsv.gz"):
    print("Error, expected 'tsv.gz' suffix for --cell_annot.")
    exit()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import pandas as pd

def transform_poolsheet(poolsheet_df):
    print("Updating poolsheet")
    transform = {
        "Pool": None,
        "Counts": {
            "CellBender": "CellBender/POOL/cellbender_remove_background_output_filtered.h5",
            "CellRanger": "CellRanger/POOL/outs/filtered_feature_bc_matrix.h5"
        },
        "Barcodes": {
            "CellBender": "CellBender/POOL/cellbender_remove_background_output_cell_barcodes.csv",
            "CellRanger": "CellRanger/POOL/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
        },
        "Bam": None
    }

    updated_data = {column: [] for column in poolsheet_df.columns}
    for _, row in poolsheet_df.iterrows():
        pool = row["Pool"]
        for column, value in row.to_dict().items():
            transform_settings = transform[column]

            updated_value = value
            if transform_settings is not None is None:
                if args.mode == "cellbender_to_cellranger":
                    updated_value = value.replace(transform_settings["CellBender"].replace("POOL", pool), transform_settings["CellRanger"].replace("POOL", pool))
                elif args.mode == "cellranger_to_cellbender":
                    updated_value = value.replace(transform_settings["CellRanger"].replace("POOL", pool), transform_settings["CellBender"].replace("POOL", pool))
                else:
                    print("Unexpected mode {}".format(args.mode))
                    exit()

            if column != "Pool":
                if not os.path.exists(updated_value):
                    print("Error, updated value '{}' does not exist".format(updated_value))
                    exit()

            updated_data[column].append(updated_value)
    return pd.DataFrame(updated_data)

def load_barcodes(poolsheet_df):
    data = []
    n_rows = 0
    for _, row in poolsheet_df.iterrows():
        df = pd.read_csv(row["Barcodes"], sep="\t", header=None, index_col=None)
        n_rows += df.shape[0]
        df.columns = ["Barcode"]
        df["Pool"] = row["Pool"]
        df["Barcode"] = df["Barcode"].str.split("-", n=1, expand=True)[0] + "_" + df["Pool"]
        data.append(df)
    print("  loaded {:,} cells".format(n_rows))
    return pd.concat(data, axis=0)

def transform_metadata(metadata_df, barcodes_df):
    overlapping_barcodes = len(set(metadata_df["Barcode"]).intersection(set(barcodes_df["Barcode"])))
    print("\tinfo for {:,} / {:,} barcodes found ({:.2f}%)".format(overlapping_barcodes, barcodes_df.shape[0], (100 / barcodes_df.shape[0]) * overlapping_barcodes))
    return barcodes_df.merge(metadata_df, how="left").fillna("NA")

####################################################################################

suffix = ""
if args.mode == "cellbender_to_cellranger":
    suffix = "CellRanger"
elif args.mode == "cellranger_to_cellbender":
    suffix = "CellBender"
else:
    print("Unexpected mode {}".format(args.mode))
    exit()

print("Transforming poolsheet")
poolsheet_df = pd.read_csv(args.poolsheet, sep="\t", header=0, index_col=None)
transformed_poolsheet = transform_poolsheet(poolsheet_df=poolsheet_df)
print(transformed_poolsheet)
del poolsheet_df
#print(args.poolsheet.replace(".tsv", "_" + suffix + ".tsv"))
transformed_poolsheet.to_csv(args.poolsheet.replace(".tsv", "_" + suffix + ".tsv"), sep="\t", header=True, index=False)

print("\nLoading barcodes")
barcodes_df = load_barcodes(poolsheet_df=transformed_poolsheet)
del transformed_poolsheet

print("\nUpdating workgroup1 metadata")
wg1_metadata_df = pd.read_csv(args.wg1, sep="\t", header=0, index_col=None)
transformed_wg1_df = transform_metadata(metadata_df=wg1_metadata_df, barcodes_df=barcodes_df)
#print(args.wg1.replace(".tsv.gz", "_" + suffix + ".tsv.gz"))
transformed_wg1_df.to_csv(args.wg1.replace(".tsv.gz", "_" + suffix + ".tsv.gz"), sep="\t", header=True, index=False, compression="gzip")
del wg1_metadata_df, transformed_wg1_df

print("\nUpdating workgroup2 metadata")
wg2_metadata_df = pd.read_csv(args.wg2, sep="\t", header=0, index_col=None)
transformed_wg2_df = transform_metadata(metadata_df=wg2_metadata_df, barcodes_df=barcodes_df)
#print(args.wg2.replace(".tsv.gz", "_" + suffix + ".tsv.gz"))
transformed_wg2_df.to_csv(args.wg2.replace(".tsv.gz", "_" + suffix + ".tsv.gz"), sep="\t", header=True, index=False, compression="gzip")
del wg2_metadata_df, transformed_wg2_df

print("\nUpdating cell annotations")
cell_annot_df = pd.read_csv(args.cell_annot, sep="\t", header=0, index_col=None)
transformed_cell_annot_df = transform_metadata(metadata_df=cell_annot_df, barcodes_df=barcodes_df)
#print(args.cell_annot.replace(".tsv.gz", "_" + suffix + ".tsv.gz"))
transformed_cell_annot_df.to_csv(args.cell_annot.replace(".tsv.gz", "_" + suffix + ".tsv.gz"), sep="\t", header=True, index=False, compression="gzip")
del cell_annot_df, transformed_cell_annot_df

print("\nDone")