#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse


"""
Syntax: 
./check_cellbender_correction.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--work_dir", type=str, required=True, help="")
parser.add_argument("--cellranger_dir", type=str, required=True, help="")
parser.add_argument("--cellbender_dir", type=str, required=True, help="")
parser.add_argument("--wg2_metadata_path", type=str, required=True, help="")
parser.add_argument("--wg2_pairing_path", type=str, required=True, help="")
parser.add_argument("--pool", type=str, required=True, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")


import os
import scanpy
import pandas as pd

cellranger_path = os.path.join(args.cellranger_dir, args.pool, "outs", "filtered_feature_bc_matrix.h5")
if not os.path.exists(cellranger_path):
    print("Error, CellRanger path '{}' does not exist.".format(cellranger_path))
    exit()

cellbender_path = os.path.join(args.cellbender_dir, args.pool, "cellbender_remove_background_output_filtered.h5")
if not os.path.exists(cellbender_path):
    print("Error, CellBender path '{}' does not exist.".format(cellbender_path))
    exit()

if not os.path.exists(args.wg2_metadata_path):
    print("Error, WG2 metadata path '{}' does not exist.".format(args.wg2_metadata_path))
    exit()

if not os.path.exists(args.wg2_pairing_path):
    print("Error, WG2 pairing path '{}' does not exist.".format(args.wg2_pairing_path))
    exit()

def load_h5_file(inpath):
    print("\tLoading '{}'".format(inpath))
    data = scanpy.read_10x_h5(inpath)
    df = pd.DataFrame(data.X.toarray(), index=data.obs_names, columns=data.var_names)
    print("\tLoaded data frame with shape: {}".format(df.shape))
    return df

print("Loading CellRanger count matrix")
cellranger_df = load_h5_file(inpath=cellranger_path)

print("Loading CellBender count matrix")
cellbender_df = load_h5_file(inpath=cellbender_path)

print("Loading cell annotation")
wg2_metadata = pd.read_csv(args.wg2_metadata_path, sep="\t", header=0, index_col=None)
wg2_metadata = wg2_metadata.loc[wg2_metadata["Pool"] == args.pool, :]
wg2_pairing = pd.read_csv(args.wg2_pairing_path, sep=";", header=0, index_col=None)
wg2_metadata = wg2_metadata.merge(wg2_pairing, how="left")
wg2_metadata.index = [barcode.split("_")[0] + "-1" for barcode in wg2_metadata["Barcode"]]
del wg2_pairing

print("Overlap")
overlap_barcodes = list(set(cellbender_df.index).intersection(set(cellranger_df.index)).intersection(set(wg2_metadata.index)))
print("\tN-barcodes overlap: {:,}".format(len(overlap_barcodes)))
overlap_genes = list(set(cellbender_df.columns).intersection(set(cellranger_df.columns)))
print("\tN-genes overlap: {:,}".format(len(overlap_genes)))
cellranger_df = cellranger_df.loc[overlap_barcodes, overlap_genes]
cellbender_df = cellbender_df.loc[overlap_barcodes, overlap_genes]
wg2_metadata = wg2_metadata.loc[overlap_barcodes, :]

print(cellranger_df)
print(cellbender_df)
print(wg2_metadata)

delta_df = (cellranger_df - cellbender_df).sum(axis=1).to_frame()
delta_df.columns = ["delta.counts"]
print(delta_df)
delta_df = delta_df.merge(wg2_metadata, left_index=True, right_index=True)
# cell_types = ["predicted.class", "predicted.subclass", "L1"]
cell_types = ["L1"]
for cell_type in cell_types:
    print("Sum:")
    print(delta_df.groupby(cell_type)["delta.counts"].sum())
    print("Sum %:")
    print(delta_df.groupby(cell_type)["delta.counts"].sum() / delta_df["delta.counts"].sum())
    print("Min:")
    print(delta_df.groupby(cell_type)["delta.counts"].min())
    print("Max:")
    print(delta_df.groupby(cell_type)["delta.counts"].max())
    print("Mean:")
    print(delta_df.groupby(cell_type)["delta.counts"].mean())

