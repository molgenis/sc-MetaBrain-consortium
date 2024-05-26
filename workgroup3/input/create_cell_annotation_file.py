#!/usr/bin/env python
# Author: M. Vochteloo
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--poolsheet", required=True, type=str, help="")
parser.add_argument("--outdir", required=True, type=str, help="")
parser.add_argument("--outfile", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import pandas as pd
import os

poolsheet_df = pd.read_csv(args.poolsheet, sep="\t", header=0, index_col=None)
print(poolsheet_df)

output_columns = ["Barcode", "sequencing_platform", "sequencing_run", "sequencing_lane", "scrna_platform", "plate_based", "umi_based", "biomaterial", "sorting", "cell_treatment", "sample_condition"]

cell_annot_df_list = []
for _, row in poolsheet_df.iterrows():
    cell_annot_df = pd.read_csv(row["Barcodes"], sep="\t", header=None, index_col=None)
    cell_annot_df.columns = ["Barcode"]
    cell_annot_df["Barcode"] = [barcode.split("-")[0] + "_" + str(row["Pool"]) for barcode in cell_annot_df["Barcode"]]
    cell_annot_df["sequencing_run"] = row["Pool"]
    cell_annot_df_list.append(cell_annot_df)
cell_annot_df = pd.concat(cell_annot_df_list, axis=0)
cell_annot_df["cell_treatment"] = "UT"
for column in output_columns:
    if column in cell_annot_df.columns:
        continue
    cell_annot_df[column] = "NONE"
print(cell_annot_df)

cell_annot_df[output_columns].to_csv(os.path.join(args.outdir, args.outfile + "_empty.tsv.gz"), sep="\t", header=True, index=False)
print("Done")