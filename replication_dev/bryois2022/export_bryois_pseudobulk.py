#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse

"""
Syntax: 
./export_bryois_pseudobulk.py -h
"""

parser = argparse.ArgumentParser(description="")
parser.add_argument("--bryois_indir", type=str, required=True, help="")
parser.add_argument("--outdir", type=str, required=True, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import os
import pyreadr
import pandas as pd

print("Loading data")
ad_df = pyreadr.read_r(os.path.join(args.bryois_indir, "Brain_single_cell_eqtl", "data_sensitive", "expression", "ad_sum_expression.individual_id.rds"))[None]
ad_df["dataset"] = "ad"
ms_df = pyreadr.read_r(os.path.join(args.bryois_indir, "Brain_single_cell_eqtl", "data_sensitive", "expression", "ms_sum_expression.individual_id.rds"))[None]
ms_df["dataset"] = "ms"
df = pd.concat([ad_df, ms_df], axis=0)
del ad_df, ms_df
df.to_csv(os.path.join(args.outdir, "sum_expression.individual_id.tsv.gz"), sep="\t", header=True, index=False, compression="gzip")
print(df)

print("Unique individual_id")
ct_ind_df = df[["cell_type", "individual_id", "n_cells", "dataset"]].drop_duplicates()
ct_ind_df.to_csv(os.path.join(args.outdir, "cells_per_ct_per_individual.tsv.gz"), sep="\t", header=True, index=False, compression="gzip")
print(ct_ind_df)

for cell_type in ct_ind_df["cell_type"].unique():
    ct = cell_type.replace(" / ", "...").replace(" ", ".")
    df.loc[df["cell_type"] == cell_type, :].to_csv(os.path.join(args.outdir, ct + ".sum_expression.individual_id.tsv.gz"), sep="\t", header=True, index=False, compression="gzip")
    ct_ind_df.loc[ct_ind_df["cell_type"] == cell_type, :].to_csv(os.path.join(args.outdir, ct + ".cells_per_individual.tsv.gz"), sep="\t", header=True, index=False, compression="gzip")

ind_df = ct_ind_df[["individual_id", "n_cells"]].groupby(["individual_id"]).sum()
ind_df.to_csv(os.path.join(args.outdir, "cells_per_individual.tsv.gz"), sep="\t", header=True, index=False, compression="gzip")

for cell_type in ct_ind_df["cell_type"].unique():
    ct = cell_type.replace(" / ", "...").replace(" ", ".")
    if ct != "Excitatory.neurons":
        continue
    for individual_id in ct_ind_df["individual_id"].unique():
        if individual_id != "SM-CTECO":
            continue
        individual_df = df.loc[(df["cell_type"] == cell_type) & (df["individual_id"] == individual_id), ["symbol", "counts"]]
        individual_df.set_index("symbol", inplace=True)
        individual_df.index.name = None
        individual_df.columns = [individual_id]
        individual_df.to_csv(os.path.join(args.outdir, ct + "." + individual_id + ".pseudobulk.tsv.gz"), sep="\t", header=True, index=True, compression="gzip")

