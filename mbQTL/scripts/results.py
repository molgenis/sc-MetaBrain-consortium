#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import pandas as pd
from scipy import stats
import gzip
import glob
import os
import re

parser = argparse.ArgumentParser(description="")
parser.add_argument("--indir", required=True, type=str, help="")
parser.add_argument("--alpha", required=False, type=float, default=0.05, help="")
parser.add_argument("--outfile", required=True, type=str, help="The output file where results will be saved.")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

dir = os.path.dirname(args.outfile)
if dir != "":
    os.makedirs(dir, exist_ok=True)

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)

# Define the columns we want to keep track of and how we want to summarise them.
columns = [
    ("NrTestedSNPs", "sum"),
    ("SNPAlleles", "snp"),
    ("MetaPN", "average"),
    ("NrDatasets", "average"),
    ("MetaP", "signif"),
    ("MetaP-Random", "signif"),
    ("BetaAdjustedMetaP", "signif"),
    ("BonfAdjustedMetaP", "signif"),
    ("BonfBHAdjustedMetaP", "signif"),
    ("DatasetZScores(DATASET)", "dataset_z_signif"),
    ("DatasetFisherZ(DATASET)", "dataset_z_signif"),
    ("DatasetFisherZ(DATASET)-Random", "dataset_z_signif"),
    ("bh_fdr", "signif"),
    ("qval", "signif"),
]
rename_columns = {"SNPAlleles": "NrSNPVariants", "MetaPN": "AvgMetaPN", "NrDatasets": "AvgNrDatasets"}
columns_order = {column[0]: index for index, column in enumerate(columns)}

print("Counting number of eQTLs ...")
row_index = 0
data = {}
colnames_order = {}
for i, fpath in enumerate(glob.glob(os.path.join(args.indir, "*/*-TopEffectsWithMultTest.txt"))):
    # Parse the file path.
    cov = fpath.split(os.sep)[-2]
    prefix = fpath.split(os.sep)[-1].replace("-TopEffectsWithMultTest.txt", "")
    if prefix != "Excitatory.neurons.BryoisEffects.QTLMeta.FISHERZMETA":
        continue

    # Prepare the counters.
    row = {"Cov": cov, "Prefix": prefix, "NrTestedGenes": 0}

    # Loop through the file.
    indices = {}
    sub_indices = {}
    with gzopen(fpath, 'r') as f:
        for j, line in enumerate(f):
            values = line.rstrip("\n").split("\t")
            if j == 0:
                # Process the header.
                found = False
                for index, value in enumerate(values):
                    match = re.match("[A-Za-z]+\(([^()]+)\)[-A-Za-z]*", value)
                    if match:
                        colname = value.replace(match.group(1), "DATASET")
                        ds_indices = {}
                        for ds_index, dataset in enumerate(match.group(1).split(";")):
                            ds_indices[ds_index] = dataset
                            ds_colname = colname.replace("DATASET", dataset)
                            if colname in columns_order:
                                row[ds_colname] = 0
                                if ds_colname not in colnames_order:
                                    colnames_order[ds_colname] = columns_order[colname]
                        sub_indices[colname] = ds_indices
                    else:
                        colname = value
                        if colname in columns_order:
                            row[colname] = 0
                            if colname not in colnames_order:
                                colnames_order[colname] = columns_order[colname]

                    indices[colname] = index
                continue

            # Count the row.
            row["NrTestedGenes"] += 1

            # Process the columns of interest based on their mode.
            for colname, mode in columns:
                if colname in indices:
                    value = values[indices[colname]]
                    if mode == "snp":
                        if len(value) == 3:
                            row[colname] += 1
                    elif mode == "sum" or mode == "average":
                        row[colname] += int(value)
                    elif mode == "signif":
                        if float(value) < args.alpha:
                            row[colname] += 1
                    elif mode == "dataset_z_signif":
                        for ds_index, ds_zscore in enumerate(value.split(";")):
                            if ds_zscore == "-":
                                continue
                            ds_colname = colname.replace("DATASET", sub_indices[colname][ds_index])
                            ds_pvalue = stats.norm.cdf(-abs(float(ds_zscore))) * 2
                            if ds_pvalue < args.alpha:
                                row[ds_colname] += 1
                    else:
                        print("Error, unexpected mode '{}'.".format(mode))
                        exit()
    f.close()

    # Post-process the averge columns.
    for colname, mode in columns:
        if mode == "average":
            row[colname] = round(row[colname] / row["NrTestedGenes"], 1)

    # Save.
    data[row_index] = row
    row_index += 1

if len(data) == 0:
    print("Error, no Data")
    exit()

# Convert to Pandas.
df = pd.DataFrame(data).T

# Sort.
pvalue_sort_col = "MetaP"
if "MetaP-Random" in df:
    df["tmpP"] = df[["MetaP", "MetaP-Random"]].min(axis=1)
    pvalue_sort_col = "tmpP"
df.sort_values(by=[pvalue_sort_col, "NrTestedGenes", "Cov", "Prefix"], ascending=[False, True, True, True])
if "tmpP" in df.columns:
    df.drop(["tmpP"], axis=1, inplace=True)

# Reorder and rename.
colnames_order = list(colnames_order.items())
colnames_order.sort(key=lambda x: (x[1], x[0]))
order = [co[0] for co in colnames_order]
df = df.loc[:, ["Cov", "Prefix", "NrTestedGenes"] + order].rename(columns=rename_columns).T

print("\nResults:")
print(df)

print("\nWriting output ...")
df.to_csv(args.outfile, sep="\t", header=False, index=True)

print("\nEND")