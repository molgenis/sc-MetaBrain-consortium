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

if args.alpha < 0 or args.alpha > 1:
    print("Error in --alpha outside range.")
    exit()

dir = os.path.dirname(args.outfile)
if dir != "":
    os.makedirs(dir, exist_ok=True)

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)

def pvalue_to_zscore(pvalue):
    if pvalue >= 0.9999999999999999:
        return -1.3914582123358836e-16
    elif pvalue <= 1e-323:
        return -38.467405617144344

    return stats.norm.ppf(pvalue / 2)

######################################################

# Convert the alpha to a z-score for speed-up later.
alpha_zscore = abs(pvalue_to_zscore(pvalue=args.alpha))

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

    # Prepare the counters.
    row = {"Cov": cov, "Prefix": prefix, "NrTestedGenes": 0}

    # Loop through the file.
    indices = {}
    sub_indices = {}
    warnings = {}
    with gzopen(fpath, 'r') as f:
        for j, line in enumerate(f):
            values = line.rstrip("\n").split("\t")

            # Process the header.
            if j == 0:
                found = False
                for index, value in enumerate(values):
                    # Check if it is a column with brackets in it (e.g. per dataset column)
                    match = re.match("[A-Za-z]+\(([^()]+)\)[-A-Za-z]*", value)
                    if match:
                        # Replace the dataset string with a placeholder to standardize the
                        # column name.
                        colname = value.replace(match.group(1), "DATASET")

                        # Save which dataset is stored at which index in the dataset
                        # string.
                        ds_indices = {}
                        for ds_index, dataset in enumerate(match.group(1).split(";")):
                            ds_indices[ds_index] = dataset

                            # Replace the placeholder with the actual dataset and add this
                            # colname to the row and colnames_order if we wish to keep
                            # it.
                            ds_colname = colname.replace("DATASET", dataset)
                            if colname in columns_order:
                                row[ds_colname] = 0
                                if ds_colname not in colnames_order:
                                    colnames_order[ds_colname] = columns_order[colname]

                        sub_indices[colname] = ds_indices
                    else:
                        # Else, simply store the colname if we wish to save it.
                        colname = value
                        if colname in columns_order:
                            row[colname] = 0
                            if colname not in colnames_order:
                                colnames_order[colname] = columns_order[colname]

                    # Store the location of the standardized column name to use
                    # as a selection dict for the next lines.
                    indices[colname] = index
                continue

            # Count the row.
            row["NrTestedGenes"] += 1

            # Process the columns of interest based on their mode.
            for colname, mode in columns:
                if colname in indices:
                    # Extract the value of the current column using the index dict
                    # we created from the header.
                    value = values[indices[colname]]

                    # Process the different modes.
                    if mode == "snp":
                        # Count SNPs (e.g. alleles string has length 3)
                        if len(value) == 3:
                            row[colname] += 1
                    elif mode == "sum" or mode == "average":
                        # Keep a total count.
                        row[colname] += int(value)
                    elif mode == "signif":
                        # Check if a value is below a certain threshold.
                        try:
                            value = float(value)
                        except ValueError:
                            if colname not in warnings:
                                warnings[colname] = {}
                            if "ValueError" not in warnings[colname]:
                                warnings[colname]["ValueError"] = 0
                            warnings[colname]["ValueError"] += 1
                            continue

                        if value < args.alpha:
                            row[colname] += 1
                    elif mode == "dataset_z_signif":
                        # Per dataset, convert the z-score to p-values and
                        # check if it is below a certain threshold.
                        for ds_index, ds_zscore in enumerate(value.split(";")):
                            ds_colname = colname.replace("DATASET", sub_indices[colname][ds_index])
                            try:
                                ds_zscore = float(ds_zscore)
                            except ValueError:
                                if ds_colname not in warnings:
                                    warnings[ds_colname] = {}
                                if "ValueError" not in warnings[ds_colname]:
                                    warnings[ds_colname]["ValueError"] = 0
                                warnings[ds_colname]["ValueError"] += 1
                                continue

                            if abs(ds_zscore) >= alpha_zscore:
                                row[ds_colname] += 1
                    else:
                        print("Error, unexpected mode '{}'.".format(mode))
                        exit()
    f.close()

    # Print warnings.
    for colname, colname_warnings in warnings.items():
        warnings_str = ""
        total_warning_count = 0
        for warning_label, warning_count in colname_warnings.items():
            warnings_str += "{} {}s".format(warning_count, warning_label)
            total_warning_count += warning_count
        print("  File {} column '{}' had {}".format(os.sep.join(fpath.split(os.sep)[-2:]), colname, warnings_str))

        if total_warning_count == row["NrTestedGenes"]:
            row[colname] = float("nan")
    if warnings:
        print("")

    # Post-process the average columns.
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

# Sort. In the case of fisherzmetarandom the MetaP is replaced with MetaP-Random and
# as such we cannot assume that MetaP exists. To combat this I will take the lowest value
# of either column, sort on that, and then remove the tmp column.
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
df = df.loc[:, ["Cov", "Prefix", "NrTestedGenes"] + order].rename(columns=rename_columns)

print("\nResults:")
print(df)

print("\nWriting output ...")
df.to_csv(args.outfile, sep="\t", header=True, index=False)

print("\nEND")