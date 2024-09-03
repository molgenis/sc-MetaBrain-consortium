#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import numpy as np
import pandas as pd
from scipy import stats
import gzip
import glob
import os
import re

parser = argparse.ArgumentParser(description="")
parser.add_argument("--indir", required=True, type=str, help="")
parser.add_argument("--n_tests_col", required=False, type=str, default="ntests", help="")
parser.add_argument("--nom_pvalue_col", required=False, type=str, default="pvalue", help="")
parser.add_argument("--beta_adj_pvalue_col", required=False, type=str, default="perm-pvalue", help="")
parser.add_argument("--bonf_pvalue_col", required=False, type=str, default="bonf_pvalue", help="")
parser.add_argument("--bonf_bh_fdr_col", required=False, type=str, default="bonf_bh_fdr", help="")
parser.add_argument("--bh_fdr_col", required=False, type=str, default="bh_fdr", help="")
parser.add_argument("--qvalue_col", required=False, type=str, default="qvalue_col", help="")
parser.add_argument("--alpha", required=False, type=float, default=0.05, help="")
parser.add_argument("--out", required=True, type=str, help="The output file where results will be saved.")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

os.makedirs(os.path.dirname(args.out), exist_ok=True)

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)

print(os.path.join(args.indir, "*/*-TopEffectsWithMulTest.txt"))


print("Counting number of eQTLs ...")
columns = [args.nom_pvalue_col, args.beta_adj_pvalue_col, args.bonf_pvalue_col, args.bonf_bh_fdr_col, args.bh_fdr_col, args.qvalue_col]
lookup = {}
row_index = 0
data = {}
dataset_nom_pvalue_columns = []
for i, fpath in enumerate(glob.glob(os.path.join(args.indir, "*/*-TopEffectsWithMultTest.txt"))):
    cov = fpath.split(os.sep)[-2]
    prefix = fpath.split(os.sep)[-1].replace("-TopEffectsWithMultTest.txt", "")

    if cov not in lookup:
        lookup[cov] = {}
    if prefix not in lookup[cov]:
        lookup[cov][prefix] = row_index
        row_index += 1

    current_row_index = lookup[cov][prefix]
    print("\tProcessing covariate {} with prefix {} to row index: {}".format(cov, prefix, current_row_index))

    row = {"Cov": cov, "Prefix": prefix, "NrTestedGenes": 0, args.n_tests_col: 0}
    for column in columns:
        row[column] = 0

    dataset_labels = {}
    dataset_stats = {}
    n_datasets = 0
    indices = {}
    with gzopen(fpath, 'r') as f:
        for j, line in enumerate(f):
            values = line.rstrip("\n").split("\t")
            if j == 0:
                for index, colname in enumerate(values):
                    match = re.match("([A-Za-z]+)\(([^()]+)\)", colname)
                    if match:
                        colname = match.group(1)

                        if colname == "DatasetZScores":
                            for ds_index, dataset in enumerate(match.group(2).split(";")):
                                dataset_labels[ds_index] = dataset
                                dataset_stats[ds_index] = {}
                                dataset_stats[ds_index]["NrTestedGenes"] = 0
                                dataset_stats[ds_index][args.nom_pvalue_col] = 0
                                n_datasets += 1

                    indices[colname] = index
                continue

            for column in columns:
                if column not in indices:
                    row[column] = np.nan
                else:
                    if float(values[indices[column]]) < args.alpha:
                        row[column] += 1

            if args.n_tests_col not in indices:
                row[args.n_tests_col] = np.nan
            else:
                row[args.n_tests_col] += int(values[indices[args.n_tests_col]])

            if n_datasets > 1:
                for ds_index, zscore in enumerate(values[indices["DatasetZScores"]].split(";")):
                    if zscore == "-":
                        continue
                    pvalue = stats.norm.cdf(-abs(float(zscore))) * 2
                    dataset_stats[ds_index]["NrTestedGenes"] += 1
                    if pvalue < args.alpha:
                        dataset_stats[ds_index][args.nom_pvalue_col] += 1
            row["NrTestedGenes"] += 1
    f.close()

    if n_datasets > 1:
        for ds_index in range(n_datasets):
            # Extract and print the stats per dataset.
            dataset_label = dataset_labels[ds_index]
            dataset_n = dataset_stats[ds_index]["NrTestedGenes"]
            dataset_n_na = row["NrTestedGenes"] - dataset_n
            dataset_n_nom = dataset_stats[ds_index][args.nom_pvalue_col]
            print("\t  Dataset: {}\tNrTestedGenes: {:,} \tN-missing: {:,}\tN-nom signif. {:,} [{:.2f}%] (<{})".format(dataset_label, dataset_n, dataset_n_na, dataset_n_nom, (100 / dataset_n) * dataset_n_nom, args.minimimal_reporting_p))

            # Add number of nominal significant effects to the summary stats.
            dataset_nom_pvalue_col = dataset_label + "P"
            if dataset_nom_pvalue_col not in dataset_nom_pvalue_columns:
                dataset_nom_pvalue_columns.append(dataset_nom_pvalue_col)
            row[dataset_nom_pvalue_col] = dataset_n_nom

    data[current_row_index] = row

if len(data) == 0:
    print("Error, no Data")
    exit()

print("\nResults:")
order = ["Cov", "Prefix", "NrTestedGenes", args.n_tests_col] + dataset_nom_pvalue_columns + columns
df = pd.DataFrame(data).T.loc[:, order].sort_values(by=[args.nom_pvalue_col, "NrTestedGenes", "Cov", "Prefix"], ascending=[False, True, True, True])
print(df)

print("\nWriting output ...")
df.to_csv(args.out, sep="\t", header=True, index=False)

print("\nEND")