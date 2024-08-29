#!/usr/bin/env python
# Author: M. Vochteloo
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--data", required=True, type=str, help="")
parser.add_argument("--nom_pvalue_column", required=False, type=str, default="pvalue", help="")
parser.add_argument("--perm_pvalue_column", required=False, type=str, default="perm-pvalue", help="")
parser.add_argument("--qvalue_column", required=False, type=str, default="qvalue", help="")
parser.add_argument("--minimimal_reporting_p", required=False, type=float, default=0.05, help="")
parser.add_argument("--out", required=True, type=str, help="The output file where results will be saved.")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import numpy as np
import pandas as pd
from scipy import stats
import gzip
import glob
import os
import re

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


print("Counting number of eQTLs ...")
data = {}
dataset_nom_pvalue_columns = []
for i, fpath in enumerate(glob.glob(args.data)):
    cov = fpath.split(os.sep)[-2]
    n = 0
    n_nom = 0
    n_datasets = 0
    dataset_labels = {}
    dataset_stats = {}
    n_perm = 0
    n_qval = 0
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
                                dataset_stats[ds_index]["N-effects"] = 0
                                dataset_stats[ds_index][args.nom_pvalue_column] = 0
                                n_datasets += 1

                    indices[colname] = index

                if args.nom_pvalue_column not in indices:
                    n_nom = np.nan
                if args.perm_pvalue_column not in indices:
                    n_perm = np.nan
                if args.qvalue_column not in indices:
                    n_qval = np.nan
                continue

            if args.nom_pvalue_column in indices and float(values[indices[args.nom_pvalue_column]]) < args.minimimal_reporting_p:
                n_nom += 1
            if args.perm_pvalue_column in indices and float(values[indices[args.perm_pvalue_column]]) < args.minimimal_reporting_p:
                n_perm += 1
            if args.qvalue_column in indices and float(values[indices[args.qvalue_column]]) < args.minimimal_reporting_p:
                n_qval += 1
            if n_datasets > 1:
                for ds_index, zscore in enumerate(values[indices["DatasetZScores"]].split(";")):
                    if zscore == "-":
                        continue
                    pvalue = stats.norm.cdf(-abs(float(zscore))) * 2
                    dataset_stats[ds_index]["N-effects"] += 1
                    if pvalue < args.minimimal_reporting_p:
                        dataset_stats[ds_index][args.nom_pvalue_column] += 1
            n += 1
    f.close()

    row = {
        "Cov": cov,
        "N-effects": n,
        args.nom_pvalue_column: n_nom,
        args.perm_pvalue_column: n_perm,
        args.qvalue_column: n_qval
    }
    if n_datasets > 1:
        for ds_index in range(n_datasets):
            # Extract and print the stats per dataset.
            dataset_label = dataset_labels[ds_index]
            dataset_n = dataset_stats[ds_index]["N-effects"]
            dataset_n_na = n - dataset_n
            dataset_n_nom = dataset_stats[ds_index][args.nom_pvalue_column]
            print("\tDataset: {}\tN-effects: {:,} \tN-missing: {:,}\tN-nom signif. {:,} [{:.2f}%] (<{})".format(dataset_label, dataset_n, dataset_n_na, dataset_n_nom, (100 / dataset_n) * dataset_n_nom, args.minimimal_reporting_p))

            # Add number of nominal significant effects to the summary stats.
            dataset_nom_pvalue_column = dataset_label + "P"
            if dataset_nom_pvalue_column not in dataset_nom_pvalue_columns:
                dataset_nom_pvalue_columns.append(dataset_nom_pvalue_column)
            row[dataset_nom_pvalue_column] = dataset_n_nom
    data[i] = row

print("\nResults:")
order = ["Cov", "N-effects"] + dataset_nom_pvalue_columns + [args.nom_pvalue_column, args.perm_pvalue_column, args.qvalue_column]
df = pd.DataFrame(data).T.loc[:, order].sort_values(by=[args.nom_pvalue_column, "N-effects", "Cov"], ascending=[False, True, True])
print(df)

print("\nWriting output ...")
df.to_csv(args.out, sep="\t", header=True, index=False)

print("\nEND")