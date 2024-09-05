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
parser.add_argument("--snp_alleles_col", required=False, type=str, default="SNPAlleles", help="")
parser.add_argument("--nom_pvalue_col", required=False, type=str, default="MetaP", help="")
parser.add_argument("--n_samples_col", required=False, type=str, default="MetaPN", help="")
parser.add_argument("--n_datasets_col", required=False, type=str, default="NrDatasets", help="")
parser.add_argument("--dataset_zcore", required=False, type=str, default="DatasetZScores", help="")
parser.add_argument("--n_tests_col", required=False, type=str, default="NrTestedSNPs", help="")
parser.add_argument("--beta_adj_pvalue_col", required=False, type=str, default="BetaAdjustedMetaP", help="")
parser.add_argument("--bonf_pvalue_col", required=False, type=str, default="BonfAdjustedMetaP", help="")
parser.add_argument("--bonf_bh_fdr_col", required=False, type=str, default="BonfBHAdjustedMetaP", help="")
parser.add_argument("--bh_fdr_col", required=False, type=str, default="bh_fdr", help="")
parser.add_argument("--qvalue_col", required=False, type=str, default="qval", help="")
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

print(os.path.join(args.indir, "*/*-TopEffectsWithMulTest.txt"))

# Create groups of columns that are handled the same way.
counter_columns = [args.n_tests_col]
average_columns = [args.n_samples_col, args.n_datasets_col]
signif_columns = [args.nom_pvalue_col, args.beta_adj_pvalue_col, args.bonf_pvalue_col, args.bonf_bh_fdr_col, args.bh_fdr_col, args.qvalue_col]
all_columns = [args.snp_alleles_col] + counter_columns + average_columns + signif_columns

print("Counting number of eQTLs ...")
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

    row = {"Cov": cov, "Prefix": prefix, "NrTestedGenes": 0, "NrTopINDELS": 0}
    for column in all_columns:
        row[column] = 0

    dataset_labels = {}
    dataset_stats = {}
    n_datasets = 0
    indices = {}
    with gzopen(fpath, 'r') as f:
        for j, line in enumerate(f):
            values = line.rstrip("\n").split("\t")
            if j == 0:
                found = False
                for index, colname in enumerate(values):
                    match = re.match("([A-Za-z]+)\(([^()]+)\)([-A-Za-z]*)", colname)
                    if match:
                        colname = match.group(1) + match.group(3)

                        if colname.startswith(args.dataset_zcore):
                            if found:
                                print("Warning, multiple columns start with '{}', using first.".format(args.dataset_zcore))
                                continue
                            found = True
                            for ds_index, dataset in enumerate(match.group(2).split(";")):
                                dataset_labels[ds_index] = dataset
                                dataset_stats[ds_index] = {}
                                dataset_stats[ds_index]["NrTestedGenes"] = 0
                                dataset_stats[ds_index][args.nom_pvalue_col] = 0
                                n_datasets += 1

                    indices[colname] = index
                continue

            if args.snp_alleles_col not in indices:
                row["NrTopINDELS"] = np.nan
            else:
                snp_allele_len = len(values[indices[args.snp_alleles_col]])
                if snp_allele_len < 3:
                    print("Warning, unexpected SNP alleles length.")
                    row["NrTopINDELS"] = np.nan
                elif snp_allele_len > 3:
                    row["NrTopINDELS"] += 1

            for column in counter_columns + average_columns:
                if column not in indices:
                    row[column] = np.nan
                else:
                    row[column] += int(values[indices[column]])

            for column in signif_columns:
                if column not in indices:
                    row[column] = np.nan
                else:
                    if float(values[indices[column]]) < args.alpha:
                        row[column] += 1

            if n_datasets > 1:
                for ds_index, zscore in enumerate(values[indices[args.dataset_zcore]].split(";")):
                    if zscore == "-":
                        continue
                    pvalue = stats.norm.cdf(-abs(float(zscore))) * 2
                    dataset_stats[ds_index]["NrTestedGenes"] += 1
                    if pvalue < args.alpha:
                        dataset_stats[ds_index][args.nom_pvalue_col] += 1
            row["NrTestedGenes"] += 1
    f.close()

    for column in average_columns:
        row[column] = row[column] / row["NrTestedGenes"]

    if n_datasets > 1:
        for ds_index in range(n_datasets):
            # Extract and print the stats per dataset.
            dataset_label = dataset_labels[ds_index]
            dataset_n = dataset_stats[ds_index]["NrTestedGenes"]
            dataset_n_na = row["NrTestedGenes"] - dataset_n
            dataset_n_nom = dataset_stats[ds_index][args.nom_pvalue_col]
            print("\t  Dataset: {}\tNrTestedGenes: {:,} \tN-missing: {:,}\tN-nom signif. {:,} [{:.2f}%] (<{})".format(dataset_label, dataset_n, dataset_n_na, dataset_n_nom, (100 / dataset_n) * dataset_n_nom, args.alpha))

            # Add number of nominal significant effects to the summary stats.
            dataset_nom_pvalue_col = args.dataset_zcore + "ToP"
            if dataset_nom_pvalue_col not in dataset_nom_pvalue_columns:
                dataset_nom_pvalue_columns.append(dataset_nom_pvalue_col)
            row[dataset_nom_pvalue_col] = dataset_n_nom

    data[current_row_index] = row

if len(data) == 0:
    print("Error, no Data")
    exit()

print("\nResults:")
order = ["Cov", "Prefix", "NrTestedGenes"] + counter_columns + ["NrTopINDELS"] + average_columns + dataset_nom_pvalue_columns + signif_columns
df = pd.DataFrame(data).T.loc[:, order].sort_values(by=[args.nom_pvalue_col, "NrTestedGenes", "Cov", "Prefix"], ascending=[False, True, True, True])
print(df)

print("\nWriting output ...")
df.to_csv(args.outfile, sep="\t", header=True, index=False)

print("\nEND")