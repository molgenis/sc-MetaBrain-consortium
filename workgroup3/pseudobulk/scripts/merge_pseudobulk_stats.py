#!/usr/bin/env python
# Author: M. Vochteloo

import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--poolsheet", required=True, type=str, help="")
parser.add_argument("--indir", required=True, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

os.makedirs(args.out, exist_ok=True)

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

print("\nLoading poolsheet ...")
poolsheet = pd.read_csv(args.poolsheet, sep="\t")
has_dataset = "Dataset" in poolsheet.columns

print("\nLoading cell stats ...")
data = []
for _, row in poolsheet.iterrows():
    fpath = os.path.join(args.indir, row["Pool"] + ".pseudobulk.stats.tsv.gz")
    if not os.path.exists(fpath):
        print("\tWarning, {} does not exist".format(fpath))
        continue

    try:
        stats = pd.read_csv(fpath, sep="\t")
        print("\tLoaded {} with shape: {}".format(os.path.basename(fpath), stats.shape))
    except pd.errors.EmptyDataError:
        print("\tFailed to load {}: no data".format(os.path.basename(fpath)))
        continue

    stats.insert(0, "pool", str(row["Pool"]))
    if has_dataset:
        stats.insert(0, "dataset", str(row["Dataset"]))
    data.append(stats)
if len(data) == 0:
    exit()
stats = pd.concat(data, axis=0)
print("\tLoaded cell stats with shape: {}".format(stats.shape))

print("\nBarcode selection stats output:")
print(stats)

print("\nBarcode selection summary stats:")
groupby = ["cell type"]
sortby = ["ncells"]
sortorder = [False]
if has_dataset:
    groupby += ["dataset"]
    sortby += ["dataset"]
    sortorder += [False, True]
sumstats = stats.loc[:, [col for col in stats.columns if col not in ["pool", "sample", "ncells_pool", "ncells_sample"]]].groupby(groupby).sum()
sumstats.sort_values(by=sortby, ascending=sortorder, inplace=True)
sumstats = sumstats.T
sumstats["Total"] = sumstats.sum(axis=1)
print(sumstats)

print("\nSaving files")
stats.to_csv(os.path.join(args.out, "pseudobulk.stats.tsv.gz"), sep="\t", header=True, index=True, compression="gzip")
sumstats.to_csv(os.path.join(args.out, "pseudobulk.sumstats.tsv.gz"), sep="\t", header=True, index=True, compression="gzip")

print("Done")