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

if not os.path.isdir(args.out):
    os.makedirs(args.out, exist_ok=True)

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

print("\nLoading poolsheet ...")
poolsheet = pd.read_csv(args.poolsheet, sep="\t")

print("\nLoading cell stats ...")
data = []
for pool in poolsheet["Pool"]:
    fpath = os.path.join(args.indir, pool + ".pseudobulk.stats.tsv.gz")
    if not os.path.exists(fpath):
        print("\tWarning, {} does not exist".format(fpath))
        continue
    df = pd.read_csv(fpath, sep="\t")
    df.insert(0, "pool", str(pool))
    data.append(df)
if len(data) == 0:
    exit()
df = pd.concat(data, axis=0)
print("\tLoaded cell stats with shape: {}".format(df.shape))

print("\nSummary stats per cell type:")
sumstats_df = df.loc[:, [col for col in df.columns if col not in ["pool", "sample", "ncells_pool", "ncells_sample"]]].groupby("cell type").sum()
sumstats_df.sort_values(by="ncells", ascending=False, inplace=True)
print(sumstats_df)

print("\tSaving files")
df.to_csv(os.path.join(args.out, "pseudobulk.stats.tsv.gz"), sep="\t", header=True, index=True, compression="gzip")

print("Done")