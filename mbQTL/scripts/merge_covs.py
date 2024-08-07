#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--covs", nargs="+", required=True, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import pandas as pd

print("Loading data...")
covs_list = []
for cov_fpath in args.covs:
    covs_list.append(pd.read_csv(cov_fpath, sep="\t", header=0, index_col=0))
df = pd.concat(covs_list).dropna(axis=1, how="any")

print("Saving results...")
df.to_csv(args.out + ".txt.gz", sep="\t", header=True, index=True, compression="gzip")

print("Done")
