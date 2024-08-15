#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--data", nargs="+", required=True, type=str, help="")
parser.add_argument("--axis", required=False, type=int, default=0, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import pandas as pd

print("Loading data...")
covs_list = []
for cov_fpath in args.data:
    covs_list.append(pd.read_csv(cov_fpath, sep="\t", header=0, index_col=0))

if len(covs_list) == 1:
    df = covs_list[0]
else:
    df = pd.concat(covs_list, axis=args.axis).dropna(axis=1, how="any")

print("Saving results...")
df.to_csv(args.out + ".txt.gz", sep="\t", header=True, index=True, compression="gzip")

print("Done")
