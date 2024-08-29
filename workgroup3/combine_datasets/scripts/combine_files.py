#!/usr/bin/env python
# Author: M. Vochteloo
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--files", required=True, type=str, nargs="*", help="")
parser.add_argument("--datasets", required=True, type=str, nargs="*", help="")
parser.add_argument("--outfile", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

if len(args.files) != len(args.datasets):
    print("Error, --files and --datasets must have the same length.")
    exit()

import pandas as pd

print("Loading data ...")
data = []
for file, dataset in zip(args.files, args.datasets):
    print("\t{}: {}".format(dataset, file))
    df = pd.read_csv(file, sep="\t", header=0, index_col=None)
    df["Dataset"] = dataset
    data.append(df)

df = pd.concat(data, axis=0)
print(df)

print("Saving data ...")
df.to_csv(args.outfile, sep="\t", header=True, index=False)

print("Done")