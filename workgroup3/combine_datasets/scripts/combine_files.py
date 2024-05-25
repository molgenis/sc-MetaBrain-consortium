#!/usr/bin/env python
# Author: M. Vochteloo
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--files", required=True, type=str, nargs="*", help="")
parser.add_argument("--outfile", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import pandas as pd

data = []
for file in args.files:
    data.append(pd.read_csv(file, sep="\t", header=0, index_col=None))

df = pd.concat(data, axis=0)
print(df)
df.to_csv(args.outfile, sep="\t", header=True, index=False)