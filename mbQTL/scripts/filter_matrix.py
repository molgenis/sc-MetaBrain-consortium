#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import pandas as pd
import gzip
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--data", required=True, type=str, help="")
parser.add_argument("--indices", required=False, type=str, default=None, help="")
parser.add_argument("--columns", required=False, type=str, default=None, help="")
parser.add_argument("--head", required=False, type=int, default=None, help="")
parser.add_argument("--tail", required=False, type=int, default=None, help="")
parser.add_argument("--outfile", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

os.makedirs(os.path.dirname(args.outfile), exist_ok=True)

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)

def load_file(fpath):
    values = []
    with gzopen(fpath, 'r') as f:
        for line in f:
            values.append(line.rstrip("\n"))
    f.close()
    values = list(set(values))
    values.sort()

    return values

print("Filter data...")
df = pd.read_csv(args.data, sep="\t", header=0, index_col=0)
print("\tLoaded dataframe: {} with shape: {}".format(os.path.basename(args.data), df.shape))

if args.indices is not None:
    indices = load_file(args.indices)
    df = df.loc[[index for index in indices if index in df.index], :]

if args.columns is not None:
    columns = load_file(args.columns)
    df = df.loc[:, [column for column in columns if column in df.columns]]

if args.head is not None:
    if args.head >= df.shape[0]:
        print("\tWarning, requesting more or all rows than are present in the input file.")
    df = df.head(args.head)

if args.tail is not None:
    if args.tail >= df.shape[0]:
        print("\tWarning, requesting more or all rows than are present in the input file.")
    df = df.tail(args.tail)

df.to_csv(args.outfile, sep="\t", header=True, index=True)
print("\tSaved dataframe: {} with shape: {}".format(os.path.basename(args.outfile), df.shape))

print("\nDone")
