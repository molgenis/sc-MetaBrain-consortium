#!/usr/bin/env python
# Author: M. Vochteloo

import pandas as pd
import argparse
import gzip
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--poolsheet", required=True, type=str, help="")
parser.add_argument("--indir", required=True, type=str, help="")
parser.add_argument("--cell_type", required=True, type=str, help="")
parser.add_argument("--aggregate_fun", required=False, type=str, default="sum", choices=["sum", "mean"], help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

if not os.path.isdir(args.out):
    os.makedirs(args.out, exist_ok=True)

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


def load_file(fpath, sep="\t", header=0, index_col=None, must_contain=None):
    data = []
    columns = None
    with gzopen(fpath, mode="r") as fh:
        for index, line in enumerate(fh):
            if index == 0 or (index + 1) % 1e5 == 0:
                print("\tParsed {:,} lines".format(index + 1), end='\r')

            if index != header and must_contain is not None and must_contain not in line:
                continue

            values = line.rstrip("\n").split(sep)
            if index == header:
                columns = values
                continue
            data.append(values)
    print("\tParsed {:,} lines".format(index + 1), end='\n')

    df = pd.DataFrame(data, columns=columns)
    if index_col is not None:
        if len(set(df.columns)) != df.shape[1]:
            print("Error, not all columns are unique.")
            exit()

        column = df.columns[index_col]
        df.index = df[column]
        df.drop(column, axis=1, inplace=True)

    print("\tLoaded {} with shape: {}".format(os.path.basename(fpath), df.shape))
    return df

print("\nLoading poolsheet ...")
poolsheet = load_file(args.poolsheet)

print("\nLoading expression data ...")
data = []
for pool in poolsheet["Pool"]:
    fpath = os.path.join(args.indir, pool + ".pseudobulk.tsv.gz")
    if not os.path.exists(fpath):
        print("\tWarning, {} does not exist".format(fpath))
        continue
    df = load_file(fpath, index_col=0, must_contain=args.cell_type)
    df.index = str(pool) + "_" + df.index
    data.append(df)
if len(data) == 0:
    exit()
df = pd.concat(data, axis=0).astype(float)
print("\tLoaded expression with shape: {}".format(df.shape))

print("\nValidating expression data ...")
# Extra filter on cell type
indices = df.index.str.split("_", expand=True).to_frame().reset_index(drop=True)
mask = (indices.iloc[:, -1] == args.cell_type).to_numpy()
df = df.loc[mask, :]
df.index = indices.loc[mask, :].iloc[:, -2]

if len(set(df.index)) != df.shape[0]:
    print("\nAggregating expression data ...")
    if args.aggregate_fun == "sum":
        df = df.groupby(df.index).sum()
    elif args.aggregate_fun == "mean":
        df = df.groupby(df.index).mean()
    else:
        print("Error, unexpected aggregate_fun type {}".format(args.aggregate_fun))
        exit()
    print("\tAggregated expression to shape: {}".format(df.shape))

print("\tSaving files")
df.T.to_csv(os.path.join(args.out, args.cell_type + ".pseudobulk.tsv.gz"), sep="\t", header=True, index=True, compression="gzip")
pd.DataFrame({0: df.columns, 1: df.columns}).to_csv(os.path.join(args.out, args.cell_type + ".pseudobulk.smf"), sep="\t", header=False, index=False)

print("Done")