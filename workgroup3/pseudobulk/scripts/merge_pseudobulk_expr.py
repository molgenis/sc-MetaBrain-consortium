#!/usr/bin/env python
# Author: M. Vochteloo

import numpy as np
import pandas as pd
import argparse
import gzip
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--poolsheet", required=True, type=str, help="")
parser.add_argument("--indir", required=True, type=str, help="")
parser.add_argument("--cell_type", required=True, type=str, help="")
parser.add_argument("--min_cells", required=False, type=int, default=0, help="")
parser.add_argument("--aggregate_fun", required=False, type=str, default="sum", choices=["sum", "mean"], help="")
parser.add_argument("--split_per_dataset", action="store_true", default=False, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

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
            if index != header and must_contain is not None and must_contain not in line:
                continue

            values = line.rstrip("\n").split(sep)
            if index == header:
                columns = values
                continue
            data.append(values)

    if len(data) == 0:
        print("\tFailed to load {}: no data".format(os.path.basename(fpath)))
        return None

    df = pd.DataFrame(data, columns=columns)
    if index_col is not None:
        if len(set(df.columns)) != df.shape[1]:
            print("Error, not all columns are unique.")
            exit()

        column = df.columns[index_col]
        df.index = df[column]
        df.drop(column, axis=1, inplace=True)
    return df


def load_file_full(inpath, header, index_col, sep="\t", low_memory=True,
              nrows=None, skiprows=None):
    df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                     low_memory=low_memory, nrows=nrows, skiprows=skiprows)
    return df


print("\nLoading poolsheet ...")
poolsheet = load_file(args.poolsheet)
has_dataset = "Dataset" in poolsheet.columns

print("\nLoading expression data ...")
expr_data = []
stats_data = []
pool_index = 0
n_pools = poolsheet.shape[0]
for _, row in poolsheet.iterrows():
    print("\tParsed {:,} / {:,} pools".format(pool_index, n_pools), end='\r')

    # Define input filepaths and check if they exist.
    expr_fpath = os.path.join(args.indir, row["Pool"] + ".pseudobulk.tsv.gz")
    stats_fpath = os.path.join(args.indir, row["Pool"] + ".pseudobulk.stats.tsv.gz")
    all_exist = True
    for fpath in [expr_fpath, stats_fpath]:
        if not os.path.exists(fpath):
            print("\tWarning, {} does not exist".format(fpath))
            all_exist = False
    if not all_exist:
        continue

    # Load data and pre-process data.
    ct_suffix = "_" + args.cell_type
    expr_df = load_file(expr_fpath, index_col=0, must_contain=ct_suffix)
    if expr_df is None:
        continue
    if len(set(expr_df.index)) != expr_df.shape[0]:
        print("\tError, sample indices are not unique.")
        exit()
    mask = expr_df.index.str.endswith(ct_suffix)
    expr_df = expr_df.loc[mask, :]
    expr_df.index = expr_df.index.str.removesuffix(ct_suffix)

    stats_df = load_file_full(stats_fpath, header=0, index_col=None)
    mask = stats_df["cell type"] == args.cell_type
    stats_df = stats_df.loc[mask, :]
    stats_df.index = stats_df["sample"]
    if len(set(stats_df.index)) != stats_df.shape[0]:
        print("\tError, expression indices are not unique.")
        exit()
    stats_df.drop(["sample", "cell type"], axis=1, inplace=True)
    stats_df = stats_df.loc[expr_df.index, :]
    if not expr_df.index.equals(stats_df.index):
        print("\tError, indices are not identical.")
        exit()

    # Post-process.
    dataset = "Dataset"
    if has_dataset:
        dataset = row["Dataset"]
    multiindex_df = pd.DataFrame({"Dataset": dataset, "Pool": row["Pool"], "Sample": expr_df.index})
    expr_df.index = pd.MultiIndex.from_frame(multiindex_df)
    stats_df.index = pd.MultiIndex.from_frame(multiindex_df)

    # Save
    expr_data.append(expr_df)
    stats_data.append(stats_df)
    pool_index += 1

print("\tParsed {:,} / {:,} genes".format(pool_index, n_pools))

if len(expr_data) == 0:
    print("\tError, no data.")
    exit()
expr_df = pd.concat(expr_data, axis=0).astype(float)
stats_df = pd.concat(stats_data, axis=0).astype(float)
print("\tLoaded expression with {:,} cells and shape: {}".format(stats_df["ncells"].sum(), expr_df.shape))
if not expr_df.index.equals(stats_df.index):
    print("\tError, indices are not identical.")
    exit()

print("\nAggregate on sample level ...")
expr_df = expr_df.groupby(level=["Dataset", "Sample"]).sum()
stats_df = stats_df.groupby(level=["Dataset", "Sample"]).sum()

print("\nExcluding samples ...")
dataset_pre_filter_n = dict(zip(*np.unique(stats_df.index.get_level_values("Dataset"), return_counts=True)))

# Filter on minimum number of cells
pre_filter_n = stats_df.shape[0]
stats_df = stats_df.loc[stats_df["ncells"] > args.min_cells, :]
print("\tRemoved {:,} samples with <{} cells".format(pre_filter_n - stats_df.shape[0], args.min_cells))

# Resolve duplicate samples between datasets; pick the one with the most cells.
pre_filter_n = stats_df.shape[0]
keep_samples = stats_df.groupby("Sample")["ncells"].idxmax()
expr_df = expr_df.loc[keep_samples, :]
stats_df = stats_df.loc[keep_samples, :]
print("\tRemoved {:,} duplicated samples\n".format(pre_filter_n - stats_df.shape[0]))

dataset_post_filter_n = dict(zip(*np.unique(stats_df.index.get_level_values("Dataset"), return_counts=True)))
for dataset, pre_filter_n in dataset_pre_filter_n.items():
    post_filter_n = dataset_post_filter_n[dataset]
    print("\tRemoved {:,} / {:,} samples from {}".format(pre_filter_n - post_filter_n, pre_filter_n, dataset))

# Make sure the indices still match and that the samples are now distinct.
if not expr_df.index.equals(stats_df.index):
    print("\tError, indices are not identical.")
    exit()
samples = expr_df.index.get_level_values("Sample")
if len(samples) != len(samples.unique()):
    print("\tError, sample indices are not unique.")
    exit()

print("\nExpression file:")
expr_df = expr_df.droplevel(["Dataset"]).T
expr_df.columns.name = None
print(expr_df)

print("\nSummary stats file:")
stats_df = stats_df.reset_index()
stats_df.insert(3, "Cell type", args.cell_type)
print(stats_df)

print("\nGTE file:")
gte_df = stats_df[["Sample", "Sample", "Dataset"]].copy()
print(gte_df)

print("\nSaving files")
expr_df.to_csv(os.path.join(args.out, args.cell_type + ".pseudobulk.tsv.gz"), sep="\t", header=True, index=True, compression="gzip")
stats_df.to_csv(os.path.join(args.out, args.cell_type + ".pseudobulk.stats.tsv"), sep="\t", header=True, index=True)
gte_df.to_csv(os.path.join(args.out, args.cell_type + ".pseudobulk.gte.txt"), sep="\t", header=False, index=False)

if has_dataset:
    for dataset in gte_df.iloc[:, 2].unique():
        gte_df.loc[gte_df.iloc[:, 2] == dataset, :].to_csv(os.path.join(args.out, args.cell_type + ".pseudobulk.gte." + dataset + ".txt"), sep="\t", header=False, index=False)

    gte_df.iloc[:, 2] = "Dataset"
    gte_df.to_csv(os.path.join(args.out, args.cell_type + ".pseudobulk.gte.NoDataset.txt"), sep="\t", header=False, index=False)

print("Done")