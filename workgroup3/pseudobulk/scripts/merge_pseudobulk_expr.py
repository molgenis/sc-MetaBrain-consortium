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
    index = 0
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

    print("\tLoaded {} with shape: {}".format(os.path.basename(fpath), df.shape))
    return df


def load_file_full(inpath, header, index_col, sep="\t", low_memory=True,
              nrows=None, skiprows=None):
    df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                     low_memory=low_memory, nrows=nrows, skiprows=skiprows)
    print("\tLoaded dataframe: {} "
          "with shape: {}".format(os.path.basename(inpath),
                                  df.shape))
    return df


print("\nLoading poolsheet ...")
poolsheet = load_file(args.poolsheet)
has_dataset = "Dataset" in poolsheet.columns and args.split_per_dataset

print("\nLoading expression data ...")
gte_data = []
ncells_data = {}
ncells_total = 0
expr_data = []
for _, row in poolsheet.iterrows():
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
    expr_df = load_file(expr_fpath, index_col=0, must_contain=args.cell_type)
    if expr_df is None:
        continue
    if len(set(expr_df.index)) != expr_df.shape[0]:
        print("\tError, expression indices are not unique.")
        exit()
    mask = expr_df.index.str.endswith("_" + args.cell_type)
    expr_df = expr_df.loc[mask, :]
    expr_df.index = expr_df.index.str.removesuffix("_" + args.cell_type)

    stats_df = load_file_full(stats_fpath, header=0, index_col=None)
    mask = stats_df["cell type"] == args.cell_type
    stats_df = stats_df.loc[mask, :]
    stats_df.index = stats_df["sample"]
    if len(set(stats_df.index)) != stats_df.shape[0]:
        print("\tError, expression indices are not unique.")
        exit()
    stats_df = stats_df.loc[expr_df.index, :]

    # Remove samples with too few cells.
    mask = stats_df["ncells"] >= args.min_cells
    if mask.sum() != stats_df.shape[0]:
        for sample, ncells in zip(expr_df.index[~mask], stats_df["ncells"][~mask]):
            print("\tExcluding sample '{}' with {:,} cells".format(sample, ncells))
    if mask.sum() == 0:
        continue
    expr_df = expr_df.loc[mask, :]
    stats_df = stats_df.loc[mask, :]

    # Post-process.
    if has_dataset:
        updated_indices = []
        for sample in expr_df.index:
            updated_sample = str(row["Dataset"]) + "_" + sample
            updated_indices.append(updated_sample)
            gte_data.append([sample, updated_sample, row["Dataset"]])
        expr_df.index = updated_indices
        del updated_indices

        if not row["Dataset"] in ncells_data:
            ncells_data[row["Dataset"]] = 0
        ncells_data[row["Dataset"]] += stats_df["ncells"].sum()
    else:
        for sample in expr_df.index:
            gte_data.append([sample, sample, "Dataset"])

        if not "Dataset" in ncells_data:
            ncells_data["Dataset"] = 0
        ncells_data["Dataset"] += stats_df["ncells"].sum()
    ncells_total += stats_df["ncells"].sum()
    expr_data.append(expr_df)

if len(expr_data) == 0:
    print("\tError, no data.")
    exit()
df = pd.concat(expr_data, axis=0).astype(float)
print("\tLoaded expression with {:,} cells and shape: {}".format(ncells_total, df.shape))

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

print("\nPostprocessing data")
df = df.T
gte = pd.DataFrame(gte_data).drop_duplicates()
stats = pd.DataFrame(ncells_data, index=[args.cell_type])

print("\nSaving files")
df.to_csv(os.path.join(args.out, args.cell_type + ".pseudobulk.tsv.gz"), sep="\t", header=True, index=True, compression="gzip")
gte.to_csv(os.path.join(args.out, args.cell_type + ".pseudobulk.gte.txt"), sep="\t", header=False, index=False)
stats.to_csv(os.path.join(args.out, args.cell_type + ".pseudobulk.stats.tsv"), sep="\t", header=True, index=True)

if has_dataset:
    for dataset in gte.iloc[:, 2].unique():
        gte.loc[gte.iloc[:, 2] == dataset, :].to_csv(os.path.join(args.out, args.cell_type + ".pseudobulk.gte." + dataset + ".txt"), sep="\t", header=False, index=False)

    gte.iloc[:, 2] = "Dataset"
    gte.to_csv(os.path.join(args.out, args.cell_type + ".pseudobulk.gte.NoDataset.txt"), sep="\t", header=False, index=False)

print("Done")