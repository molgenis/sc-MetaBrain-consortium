#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse

"""
Syntax: 
./compare_cell_annotation.py -h
"""

CELL_LEVELS = ["class", "subclass", "cluster", "cross_species_cluster"]

parser = argparse.ArgumentParser(description="")
parser.add_argument("--indir1", type=str, required=True, help="")
parser.add_argument("--label1", type=str, required=True, help="")
parser.add_argument("--indir2", type=str, required=True, help="")
parser.add_argument("--label2", type=str, required=True, help="")
parser.add_argument("--cell_level", type=str, required=False, choices=CELL_LEVELS, default="subclass", help="")
parser.add_argument("--cell_mapping", type=str, required=False, help="")
parser.add_argument("--palette", type=str, required=False, help="")
parser.add_argument("--outdir", type=str, required=True, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import os
import glob
import json
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)
if not os.path.exists(os.path.join(args.outdir, "plot")):
    os.makedirs(os.path.join(args.outdir, "plot"))

palette = None
if args.palette is not None:
    with open(args.palette) as f:
        palette = json.load(f)
    f.close()
    palette["all"] = "#000000"
    palette["DATASET"] = "#000000"

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


###############################################################

def load_data(indir, column):
    columns = ["Pool", "Barcode", column, column + ".score"]

    full_fpath = None
    single_pool_fpath = os.path.join(indir, "map", "azimuth_all.metadata.tsv.gz")
    if os.path.exists(single_pool_fpath):
        full_fpath = single_pool_fpath

    merged_fpath = os.path.join(indir, "map", "azimuth.metadata.tsv.gz")
    if os.path.exists(merged_fpath):
        full_fpath = merged_fpath

    if full_fpath is not None:
        df = pd.read_csv(merged_fpath, sep="\t", header=0, index_col=None)
        # print(df.columns.tolist())
        print("\tLoaded {:,} pools".format(len(df["Pool"].unique())))
        return df[columns]

    fpaths = glob.glob(os.path.join(indir, "map", "azimuth_*.metadata.tsv.gz"))
    print(os.path.join(indir, "map", "azimuth_*.metadata.tsv.gz"))
    print("\tFound {:,} pools".format(len(fpaths)))

    df_list = []
    for fpath in fpaths:
        df = pd.read_csv(fpath, sep="\t", header=0, index_col=None)
        df.insert(0, "Pool", os.path.basename(fpath).replace("azimuth_", "").replace(".metadata.tsv.gz", ""))
        df_list.append(df[columns])

    return pd.concat(df_list, axis=0)

def create_confusion_matrix(df, x, y):
    row_counts = list(zip(*np.unique(df[y], return_counts=True)))
    row_counts.sort(key=lambda x: x[0])
    row_labels = ["{}\n[n={:,.0f}]".format(label, size) for label, size in row_counts]

    col_counts = list(zip(*np.unique(df[x], return_counts=True)))
    col_counts.sort(key=lambda x: x[0])
    col_labels = ["{}\n[n={:,.0f}]".format(label, size) for label, size in col_counts]

    confusion_df = pd.DataFrame(np.nan, index=row_labels, columns=col_labels)
    annotation_df = pd.DataFrame("", index=row_labels, columns=col_labels)
    tp_count = 0
    tp_total = 0
    for row_label, (row_value, row_count) in zip(row_labels, row_counts):
        for col_label, (col_value, _) in zip(col_labels, col_counts):
            n_overlap = df.loc[(df[y] == row_value) & (df[x] == col_value), :].shape[0]
            frac_overlap = np.nan
            if n_overlap > 0:
                frac_overlap = n_overlap / row_count
                annotation_df.loc[row_label, col_label] = "{:.2f}\nn={:,.0f}".format(frac_overlap, n_overlap)

            if row_value == col_value:
                tp_count += n_overlap
                tp_total += row_count
            confusion_df.loc[row_label, col_label] = frac_overlap

    tpr = np.nan
    if tp_total > 0:
        tpr = tp_count / tp_total

    return confusion_df, annotation_df, tpr


def plot_heatmap(df, annot_df, xlabel="", ylabel="", title="", outfile="plot"):
    cmap = sns.diverging_palette(246, 24, as_cmap=True)

    fig, axes = plt.subplots(nrows=2,
                             ncols=2,
                             figsize=(1 * df.shape[1] + 5, 1 * df.shape[0] + 5),
                             gridspec_kw={"width_ratios": [0.2, 0.8],
                                          "height_ratios": [0.8, 0.2]})
    sns.set(color_codes=True)

    row_index = 0
    col_index = 0
    for _ in range(4):
        ax = axes[row_index, col_index]
        if row_index == 0 and col_index == 1:

            sns.heatmap(df, cmap=cmap, vmin=-1, vmax=1, center=0,
                        square=True, annot=annot_df, fmt='',
                        cbar=False, annot_kws={"size": 12},
                        ax=ax)

            plt.setp(ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=20, rotation=0))
            plt.setp(ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=20, rotation=90))

            ax.set_xlabel(xlabel, fontsize=14)
            ax.xaxis.set_label_position('bottom')

            ax.set_ylabel(ylabel, fontsize=14)
            ax.yaxis.set_label_position('left')

            ax.set_title(title, fontsize=20)
        else:
            ax.set_axis_off()

        col_index += 1
        if col_index > 1:
            col_index = 0
            row_index += 1

    fig.savefig(os.path.join(args.outdir, "plot", "{}.png".format(outfile)))
    plt.close()

###############################################################

print("Loading data ...")
extract_column = args.cell_level
plot_column = args.cell_level
cell_mapping = None
if args.cell_mapping is not None:
    cell_mapping = pd.read_csv(args.cell_mapping, sep=";", header=0, index_col=None)
    print(cell_mapping)
    if cell_mapping.shape[1] != 2:
        print("Error, expected data frame with 2 columns.")
        exit()

    if cell_mapping.columns[0].replace("predicted.", "") in CELL_LEVELS:
        extract_column = cell_mapping.columns[0]
        plot_column = cell_mapping.columns[1]
    else:
        extract_column = cell_mapping.columns[1]
        plot_column = cell_mapping.columns[0]

df1 = load_data(indir=args.indir1, column=extract_column)
print(df1)
print("  {} data frame has shape: {}".format(args.label1, df1.shape))
if cell_mapping is not None:
    df1 = df1.merge(cell_mapping, how="left")

df2 = load_data(indir=args.indir2, column=extract_column)
print(df2)
print("  {} data frame has shape: {}".format(args.label2, df2.shape))
if cell_mapping is not None:
    df2 = df2.merge(cell_mapping, how="left")

print("Overlapping data")
df = df1.merge(df2, how="inner", on=["Pool", "Barcode"], suffixes=("_" + args.label1, "_" + args.label2))
print(df)
print("  Merged data frame has shape: {}".format(df.shape))
n_agree = sum(df[plot_column + "_" + args.label1] == df[plot_column + "_" + args.label2])
print("    {:,} / {:,} times they agree about cell type assignment ({:.2f}%)".format(n_agree, df2.shape[0], (100 / df.shape[0]) * n_agree))

confusion_df, annotation_df, tpr = create_confusion_matrix(df=df, x=plot_column + "_" + args.label1, y=plot_column + "_" + args.label2)
# print(confusion_df)

plot_heatmap(
    df=confusion_df,
    annot_df=annotation_df,
    xlabel=args.label1,
    ylabel=args.label2,
    title="Azimuth\nTPR: {:.3f}  N={:,}".format(tpr, df.shape[0]),
    outfile="azimuth_{}_vs_{}_confusion_matrix".format(args.label1, args.label2))