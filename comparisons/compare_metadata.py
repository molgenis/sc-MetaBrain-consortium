#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse

"""
Syntax: 
./compare_metadata.py -h
"""

parser = argparse.ArgumentParser(description="")
parser.add_argument("--meta1", type=str, required=True, help="")
parser.add_argument("--label1", type=str, required=True, help="")
parser.add_argument("--link1", type=str, required=False, default=None, help="")
parser.add_argument("--keys1", type=str, nargs="+", required=False, default=["Barcode"], help="")
parser.add_argument("--meta2", type=str, required=True, help="")
parser.add_argument("--label2", type=str, required=True, help="")
parser.add_argument("--link2", type=str, required=False, default=None, help="")
parser.add_argument("--keys2", type=str, nargs="+", required=False, default=["Barcode"], help="")
parser.add_argument("--values", type=str, nargs="+", required=False, default=None, help="")
parser.add_argument("--outdir", type=str, required=True, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import os
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

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def load_metadata(data_fpath, link_fpath=None, keys=None):
    meta_df = pd.read_csv(data_fpath, sep="\t", header=0, index_col=None, nrows=None)
    print("\tMetadata file: {}".format(meta_df.shape))

    if link_fpath is not None:
        print("Adding link data.")
        meta_df = meta_df.merge(pd.read_csv(args.link2, sep=";", header=0, index_col=None))

    if keys is not None:
        for key_index, key in enumerate(keys):
            if key_index == 0:
                meta_df.index = meta_df[key]
            else:
                meta_df.index = meta_df.index + "_" + meta_df[key]

    if len(set(meta_df.index)) != len(meta_df.index):
        print("Error, indices in data frame are not unique.")
        exit()
    if len(set(meta_df.columns)) != len(meta_df.columns):
        print("Error, columns in data frame are not unique.")
        exit()

    return meta_df.loc[:, [column for column in meta_df.columns if column not in keys]]

def calc_confusion_matrix(ref_matrix, alt_matrix):
    ref_counts = list(zip(*np.unique(ref_matrix[:, value_index].astype(str), return_counts=True)))
    ref_counts.sort(key=lambda x: x[0])

    alt_counts = list(zip(*np.unique(alt_matrix[:, value_index].astype(str), return_counts=True)))
    alt_counts.sort(key=lambda x: x[0])

    confusion_m = np.empty((len(ref_counts), len(alt_counts)), dtype=np.float64)
    annotation_m = np.empty((len(ref_counts), len(alt_counts)), dtype=object)
    total_n_same = 0
    total_n_total = 0
    for index_index, (index_label, _) in enumerate(ref_counts):
        ref_mask = ref_matrix[:, value_index] == index_label
        n_total = np.sum(ref_mask)
        for column_index, (column_label, _) in enumerate(alt_counts):
            alt_mask = alt_matrix[:, value_index] == column_label
            n_same = np.sum(np.logical_and(ref_mask, alt_mask))
            prop_same = n_same / n_total

            confusion_m[index_index, column_index] = prop_same
            annotation_m[index_index, column_index] = "{:.2f}%\nN={:,}".format(prop_same * 100, n_same)

            total_n_same += n_same
            total_n_total += n_total

    indices = ["{}\n[n={:,.0f}]".format(label, size) for label, size in ref_counts]
    columns = ["{}\n[n={:,.0f}]".format(label, size) for label, size in alt_counts]

    confusion_df = pd.DataFrame(confusion_m, index=indices, columns=columns)
    annotation_df = pd.DataFrame(annotation_m, index=indices, columns=columns)
    return confusion_df, annotation_df


def plot_heatmap(df, annot_df, xlabel="", ylabel="", title="", outfile="plot"):
    cmap = sns.diverging_palette(246, 24, as_cmap=True)

    sns.set_style("ticks")
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
                        cbar=False, annot_kws={"size": 10},
                        ax=ax)

            plt.setp(ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=20, rotation=0))
            plt.setp(ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=20, rotation=90))

            ax.set_xlabel(xlabel, fontsize=14)
            ax.xaxis.set_label_position('top')

            ax.set_ylabel(ylabel, fontsize=14)
            ax.yaxis.set_label_position('right')

            ax.set_title(title, fontsize=40)
        else:
            ax.set_axis_off()

        col_index += 1
        if col_index > 1:
            col_index = 0
            row_index += 1

    fig.savefig(os.path.join(args.outdir, "plot", outfile + ".png"))
    plt.close()

###############################################################

print("Loading data")
meta_df1 = load_metadata(data_fpath=args.meta1, link_fpath=args.link1, keys=args.keys1)
print(meta_df1)
meta_df2 = load_metadata(data_fpath=args.meta2, link_fpath=args.link2, keys=args.keys2)
print(meta_df2)

print("Overlapping data")
overlapping_indices = list(set(meta_df1.index).intersection(set(meta_df2.index)))
print("\t{:,} indices overlapping".format(len(overlapping_indices)))
if len(overlapping_indices) == 0:
    exit()
value_columns = list(set(meta_df1.columns).intersection(set(meta_df2.columns)))
value_columns1 = value_columns
value_columns2 = value_columns
del value_columns

if args.values is not None:
    value_columns1 = []
    value_columns2 = []
    for value in args.values:
        value1 = value
        value2 = value
        if ";" in value:
            value1, value2 = value.split(";")

        if value1 not in meta_df1.columns or value2 not in meta_df2.columns:
            continue
        value_columns1.append(value1)
        value_columns2.append(value2)

if len(value_columns1) != len(value_columns2):
    exit()
n_value_columns = len(value_columns1)
if n_value_columns == 0:
    exit()
print("\tcomparing {} value columns".format(n_value_columns))
print(value_columns1)
print(value_columns2)

meta_m1 = meta_df1.loc[overlapping_indices, value_columns1].to_numpy()
meta_m2 = meta_df2.loc[overlapping_indices, value_columns2].to_numpy()
del meta_df1, meta_df2

print("Comparing data")
value_columns = [value_columns1, value_columns2]
meta_m = [meta_m1, meta_m2]
labels = [args.label1, args.label2]
for value_index in range(n_value_columns):
    for ref_index, alt_index in [(0, 1), (1, 0)]:
        ref_value = value_columns[ref_index][value_index]
        ref_matrix = meta_m[ref_index]
        ref_label = labels[ref_index]

        alt_value = value_columns[alt_index][value_index]
        alt_matrix = meta_m[alt_index]
        alt_label = labels[alt_index]

        print("\treference: {} - {}\talternative: {} - {}".format(ref_label, ref_value, alt_label, alt_value))
        confusion_df, annotation_df = calc_confusion_matrix(ref_matrix=ref_matrix, alt_matrix=alt_matrix)

        title = "{} - {}".format(ref_value, alt_value)
        if ref_value == alt_value:
            title = ref_value

        plot_heatmap(
            df=confusion_df,
            annot_df=annotation_df,
            xlabel="{} [N = {:,.0f}]".format(alt_label, alt_matrix.shape[0]),
            ylabel="{} [N = {:,.0f}]".format(ref_label, ref_matrix.shape[0]),
            title=title,
            outfile="ref_{}_{}_vs_alt_{}_{}_confusion_matrix".format(ref_label, ref_value, alt_label, alt_value)
        )

print("END")