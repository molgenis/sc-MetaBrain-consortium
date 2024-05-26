#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse

"""
Syntax: 
./dd_confusion_matrix.py -h
"""

parser = argparse.ArgumentParser(description="")
parser.add_argument("--indir", type=str, required=True, help="")
parser.add_argument("--outdir", type=str, required=True, help="")
parser.add_argument("--extensions", nargs="+", type=str, choices=["png", "pdf", "eps"], default=["png"], help="The figure file extension.. Default: 'png'.")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import os
import glob
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

###############################################################


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

    for extension in args.extensions:
        fig.savefig(os.path.join(args.outdir, "{}.{}".format(outfile, extension)))
    plt.close()

###############################################################

def load_data(folder):
    df_list = []
    for fpath in glob.glob("/groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/*/Step2-DemultiplexingAndDoubletRemoval/*/" + folder + "Run1/" + folder + "_doublets_singlets.tsv.gz"):
        df = pd.read_csv(fpath, sep="\t", header=0, index_col=None)
        df["study"] = fpath.split("/")[7]
        df["pool"] = fpath.split("/")[9]
        df_list.append(df)

    return pd.concat(df_list, axis=0)

print("Loading data ...")
doubletfinder_df = load_data(folder="DoubletFinder")
scblfinder_df = load_data(folder="scDblFinder")
df = doubletfinder_df.merge(scblfinder_df)
print(df)
print(df["study"].value_counts())

confusion_df, annotation_df, tpr = create_confusion_matrix(df=df, x="DoubletFinder_DropletType", y="scDblFinder_DropletType")
print(confusion_df)

plot_heatmap(
    df=confusion_df,
    annot_df=annotation_df,
    xlabel="DoubletFinder",
    ylabel="scDblFinder",
    title="TPR: {:.3f}".format(tpr),
    outfile="DoubletFinder_vs_scDblFinder_confusion_matrix")
