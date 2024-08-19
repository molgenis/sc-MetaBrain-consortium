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
    row_counts = list(zip(*np.unique(df[y].astype(str), return_counts=True)))
    row_counts.sort(key=lambda x: x[0])
    row_labels = ["{}\n[n={:,.0f}]".format(label, size) for label, size in row_counts]

    col_counts = list(zip(*np.unique(df[x].astype(str), return_counts=True)))
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

def load_data(method):
    if method in ["DoubletFinder", "scDblFinder"]:
        fpath = os.path.join(method + "Run1", method + "_doublets_singlets.tsv.gz")
        load_function = load_default
    elif method in ["demuxlet"]:
        fpath = os.path.join("popscle", "demuxlet", "demuxletOUT.best")
        load_function = load_demuxlet
    elif method in ["Souporcell"]:
        fpath = os.path.join("souporcell", "clusters.tsv.gz")
        load_function = load_souporcell
    else:
        print("Error in load_data")
        exit()

    inpaths = glob.glob(os.path.join(args.indir + "Step2-DemultiplexingAndDoubletRemoval", "*", fpath))
    if len(inpaths) == 0:
        return None

    df_list = []
    for i, fpath in enumerate(inpaths):
        # if i > 10:
        #     break
        df_list.append(load_function(fpath=fpath))

    return pd.concat(df_list, axis=0)

def load_default(fpath):
    df = pd.read_csv(fpath, sep="\t", header=0, index_col=None)
    df["pool"] = fpath.split("/")[-3]
    return df

def load_demuxlet(fpath):
    df = pd.read_csv(fpath, sep="\t", header=0, index_col=None)
    df["demuxlet_DropletType"] = df["DROPLET.TYPE"].map({"SNG": "singlet", "DBL": "doublet", "AMB": "unassigned"})
    df = df[["BARCODE", "demuxlet_DropletType"]].rename(columns={"BARCODE": "Barcode"})
    df["pool"] = fpath.split("/")[-4]
    return df

def load_souporcell(fpath):
    df = pd.read_csv(fpath, sep="\t", header=0, index_col=None)
    df = df[["barcode", "status", "singlet_posterior", "doublet_posterior"]].rename(columns={"barcode": "Barcode", "status": "Souporcell_DropletType", "singlet_posterior": "souporcell_SingletScore", "doublet_posterior": "souporcell_DoubletScore"})
    df["pool"] = fpath.split("/")[-3]
    return df


print("Loading data ...")
df = None
methods = []
for method in ["demuxlet", "Souporcell", "DoubletFinder", "scDblFinder"]:
    print("\tLoading {}".format(method))
    method_df = load_data(method=method)
    print(method_df[method + "_DropletType"].value_counts())
    if method_df is None:
        continue
    methods.append(method)
    if df is None:
        df = method_df
    else:
        df = df.merge(method_df, on=["Barcode", "pool"], how="outer")
df = df.fillna("unknown")
print(df)
print(df.loc[df["demuxlet_DropletType"] == "unknown", "pool"].value_counts())

print("\nPlotting data ...")
for i, method1 in enumerate(methods):
    for j, method2 in enumerate(methods):
        if i >= j:
            continue
        confusion_df, annotation_df, tpr = create_confusion_matrix(df=df, x=method1 + "_DropletType", y=method2 + "_DropletType")

        plot_heatmap(
            df=confusion_df,
            annot_df=annotation_df,
            xlabel=method1,
            ylabel=method2,
            title="TPR: {:.3f}".format(tpr),
            outfile=method1 + "_vs_" + method2 + "_confusion_matrix")
