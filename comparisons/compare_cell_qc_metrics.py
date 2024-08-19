#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse

import numpy as np

"""
Syntax: 
./compare_cell_qc_metrics.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--wg3_dir1", type=str, required=True, help="")
parser.add_argument("--label1", type=str, required=True, help="")
parser.add_argument("--wg3_dir2", type=str, required=True, help="")
parser.add_argument("--label2", type=str, required=True, help="")
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
from scipy import stats
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)
if not os.path.exists(os.path.join(args.outdir, "plot")):
    os.makedirs(os.path.join(args.outdir, "plot"))

palette = None
if args.palette is not None:
    with open(args.palette) as f:
        palette = json.load(f)
    f.close()

palette.update({
    "neither": "#808080",
    args.label1: "#0072B2",
    args.label2: "#D55E00",
    "both": "#009E73"
})

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def load_metadata(in_dir, value_name="variable"):
    qc_metrics = ["nCount_RNA", "nFeature_RNA", "complexity", "percent.mt", "percent.rb"]

    pool_meta_fpaths = glob.glob(os.path.join(in_dir, "expression_input", "pools", "*.metadata.tsv.gz"))
    print("\tFound {:,} pools".format(len(pool_meta_fpaths)))

    data = []
    n_rows = 0
    for pool_meta_fpath in pool_meta_fpaths:
        df = pd.read_csv(pool_meta_fpath, sep="\t", header=0, index_col=None)
        n_rows += df.shape[0]
        # if df["Pool"][0] != "EGAN00003566084":
        #     continue
        if len(set(df["Barcode"])) != len(df["Barcode"]):
            print("Error, barcodes are not unique.")
            exit()
        pivot_df = df.melt(id_vars=["Pool", "IID", "Barcode"], value_vars=qc_metrics, value_name=value_name)
        data.append(pivot_df)
        # break
    print("  loaded {:,} cells".format(n_rows))
    return pd.concat(data, axis=0)



def plot_regplot(df, x="x", y="y", hue=None, palette=None, lines=None, xlabel="", ylabel="", title="", filename="plot"):
    sns.set(rc={'figure.figsize': (12, 9)})
    sns.set_style("ticks")
    fig, (ax1, ax2) = plt.subplots(nrows=1,
                                   ncols=2,
                                   gridspec_kw={"width_ratios": [0.99, 0.01]})
    sns.despine(fig=fig, ax=ax1)

    if hue is not None and palette is not None:
        df["color"] = [palette[hue_value] for hue_value in df[hue]]

    sns.regplot(x=x, y=y, data=df, ci=None,
                scatter_kws={'facecolors': df["color"],
                             'linewidths': 0},
                line_kws={"color": "#000000"},
                ax=ax1)

    if lines is not None:
        for threshold in lines:
            if threshold is None:
                continue
            ax1.axhline(threshold, ls='--', color="#808080", alpha=0.3, zorder=-1)
            ax1.axvline(threshold, ls='--', color="#808080", alpha=0.3, zorder=-1)

    xmin, xmax = ax1.get_xlim()
    ymin, ymax = ax1.get_ylim()

    min_value = max(xmin, ymin)
    max_value = min(xmax, ymax)
    ax1.plot([min_value, max_value], [min_value, max_value], ls="--", c=".3")

    ax1.set_xlabel(xlabel,
                  fontsize=14,
                  fontweight='bold')
    ax1.set_ylabel(ylabel,
                  fontsize=14,
                  fontweight='bold')
    ax1.set_title(title,
                 fontsize=18,
                 fontweight='bold')

    coef, _ = stats.pearsonr(df[x], df[y])
    ax1.annotate(
        'n = {:,}'.format(df.shape[0]),
        xy=(0.03, 0.94),
        xycoords=ax1.transAxes,
        color="#000000",
        fontsize=14,
        fontweight='bold')
    ax1.annotate(
        'r = {:.2f}'.format(coef),
        xy=(0.03, 0.90),
        xycoords=ax1.transAxes,
        color="#000000",
        fontsize=14,
        fontweight='bold')

    if hue is not None and palette is not None:
        handles = []
        for hue_value in df["hue"].unique():
            color = palette[hue_value]
            n = df.loc[df["hue"] == hue_value, :].shape[0]
            label = "{} [n={:,}]".format(hue_value, n)
            handles.append(mpatches.Patch(color=color, label=label))
        ax2.legend(handles=handles, loc="center")
    ax2.set_axis_off()

    plt.tight_layout()
    fpath = os.path.join(args.outdir, "plot", "{}_{}_{}.png".format(x, y, filename))
    print("Saving {}".format(os.path.basename(fpath)))
    fig.savefig(fpath)
    plt.close()


################################################

print("Loading data")
cell_qc_df1 = load_metadata(in_dir=args.wg3_dir1, value_name=args.label1)
print("\tCell QC metrics matrix 1: {}".format(cell_qc_df1.shape))
print(cell_qc_df1)

cell_qc_df2 = load_metadata(in_dir=args.wg3_dir2, value_name=args.label2)
print("\tCell QC metrics matrix 2: {}".format(cell_qc_df2.shape))
print(cell_qc_df2)

cell_qc_df = cell_qc_df1.merge(cell_qc_df2)
del cell_qc_df1, cell_qc_df2
print("\tOverlapping cell QC metrics matrix 2: {}".format(cell_qc_df.shape))
print(cell_qc_df)

thresholds = {"nCount_RNA": (500, None), "nFeature_RNA": (300, None), "complexity": (None, None), "percent.mt": (None, 5), "percent.rb": (None, None)}
for variable in cell_qc_df["variable"].unique():
    print("Plotting {}".format(variable))

    variable_df = cell_qc_df.loc[cell_qc_df["variable"] == variable, [args.label1, args.label2]].copy()
    lower, upper = thresholds[variable]
    if lower is None:
        lower = -np.inf
    if upper is None:
        upper = np.inf
    mask1 = (variable_df[args.label1] > lower) & (variable_df[args.label1] < upper)
    mask2 = (variable_df[args.label2] > lower) & (variable_df[args.label2] < upper)
    variable_df["hue"] = "neither"
    variable_df.loc[mask1 & mask2, "hue"] = "both"
    variable_df.loc[mask1 & ~mask2, "hue"] = args.label1
    variable_df.loc[~mask1 & mask2, "hue"] = args.label2
    # print(variable_df)
    # print(variable_df["hue"].value_counts())

    plot_regplot(
        df=variable_df,
        x=args.label1,
        y=args.label2,
        hue="hue",
        palette=palette,
        lines=thresholds[variable],
        xlabel=args.label1,
        ylabel=args.label2,
        title=variable,
        filename=variable + "_cell_qc"
    )

print("END")