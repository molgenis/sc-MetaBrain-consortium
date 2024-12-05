#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse


"""
Syntax: 
./plot_pc_titration_curve.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--work_dir", type=str, required=True, help="")
parser.add_argument("--wg3_dir", type=str, required=True, help="")
parser.add_argument("--ancestry", type=str, required=True, help="")
parser.add_argument("--cell_level", type=str, required=True, help="")
parser.add_argument("--palette", type=str, required=False, default=None, help="A color palette file.")
parser.add_argument("--extensions", nargs="+", type=str, choices=["png", "pdf", "eps"], default=["png"], help="The figure file extension.. Default: 'png'.")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import json
import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

outdir = os.path.join(args.work_dir, 'plots', 'pc_titration_curve')
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Load  the color palette.
palette = None
if args.palette is not None:
    with open(args.palette) as f:
        palette = json.load(f)
    f.close()

# Set the right pdf font for exporting.
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def lineplot(df, x="x", y="y", units=None, hue=None, palette=None, title="",
             xlabel="", ylabel="", filename="plot", info=None, outdir=None):
    sns.set(rc={'figure.figsize': (12, 9)})
    sns.set_style("ticks")
    fig, (ax1, ax2) = plt.subplots(nrows=1,
                                   ncols=2,
                                   gridspec_kw={"width_ratios": [0.99, 0.01]})
    sns.despine(fig=fig, ax=ax1)

    g = sns.lineplot(data=df,
                     x=x,
                     y=y,
                     units=units,
                     hue=hue,
                     palette=palette,
                     estimator=None,
                     legend=None,
                     ax=ax1)

    ax1.set_title(title,
                  fontsize=14,
                  fontweight='bold')
    ax1.set_xlabel(xlabel,
                   fontsize=10,
                   fontweight='bold')
    ax1.set_ylabel(ylabel,
                   fontsize=10,
                   fontweight='bold')

    if palette is not None:
        handles = []
        for key, color in palette.items():
            if key in df[hue].values.tolist():
                label = key
                if info is not None and key in info:
                    label = "{} [{}]".format(key, info[key])
                handles.append(mpatches.Patch(color=color, label=label))
        ax2.legend(handles=handles, loc="center")
    ax2.set_axis_off()

    plt.tight_layout()
    outpath = "{}.png".format(filename)
    if outdir is not None:
        outpath = os.path.join(outdir, outpath)
    fig.savefig(outpath)
    plt.close()

indir = os.path.join(args.wg3_dir, "expression_input", args.ancestry, args.cell_level)

print("Loading data")
df_list = []
for subfolder in glob.glob(os.path.join(indir, "*")):
    cell_type = os.path.basename(subfolder)
    fpath = os.path.join(indir, cell_type, "manual_selection", cell_type + "_mbqtl_n_pc_selection.txt")
    if not os.path.isfile(fpath):
        print("Error, skipping '{}'".format(fpath))
        continue

    df = pd.read_csv(fpath, sep="\t", header=0, index_col=None)
    df["ancestry"] = args.ancestry
    df["cell_level"] = args.cell_level
    df["cell_type"] = cell_type
    df["%repl"] = (100 / df["N Genes"]) * df["N eQTLs"]
    df_list.append(df)
df = pd.concat(df_list, axis=0)
df["VAR removed"] = df["VAR removed"] * 100
print(df)

print("Plotting data")
lineplot(
    df=df,
    x="N PCs",
    y="%repl",
    hue="cell_type",
    palette=palette,
    xlabel="#PCs removed",
    ylabel="% repl. eQTLs",
    filename=os.path.basename(os.path.normpath(args.wg3_dir)) + "_pc_titration_curve",
    outdir=outdir)

print("Loading data")
df_list = []
for subfolder in glob.glob(os.path.join(indir, "*")):
    cell_type = os.path.basename(subfolder)
    fpath = os.path.join(indir, cell_type, "PostQC", cell_type + ".qtlInput.Pcs.var.txt")
    if not os.path.isfile(fpath):
        print("Error, skipping '{}'".format(fpath))
        continue

    df = pd.read_csv(fpath, sep="\t", header=0, index_col=None)
    df.columns = ["N PCs", "VAR"]
    df["ancestry"] = args.ancestry
    df["cell_level"] = args.cell_level
    df["cell_type"] = cell_type
    df["VAR removed"] = df["VAR"].cumsum()
    df_list.append(df)
df = pd.concat(df_list, axis=0)
df["VAR removed"] = df["VAR removed"] * 100
print(df)

lineplot(
    df=df,
    x="N PCs",
    y="VAR removed",
    hue="cell_type",
    palette=palette,
    xlabel="#PCs removed",
    ylabel="%VAR removed",
    filename=os.path.basename(os.path.normpath(args.wg3_dir)) + "_var_removed_curve",
    outdir=outdir)

print("Done")