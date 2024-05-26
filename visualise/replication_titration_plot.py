#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse


"""
Syntax: 
./replication_titration_plot.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--work_dir", type=str, required=True, help="")
parser.add_argument("--replication_dir", type=str, required=True, help="")
parser.add_argument("--discovery_name", type=str, required=True, help="")
parser.add_argument("--replication_name", type=str, required=True, help="")
parser.add_argument("--palette", type=str, required=False, default=None, help="A color palette file.")
parser.add_argument("--extensions", nargs="+", type=str, choices=["png", "pdf", "eps"], default=["png"], help="The figure file extension.. Default: 'png'.")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import os
import json
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

outdir = os.path.join(args.work_dir, 'plots', 'replication_titration_plot')
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Set the right pdf font for exporting.
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Load  the color palette.
palette = None
if args.palette is not None:
    with open(args.palette) as f:
        palette = json.load(f)
    f.close()

def plot(df, type, x="x", y="y", units=None, hue=None, palette=None, title="",
         xlabel="", ylabel="", filename="plot", info=None, outdir=None):
    # print(df[[x, y]])
    sns.set(rc={'figure.figsize': (12, 9)})
    sns.set_style("ticks")
    fig, (ax1, ax2) = plt.subplots(nrows=1,
                                   ncols=2,
                                   gridspec_kw={"width_ratios": [0.99, 0.01]})
    sns.despine(fig=fig, ax=ax1)

    if type == "lineplot":
        g = sns.lineplot(data=df,
                         x=x,
                         y=y,
                         units=units,
                         hue=hue,
                         palette=palette,
                         estimator=None,
                         legend=None,
                         ax=ax1)
    elif type == "box":
        sns.violinplot(x=x,
                       y=y,
                       hue=hue,
                       data=df,
                       palette=palette,
                       color=None if palette is not None else "#808080",
                       cut=0,
                       dodge=True,
                       ax=ax1)

        plt.setp(ax1.collections, alpha=.75)

        sns.boxplot(x=x,
                    y=y,
                    hue=hue,
                    data=df,
                    whis=np.inf,
                    color="#808080",
                    dodge=False,
                    ax=ax1)

        for patch in ax1.artists:
            r, g, b, a = patch.get_facecolor()
            patch.set_facecolor((r, g, b, .75))

        if ax1.get_legend() is not None:
            ax1.get_legend().remove()

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
    print(outpath)
    fig.savefig(outpath)
    plt.close()


df_list = []
for i in range(0, 101):
    fpath = os.path.join(
        args.work_dir,
        args.replication_dir,
        "replication_data",
        "{}discovery_{}replication".format(args.discovery_name.replace("N", str(i)), args.replication_name.replace("N", str(i))),
        "{}_Disc_{}_Repl_ReplicationStats.txt.gz".format(args.discovery_name.replace("N", str(i)), args.replication_name.replace("N", str(i))))
    if not os.path.exists(fpath):
        continue

    df = pd.read_csv(fpath, sep="\t", header=0, index_col=0)
    df["NPcs"] = i
    df["hue"] = df["discovery_ct"] + "_" + df["replication_ct"]
    df.index = df["hue"]
    df.index.name = None
    df = df.loc[df["discovery_ct"] != df["replication_ct"], :]
    df_list.append(df)
df = pd.concat(df_list, axis=0)
print(df)

for statistic in ["AC", "pi1", "Rb"]:
    plot(
        df=df.loc[df["Disc significant"] & ~df["Repl significant"], :],
        type="lineplot",
        x="NPcs",
        y=statistic,
        units="hue",
        hue="discovery_ct",
        palette=palette,
        xlabel="#PCs removed",
        ylabel=statistic,
        title="Replication between cell types",
        outdir=outdir,
        filename="replication_titration_" + statistic + "_lineplot"
    )

    plot(
        df=df.loc[df["Disc significant"] & ~df["Repl significant"], :],
        type="box",
        x="NPcs",
        y=statistic,
        xlabel="#PCs removed",
        ylabel=statistic,
        title="Replication between cell types",
        outdir=outdir,
        filename="replication_titration_" + statistic + "_boxplot"
    )