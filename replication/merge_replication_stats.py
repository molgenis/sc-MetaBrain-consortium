#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse


"""
Syntax: 
./merge_replication_stats.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--workdir", type=str, required=True, help="")
parser.add_argument("--discovery_name", type=str, required=True, help="")
parser.add_argument("--replication_name", type=str, required=True, help="")
parser.add_argument("--cell_types", nargs="*", type=str, required=False, default=["AST", "END", "EX", "IN", "MIC", "OLI", "OPC", "PER"], help="")
parser.add_argument("--plotdir", type=str, required=True, help="")
parser.add_argument("--palette", type=str, required=False, default=None, help="A color palette file.")
parser.add_argument("--extensions", nargs="+", type=str, choices=["png", "pdf", "eps"], default=["png"], help="The figure file extension.. Default: 'png'.")
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

df_list = []
for discovery_ct in args.cell_types:
    for replication_ct in args.cell_types:
        fpath = os.path.join(args.workdir, args.discovery_name + "_" + discovery_ct + "_Disc_" + args.replication_name + "_" + replication_ct + "_Repl_ReplicationStats.txt.gz")
        if not os.path.exists(fpath):
            continue

        df = pd.read_csv(fpath, sep="\t", header=0, index_col=None)
        df["discovery_name"] = args.discovery_name
        df["discovery_ct"] = discovery_ct
        df["replication_name"] = args.replication_name
        df["replication_ct"] = replication_ct
        df_list.append(df)
df = pd.concat(df_list, axis=0)
print(df)
df.to_csv(os.path.join(args.workdir, args.discovery_name + "_Disc_" + args.replication_name + "_Repl_ReplicationStats.txt.gz"), sep="\t", header=True, index=False)


def create_pivot_table(df, index_col, column_col, value_col, n_digits=2):
    index = list(df[index_col].unique())
    index.sort()

    columns = list(df[column_col].unique())
    columns.sort()

    value_df = pd.DataFrame(np.nan, index=index, columns=columns)
    annot_df = pd.DataFrame("", index=index, columns=columns)
    for index, row in df.iterrows():
        value_df.loc[row[index_col], row[column_col]] = row[value_col]
        annot_df.loc[row[index_col], row[column_col]] = "{:,}\n{} = {:.{}f}".format(row["N"], value_col, row[value_col], n_digits)

    return value_df, annot_df

def plot(df, annot_df, vmin=None, vmax=None, center=None, xlabel="", ylabel="", title="", filename="heatmap"):
    sns.set_style("ticks")
    annot_df.fillna("", inplace=True)

    fig, ax = plt.subplots(figsize=(df.shape[1], df.shape[0]))
    sns.set(color_codes=True)

    sns.heatmap(df,
                vmin=vmin,
                vmax=vmax,
                cmap=sns.diverging_palette(246, 24, as_cmap=True),
                cbar=False,
                center=center,
                square=True,
                annot=annot_df,
                fmt='',
                annot_kws={"size": 8, "color": "#000000"},
                ax=ax)

    plt.setp(ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=20,
                                rotation=0))
    plt.setp(ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=20,
                                rotation=90))

    ax.set_xlabel(xlabel, fontsize=14)
    ax.xaxis.set_label_position('bottom')

    ax.set_ylabel(ylabel, fontsize=14)
    ax.yaxis.set_label_position('left')

    fig.suptitle(title,
                 fontsize=22,
                 fontweight='bold')

    plt.tight_layout()
    fig.savefig(os.path.join(args.plotdir, "{}.png".format(filename)))
    plt.close()

for metric in ["AC", "pi1", "Rb"]:
    value_df, annot_df = create_pivot_table(
        df=df.loc[df["Disc significant"] & ~df["Repl significant"], :],
        index_col="replication_ct",
        column_col="discovery_ct",
        value_col=metric,
        n_digits=0 if metric == "AC" else 2
    )

    plot(
        df=value_df,
        annot_df=annot_df,
        center=value_df.min().min(),
        xlabel=args.discovery_name + " discovery",
        ylabel=args.replication_name + " replication",
        title=metric,
        filename="{}_Disc_{}_Repl_{}".format(args.discovery_name, args.replication_name, metric)
    )


