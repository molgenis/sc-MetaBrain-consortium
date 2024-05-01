#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse


"""
Syntax: 
./plot_replication_stats.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--workdir", type=str, required=True, help="")
parser.add_argument("--discovery_name", type=str, required=True, help="")
parser.add_argument("--replication_name", type=str, required=True, help="")
parser.add_argument("--disc_cell_types", nargs="*", type=str, required=False, default=["AST", "END", "EX", "IN", "MIC", "OLI", "OPC", "PER"], help="")
parser.add_argument("--repl_cell_types", nargs="*", type=str, required=False, default=["AST", "END", "EX", "IN", "MIC", "OLI", "OPC", "PER"], help="")
parser.add_argument("--standardise", action='store_true', help="")
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

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def plot(df, title="", filename="heatmap"):
    sns.set_style("ticks")

    nrows = len(args.disc_cell_types)
    ncols = len(args.repl_cell_types)
    sns.set(rc={'figure.figsize': (3 * nrows, 1 * ncols)})
    sns.set_style("ticks")
    fig, axes = plt.subplots(nrows=1,
                             ncols=3,
                             sharex='none',
                             sharey='none')

    for col_index, (metric, vmin, center, vmax) in enumerate([("AC", -100, 0, 100),
                                                              ("pi1", -1, 0, 1),
                                                              ("Rb", -1, 0, 1)]):
        value_df, annot_df = create_pivot_table(
            df=df.loc[df["Disc significant"] & ~df["Repl significant"], :],
            index_col="replication_ct",
            column_col="discovery_ct",
            value_col=metric,
            n_digits=0 if metric == "AC" else 2
        )

        plot_heatmap(
            ax=axes[col_index],
            df=value_df,
            annot_df=annot_df,
            vmin=vmin if args.standardise else None,
            center=center if args.standardise else value_df.min().min(),
            vmax=vmax if args.standardise else None,
            xlabel=args.discovery_name + " discovery",
            ylabel=args.replication_name + " replication",
            title=metric
        )

    fig.suptitle(title,
                 fontsize=14,
                 fontweight='bold')

    plt.tight_layout()
    for extension in args.extensions:
        fig.savefig(os.path.join(args.plotdir, "{}.{}".format(filename, extension)))
    plt.close()


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

    # TODO: this might not be the most relevant info to show.
    ss_df = df.loc[df[index_col] == df[column_col], [index_col, column_col, "N"]].copy()
    if ss_df.shape[0] > 0:
        index_ss = dict(zip(ss_df[index_col], ss_df["N"]))
        columns_ss = dict(zip(ss_df[column_col], ss_df["N"]))
        for out_df in [value_df, annot_df]:
            out_df.index = ['{}\n[N={:,}]'.format(index, index_ss[index]) for index in out_df.index]
            out_df.columns = ['{}\n[N={:,}]'.format(column, columns_ss[column]) for column in out_df.columns]
    return value_df, annot_df


def plot_heatmap(ax, df, annot_df, vmin=None, vmax=None, center=None, xlabel="", ylabel="", title=""):
    n = max(df.shape[0], df.shape[1])

    sns.heatmap(df,
                vmin=vmin,
                vmax=vmax,
                cmap=sns.diverging_palette(246, 24, as_cmap=True),
                cbar=False,
                center=center,
                square=True,
                annot=annot_df,
                fmt='',
                annot_kws={"size": 1 * n, "color": "#000000"},
                ax=ax)

    plt.setp(ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=10,
                                rotation=0))
    plt.setp(ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=10,
                                rotation=90, horizontalalignment='right'))

    ax.set_xlabel("\n" + xlabel + "\n", fontsize=8)
    ax.xaxis.set_label_position('bottom')

    ax.set_ylabel("\n" + ylabel + "\n", fontsize=8)
    ax.yaxis.set_label_position('left')

    ax.set_title(title, fontsize=12)


df_list = []
for discovery_ct in args.disc_cell_types:
    for replication_ct in args.repl_cell_types:
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

filename = args.discovery_name + "_Disc_" + args.replication_name + "_Repl_ReplicationStats"
df.to_csv(os.path.join(args.workdir, filename + ".txt.gz"), sep="\t", header=True, index=False)

plot(df=df,
     title="{} vs {}".format(args.discovery_name, args.replication_name),
     filename=filename)