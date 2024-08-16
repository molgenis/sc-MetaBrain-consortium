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
parser.add_argument("--disc_name", type=str, required=True, help="")
parser.add_argument("--disc_cell_types", nargs="*", type=str, required=False, default=["AST", "END", "EX", "IN", "MIC", "OLI", "OPC", "PER"], help="")
parser.add_argument("--disc_settings", nargs="*", type=str, required=False, default=["Main"], help="")
parser.add_argument("--repl_name", type=str, required=True, help="")
parser.add_argument("--repl_cell_types", nargs="*", type=str, required=False, default=["AST", "END", "EX", "IN", "MIC", "OLI", "OPC", "PER"], help="")
parser.add_argument("--repl_settings", nargs="*", type=str, required=False, default=["Main"], help="")
parser.add_argument("--standardise", action='store_true', help="")
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

outfolder = args.disc_name + "discovery_" + args.repl_name + "replication/"
data_outdir = os.path.join(args.workdir, "replication_data", outfolder)
plot_outdir = os.path.join(args.workdir, "replication_plot", outfolder)
for outdir in [data_outdir, plot_outdir]:
    if not os.path.exists(outdir):
        os.makedirs(outdir)


def plot(df, title="", filename="heatmap"):
    sns.set_style("ticks")

    nrows = 1
    ncols = 4
    ncolumn_cells = len(df["discovery_id"].unique())
    nrow_cells = len(df["replication_id"].unique())
    sns.set(rc={'figure.figsize': (ncols * ncolumn_cells, nrows * nrow_cells if nrow_cells > 1 else 9)})
    sns.set_style("ticks")
    fig, axes = plt.subplots(nrows=nrows,
                             ncols=ncols,
                             sharex='none',
                             sharey='none')

    for col_index, (metric, vmin, center, vmax) in enumerate([("AC", -100, 0, 100),
                                                              ("Coef", -1, 0, 1),
                                                              ("pi1", -1, 0, 1),
                                                              ("Rb", -1, 0, 1)]):

        if ncolumn_cells > 1 and nrow_cells > 1:
            value_df, annot_df = create_pivot_table(
                df=df.loc[df["Disc significant"] & ~df["Repl significant"], :],
                index_col="replication_id",
                column_col="discovery_id",
                value_col=metric,
                n_digits=0 if metric == "AC" else 2
            )

            if args.repl_settings == ["Main"]:
                value_df.index = [index.split("_")[1] for index in value_df.index]
                annot_df.index = [index.split("_")[1] for index in annot_df.index]

            if args.disc_settings == ["Main"]:
                value_df.columns = [column.split("_")[1] for column in value_df.columns]
                annot_df.columns = [column.split("_")[1] for column in annot_df.columns]

            min_value = value_df.min().min()
            max_value = value_df.max().max()

            plot_heatmap(
                ax=axes[col_index],
                df=value_df,
                annot_df=annot_df,
                vmin=vmin if args.standardise else (min_value if min_value < 0 else None),
                center=center if args.standardise else (min_value if min_value > 0 else 0),
                vmax=vmax if args.standardise else max_value,
                xlabel=args.disc_name + " discovery",
                ylabel=args.repl_name + " replication",
                title=metric
            )
        elif (ncolumn_cells > 1 and nrow_cells == 1) or (ncolumn_cells == 1 and nrow_cells > 1):
            plot_df = df.loc[df["Disc significant"] & ~df["Repl significant"], :].copy()

            if ncolumn_cells > 1 and nrow_cells == 1:
                x = "discovery_ct"
                xlabel = args.disc_name + " discovery"
                ylabel = args.repl_cell_types[0]
            else:
                x = "replication_ct"
                xlabel = args.repl_name + " replication"
                ylabel = args.disc_cell_types[0]

            # TODO: this plot does not look great.
            plot_barplot(
                fig=fig,
                ax=axes[col_index],
                df=plot_df,
                x=x,
                y=metric,
                n="N",
                xlabel=xlabel,
                ylabel=ylabel if col_index == 0 else "",
                title=metric
            )

    fig.suptitle(title,
                 fontsize=14,
                 fontweight='bold')

    plt.tight_layout()
    for extension in args.extensions:
        fig.savefig(os.path.join(plot_outdir, "{}.{}".format(filename, extension)))
    plt.close()


def create_pivot_table(df, index_col, column_col, value_col, n_digits=2):
    index = [folder + "_" + ct for folder in args.repl_settings for ct in args.repl_cell_types if folder + "_" + ct in df[index_col].unique()]
    index.sort()
    columns = [folder + "_" + ct for folder in args.disc_settings for ct in args.disc_cell_types if folder + "_" + ct in df[column_col].unique()]
    columns.sort()

    value_df = pd.DataFrame(np.nan, index=index, columns=columns)
    annot_df = pd.DataFrame("", index=index, columns=columns)
    for index, row in df.iterrows():
        value_df.loc[row[index_col], row[column_col]] = row[value_col]
        annot_df.loc[row[index_col], row[column_col]] = "{:,}\n{} = {:.{}f}".format(row["N"], value_col, row[value_col], n_digits)

    # TODO: this might not be the most relevant info to show.
    # ss_df = df.loc[df[index_col] == df[column_col], [index_col, column_col, "N"]].copy()
    # if ss_df.shape[0] > 0:
    #     index_ss = list(zip(ss_df[index_col], ss_df["N"]))
    #     index_ss.sort(key = lambda x: -x[1])
    #     index_order = [ss[0] for ss in index_ss]
    #
    #     columns_ss = list(zip(ss_df[column_col], ss_df["N"]))
    #     columns_ss.sort(key = lambda x: -x[1])
    #     columns_order = [ss[0] for ss in index_ss]
    #
    #     value_df = value_df.loc[index_order, columns_order]
    #     annot_df = annot_df.loc[index_order, columns_order]
    #
    #     for out_df in [value_df, annot_df]:
    #         out_df.index = ['{}\n[N={:,}]'.format(index, n) for index, n in index_ss]
    #         out_df.columns = ['{}\n[N={:,}]'.format(column, n) for column, n in columns_ss]
    return value_df, annot_df


def plot_heatmap(ax, df, annot_df, vmin=None, vmax=None, center=None, xlabel="", ylabel="", title=""):
    n = min(df.shape[0], df.shape[1])

    sns.heatmap(df,
                vmin=vmin,
                vmax=vmax,
                cmap=sns.diverging_palette(246, 24, as_cmap=True),
                cbar=False,
                center=center,
                square=True,
                annot=annot_df,
                fmt='',
                annot_kws={"size": 1 * n},
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


def plot_barplot(fig, ax, df, x="x", y="y", n="n", xlabel="", ylabel="", title=""):
    sns.despine(fig=fig, ax=ax)

    df.sort_values(by=y, inplace=True)

    g = sns.barplot(x=x,
                    y=y,
                    data=df,
                    dodge=False,
                    order=df[x],
                    ax=ax)

    for i, (index, row) in enumerate(df.iterrows()):
        g.text(i,
               row[y],
               "{:.2f}\n[N={:,}]\n ".format(row[y], row[n]),
               fontsize=14,
               color='black',
               ha="center")

    ax.set_title(title,
                 fontsize=22,
                 fontweight='bold')
    ax.set_ylabel(ylabel,
                  fontsize=14,
                  fontweight='bold')
    ax.set_xlabel(xlabel,
                  fontsize=14,
                  fontweight='bold')

df_list = []
for disc_settings in args.disc_settings:
    for disc_ct in args.disc_cell_types:
        for repl_settings in args.repl_settings:
            for repl_ct in args.repl_cell_types:
                fpath = os.path.join(args.workdir, "replication_data", args.disc_name + disc_settings + "discovery_" + args.repl_name + repl_settings + "replication", args.disc_name + disc_settings  + "_" + disc_ct + "_Disc_" + args.repl_name + repl_settings + "_" + repl_ct + "_Repl_ReplicationStats.txt.gz")
                if not os.path.exists(fpath):
                    continue

                df = pd.read_csv(fpath, sep="\t", header=0, index_col=None)
                df["discovery_name"] = args.disc_name
                df["discovery_folder"] = disc_settings
                df["discovery_ct"] = disc_ct
                df["discovery_id"] = disc_settings + "_" + disc_ct
                df["replication_name"] = args.repl_name
                df["replication_folder"] = repl_settings
                df["replication_ct"] = repl_ct
                df["replication_id"] = repl_settings + "_" + repl_ct
                df_list.append(df)
df = pd.concat(df_list, axis=0)
print(df)

filename = args.disc_name + "_Disc_" + args.repl_name + "_Repl_ReplicationStats"
df.to_csv(os.path.join(data_outdir, filename + ".txt.gz"), sep="\t", header=True, index=False)

print("Replication stats:")
replication_df = df.loc[df["Disc significant"] & ~df["Repl significant"], ["N", "discovery_id", "replication_id"]].merge(
df.loc[df["Disc significant"] & df["Repl significant"], ["N", "discovery_id", "replication_id"]], on=["discovery_id", "replication_id"], suffixes=('_discovery', '_replicating')
)
for _, row in replication_df.iterrows():
    print("\t{} - {}:\t{:,} / {:,} ({:.0f}%)".format(row["discovery_id"], row["replication_id"], row["N_replicating"], row["N_discovery"], (100 / row["N_discovery"]) * row["N_replicating"]))
print("")

plot(df=df,
     title="{} vs {}".format(args.disc_name, args.repl_name),
     filename=filename)