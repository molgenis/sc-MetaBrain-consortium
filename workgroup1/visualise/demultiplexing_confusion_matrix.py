#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse

"""
Syntax: 
./demultiplexing_confusion_matrix.py -h
"""

parser = argparse.ArgumentParser(description="")
parser.add_argument("--ref", type=str, required=True, help="")
parser.add_argument("--ref_ind_coupling", type=str, required=True, help="")
parser.add_argument("--alt", type=str, required=True, help="")
parser.add_argument("--center", type=str, required=False, default=None, choices=["Broad", "NYGC"], help="")
parser.add_argument("--outdir", type=str, required=True, help="")
parser.add_argument("--extensions", nargs="+", type=str, choices=["png", "pdf", "eps"], default=["png"], help="The figure file extension.. Default: 'png'.")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import glob
import os
import re
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
    row_labels = ["{} [n={:,.0f}]".format(label, size) for label, size in row_counts]

    col_counts = list(zip(*np.unique(df[x].astype(str), return_counts=True)))
    col_counts.sort(key=lambda x: x[0])
    col_labels = ["{} [n={:,.0f}]".format(label, size) for label, size in col_counts]

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

print("Loading reference data ...")
ref_df = pd.read_csv(args.ref, sep=",", header=0, index_col=None)
ref_df["PartialscrnaSeqID"] = ref_df["barcode"].str.split("_", n=1, expand=True)[0]
print(ref_df)

ref_coupling_df = pd.read_csv(args.ref_ind_coupling, sep="\t", header=0, index_col=None)
if args.center is not None:
    ref_coupling_df["Centre"] = ref_coupling_df["id"].str.split("_", n=1, expand=True)[1]
    ref_coupling_df = ref_coupling_df.loc[ref_coupling_df["Centre"] == args.center, :]
duplicates = set(ref_coupling_df.loc[ref_coupling_df["scrnaSeqID"].duplicated(), "scrnaSeqID"].values)
print("\tExcluding {} duplicates".format(len(duplicates)))
ref_coupling_df = ref_coupling_df.loc[~ref_coupling_df["scrnaSeqID"].isin(duplicates), :]
print(ref_coupling_df)

if sum(ref_coupling_df[["PartialscrnaSeqID", "individualID"]].duplicated()) > 0:
    print("Error, reference coupling has duplicates.")
    exit()

ref_df = ref_df.merge(ref_coupling_df, on=["PartialscrnaSeqID", "individualID"], how="left")
del ref_coupling_df
print(ref_df)

if len(ref_df["barcode"].unique()) != ref_df.shape[0]:
    print("Error, reference contains duplicate barcodes")
    exit()

print("Loading alternative data ...")
if args.center is not None:
    print("\tLimiting to pools from '{}' ...".format(args.center))
alt_df_list = []
i = 0
for fpath in glob.glob(args.alt.replace("POOL", "*")):
    before, after = args.alt.split("POOL")
    pool, centre = fpath.replace(before, "").replace(after, "").split("_")
    if args.center is not None and centre != args.center:
        continue
    alt_df = pd.read_csv(fpath, sep="\t", header=0, index_col=None)
    alt_df["Pool"] = pool + "_" + centre
    alt_df["Pool_short"] = pool
    alt_df["Centre"] = centre
    alt_df["barcode"] = pool + "_" + alt_df["Barcode"]
    alt_df_list.append(alt_df)

    i += 1
    # if i > 2:
    #     break
alt_df = pd.concat(alt_df_list, axis=0)
print(alt_df)

for method in ["Demuxlet", "Souporcell", "AnyDoublet"]:
    counts = alt_df[method + "_Individual_Assignment"].value_counts()
    total = alt_df.shape[0]
    doublet = counts["doublet"] if "doublet" in counts else 0
    unassigned = counts["unassigned"] if "unassigned" in counts else 0
    assigned = total - doublet - unassigned
    print("{}\tassigned: {:,} [{:.2f}%]\tdoublet: {:,} [{:.2f}%]\tunassigned: {:,} [{:.2f}%]".format(method, assigned, (100 / total) * assigned, doublet, (100 / total) * doublet, unassigned, (100 / total) * unassigned))
exit()

if len(alt_df["barcode"].unique()) != alt_df.shape[0]:
    print("Error, alternative contains duplicate barcodes")
    exit()

# print("Creating alternative confusion matrix ...")
# confusion_df, annotation_df, tpr = create_confusion_matrix(df=alt_df, x="Demuxlet_Individual_Assignment", y="Souporcell_Individual_Assignment")
# print(confusion_df)
# print(tpr)
#
# plot_heatmap(
#     df=confusion_df,
#     annot_df=annotation_df,
#     xlabel="Demuxlet",
#     ylabel="Souporcell",
#     title="TPR: {:.3f}".format(tpr),
#     outfile="scMetaBrain_demuxlet_vs_souporcell_demultiplex_confusion_matrix" + ("_" + args.center if args.center is not None else "")
# )

print("Merging data ...")
df = ref_df.merge(alt_df, on="barcode", how="inner")
print(df)

for method in ["Demuxlet", "Souporcell", "AnyDoublet"]:
    print("Comparing with method '{}' ...".format(method))

    df["Match"] = df["wholeGenomeSeqID"] == df[method + "_Individual_Assignment"]
    counts = df["Match"].value_counts()
    assigned = df.loc[(df[method + "_Individual_Assignment"] != "doublet") & (df[method + "_Individual_Assignment"] != "unassigned"), :].shape[0]
    print("\tMatch: {:,}/{:,} [{:.2f}%]\tMatch (excl. unassign.): {:,}/{:,} [{:.2f}%]".format(counts[True], counts.sum(), (100 / counts.sum()) * counts[True], counts[True], assigned, (100 / assigned) * counts[True]))

    # print(df.loc[~df["Match"], "wholeGenomeSeqID"].value_counts())
    stats_per_pool = []
    for pool in df["Pool"].unique():
        pool_counts = df.loc[df["Pool"] == pool, "Match"].value_counts()
        pool_assigned = df.loc[(df["Pool"] == pool) & (df[method + "_Individual_Assignment"] != "doublet") & (df[method + "_Individual_Assignment"] != "unassigned"), :].shape[0]

        n_true = 0
        if True not in pool_counts:
            stats_per_pool.append([pool, np.nan, pool_counts.sum(), np.nan, np.nan, pool_assigned, np.nan])
            continue

        # print(pool_counts)
        # print(pool_assigned)
        stats_per_pool.append([pool, pool_counts[True], pool_counts.sum(), (100 / pool_counts.sum()) * pool_counts[True], pool_counts[True], pool_assigned, (100 / pool_assigned) * pool_counts[True]])
    stats_per_pool.sort(key=lambda x: -x[6])
    for stats in stats_per_pool:
        print("\tPool: {}\tMatch: {:,}/{:,} [{:.2f}%]\tMatch (excl. unassign.): {:,}/{:,} [{:.2f}%]".format(*stats))

    # print("\tCreating confusion matrix ...".format(method))
    # confusion_df, annotation_df, tpr = create_confusion_matrix(df=df, x="wholeGenomeSeqID", y=method + "_Individual_Assignment")
    # print(confusion_df)
    # print(tpr)
    #
    # plot_heatmap(
    #     df=confusion_df,
    #     annot_df=annotation_df,
    #     xlabel="Fujita",
    #     ylabel="scMetaBrain",
    #     title="TPR: {:.3f}".format(tpr),
    #     outfile="Fujita_vs_scMetaBrain_{}_demultiplex_confusion_matrix{}".format(method.lower(), "_" + args.center if args.center is not None else "")
    # )
