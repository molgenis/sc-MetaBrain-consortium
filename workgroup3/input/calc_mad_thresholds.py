#!/usr/bin/env python
# Author: M. Vochteloo
import argparse

"""
Syntax: 
./calc_mad_thresholds.py -h
"""


parser = argparse.ArgumentParser(description="")
parser.add_argument("--wg3_indir", required=True, type=str, help="")
parser.add_argument("--prefix", required=False, type=str, default="", help="")
parser.add_argument("--cell_level", required=False, default="L1", type=str, help="")
parser.add_argument("--individual", required=False, default="Assignment", type=str, help="")
parser.add_argument("--sample", required=False, default="Assignment_Run_Lane", type=str, help="")
parser.add_argument("--features", nargs="+", required=False, default=["nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"], type=str, help="")
parser.add_argument("--mad_thresholds", nargs="+", required=False, default=[0, 1, 2, 3, 4, 5], type=int, help="")
parser.add_argument("--feature_cutoffs", nargs="+", required=False, default=["nCount_RNA<500", "percent.mt>5"], type=str, help="")
parser.add_argument("--outdir", required=True, type=str, help="")
parser.add_argument("--palette", type=str, required=False, help="")
parser.add_argument("--extensions", nargs="+", type=str, choices=["png", "pdf", "eps"], default=["png"], help="The figure file extension.. Default: 'png'.")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import json
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import glob

palette = None
if args.palette is not None:
    with open(args.palette) as f:
        palette = json.load(f)
    f.close()
    palette["NA"] = "#000000"

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def calculate_mad(data, constant=1.4826):
    median = np.median(data)
    mad = (np.median(np.abs(data - median))) * constant
    return median, mad


def threshold_to_mad(df):
    data = []
    for feature in args.features:
        for feature_cutoff in args.feature_cutoffs:
            if feature_cutoff.startswith(feature):
                median, mad = calculate_mad(data=df[feature])

                cutoff_str = feature_cutoff.replace(feature, "")
                cutoff_direction = cutoff_str[0]
                cutoff = float(cutoff_str[1:])

                if cutoff_direction == "<":
                    direction = "lower"
                    threshold = round((median - cutoff) / mad, 4)
                elif cutoff_direction == ">":
                    direction = "upper"
                    threshold = round(-(median - cutoff) / mad, 4)
                else:
                    print("Error")
                    exit()

                data.append([feature, threshold, direction, cutoff])

    return pd.DataFrame(data, columns=["feature", "threshold", "direction", "cutoff"])

def mad_to_threshold(df):
    data = []
    for feature in args.features:
        median, mad = calculate_mad(data=df[feature])
        for threshold in args.mad_thresholds + [0.3427, 0.2767]:
            lower_cutoff = median - mad * threshold
            data.append([feature, threshold, "lower", lower_cutoff])

            upper_cutoff = median + mad * threshold
            data.append([feature, threshold, "upper", upper_cutoff])

    return pd.DataFrame(data, columns=["feature", "threshold", "direction", "cutoff"])

def plot_single_swarmplot(data, x="variable", y="value", hue=None, palette=None,
                           xlabel="", ylabel="", title="", filename=""):
    sns.set(rc={'figure.figsize': (12, 9)})
    sns.set_style("ticks")
    fig, (ax1, ax2) = plt.subplots(nrows=1,
                                   ncols=2,
                                   gridspec_kw={"width_ratios": [0.99, 0.01]})
    swarmplot(fig=fig,
             ax=ax1,
             data=data,
             x=x,
             y=y,
             hue=hue,
             palette=palette,
             xlabel=xlabel,
             ylabel=ylabel,
             title=title)

    if palette is not None and hue is not None:
        handles = []
        for key, color in palette.items():
            if key in data[hue].values.tolist():
                handles.append(mpatches.Patch(color=color, label=key))
        ax2.legend(handles=handles, loc="center")
    ax2.set_axis_off()

    plt.tight_layout()
    for extension in args.extensions:
        outpath = os.path.join(args.outdir, "{}.{}".format(filename, extension))
        fig.savefig(outpath)
    plt.close()


def plot_multiple_swarmplot(data, panels, x="variable", y="value", hue=None, palette=None,
                             xlabel="", ylabel="", title="", filename=""):
    panel_values = list(data[panels].unique())
    panel_values.sort()

    npanels = len(panel_values)
    ncols = int(np.ceil(np.sqrt((npanels + 1))))
    nrows = int(np.ceil((npanels + 1) / ncols))

    sns.set_style("ticks")
    fig, axes = plt.subplots(nrows=nrows,
                             ncols=ncols,
                             sharex="none",
                             sharey="none",
                             figsize=(12 * ncols, 9 * nrows))
    sns.set(color_codes=True)

    unique_hues = set()
    row_index = 0
    col_index = 0
    for i in range(ncols * nrows):
        if nrows == 1 and ncols == 1:
            ax = axes
        elif nrows == 1 and ncols > 1:
            ax = axes[col_index]
        elif nrows > 1 and ncols == 1:
            ax = axes[row_index]
        else:
            ax = axes[row_index, col_index]

        sns.despine(fig=fig, ax=ax)

        if i < npanels:
            panel = panel_values[i]
            panel_data = data.loc[data[panels] == panel, :]

            if palette is not None and hue is not None:
                unique_hues.update(panel_data[hue].unique())

            swarmplot(fig=fig,
                       ax=ax,
                       data=panel_data,
                       x=x,
                       y=y,
                       hue=hue,
                       palette=palette,
                       xlabel=xlabel if row_index == nrows else "",
                       ylabel=ylabel if col_index == 0 else "",
                       title=title + " - " + panel)

        elif i == npanels:
            ax.set_axis_off()
            if palette is not None and hue is not None:
                handles = []
                for key, value in palette.items():
                    if key in unique_hues:
                        handles.append(mpatches.Patch(color=value, label=key))
                ax.legend(handles=handles, loc="center")
        else:
            ax.set_axis_off()

        col_index += 1
        if col_index > (ncols - 1):
            col_index = 0
            row_index += 1

    plt.tight_layout()
    for extension in args.extensions:
        outpath = os.path.join(args.outdir, "{}.{}".format(filename, extension))
        fig.savefig(outpath)
    plt.close()


def swarmplot(fig, ax, data, x="x", y="y", hue=None, palette=None, xlabel="", ylabel="", title=""):
    sns.despine(fig=fig, ax=ax)

    sns.swarmplot(x=x,
                   y=y,
                   hue=hue,
                   data=data,
                   palette=palette,
                  legend=False,
                   ax=ax)

    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, horizontalalignment="right")

    ax.set_title(title,
                 fontsize=22,
                 fontweight='bold')
    ax.set_xlabel(xlabel,
                  fontsize=14,
                  fontweight='bold')
    ax.set_ylabel(ylabel,
                  fontsize=14,
                  fontweight='bold')


def plot_single_lineplot(data, x="variable", y="value", units=None, style=None, hue=None,
                        palette=None, xlabel="", ylabel="", title="", filename=""):
    sns.set(rc={'figure.figsize': (12, 9)})
    sns.set_style("ticks")
    fig, (ax1, ax2) = plt.subplots(nrows=1,
                                   ncols=2,
                                   gridspec_kw={"width_ratios": [0.99, 0.01]})
    lineplot(fig=fig,
             ax=ax1,
             data=data,
             x=x,
             y=y,
             units=units,
             style=style,
             hue=hue,
             palette=palette,
             xlabel=xlabel,
             ylabel=ylabel,
             title=title)

    if palette is not None and hue is not None:
        handles = []
        for key, color in palette.items():
            if key in data[hue].values.tolist():
                handles.append(mpatches.Patch(color=color, label=key))
        ax2.legend(handles=handles, loc="center")
    ax2.set_axis_off()

    plt.tight_layout()
    for extension in args.extensions:
        outpath = os.path.join(args.outdir, "{}.{}".format(filename, extension))
        fig.savefig(outpath)
    plt.close()

def plot_multiple_lineplot(data, panels, x="variable", y="value", units=None, style=None, hue=None,
                           palette=None, xlabel="", ylabel="", title="", filename=""):
    panel_values = list(data[panels].unique())
    panel_values.sort()

    npanels = len(panel_values)
    ncols = int(np.ceil(np.sqrt((npanels + 1))))
    nrows = int(np.ceil((npanels + 1) / ncols))

    sns.set_style("ticks")
    fig, axes = plt.subplots(nrows=nrows,
                             ncols=ncols,
                             sharex="none",
                             sharey="none",
                             figsize=(12 * ncols, 9 * nrows))
    sns.set(color_codes=True)

    unique_hues = set()
    row_index = 0
    col_index = 0
    for i in range(ncols * nrows):
        if nrows == 1 and ncols == 1:
            ax = axes
        elif nrows == 1 and ncols > 1:
            ax = axes[col_index]
        elif nrows > 1 and ncols == 1:
            ax = axes[row_index]
        else:
            ax = axes[row_index, col_index]

        sns.despine(fig=fig, ax=ax)

        if i < npanels:
            panel = panel_values[i]
            panel_data = data.loc[data[panels] == panel, :]

            if palette is not None and hue is not None:
                unique_hues.update(panel_data[hue].unique())

            lineplot(fig=fig,
                     ax=ax,
                     data=panel_data,
                     x=x,
                     y=y,
                     units=units,
                     style=style,
                     hue=hue,
                     palette=palette,
                     xlabel=xlabel if row_index == nrows else "",
                     ylabel=ylabel if col_index == 0 else "",
                     title=title + " - " + panel)

        elif i == npanels:
            ax.set_axis_off()
            if palette is not None and hue is not None:
                handles = []
                for key, value in palette.items():
                    if key in unique_hues:
                        handles.append(mpatches.Patch(color=value, label=key))
                ax.legend(handles=handles, loc="center")
        else:
            ax.set_axis_off()

        col_index += 1
        if col_index > (ncols - 1):
            col_index = 0
            row_index += 1

    plt.tight_layout()
    for extension in args.extensions:
        outpath = os.path.join(args.outdir, "{}.{}".format(filename, extension))
        fig.savefig(outpath)
    plt.close()


def lineplot(fig, ax, data, x="x", y="y", units=None, style=None, hue=None,
             palette=None, xlabel="", ylabel="", title=""):
    sns.despine(fig=fig, ax=ax)

    g = sns.lineplot(data=data,
                     x=x,
                     y=y,
                     units=units,
                     style=style,
                     hue=hue,
                     palette=palette,
                     estimator=None,
                     legend=None,
                     ax=ax)

    ax.set_title(title,
                 fontsize=14,
                 fontweight='bold')
    ax.set_xlabel(xlabel,
                  fontsize=10,
                  fontweight='bold')
    ax.set_ylabel(ylabel,
                  fontsize=10,
                  fontweight='bold')

# print("Loading metadata...")
# mad_df_list = []
# threshold_df_list = []
# for metadata_fpath in glob.glob(os.path.join(args.wg3_indir, args.prefix + "*", "expression_input", "metadata.tsv.gz")):
#     dataset = metadata_fpath.split(os.sep)[-3].split("-")[-1]
#     print("\tProcessing {}".format(dataset))
#
#     metadata_df = pd.read_csv(metadata_fpath, sep="\t", header=0, index_col=None, low_memory=False, nrows=None)
#     # print(metadata_df)
#
#     cell_types = ["NA"] + list(metadata_df[args.cell_level].unique())
#     for cell_type in cell_types:
#         ct_metadata_df = metadata_df
#         if cell_type != "NA":
#             ct_metadata_df = metadata_df.loc[metadata_df[args.cell_level] == cell_type, :]
#
#         mad_df = threshold_to_mad(df=ct_metadata_df)
#         mad_df["dataset"] = dataset
#         mad_df["cell type"] = cell_type
#         mad_df["N"] = ct_metadata_df.shape[0]
#         mad_df_list.append(mad_df)
#
#         threshold_df = mad_to_threshold(df=ct_metadata_df)
#         threshold_df["dataset"] = dataset
#         threshold_df["cell type"] = cell_type
#         threshold_df["N"] = ct_metadata_df.shape[0]
#         threshold_df_list.append(threshold_df)
#
#     # if dataset == "Zhou2020":
#     #     break
#
# mad_df = pd.concat(mad_df_list, axis=0)
# mad_df["direction_str"] = mad_df["direction"].map({"lower": " < ", "upper": " > "})
# mad_df["selection"] = mad_df["feature"] + mad_df["direction_str"] + mad_df["cutoff"].astype(str)
# print(mad_df)
#
# mad_df.to_csv("mad.tsv", sep="\t", header=True, index=False)
mad_df = pd.read_csv("mad.tsv", sep="\t", header=0, index_col=None)
mad_df["cell type"].fillna("NA", inplace=True)
print(mad_df)

# for ct in mad_df["cell type"].unique():
#     print("Cell type: {}".format(ct))
#     table = pd.pivot_table(mad_df.loc[mad_df["cell type"] == ct, :],
#                            values='threshold',
#                            index=['feature', 'direction', 'cutoff'],
#                            columns='dataset')
#     print(table)
#     print("")
# exit()

# threshold_df = pd.concat(threshold_df_list, axis=0)
# print(threshold_df)
#
# threshold_df.to_csv("threshold.tsv", sep="\t", header=True, index=False)
# threshold_df = pd.read_csv("threshold.tsv", sep="\t", header=0, index_col=None)
# threshold_df["cell type"].fillna("NA", inplace=True)
# print(threshold_df)

print("Plotting")
plot_multiple_swarmplot(
    data=mad_df,
    panels="selection",
    x="cell type",
    y="threshold",
    hue="dataset",
    palette=palette,
    xlabel="feature cutoff",
    ylabel="MAD threshold",
    title="Selection threshold to MAD cutoff",
    filename="selection_threshold_to_mad_cutoff"
)

for feature in args.features:
    plot_single_lineplot(
        data=mad_df.loc[(mad_df["cell type"] == "NA") & (mad_df["feature"] == feature), :],
        x="threshold",
        y="cutoff",
        hue="dataset",
        palette=palette,
        style="Direction",
        xlabel="MAD threshold",
        ylabel="feature cutoff",
        title=feature + " - ALL CELLS",
        filename="mad_thresholds_" + feature
    )

    plot_single_lineplot(
        data=mad_df.loc[mad_df["feature"] == feature, :],
        x="threshold",
        y="cutoff",
        hue="cell type",
        palette=palette,
        units="study",
        style="direction",
        xlabel="MAD threshold",
        ylabel="feature cutoff",
        title=feature,
        filename="mad_thresholds_" + feature + "_per_ct"
    )

    plot_multiple_lineplot(
        data=mad_df.loc[mad_df["feature"] == feature, :],
        panels="cell type",
        x="threshold",
        y="cutoff",
        hue="dataset",
        palette=palette,
        units="dataset",
        style="direction",
        xlabel="MAD threshold",
        ylabel="feature cutoff",
        title=feature,
        filename="mad_thresholds_" + feature + "_per_ct_per_dataset"
    )
