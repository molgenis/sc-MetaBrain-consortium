#!/usr/bin/env python
# Author: M. Vochteloo

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import json
import math
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--poolsheet", required=True, type=str, help="")
parser.add_argument("--indir", required=True, type=str, help="")
parser.add_argument("--barcode_qc_columns", required=True, type=str, nargs="*", help="")
parser.add_argument("--cell_level", required=False, type=str, default="L1", help="")
parser.add_argument("--tag", required=False, type=str, default=None, help="")
parser.add_argument("--palette", type=str, required=False, help="")
parser.add_argument("--data_out", required=True, type=str, help="")
parser.add_argument("--plot_out", required=True, type=str, help="")
args = parser.parse_args()

os.makedirs(args.data_out, exist_ok=True)
os.makedirs(args.plot_out, exist_ok=True)

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

palette = None
if args.palette is not None:
    with open(args.palette) as f:
        palette = json.load(f)
    f.close()


def load_barcode_qc():
    data = []
    coi = args.barcode_qc_columns + [args.cell_level, args.tag]
    for pool in poolsheet["Pool"]:
        fpath = os.path.join(args.indir, pool + ".full.metadata.tsv.gz")
        if not os.path.exists(fpath):
            print("\tWarning, {} does not exist".format(fpath))
            continue

        try:
            df = pd.read_csv(fpath, sep="\t")
            print("\tLoaded {} with shape: {}".format(os.path.basename(fpath), df.shape))
        except pd.errors.EmptyDataError:
            print("\tFailed to load {}: no data".format(os.path.basename(fpath)))
            continue

        if len(set(coi).intersection(set(df.columns))) != len(coi):
            print("Error, not all expected columns are present.")
            exit()

        pool_coi = coi.copy()
        if "Dataset" in df.columns:
            pool_coi.append("Dataset")

        barcode_qc = df.loc[:, pool_coi].copy()
        barcode_qc.insert(0, "pool", str(pool))
        barcode_qc[args.cell_level] = barcode_qc[args.cell_level].fillna("NA")
        data.append(barcode_qc)
        del df

    if len(data) == 0:
        exit()
    barcode_qc = pd.concat(data, axis=0)
    print("\tLoaded barcode QC metrics with shape: {}".format(barcode_qc.shape))
    return barcode_qc

def calculate_mad_thresholds(data, constant=1.4826, threshold=1.):
    median = np.median(data)
    lower_threshold = median - (np.median(np.abs(data - median))) * constant * threshold
    upper_threshold = median + (np.median(np.abs(data - median))) * constant * threshold
    return lower_threshold, upper_threshold

def add_mad_lines(ax, values, horizontal=False, alpha=0.3):
    outer_bound = ax.get_xlim()
    line_func = ax.axvline
    rotation = 90
    if horizontal:
        outer_bound = ax.get_ylim()
        line_func = ax.axhline
        rotation = 0

    median = np.median(values)

    line_func(median, ls='--', color="#000000", alpha=min(1., alpha * 2), zorder=-1)

    for mult in range(1, 6):
        lower_threshold, upper_threshold = calculate_mad_thresholds(data=values, threshold=mult)

        if upper_threshold < outer_bound[1]:
            line_func(upper_threshold, ls='--', color="#808080", alpha=min(1., alpha), zorder=-1)
        if lower_threshold > outer_bound[0]:
            line_func(lower_threshold, ls='--', color="#808080", alpha=min(1., alpha), zorder=-1)

def plot(df, x="x", y="y", panels="z", type="scatter", mask=None, palette=None, xinclude=None, yinclude=None, add_mad=True, xlabel=None, ylabel=None, title=None, filename="plot"):
    panel_values = list(df[panels].unique())
    panel_values.sort()

    if xlabel is None:
        xlabel = x
    if ylabel is None:
        ylabel = y
    if title is None:
        title = ""

    if y not in df:
        df[y] = df[x]

    total_n = df.shape[0]

    if mask is not None:
        df[x + "include"] = df[mask]
        df[y + "include"] = df[mask]
    else:
        df[x + "include"] = True
        df[y + "include"] = True
        for column, include in [(x, xinclude), (y, yinclude)]:
            if include is None:
                continue
            lower, upper = include
            df[column + "include"] = False
            if lower is None:
                lower = -np.inf
            if upper is None:
                upper = np.inf
            df.loc[(df[column] >= lower) & (df[column] <= upper), column + "include"] = True

    nplots = len(panel_values)
    ncols = math.ceil(np.sqrt(nplots))
    nrows = math.ceil(nplots / ncols)

    sns.set_style("ticks")
    fig, axes = plt.subplots(nrows=nrows,
                             ncols=ncols,
                             sharex='none',
                             sharey='none',
                             figsize=(6 * ncols, 6 * nrows))
    sns.set(color_codes=True)

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

        if i < nplots:
            data = df.loc[df[panels] == panel_values[i], [x, y, x + "include", y + "include"]].copy()
            if data.shape[0] == 0:
                continue

            color = "#808080"
            if palette is not None and panel_values[i] in palette:
                color = palette[panel_values[i]]

            data["hue"] = False
            data.loc[(data[x + "include"]) & (data[y + "include"]), "hue"] = True

            sns.despine(fig=fig, ax=ax)

            if (~data[x].isna() & ~data[y].isna()).sum() > 0:
                if type == "scatter":
                    sns.scatterplot(x=x,
                                    y=y,
                                    hue="hue",
                                    data=data,
                                    palette={False: "#000000", True: color},
                                    alpha=0.1,
                                    linewidth=0,
                                    legend=False,
                                    ax=ax)
                elif type == "hist":
                    # g = sns.histplot(x=x,
                    #                  data=data,
                    #                  hue="hue",
                    #                  palette={False: "#000000", True: color},
                    #                  stat="percent",
                    #                  multiple="stack",
                    #                  hue_order=[False, True],
                    #                  legend=False,
                    #                  ax=ax)
                    # g = sns.histplot(x=x,
                    #                  data=data.loc[~data["hue"], :],
                    #                  color="#000000",
                    #                  stat="percent",
                    #                  alpha=0.2,
                    #                  legend=False,
                    #                  ax=ax)
                    g = sns.histplot(x=x,
                                     data=data.loc[data["hue"], :],
                                     color=color,
                                     stat="percent",
                                     alpha=1,
                                     legend=False,
                                     ax=ax)
                else:
                    print("Unexpected plot type ''".format(type))
                    exit()

                if add_mad:
                    add_mad_lines(ax=ax, values=data[x], horizontal=False)
                    add_mad_lines(ax=ax, values=data[y], horizontal=True)

            if xinclude is not None:
                for threshold in xinclude:
                    if threshold is None:
                        continue
                    ax.axvline(threshold, ls='--', color="#b22222", alpha=0.5, zorder=-1)
            if yinclude is not None:
                for threshold in yinclude:
                    if threshold is None:
                        continue
                    ax.axhline(threshold, ls='--', color="#b22222", alpha=0.5, zorder=-1)

            # Set annotation.
            n = data.shape[0]
            n_include = sum(data["hue"])
            n_exclude = n - n_include
            ax.annotate(
                'total N = {:,} [{:.0f}%]'.format(n, (100 / total_n) * n),
                xy=(0.5, 0.90),
                xycoords=ax.transAxes,
                color="#000000",
                fontsize=12,
                fontweight='bold')
            ax.annotate(
                'N include = {:,} [{:.0f}%]'.format(n_include, (100 / n) * n_include),
                xy=(0.5, 0.85),
                xycoords=ax.transAxes,
                color=color,
                fontsize=12,
                fontweight='bold')
            ax.annotate(
                'N exclude = {:,} [{:.0f}%]'.format(n_exclude, (100 / n) * n_exclude),
                xy=(0.5, 0.80),
                xycoords=ax.transAxes,
                color="#000000",
                fontsize=12,
                fontweight='bold')
            if type == "hist":
                ax.annotate(
                    'average x = {:.2f}'.format(data.loc[data["hue"], x].mean()),
                    xy=(0.5, 0.75),
                    xycoords=ax.transAxes,
                    color=color,
                    fontsize=12,
                    fontweight='bold')
                ax.annotate(
                    'average x = {:.2f}'.format(data.loc[~data["hue"], x].mean()),
                    xy=(0.5, 0.70),
                    xycoords=ax.transAxes,
                    color="#000000",
                    fontsize=12,
                    fontweight='bold')
            if type == "scatter":
                ax.annotate(
                    'average x = {:.2f}'.format(data[x].mean()),
                    xy=(0.5, 0.75),
                    xycoords=ax.transAxes,
                    color="#000000",
                    fontsize=12,
                    fontweight='bold')
                ax.annotate(
                    'average y = {:.2f}'.format(data[y].mean()),
                    xy=(0.5, 0.70),
                    xycoords=ax.transAxes,
                    color="#000000",
                    fontsize=12,
                    fontweight='bold')

                pearson_coef = np.nan
                if data[x].std() > 0 and data[y].std() > 0:
                    pearson_coef, _ = stats.pearsonr(data[y], data[x])

                ax.annotate(
                    'total r = {:.2f}'.format(pearson_coef),
                    xy=(0.5, 0.65),
                    xycoords=ax.transAxes,
                    color="#000000",
                    fontsize=12,
                    fontweight='bold')

            ax.set_xlabel(xlabel if row_index == (nrows - 1) else "",
                          fontsize=10,
                          fontweight='bold')
            ax.set_ylabel(ylabel if col_index == 0 else "",
                          fontsize=10,
                          fontweight='bold')
            ax.set_title(panel_values[i],
                         fontsize=14,
                         fontweight='bold')
        else:
            ax.set_axis_off()

        col_index += 1
        if col_index > (ncols - 1):
            col_index = 0
            row_index += 1

    fig.suptitle(title,
                 fontsize=16,
                 color="#000000",
                 weight='bold')

    if not os.path.exists("plot"):
        os.makedirs("plot")

    plt.tight_layout()
    fpath = os.path.join(args.plot_out, "{}.png".format(filename))
    print("Saving {}".format(os.path.basename(fpath)))
    fig.savefig(fpath)
    plt.close()

#############################################

print("\nLoading poolsheet ...")
poolsheet = pd.read_csv(args.poolsheet, sep="\t")

print("\nLoading cell stats ...")
barcode_qc = load_barcode_qc()

print("\nBarcode QC stats:")
print(barcode_qc)

print("\nSaving file")
barcode_qc.to_csv(os.path.join(args.data_out, "qc_metrics.tsv.gz"), sep="\t", header=True, index=False, compression="gzip")

print("\nVisualising barcode QC ...")
barcode_qc["mask"] = barcode_qc[args.tag] == "NotOutlier"
plot(
    df=barcode_qc,
    x="nCount_RNA",
    y="percent.mt",
    mask="mask",
    panels=args.cell_level,
    palette=palette,
    filename="ncount_rna_vs_percent_mt_all_barcodes"
)

id_vars = ["pool", args.cell_level, args.tag, "mask"]
if "Dataset" in barcode_qc:
    id_vars.append("Dataset")
barcode_qcm = barcode_qc.melt(id_vars=id_vars, value_vars=args.barcode_qc_columns)

for variable in barcode_qcm["variable"].unique():
    plot(
        df=barcode_qcm.loc[barcode_qcm["variable"] == variable, :].copy(),
        x="value",
        panels=args.cell_level,
        mask="mask",
        type="hist",
        palette=palette,
        add_mad=False,
        xlabel=variable,
        ylabel="% of cells",
        title=variable,
        filename="barcode_qc_{}_stats".format(variable)
    )

for cell_type in barcode_qcm[args.cell_level].unique():
    if cell_type in palette:
        for variable in barcode_qcm["variable"].unique():
            palette[variable] = palette[cell_type]

    plot(
        df=barcode_qcm.loc[barcode_qcm[args.cell_level] == cell_type, :].copy(),
        x="value",
        panels="variable",
        mask="mask",
        type="hist",
        palette=palette,
        add_mad=False,
        xlabel=cell_type,
        ylabel="% of cells",
        title=cell_type,
        filename="barcode_qc_{}_stats".format(cell_type)
    )

print("Done")