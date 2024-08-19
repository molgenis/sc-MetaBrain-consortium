#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse

import numpy as np

"""
Syntax: 
./visualise_cell_qc.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--wg3_dir", type=str, required=True, help="")
parser.add_argument("--droplet_type", type=str, required=False, default="singlet", help="")
parser.add_argument("--cell_level", type=str, required=False, default="L1", help="")
parser.add_argument("--cell_type", type=str, required=False, default="EX", help="")
parser.add_argument("--ancestry", type=str, required=False, default="EUR", help="")
parser.add_argument("--min_umi", type=int, required=False, default=500, help="")
parser.add_argument("--max_mt", type=int, required=False, default=5, help="")
parser.add_argument("--wg1_metadata", type=str, required=False, default=None, help="")
parser.add_argument("--limit_cells", type=str, required=False, default=None, help="")
parser.add_argument("--wg0_dir", type=str, required=False, help="")
parser.add_argument("--gene_qc", type=str, required=False, default="MALAT1", help="")
parser.add_argument("--palette", type=str, required=False, help="")
parser.add_argument("--outdir", type=str, required=True, help="")
parser.add_argument("--label", type=str, required=False, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import warnings
warnings.simplefilter("ignore", UserWarning)

import os
import math
import glob
import json
from scipy import stats
import scanpy
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

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def load_metadata(in_dir):
    full_in_dir = os.path.join(in_dir, "expression_input", "pools")
    full_meta_fpaths = glob.glob(os.path.join(full_in_dir, "*.full.metadata.tsv.gz"))
    full_meta_pools = [os.path.basename(fpath).replace(".full.metadata.tsv.gz", "") for fpath in full_meta_fpaths]

    meta_fpaths = glob.glob(os.path.join(full_in_dir, "*.qc_metrics.tsv.gz"))
    meta_filename = ".qc_metrics.tsv.gz"
    if len(meta_fpaths) == 0:
        meta_fpaths = glob.glob(os.path.join(full_in_dir, "*.metadata.tsv.gz"))
        meta_filename = ".metadata.tsv.gz"
    meta_pools = [os.path.basename(fpath).replace(meta_filename, "") for fpath in meta_fpaths]

    filter_stats_fpaths = glob.glob(os.path.join(full_in_dir, "*.filter.stats.tsv"))
    filter_stats_pool = [os.path.basename(fpath).replace(".filter.stats.tsv", "") for fpath in filter_stats_fpaths]

    pools = list(set(full_meta_pools).intersection(set(meta_pools)).intersection(set(filter_stats_pool)))
    print("\tFound {:,} pools".format(len(pools)))

    data = []
    filter_stats_data = []
    n_rows = 0
    for pool in pools:
        # if pool != "EGAN00003566115":
        #     continue
        full_meta_df = pd.read_csv(os.path.join(full_in_dir, pool + ".full.metadata.tsv.gz"), sep="\t", header=0, index_col=None)
        meta_df = pd.read_csv(os.path.join(full_in_dir, pool + meta_filename), sep="\t", header=0, index_col=None)
        df = full_meta_df[["Pool", "Barcode", "DropletType", "L1", "Provided_Ancestry", "Assignment", "cell_treatment"]].merge(meta_df, how="left")
        if "Fujita2022" in full_in_dir:
            df["location"] = df["Pool"].str.split("_", n=1, expand=True)[1]
        n_rows += df.shape[0]
        del full_meta_df, meta_df

        if args.gene_qc not in df.columns and args.wg0_dir is not None and args.gene_qc is not None:
            print("Loading raw gene expression of '{}'".format(args.gene_qc))
            gene_expr_df = load_gene_expression(in_dir=args.wg0_dir, gene=args.gene_qc)
            print(gene_expr_df)

            print("\tMerging to data frame")
            df = df.merge(gene_expr_df, how="left")
            df[args.gene_qc] = df[args.gene_qc].fillna(0)

        if args.gene_qc in df.columns:
            df[args.gene_qc + "%"] = (100 / df["nCount_RNA"]) * df[args.gene_qc]

        if len(set(df["Barcode"])) != len(df["Barcode"]):
            print("Error, barcodes are not unique.")
            exit()

        data.append(df)

        filter_stats_df = pd.read_csv(os.path.join(full_in_dir, pool + ".filter.stats.tsv"), sep="\t", header=0, index_col=0)
        filter_stats_data.append(filter_stats_df)

        # break
    print("  loaded {:,} cells".format(n_rows))
    return pd.concat(data, axis=0), pd.concat(filter_stats_data, axis=1)

def load_gene_expression(in_dir, gene):
    folder = os.path.basename(os.path.normpath(in_dir))
    count_fpaths = []
    if folder == "CellRanger":
        count_fpaths = glob.glob(os.path.join(in_dir, "*", "outs", "filtered_feature_bc_matrix.h5"))
        pool_index = -3
    elif folder == "CellBender":
        count_fpaths = glob.glob(os.path.join(in_dir, "*", "cellbender_remove_background_output_filtered.h5"))
        pool_index = -2
    else:
        print("Error, unexepected folder '{}'".format(folder))
        return None

    n_count_fpaths = len(count_fpaths)
    print("\t\tFound {:,} count matrices".format(n_count_fpaths))
    if n_count_fpaths == 0:
        return None

    data = []
    for count_index, count_fpath in enumerate(count_fpaths):
        pool = count_fpath.split(os.sep)[pool_index]

        count_data = scanpy.read_10x_h5(count_fpath)
        print("\tProcessing: {}/{} Loaded matrix {} with {:,} barcodes and {:,} genes.".format(count_index, n_count_fpaths, pool, *count_data.X.shape))

        gene_mask = count_data.var_names.to_numpy() == gene
        if np.sum(gene_mask) != 1:
            df = pd.DataFrame({"barcode": count_data.obs_names, gene: np.nan})
        else:
            df = pd.DataFrame(np.hstack((np.expand_dims(count_data.obs_names, axis=1), count_data.X[:, gene_mask].toarray())), columns=["Barcode", gene])

        df["Pool"] = str(pool)
        df["Barcode"] = df["Barcode"].str.split("-", n=1, expand=True)[0] + "_" + df["Pool"]
        data.append(df)

        # if count_index > 2:
        #     break

    return pd.concat(data, axis=0)


def load_cell_qc(in_dir):
    df = pd.read_csv(os.path.join(in_dir, "expression_input", "metadata.tagged.tsv.gz"), sep="\t", header=0, index_col=None)
    return df[["Pool", "Barcode"] + [column for column in df.columns if column.endswith("tag")]]

def filter_cells(df, min_umi=None, max_mt=None, filter_stats_df=None, print_cells=True):
    if min_umi is None:
        min_umi = args.min_umi
    if max_mt is None:
        max_mt = args.max_mt

    n_barcodes = df.shape[0]
    n_seurat = np.nan
    if filter_stats_df is not None:
        n_seurat = filter_stats_df["cells >1 - features >1"]
    droplettype_mask = df["DropletType"] == "singlet"
    ncount_mask = df["nCount_RNA"].astype(float) >= min_umi
    percentmt_mask = df["percent.mt"].astype(float) <= max_mt
    tag_mask = ncount_mask & percentmt_mask
    celltype_mask = df[str(args.cell_level)] == str(args.cell_type)
    ancestry_mask = df["Provided_Ancestry"] == str(args.ancestry)
    celltreatment_mask = df["cell_treatment"] == "UT"
    mask = droplettype_mask & tag_mask & celltype_mask & ancestry_mask & celltreatment_mask

    if print_cells:
        print("Filtering barcodes:")
        print("  Input barcodes: N = {:,}".format(n_barcodes))
        print("  DropletType:")
        if "DoubletFinder_DropletType" in df.columns:
            doubletfinder_mask = df["DoubletFinder_DropletType"] == "singlet"
            print("    DoubletFinder - DropletType - Singlet: N = {:,} [{:.2f}%]".format(sum(doubletfinder_mask), (100 / n_barcodes) * sum(doubletfinder_mask)))
            scdblfinder_mask = df["scDblFinder_DropletType"] == "singlet"
            print("    scDblFinder - DropletType - Singlet: N = {:,} [{:.2f}%]".format(sum(scdblfinder_mask), (100 / n_barcodes) * sum(scdblfinder_mask)))
            print("    Agreement: N = {:,} [{:.2f}%]".format(sum(droplettype_mask), (100 / n_barcodes) * sum(droplettype_mask)))
        mask1 = droplettype_mask & celltype_mask
        print("  {} - {}: N = {:,} [{:.2f}%]".format(args.cell_level, args.cell_type, sum(mask1), (100 / n_barcodes) * sum(mask1)))

        mask2a = droplettype_mask & celltype_mask & ncount_mask
        mask2b = droplettype_mask & celltype_mask & percentmt_mask
        print("  Cell QC:")
        print("    nCount_RNA.tag - NotOutlier: N = {:,} [{:.2f}%]".format(sum(mask2a), (100 / n_barcodes) * sum(mask2a)))
        print("    percent.mt.tag - NotOutlier: N = {:,} [{:.2f}%]".format(sum(mask2b), (100 / n_barcodes) * sum(mask2b)))
        mask2 = droplettype_mask & celltype_mask & tag_mask
        print("    pass: N = {:,} [{:.2f}%]".format(sum(mask2), (100 / n_barcodes) * sum(mask2)))
        mask3 = droplettype_mask & celltype_mask & tag_mask & ancestry_mask
        print("  has genotype information: N = {:,} [{:.2f}%]".format(sum(mask3), (100 / n_barcodes) * sum(mask3)))
        print("")

        print(df.loc[mask, :])
    return mask.sum()

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

def plot(df, x="x", y="y", panels="z", type="scatter", mask=None, palette=None, xinclude=None, yinclude=None, add_mad=True, xlabel="", ylabel="", title="", filename="plot"):
    panel_values = list(df[panels].unique())
    panel_values.sort()

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
                g = sns.histplot(x=x,
                                 data=data,
                                 hue="hue",
                                 palette={False: "#000000", True: color},
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
            ax.annotate(
                'average x = {:.2f}'.format(data[x].mean()),
                xy=(0.5, 0.75),
                xycoords=ax.transAxes,
                color="#000000",
                fontsize=12,
                fontweight='bold')
            if type == "scatter":
                ax.annotate(
                    'average y = {:.2f}'.format(data[y].mean()),
                    xy=(0.5, 0.70),
                    xycoords=ax.transAxes,
                    color="#000000",
                    fontsize=12,
                    fontweight='bold')

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
    fpath = os.path.join(args.outdir, "plot", "{}.png".format(filename))
    print("Saving {}".format(os.path.basename(fpath)))
    fig.savefig(fpath)
    plt.close()

def plot_heatmap(df, annot_df, xlabel="", ylabel="", title="", filename="plot"):
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

            sns.heatmap(df, cmap=cmap, center=df.min().min(),
                        square=True, annot=annot_df, fmt='',
                        cbar=False, annot_kws={"size": 12},
                        ax=ax)

            plt.setp(ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=20, rotation=0))
            plt.setp(ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=20, rotation=0))

            ax.set_xlabel(xlabel, fontsize=14)
            ax.xaxis.set_label_position('bottom')

            ax.set_ylabel(ylabel, fontsize=14)
            ax.yaxis.set_label_position('left')

            ax.set_title(title, fontsize=20)
        else:
            ax.set_axis_off()

        col_index += 1
        if col_index > 1:
            col_index = 0
            row_index += 1

    plt.tight_layout()
    fpath = os.path.join(args.outdir, "plot", "{}.png".format(filename))
    print("Saving {}".format(os.path.basename(fpath)))
    fig.savefig(fpath)
    plt.close()

################################################

print("Loading data")
metadata_df, filter_stats_df = load_metadata(in_dir=args.wg3_dir)
print("\tMetadata matrix: {}".format(metadata_df.shape))
print(metadata_df)
print("\tFilter stats matrix: {}".format(filter_stats_df.shape))
print(filter_stats_df)

metadata_df.to_csv("visualise_cell_qc_metadata.tsv.gz", sep="\t", header=True, index=False, compression="gzip")
filter_stats_df.to_csv("visualise_cell_qc_filter_stats.tsv.gz", sep="\t", header=True, index=True, compression="gzip")
# exit()
# metadata_df = pd.read_csv("visualise_cell_qc_metadata.tsv.gz", sep="\t", header=0, index_col=None)
# filter_stats_df = pd.read_csv("visualise_cell_qc_filter_stats.tsv.gz", sep="\t", header=0, index_col=0)

if args.limit_cells is not None:
    print("Limiting cells to {}".format(args.limit_cells))
    limit_df = pd.read_csv(args.limit_cells, sep="\t", header=0, index_col=None)

    limit_df["include"] = True
    if "DropletType" in limit_df.columns:
        limit_df["include"] = False
        limit_df.loc[limit_df["DropletType"] == "singlet", "include"] = True

    metadata_df = metadata_df.merge(limit_df[["Pool", "Barcode", "include"]], how="left")
    metadata_df["include"] = metadata_df["include"].fillna(False)

    metadata_df = metadata_df.loc[metadata_df["include"], :]
    print("\tKept {:,} cells.".format(metadata_df.shape[0]))

filter_stats_df = filter_stats_df.sum(axis=1)
_ = filter_cells(df=metadata_df, filter_stats_df=filter_stats_df)

plot(
    df=metadata_df[["nCount_RNA", "percent.mt", "L1"]].dropna().copy(),
    x="nCount_RNA",
    y="percent.mt",
    panels="L1",
    palette=palette,
    xinclude=(args.min_umi, None),
    yinclude=(None, args.max_mt),
    xlabel="nCount_RNA",
    ylabel="percent.mt",
    title="",
    filename=args.label + "_ncountrna_vs_percentmt_cell_qc"
)

# umi_thresholds = range(75, 325, 25)
# mt_thresholds = range(5, 16, 1)
umi_thresholds = range(50, 550, 50)
mt_thresholds = [x / 10.0 for x in range(25, 275, 25)]
df = pd.DataFrame(np.nan, index=umi_thresholds, columns=mt_thresholds)
annot_df = pd.DataFrame(np.nan, index=umi_thresholds, columns=mt_thresholds)
for min_umi in umi_thresholds:
    for max_mt in mt_thresholds:
        n_cells = filter_cells(df=metadata_df, min_umi=min_umi, max_mt=max_mt, print_cells=False)
        df.loc[min_umi, max_mt] = n_cells
        annot_df.loc[min_umi, max_mt] = "n={:,}\n[{:.2f}%]".format(n_cells, (100 / metadata_df.shape[0]) * n_cells)

plot_heatmap(
    df=df,
    annot_df=annot_df,
    xlabel="precent.mt",
    ylabel="nCount_RNA",
    title="Barcode QC - {} - {} - {}".format(args.ancestry, args.cell_level, args.cell_type),
    filename=args.label + "_ncountrna_vs_percentmt_thresholds"
)

if args.gene_qc in metadata_df.columns:
    metadata_df["DropletTypeInclude"] = False
    metadata_df.loc[metadata_df["DropletType"] == "singlet", "DropletTypeInclude"] = True

    metadata_df["NonZeroGeneQC"] = False
    metadata_df.loc[metadata_df[args.gene_qc] != 0, "NonZeroGeneQC"] = True

    metadata_df["CellQCInclude"] = False
    metadata_df.loc[(metadata_df["nCount_RNA"] >= args.min_umi) & (metadata_df["percent.mt"] <= args.max_mt), "CellQCInclude"] = True

    plot(
        df=metadata_df[["nCount_RNA", "percent.mt", "L1", "NonZeroGeneQC"]].dropna().copy(),
        x="nCount_RNA",
        y="percent.mt",
        panels="L1",
        mask="NonZeroGeneQC",
        palette=palette,
        xinclude=(args.min_umi, None),
        yinclude=(None, args.max_mt),
        xlabel="nCount_RNA",
        ylabel="percent.mt",
        title="",
        filename=args.label + "_ncountrna_vs_percentmt_cell_qc_nonzero_" + args.gene_qc
    )

    for gene_qc_column, mask in [(args.gene_qc, "DropletTypeInclude"), (args.gene_qc + "%", "NonZeroGeneQC"), (args.gene_qc + "%", "CellQCInclude")]:
        plot(
            df=metadata_df[[gene_qc_column, "L1", mask]].dropna().copy(),
            x=gene_qc_column,
            panels="L1",
            mask=mask,
            type="hist",
            palette=palette,
            add_mad=False,
            xlabel=gene_qc_column,
            ylabel="counts",
            title="",
            filename=args.label + "_" + gene_qc_column + "_distributions_" + mask
        )

        for qc_column, include_low, include_high in [("nCount_RNA", args.min_umi, None), ("percent.mt", None, args.max_mt)]:
            plot(
                df=metadata_df[[qc_column, gene_qc_column, "L1"]].dropna().copy(),
                x=qc_column,
                y=gene_qc_column,
                panels="L1",
                palette=palette,
                xinclude=(include_low, include_high),
                yinclude=(None, None),
                add_mad=False,
                xlabel=qc_column,
                ylabel=gene_qc_column,
                title="",
                filename="{}_{}_vs_{}_cell_qc".format(args.label, qc_column, gene_qc_column)
            )


    for ct in metadata_df[args.cell_level].unique():
        if not isinstance(ct, str):
            continue
        df = pd.DataFrame(np.nan, index=umi_thresholds, columns=mt_thresholds)
        annot_df = pd.DataFrame(np.nan, index=umi_thresholds, columns=mt_thresholds)
        for min_umi in umi_thresholds:
            for max_mt in mt_thresholds:
                metadata_df["CellQCInclude"] = False
                metadata_df.loc[(metadata_df["nCount_RNA"] >= min_umi) & (metadata_df["percent.mt"] <= max_mt), "CellQCInclude"] = True

                n_cells = metadata_df.loc[(metadata_df["CellQCInclude"]) & (metadata_df[args.cell_level] == ct), :].shape[0]
                nonzero_cells = metadata_df.loc[(metadata_df["CellQCInclude"]) & (metadata_df[args.cell_level] == ct) & (metadata_df["NonZeroGeneQC"]), :].shape[0]
                prop_zero = np.nan
                if nonzero_cells > 0:
                    prop_zero = (100 / n_cells) * nonzero_cells

                df.loc[min_umi, max_mt] = prop_zero
                annot_df.loc[min_umi, max_mt] = "n={:,}\n[{:.2f}%]".format(nonzero_cells, prop_zero)

        plot_heatmap(
            df=df,
            annot_df=annot_df,
            xlabel="precent.mt",
            ylabel="nCount_RNA",
            title=">0 {} counts - {} - {} - {}".format(args.gene_qc, args.ancestry, args.cell_level, ct),
            filename=args.label + "_" + ct + "_ncountrna_vs_percentmt_nonzero_" + args.gene_qc + "_counts"
        )


print("END")