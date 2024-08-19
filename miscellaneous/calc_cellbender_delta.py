#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse

"""
Syntax: 
./calc_cellbender_delta.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--wg0_dir", type=str, required=True, help="")
parser.add_argument("--wg1_dir", type=str, required=True, help="")
parser.add_argument("--wg2_dir", type=str, required=True, help="")
parser.add_argument("--datasets", nargs="+", type=str, required=False, default=["Cain2023", "Mathys2019", "RocheAD2022", "RocheColumbia2022", "RocheMS2022", "Zhou2020"], help="")
parser.add_argument("--wg2_pairing", type=str, required=False, help="")
parser.add_argument("--cell_level", type=str, required=False, default="L1", help="")
parser.add_argument("--cell_types", nargs="+", type=str, required=False, default=["AST", "END", "EX", "IN", "MIC", "OLI", "OPC", "PER"], help="")
parser.add_argument("--genes", type=str, required=False, help="")
parser.add_argument("--palette", type=str, required=False, help="")
parser.add_argument("--out_dir", type=str, required=True, help="")
parser.add_argument("--force", action='store_true', help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import warnings
warnings.simplefilter("ignore", UserWarning)

import os
import glob
import math
import json
import numpy as np
import pandas as pd
import scanpy
import gzip
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adjustText import adjust_text
import matplotlib.patches as mpatches

palette = None
if args.palette is not None:
    with open(args.palette) as f:
        palette = json.load(f)
    f.close()

palette.update({"True": "#009E73", "False": "#D55E00", True: "#009E73", False: "#D55E00"})

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + "t")
    else:
        return open(file, mode)

def load_coo_df(fpath):
    data = scanpy.read_10x_h5(fpath)
    barcodes = {index: barcode for index, barcode in enumerate(data.obs_names)}

    var_names_df = data.var['gene_ids'].to_frame().reset_index(names="gene_symbol")
    var_names_df["gene"] = var_names_df["gene_ids"] + "_" + var_names_df["gene_symbol"]
    if (data.var_names.to_numpy() != var_names_df["gene_symbol"]).all():
        print("Error in merging genes metadata.")
        exit()
    genes = {index: barcode for index, barcode in enumerate(var_names_df["gene"])}
    del var_names_df

    coo = data.X.tocoo(copy=False)
    df = pd.DataFrame({'index': coo.row, 'col': coo.col, 'data': coo.data})
    df["barcode"] = df["index"].map(barcodes)
    df["gene"] = df["col"].map(genes)
    return df[["barcode", "gene", "data"]]


def plot(df, panels=None, x="x", y="y", facecolors=None, label=None, max_labels=10, diag=True,
                xlabel=None, ylabel=None, title="", palette=None, ci=95, y_pos=0.95, y_pos_delta=-0.05,
                outdir=None, filename="plot"):
    fig, (ax, legend_ax) = plt.subplots(nrows=1, ncols=2, figsize=(12, 9), gridspec_kw={"width_ratios": [0.99, 0.01]})
    sns.set(color_codes=True)
    sns.set_style("ticks")

    sns.despine(fig=fig, ax=ax)

    if xlabel is None:
        xlabel = x
    if ylabel is None:
        ylabel = y

    if panels is None:
        df["panels"] = "ALL"
        panels = "panels"

    panel_values = list(df[panels].unique())
    panel_values.sort()
    nplots = len(panel_values)
    ncols = math.ceil(np.sqrt(nplots))
    nrows = math.ceil(nplots / ncols)

    sns.set_style("ticks")
    fig, axes = plt.subplots(nrows=nrows,
                             ncols=ncols,
                             sharex='none',
                             sharey='none',
                             figsize=(12 * ncols, 9 * nrows))
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
            print(panel_values[i])
            data = df.loc[df[panels] == panel_values[i], :].copy()
            if data.shape[0] == 0:
                continue
            data.sort_values(by=y, inplace=True)

            color = "#000000"
            if panel_values[i] in palette:
                color = palette[panel_values[i]]

            data["colors"] = "#808080"
            if facecolors is not None:
                data["colors"] = data[facecolors].map(palette)

            sns.despine(fig=fig, ax=ax)
            sns.regplot(x=x, y=y, data=data, ci=ci,
                        scatter_kws={'facecolors': data["colors"],
                                     'edgecolors': data["colors"],
                                    'alpha': 0.3},
                        line_kws={"color": color},
                        ax=ax
                        )

            if label is not None:
                texts = []
                for j, (_, point) in enumerate(data.iterrows()):
                    if j > max_labels:
                        continue
                    texts.append(ax.text(point[x],
                                         point[y],
                                         str(point[label]),
                                         color=point["colors"]))

                    # adjust_text(texts,
                    #             ax=ax,
                    #             arrowprops=dict(arrowstyle='-', color='#808080'))

            ax.axhline(0, ls='--', color="#808080", alpha=0.3, zorder=-1)
            ax.axvline(0, ls='--', color="#808080", alpha=0.3, zorder=-1)

            if diag:
                low_x, high_x = ax.get_xlim()
                low_y, high_y = ax.get_ylim()
                low = max(low_x, low_y)
                high = min(high_x, high_y)
                ax.plot([low, high], [low, high], ls='--', color="#808080", alpha=0.3, zorder=-1)

            annotations = ["All"]
            if facecolors is not None:
                annotations.extend(list(set(data[facecolors].values)))
            tmp_y_pos = y_pos
            for annotation in annotations:
                subset = data[[x, y]].copy()
                color = "#000000"
                if annotation != "All":
                    subset = data.loc[data[facecolors] == annotation, [x, y]].copy()
                    color = palette[annotation]

                xmedian = np.median(subset[x])
                xstd = np.std(subset[x])
                ymedian = np.median(subset[y])
                ystd = np.std(subset[y])
                ax.plot([xmedian, xmedian], [ymedian + ystd, ymedian - ystd], ls='--', color=color, linewidth=2, alpha=0.3, zorder=-1)
                ax.plot([xmedian + xstd, xmedian - xstd], [ymedian, ymedian], ls='--', color=color, linewidth=2, alpha=0.3, zorder=-1)
                ax.annotate(
                    "{} [N={:,}, x={:.2f}, y={:.2f}]".format(annotation, subset.shape[0], xmedian, ymedian),
                    xy=(0.03, tmp_y_pos),
                    xycoords=ax.transAxes,
                    color=color,
                    fontsize=14,
                    fontweight='bold'
                )
                tmp_y_pos += y_pos_delta

            ax.set_title(panel_values[i],
                         fontsize=18,
                         color=color,
                         weight='bold')
            ax.set_ylabel(ylabel,
                          fontsize=14,
                          fontweight='bold')
            ax.set_xlabel(xlabel,
                          fontsize=14,
                          fontweight='bold')

            legend_ax.set_axis_off()
        else:
            ax.set_axis_off()

        col_index += 1
        if col_index > (ncols - 1):
            col_index = 0
            row_index += 1

    plt.tight_layout()
    fpath = filename + ".png"
    if outdir is not None:
        fpath = os.path.join(outdir, fpath)
    fig.savefig(fpath)
    plt.close()

dataset_dirs = {
    "wg0": {
        "Cain2023": "2023-09-10-Cain2023",
        "Mathys2019": "2023-09-10-Mathys2019",
        "RocheAD2022": "2023-09-10-RocheAD2022",
        "RocheColumbia2022": "2023-09-10-RocheColumbia2022",
        "RocheMS2022": "2023-09-10-RocheMS2022",
        "Zhou2020": "2023-09-07-Zhou2020"
    }, "wg1": {
        "Cain2023": "2024-03-20-Cain2023",
        "Mathys2019": "2024-03-20-Mathys2019",
        "RocheAD2022": "2024-03-22-RocheAD2022",
        "RocheColumbia2022": "2024-03-24-RocheColumbia2022",
        "RocheMS2022": "2024-03-24-RocheMS2022",
        "Zhou2020": "2024-03-20-Zhou2020"
    }, "wg2": {
        "Cain2023": "2024-03-21-Cain2023",
        "Mathys2019": "2024-03-21-Mathys2019",
        "RocheAD2022": "2024-03-22-RocheAD2022",
        "RocheColumbia2022": "2024-03-23-RocheColumbia2022",
        "RocheMS2022": "2024-03-22-RocheMS2022",
        "Zhou2020": "2024-03-21-Zhou2020"
    }
}

# Load the metadata files.
genes_of_interest = None
if args.genes is not None:
    print("Loading genes.")
    genes_of_interest_df = pd.read_csv(args.genes, sep="\t")
    genes_of_interest = set(genes_of_interest_df["ENSG"] + "_" + genes_of_interest_df["GeneID"])
    print("\tLoaded {:,} genes.".format(len(genes_of_interest)))
    del genes_of_interest_df

for dataset in args.datasets:
    print("Processing '{}'".format(dataset))

    # Create the output file.
    out_dir = os.path.join(args.out_dir, str(dataset))
    os.makedirs(out_dir, exist_ok=True)

    genes_results_fpath = os.path.join(out_dir, dataset + "cellbender_delta_genes.tsv.gz")
    barcodes_results_fpath = os.path.join(out_dir, dataset + "cellbender_delta_barcodes.tsv.gz")
    if not os.path.exists(genes_results_fpath) or not os.path.exists(barcodes_results_fpath) or args.force:
        # Load the metadata files.
        wg1_metadata_fpath = os.path.join(args.wg1_dir, dataset_dirs["wg1"][dataset], "Step2-DemultiplexingAndDoubletRemoval", "CombinedResults", "Final_Assignments_demultiplexing_doublets.tsv.gz")
        if not os.path.exists(wg1_metadata_fpath):
            print("Error, WG1 metadata file is missing.")
            exit()
        wg1_df = pd.read_csv(wg1_metadata_fpath, sep="\t", header=0, index_col=None)
        print(wg1_df)

        wg2_metadata_fpath = os.path.join(args.wg2_dir, dataset_dirs["wg2"][dataset], "map", "azimuth_all.metadata.tsv.gz")
        if not os.path.exists(wg2_metadata_fpath):
            wg2_metadata_fpath = os.path.join(args.wg2_dir, dataset_dirs["wg2"][dataset], "map", "azimuth.metadata.tsv.gz")
            if not os.path.exists(wg2_metadata_fpath):
                print("Error, WG2 metadata file is missing.")
                exit()
        wg2_df = pd.read_csv(wg2_metadata_fpath, sep="\t", header=0, index_col=None)

        if args.wg2_pairing is not None:
            pair_df = pd.read_csv(args.wg2_pairing, sep=";", header=0, index_col=None)
            wg2_df = wg2_df.merge(pair_df, how="left")
            del pair_df
        print(wg2_df)

        # Find the input files.
        cr_count_fpaths = glob.glob(os.path.join(args.wg0_dir, dataset_dirs["wg0"][dataset], "CellRanger", "*", "outs", "raw_feature_bc_matrix.h5"))
        cr_count_fpaths = dict(zip([fpath.split(os.sep)[-3] for fpath in cr_count_fpaths], cr_count_fpaths))
        cb_count_fpaths = glob.glob(os.path.join(args.wg0_dir, dataset_dirs["wg0"][dataset], "CellBender", "*", "cellbender_remove_background_output_filtered.h5"))
        cb_count_fpaths = dict(zip([fpath.split(os.sep)[-2] for fpath in cb_count_fpaths], cb_count_fpaths))
        pools = [pool for pool in cr_count_fpaths.keys() if pool in cb_count_fpaths.keys()]
        n_pools = len(pools)
        print("  Found {:,} pools".format(n_pools))
        if n_pools == 0:
            continue

        ################################################

        barcode_df_list = []
        genes_df_list = []
        for pool_index, pool in enumerate(pools):
            print("\tProcessing: {}/{}\t pool '{}'.".format(pool_index, n_pools - 1, pool))

            # Read in the count data.
            cb_df = load_coo_df(cb_count_fpaths[pool])
            cr_df = load_coo_df(cr_count_fpaths[pool])

            # Merge the sparse matrix coordinates.
            df = cb_df.merge(cr_df, how="left", on=["barcode", "gene"], suffixes=["_cb", "_cr"]).fillna(0)
            df["delta"] = df["data_cb"] - df["data_cr"]
            df.sort_values(by="delta", ascending=False, inplace=True)
            del cb_df, cr_df
            print(df)

            # Add a column if the gene is of interest or not.
            df["include"] = "False"
            if genes_of_interest is not None:
                df.loc[df["gene"].isin(genes_of_interest), "include"] = "True"

            # Add the metadata.
            df["Barcode"] = df["barcode"].str.split("-", n=1, expand=True)[0] + "_" + str(pool)
            df["Pool"] = str(pool)
            df["Dataset"] = str(dataset)
            df = df.merge(wg1_df[["Pool", "Barcode", "DropletType"]], how="left").merge(wg2_df[["Pool", "Barcode", args.cell_level]], how="left")
            print(df)

            # Calculate the stats per gene.
            grouped_genes_df = df.loc[df["DropletType"] == "singlet", :].groupby(["gene", "include", args.cell_level]).sum(["data_cb", "data_cr", "delta"]).reset_index(drop=False)
            grouped_genes_df["prop"] = (100 / grouped_genes_df["data_cr"]) * grouped_genes_df["delta"].abs()
            # grouped_genes_df.to_csv("test.tsv.gz", sep="\t", header=True, index=False, compression="gzip")
            # grouped_genes_df = pd.read_csv("test.tsv.gz", sep="\t", header=0, index_col=None)
            grouped_genes_df["label"] = grouped_genes_df["gene"].str.split("_", n=1, expand=True)[1]
            print(grouped_genes_df)
            plot(
                df=grouped_genes_df,
                panels=args.cell_level,
                x="data_cb",
                y="delta",
                facecolors="include",
                label="label",
                palette=palette,
                xlabel="CellBender counts",
                ylabel="(CellBender - CellRanger) counts",
                title="Gene level summary",
                y_pos=0.05,
                y_pos_delta=0.05,
                ci=None,
                diag=False,
                outdir=out_dir,
                filename="genes_delta"
            )
            plot(
                df=grouped_genes_df,
                panels=args.cell_level,
                x="prop",
                y="delta",
                facecolors="include",
                label="label",
                palette=palette,
                xlabel="% counts removed",
                ylabel="(CellBender - CellRanger) counts",
                title="Gene level summary",
                y_pos=0.05,
                y_pos_delta=0.05,
                ci=None,
                diag=False,
                outdir=out_dir,
                filename="genes_prop_delta"
            )
            genes_df_list.append(grouped_genes_df)
            del grouped_genes_df

            # Calculate the stats per cell type.
            grouped_barcodes_df = df.loc[df["DropletType"] == "singlet", :].groupby(["barcode", "include"]).sum(["data_cb", "data_cr", "delta"]).reset_index(drop=False)
            values = ["data_cb", "data_cr", "delta"]
            grouped_barcodes_df = pd.pivot_table(grouped_barcodes_df, values=values, index="barcode", columns="include").reset_index(drop=False)
            grouped_barcodes_df.columns = grouped_barcodes_df.columns.map('_'.join).str.strip('_')
            for value in values:
                grouped_barcodes_df[value + "_total"] = grouped_barcodes_df[value + "_True"] + grouped_barcodes_df[value + "_False"]
            print(grouped_barcodes_df)
            print(grouped_barcodes_df["delta_total"].mean(), grouped_barcodes_df["delta_total"].median())

            # print(grouped_barcodes_df)
            barcode_df_list.append(grouped_barcodes_df)
            del grouped_barcodes_df, values

            # TODO: remove
            if pool_index == 0:
                break

        # Merge and save the genes.
        genes_df = pd.concat(genes_df_list, axis=0)
        genes_df.sort_values(by=["delta"], ascending=True, inplace=True)
        # print(genes_df)
        genes_df.to_csv(genes_results_fpath, sep="\t", header=True, index=False, compression="gzip")

        # Merge and save the genes.
        barcodes_df = pd.concat(barcode_df_list, axis=0)
        barcodes_df.sort_values(by=["delta_total"], ascending=True, inplace=True)
        # print(barcodes_df)
        barcodes_df.to_csv(barcodes_results_fpath, sep="\t", header=True, index=False, compression="gzip")
    else:
        genes_df = pd.read_csv(genes_results_fpath, sep="\t", header=0, index_col=None)
        barcodes_df = pd.read_csv(barcodes_results_fpath, sep="\t", header=0, index_col=None)

    ################################################

    print(genes_df)

    # genes_df["label"] = grouped_genes_df["gene"].str.split("_", n=1, expand=True)[1]
    # plot(
    #     df=genes_df,
    #     x="data_cb",
    #     y="data_cr",
    #     facecolors="include",
    #     label="label",
    #     palette={"True": "#009E73", "False": "#D55E00", True: "#009E73", False: "#D55E00"},
    #     xlabel="CellBender counts",
    #     ylabel="CellRanger counts",
    #     title="Gene level summary",
    #     outdir=out_dir,
    #     filename="genes"
    # )
    # plot(
    #     df=genes_df,
    #     x="data_cb",
    #     y="delta",
    #     facecolors="include",
    #     label="label",
    #     palette={"True": "#009E73", "False": "#D55E00", True: "#009E73", False: "#D55E00"},
    #     xlabel="CellBender counts",
    #     ylabel="(CellBender - CellRanger) counts",
    #     title="Gene level summary",
    #     y_pos=0.05,
    #     y_pos_delta=0.05,
    #     outdir=out_dir,
    #     filename="genes_delta"
    # )
    # exit()

    ################################################

    print(barcodes_df)

    # # Load the metadata files.
    # wg1_metadata_fpath = os.path.join(args.wg1_dir, dataset_dirs["wg1"][dataset], "Step2-DemultiplexingAndDoubletRemoval", "CombinedResults", "Final_Assignments_demultiplexing_doublets.tsv.gz")
    # if not os.path.exists(wg1_metadata_fpath):
    #     print("Error, WG1 metadata file is missing.")
    #     exit()
    # wg1_df = pd.read_csv(wg1_metadata_fpath, sep="\t", header=0, index_col=None)
    # print(wg1_df)
    #
    # wg2_metadata_fpath = os.path.join(args.wg2_dir, dataset_dirs["wg2"][dataset], "map", "azimuth_all.metadata.tsv.gz")
    # if not os.path.exists(wg2_metadata_fpath):
    #     wg2_metadata_fpath = os.path.join(args.wg2_dir, dataset_dirs["wg2"][dataset], "map", "azimuth.metadata.tsv.gz")
    #     if not os.path.exists(wg2_metadata_fpath):
    #         print("Error, WG2 metadata file is missing.")
    #         exit()
    # wg2_df = pd.read_csv(wg2_metadata_fpath, sep="\t", header=0, index_col=None)
    #
    # if args.wg2_pairing is not None:
    #     pair_df = pd.read_csv(args.wg2_pairing, sep=";", header=0, index_col=None)
    #     wg2_df = wg2_df.merge(pair_df, how="left")
    #     del pair_df
    # print(wg2_df)
    #
    # # Merge all the data together.
    # barcodes_df = barcodes_df.merge(wg1_df[["Pool", "Barcode", "DropletType"]], how="left").merge(wg2_df[["Pool", "Barcode", args.cell_level]], how="left")
    # print(barcodes_df)
    #
    # plot(
    #     df=barcodes_df,
    #     x="delta_True",
    #     y="delta_False",
    #     facecolors=args.cell_level,
    #     palette=palette,
    #     xlabel="delta count selected genes",
    #     ylabel="delta count other genes",
    #     title="Barcode level summary",
    #     filename="barcodes"
    # )

print("\nDone")