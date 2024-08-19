#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse


"""
Syntax: 
./compare_pseudobulk.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--bryois_indir", type=str, required=True, help="")
parser.add_argument("--label1", type=str, required=True, help="")
parser.add_argument("--pseudobulk_indir", type=str, required=True, help="")
parser.add_argument("--label2", type=str, required=True, help="")
parser.add_argument("--datasets", nargs="+", type=str, required=False, default=["2023-09-10-Cain2023", "2023-09-10-Mathys2019", "2023-09-10-RocheAD2022", "2023-09-10-RocheColumbia2022", "2023-09-10-RocheMS2022", "2023-09-07-Zhou2020"], help="")
parser.add_argument("--cell_types", nargs="+", type=str, required=False, default=["Astrocytes", "Endothelial.cells", "Excitatory.neurons", "Inhibitory.neurons", "Microglia", "Oligodendrocytes", "OPCs...COPs", "Pericytes"], help="")
parser.add_argument("--value_columns", nargs="+", type=str, required=False, default=["counts", "n_expressed", "perc_expressed"], help="")
parser.add_argument("--limitgenes", type=str, required=False, help="")
parser.add_argument("--palette", type=str, required=False, help="")
parser.add_argument("--no_demean", dest="zero_mean", action='store_false', help="")
parser.add_argument("--force", action='store_true', help="")
parser.add_argument("--outdir", type=str, required=True, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import os
import json
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

palette = None
if args.palette is not None:
    with open(args.palette) as f:
        palette = json.load(f)
    f.close()
    palette["all"] = "#000000"

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def plot_histplot(df, x="x", hue=None, palette=None, xlabel="", ylabel="", title="", filename="plot"):
    sns.set_style("ticks")

    sns.set(rc={'figure.figsize': (18, 12)})
    sns.set_style("ticks")
    fig, (ax1, ax2) = plt.subplots(nrows=1,
                                   ncols=2,
                                   gridspec_kw={"width_ratios": [0.99, 0.01]})
    sns.despine(fig=fig, ax=ax1)

    g = sns.histplot(x=x,
                     hue=hue,
                     data=df,
                     color="#808080",
                     palette=palette,
                     legend=False,
                     ax=ax1)

    ax1.set_title(title,
                  fontsize=16,
                  fontweight='bold')
    ax1.set_ylabel(ylabel,
                   fontsize=14,
                   fontweight='bold')
    ax1.set_xlabel(xlabel,
                   fontsize=14,
                   fontweight='bold')

    if palette is not None:
        handles = []
        for key, color in palette.items():
            if key in df[hue].values.tolist():
                subset = df.loc[df[hue] == key, :]
                avg = subset[x].mean()
                label = "{}\n[n={:,}, avg={:.2f}]".format(key, subset.shape[0], avg)
                handles.append(mpatches.Patch(color=color, label=label))
        ax2.legend(handles=handles, loc="center")
    ax2.set_axis_off()

    plt.tight_layout()
    fpath = os.path.join(args.outdir, "plot", "{}_{}_{}.png".format(args.label1, args.label2, filename))
    print("Saving {}".format(os.path.basename(fpath)))
    fig.savefig(fpath)
    plt.close()


def plot_heatmap(df, x="x", xlabel="", ylabel="", title="", filename="plot", scale=1):
    corr_df = pd.pivot_table(df, values=x, index='individual_id1', columns='individual_id2')

    sns.set_style("ticks")
    cmap = sns.diverging_palette(246, 24, as_cmap=True)

    fig, axes = plt.subplots(nrows=2,
                             ncols=2,
                             figsize=(scale * corr_df.shape[1] + 10, scale * corr_df.shape[0] + 10),
                             gridspec_kw={"width_ratios": [0.2, 0.8],
                                          "height_ratios": [0.8, 0.2]})
    sns.set(color_codes=True)

    row_index = 0
    col_index = 0
    for _ in range(4):
        ax = axes[row_index, col_index]
        if row_index == 0 and col_index == 1:
            sns.heatmap(corr_df, cmap=cmap, vmin=-1, vmax=1, center=0,
                        square=True, annot=corr_df.round(2), fmt='',
                        cbar=False, annot_kws={"size": 14, "color": "#000000"},
                        ax=ax)

            plt.setp(ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=20, rotation=0))
            plt.setp(ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=20, rotation=90))

            ax.set_title(title, fontsize=16)

            ax.set_xlabel(xlabel, fontsize=14)
            ax.xaxis.set_label_position('top')

            ax.set_ylabel(ylabel, fontsize=14)
            ax.yaxis.set_label_position('right')
        else:
            ax.set_axis_off()

        col_index += 1
        if col_index > 1:
            col_index = 0
            row_index += 1

    plt.tight_layout()
    fpath = os.path.join(args.outdir, "plot", "{}_{}_{}.png".format(args.label1, args.label2, filename))
    print("Saving {}".format(os.path.basename(fpath)))
    fig.savefig(fpath)
    plt.close()

def pivot_table(df, ct, key):
    df = pd.pivot_table(df.loc[df['cell_type'] == ct, :], values=key, index=['cell_type', 'ensembl', 'symbol'], columns='individual_id')
    df.reset_index(drop=False, inplace=True)
    df.index = df["symbol"] + "_" + df["ensembl"]
    df.drop(["cell_type", "ensembl", "symbol"], axis=1, inplace=True)

    if len(set(df.columns)) != len(df.columns):
        print("Error, columns are not unique.")
        exit()
    if len(set(df.index)) != len(df.index):
        print("Error, indices are not unique.")
        exit()
    return df


################################################

sample_corr_fpath = os.path.join(args.outdir, "sample_corr.tsv")
gene_corr_fpath = os.path.join(args.outdir, "gene_corr.tsv")
if (not os.path.exists(sample_corr_fpath) and not os.path.exists(gene_corr_fpath)) or args.force:
    print("Loading data")
    bryois_pd_df = pd.read_csv(os.path.join(args.bryois_indir, "Brain_single_cell_eqtl_MV", "Excitatory.neurons.sum_expression.individual_id.tsv.gz"), sep="\t", header=0, index_col=None)
    print("\tBryois pseudobulk data frame: {}".format(bryois_pd_df.shape))
    print(bryois_pd_df)

    scmd_pd_list = []
    individual_id_to_dataset = {}
    for dataset in args.datasets:
        pb_fpath = os.path.join(args.pseudobulk_indir, dataset, "sum_expression.individual_id.tsv.gz")
        if not os.path.exists(pb_fpath):
            print("\tWarning, {} has no pseodobulk.".format(dataset))
            continue

        scmd_pd_df = pd.read_csv(pb_fpath, sep="\t", header=0, index_col=None)
        for individual_id in scmd_pd_df["individual_id"]:
            individual_id_to_dataset[individual_id] = dataset
        scmd_pd_list.append(scmd_pd_df)
        del scmd_pd_df

    scmd_pd_df = pd.concat(scmd_pd_list, axis=0)
    del scmd_pd_list
    print("\tscMetaBrain pseudobulk data frame: {}".format(scmd_pd_df.shape))

    gene_corr_data = []
    sample_corr_data = []
    n_datasets = set()
    for cell_type in args.cell_types:
        print("Processing '{}'.".format(cell_type))
        ct = cell_type.replace("...", " / ").replace(".", " ")

        # Loop through the numerical variables we can correlate.
        for value_column in args.value_columns:
            print("\tExtracting '{}' data".format(value_column))

            bryois_df = pivot_table(df=bryois_pd_df, ct=ct, key=value_column)
            print(bryois_df)

            scmd_df = pivot_table(df=scmd_pd_df, ct=ct, key=value_column)
            print(scmd_df)

            print("\tFiltering on overlap")

            # Filter on overlapping
            overlapping_genes = list(set(bryois_df.index).intersection(set(scmd_df.index)))
            overlapping_samples = list(set(bryois_df.columns).intersection(set(scmd_df.columns)))
            n_genes = len(overlapping_genes)
            n_samples = len(overlapping_samples)
            print("\t\t{:,} genes overlapping and {:,} samples overlapping".format(n_genes, n_samples))
            # print(bryois_df.loc[overlapping_genes, overlapping_samples])
            # print(scmd_df.loc[overlapping_genes, overlapping_samples])
            bryois_m = bryois_df.loc[overlapping_genes, overlapping_samples].to_numpy()
            scmb_m = scmd_df.loc[overlapping_genes, overlapping_samples].to_numpy()
            del bryois_df, scmd_df

            # Set the mean to 0.
            if args.zero_mean:
                print("\t\tSetting row average to zero.")
                bryois_m = bryois_m - bryois_m.mean(axis=1, keepdims=True)
                scmb_m = scmb_m - scmb_m.mean(axis=1, keepdims=True)

            print("\tCorrelating per sample")

            # Correlate the samples.
            n_correlations = n_samples * n_samples
            corr_index = 0
            for individual_index1, individual_id1 in enumerate(overlapping_samples):
                for individual_index2, individual_id2 in enumerate(overlapping_samples):
                    if corr_index == 0 or corr_index % 1000 == 0:
                        print("\t\tCalculated {:,} / {:,} correlations".format(corr_index + 1, n_correlations), end='\r')

                    pearson_coef = np.nan
                    if np.std(bryois_m[:, individual_index1]) > 1e-7 and np.std(scmb_m[:, individual_index2]) > 1e-7:
                        pearson_coef, _ = stats.pearsonr(x=bryois_m[:, individual_index1], y=scmb_m[:, individual_index2])

                    sample_corr_data.append([individual_id_to_dataset[individual_id1].split("-")[-1], individual_id1, individual_id_to_dataset[individual_id2].split("-")[-1], individual_id2, n_genes, value_column, pearson_coef])
                    corr_index += 1
            print("\t\tCalculated {:,} / {:,} correlations".format(corr_index + 1, n_correlations))

            print("\tCorrelating per gene")

            # Correlate the genes per dataset.
            datasets_a = np.array([individual_id_to_dataset[individual_id].split("-")[-1] for individual_id in overlapping_samples])
            n_correlations = (len(np.unique(datasets_a)) + 1) * n_genes
            corr_index = 0
            for dataset in list(np.unique(datasets_a)) + ["all"]:
                # Create a mask for the samples.
                sample_mask = np.ones(len(overlapping_samples), dtype=bool)
                if dataset != "all":
                    sample_mask = datasets_a == dataset
                n_samples = len(sample_mask)

                # Correlate the genes.
                for gene_index, gene in enumerate(overlapping_genes):
                    if corr_index == 0 or corr_index % 1000 == 0:
                        print("\t\tCalculated {:,} / {:,} correlations".format(corr_index, n_correlations), end='\r')

                    pearson_coef = np.nan
                    if np.std(bryois_m[gene_index, :][sample_mask]) > 1e-7 and np.std(scmb_m[gene_index, :][sample_mask]) > 1e-7:
                        pearson_coef, _ = stats.pearsonr(x=bryois_m[gene_index, :][sample_mask], y=scmb_m[gene_index, :][sample_mask])

                    gene_corr_data.append([dataset, gene, n_samples, value_column, pearson_coef])
                    corr_index += 1
            print("\t\tCalculated {:,} / {:,} correlations".format(corr_index, n_correlations))

        sample_corr_df = pd.DataFrame(sample_corr_data, columns=["dataset1", "individual_id1", "dataset2", "individual_id2", "n_genes", "variable", "pearsonr"])
        sample_corr_df.sort_values(by="pearsonr", inplace=True)
        print(sample_corr_df)
        sample_corr_df.to_csv(sample_corr_fpath, sep="\t", header=True, index=False)

        gene_corr_df = pd.DataFrame(gene_corr_data, columns=["dataset", "gene", "n_samples", "variable", "pearsonr"])
        gene_corr_df.dropna(inplace=True)
        gene_corr_df.sort_values(by="pearsonr", inplace=True)
        print(gene_corr_df)
        gene_corr_df.to_csv(gene_corr_fpath, sep="\t", header=True, index=False)

sample_corr_df = pd.read_csv(sample_corr_fpath, sep="\t", header=0, index_col=None)
print(sample_corr_df)
gene_corr_df = pd.read_csv(gene_corr_fpath, sep="\t", header=0, index_col=None)
print(gene_corr_df)

################################################

for value_column in args.value_columns:
    plot_histplot(
        df=sample_corr_df.loc[(sample_corr_df["individual_id1"] == sample_corr_df["individual_id2"]) & (sample_corr_df["variable"] == value_column), :],
        x="pearsonr",
        hue="dataset1",
        palette=palette,
        xlabel=value_column + " pearson r",
        ylabel="count",
        title=args.label1 + " vs " + args.label2 + " - " + value_column,
        filename="sample_" + value_column + "_pearsonr"
    )
    for dataset in sample_corr_df["dataset1"].unique():
        if not "Columbia" in dataset:
            continue
        plot_heatmap(
            df=sample_corr_df.loc[(sample_corr_df["dataset1"] == dataset) & (sample_corr_df["dataset2"] == dataset) & (sample_corr_df["variable"] == value_column), :],
            x="pearsonr",
            xlabel=value_column + " pearson r",
            ylabel="count",
            title=args.label1 + " vs " + args.label2 + " - " + value_column,
            filename="sample_" + value_column + "_" + dataset + "_pearsonr"
        )

    plot_histplot(
        df=gene_corr_df.loc[gene_corr_df["variable"] == value_column, :],
        x="pearsonr",
        hue="dataset",
        palette=palette,
        xlabel=value_column + " pearson r",
        ylabel="count",
        title=args.label1 + " vs " + args.label2 + " - " + value_column,
        filename="gene_" + value_column + "_pearsonr"
    )

print("END")