#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse

"""
Syntax: 
./compare_expression.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--expr1", type=str, required=True, help="")
parser.add_argument("--label1", type=str, required=True, help="")
parser.add_argument("--long1", action='store_true', help="")
parser.add_argument("--filter1", type=str, required=False, default=None, help="")
parser.add_argument("--datasets1", nargs="+", type=str, required=False, default=None, help="")
parser.add_argument("--expr2", type=str, required=True, help="")
parser.add_argument("--label2", type=str, required=True, help="")
parser.add_argument("--long2", action='store_true', help="")
parser.add_argument("--filter2", type=str, required=False, default=None, help="")
parser.add_argument("--datasets2", nargs="+", type=str, required=False, default=None, help="")
parser.add_argument("--keep_dupl", action='store_true', help="")
parser.add_argument("--gte", type=str, required=False, help="")
parser.add_argument("--sample_to_ind", nargs="+", type=int, required=False, choices=[1, 2], default=None, help="")
parser.add_argument("--limitgenes", type=str, required=False, help="")
parser.add_argument("--force", action='store_true', help="")
parser.add_argument("--method", type=str, required=False, choices=["pearson", "spearman"], default="pearson", help="")
parser.add_argument("--no_demean", dest="zero_mean", action='store_false', help="")
parser.add_argument("--corr_all_samples", action='store_true', help="")
parser.add_argument("--palette", type=str, required=False, help="")
parser.add_argument("--outdir", type=str, required=True, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

if args.sample_to_ind is not None and args.gte is None:
    print("Error, --gte need to be set if --sample_to_ind is not None")
    exit()

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

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)
if not os.path.exists(os.path.join(args.outdir, "plot")):
    os.makedirs(os.path.join(args.outdir, "plot"))

palette = None
if args.palette is not None:
    with open(args.palette) as f:
        palette = json.load(f)
    f.close()
    palette["all"] = "#000000"
    palette["DATASET"] = "#000000"

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def load_expr_data(fpath, sep="\t", header=0, index_col=0, nrows=None, long=False, filter=None, datasets=None):
    if datasets is None:
        datasets = ["DATASET"]

    if long:
        index_col=None

    df_list = []
    sample_to_dataset_dict = {}
    for dataset in datasets:
        dataset_fpath = fpath.replace("DATASET", dataset)
        print("  Loading {}".format(dataset_fpath))
        if dataset_fpath.endswith(".bed") or dataset_fpath.endswith(".bed.gz"):
            df = pd.read_csv(dataset_fpath, sep=sep, header=header, index_col=None, nrows=nrows)
            df.index = df["ID"].str.split("_", n=1, expand=True)[0]
            df.index.name = None
            df = df.iloc[:, 4:]
        else:
            df = pd.read_csv(dataset_fpath, sep=sep, header=header, index_col=index_col, nrows=nrows)

        if filter:
            column, value = filter.split(":")
            df = df.loc[df[column] == value, :]

        if long:
            df = pd.pivot_table(df, values='counts', index=['ensembl', 'symbol'], columns='individual_id')
            df = df.droplevel('ensembl')

        if len(set(df.columns)) != len(df.columns):
            print("Error, columns are not unique.")
            exit()

        # Check for duplicates in the genes or barcodes.
        n_genes = df.shape[0]
        index_mask = np.ones(n_genes, dtype=bool)
        if len(set(df.index)) != len(df.index):
            print("\tWarning, genes are not unique.")
            if not args.keep_dupl:
                seen = set()
                duplicates = []
                for index in df.index:
                    if index in seen:
                        duplicates.append(index)
                    seen.add(index)

                index_mask = np.array([gene not in duplicates for gene in df.index], dtype=bool)
                print("\tRemoving '{:,}' duplicate genes having {:,} rows in total.".format(len(duplicates), np.size(index_mask) - np.sum(index_mask)))

        sample_to_dataset_dict.update({column: dataset.split("-")[-1] for column in df.columns})

        df_list.append(df.loc[index_mask, :])

    df = pd.concat(df_list, axis=1).dropna()
    if len(set(df.columns)) != len(df.columns):
        print("Error, columns are not unqiue.")
        exit()
    if len(set(df.index)) != len(df.index):
        print("Error, indices are not unqiue.")
        exit()

    return df, sample_to_dataset_dict

def correlate(x, y, data=None):
    if data is not None:
        x = data[x]
        y = data[y]
    if args.method == "pearson":
        return stats.pearsonr(y, x)
    elif args.method == "spearman":
        return stats.spearmanr(y, x)
    else:
        print("Error, unexpected method {}".format(args.method))
        exit()

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
    fpath = os.path.join(args.outdir, "plot", "{}_{}_{}{}.png".format(args.label1, args.label2, filename, "_limitgenes" if args.limitgenes is not None else ""))
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
    fpath = os.path.join(args.outdir, "plot", "{}_{}_{}{}.png".format(args.label1, args.label2, filename, "_limitgenes" if args.limitgenes is not None else ""))
    print("Saving {}".format(os.path.basename(fpath)))
    fig.savefig(fpath)
    plt.close()

def plot_regplot(df, x="x", y="y", hue=None, palette=None, xlabel="", ylabel="", title="", filename="plot"):
    sns.set(rc={'figure.figsize': (12, 9)})
    sns.set_style("ticks")
    fig, (ax1, ax2) = plt.subplots(nrows=1,
                                   ncols=2,
                                   gridspec_kw={"width_ratios": [0.99, 0.01]})
    sns.despine(fig=fig, ax=ax1)

    if hue is None:
        df["hue"] = "all"
        hue = "hue"

    info = {}
    hue_values = ["all"] + list(df[hue].unique())
    for hue_value in hue_values:
        data = df
        if hue_value != "all":
            data = df.loc[df[hue] == hue_value, :]

        facecolors = "#808080"
        if palette is not None and hue_value in palette:
            facecolors = palette[hue_value]

        coef, _ = correlate(data=data, x=x, y=y)
        info[hue_value] = (data.shape[0], coef)

        sns.regplot(x=x, y=y, data=data, ci=None,
                    scatter_kws={'color': facecolors},
                    line_kws={"color": facecolors},
                    ax=ax1)

    xmin, xmax = ax1.get_xlim()
    ymin, ymax = ax1.get_ylim()

    min_value = max(xmin, ymin)
    max_value = min(xmax, ymax)
    ax1.plot([min_value, max_value], [min_value, max_value], ls="--", c=".3")

    ax1.set_xlabel(xlabel,
                  fontsize=14,
                  fontweight='bold')
    ax1.set_ylabel(ylabel,
                  fontsize=14,
                  fontweight='bold')
    ax1.set_title(title,
                 fontsize=18,
                 fontweight='bold')

    if hue is not None and palette is not None:
        handles = []
        for hue_value in hue_values:
            color = palette[hue_value]
            label = "{}\n[n={:,}, avg={:.2f}".format(hue_value, *info[hue_value])
            handles.append(mpatches.Patch(color=color, label=label))
        ax2.legend(handles=handles, loc="center")
    ax2.set_axis_off()

    plt.tight_layout()
    fpath = os.path.join(args.outdir, "plot", "{}_{}_{}.png".format(args.label1, args.label2, filename))
    print("Saving {}".format(os.path.basename(fpath)))
    fig.savefig(fpath)
    plt.close()


################################################

sample_corr_fpath = os.path.join(args.outdir, "sample_corr.tsv")
gene_corr_fpath = os.path.join(args.outdir, "gene_corr.tsv")
if (not os.path.exists(sample_corr_fpath) and not os.path.exists(gene_corr_fpath)) or args.force:
    print("Loading data")
    expr_df1, sample_to_dataset_dict1 = load_expr_data(fpath=args.expr1, long=args.long1, filter=args.filter1, datasets=args.datasets1)
    print("\tExpression matrix 1: {}".format(expr_df1.shape))

    expr_df2, sample_to_dataset_dict2 = load_expr_data(fpath=args.expr2, long=args.long2, filter=args.filter2, datasets=args.datasets2)
    print("\tExpression matrix 2: {}".format(expr_df2.shape))

    # Deal with the WG3 format.
    # expr_df1.columns = [column.split(";;")[1].split("_")[0] for column in expr_df1.columns]
    # expr_df2.columns = [column.split(";;")[1].split("_")[0] for column in expr_df2.columns]

    overlapping_samples = list(set(expr_df1.columns).intersection(set(expr_df2.columns)))
    if args.datasets1 is not None and args.datasets2 is None:
        sample_to_dataset_dict = sample_to_dataset_dict1
    elif args.datasets1 is None and args.datasets2 is not None:
        sample_to_dataset_dict = sample_to_dataset_dict2
    else:
        warning_printed = False
        for sample in overlapping_samples:
            if sample_to_dataset_dict1[sample] != sample_to_dataset_dict2[sample] and not warning_printed:
                print("Warning, samples come from different datasets, using dict1.")
                warning_printed = True
        sample_to_dataset_dict = sample_to_dataset_dict1

    if args.gte is not None:
        gte_df = pd.read_csv(args.gte, sep="\t", header=0, index_col=None)
        print("Gene-to-expression: {}".format(gte_df.shape))

        # Check if all genotype_ids and expression_ids are unique.
        for gte_index in [0, 1]:
            if len(set(gte_df.iloc[:, gte_index])) != len(gte_df.iloc[:, gte_index]):
                print("Error, GTE file column {} contains duplicates".format(gte_index))
                exit()
        individual_ids = set(gte_df.iloc[:, 0])
        sample_ids = set(gte_df.iloc[:, 1])
        sample_to_individual_dict = dict(zip(gte_df.iloc[:, 1], gte_df.iloc[:, 0]))
        sample_to_dataset_dict = dict(zip(gte_df.iloc[:, 0], gte_df.iloc[:, 2]))

        # Translate genotype_id to expression_id for the column names.
        if args.sample_to_ind is not None:
            if 1 in args.sample_to_ind:
                expr_df1.columns = [column if column in individual_ids else (sample_to_individual_dict[column] if column in sample_to_individual_dict else column) for column in expr_df1.columns]
            if 2 in args.sample_to_ind:
                expr_df2.columns = [column if column in individual_ids else (sample_to_individual_dict[column] if column in sample_to_individual_dict else column) for column in expr_df2.columns]

        del gte_df, sample_to_individual_dict

    genes_select = set(expr_df1.index).intersection(set(expr_df2.index))
    if args.limitgenes is not None:
        limitgenes_df = pd.read_csv(args.limitgenes, sep="\t", header=0, index_col=None)
        print("Limit genes matrix: {}".format(limitgenes_df.shape))
        genes_select = set(limitgenes_df["symbol"])

    if len(set(expr_df1.columns)) != len(expr_df1.columns):
        print("Error, columns in data frame 1 are not unique.")
        exit()
    if len(set(expr_df2.columns)) != len(expr_df2.columns):
        print("Error, columns in data frame 2 are not unique.")
        exit()

    print("Filtering on overlap")

    # Filter on overlapping
    overlapping_genes = list(set(expr_df1.index).intersection(set(expr_df2.index)).intersection(genes_select))
    overlapping_samples = list(set(expr_df1.columns).intersection(set(expr_df2.columns)).intersection(sample_to_dataset_dict.keys()))
    n_genes = len(overlapping_genes)
    n_samples = len(overlapping_samples)
    print("\t{:,} genes overlapping and {:,} samples overlapping".format(n_genes, n_samples))
    if n_genes == 0 or n_samples == 0:
        exit()
    for sample in overlapping_samples:
        if sample not in sample_to_dataset_dict:
            print("Error, sample '{}' has no dataset.".format(sample))
            exit()
    datasets_s = pd.Series([sample_to_dataset_dict[sample] for sample in overlapping_samples])
    print(datasets_s.value_counts())
    print(expr_df1.loc[overlapping_genes, overlapping_samples])
    print(expr_df2.loc[overlapping_genes, overlapping_samples])

    expr_m1 = expr_df1.loc[overlapping_genes, overlapping_samples].to_numpy()
    expr_m2 = expr_df2.loc[overlapping_genes, overlapping_samples].to_numpy()
    del expr_df1, expr_df2

    # print("Gene delta:")
    # delta_m = expr_m1 - expr_m2
    # print(pd.DataFrame(delta_m, index=overlapping_genes, columns=overlapping_samples))
    # n_values = delta_m.shape[0] * delta_m.shape[1]
    # n_no_change = np.sum(delta_m == 0)
    # print("N no change: {} {:.2f}%".format(n_no_change, (100 / n_values) * n_no_change))
    # print("Min change: {}".format(np.min(delta_m)))
    # print("Max change: {}".format(np.max(delta_m)))
    # print("Average change: {}".format(np.mean(delta_m)))
    # print("Standard deviation change: {}".format(np.std(delta_m)))
    # print("Median change: {}".format(np.median(delta_m)))
    # print("Average change (not 0): {}".format(np.mean(delta_m[delta_m != 0])))
    # print("Standard deviation change (not 0): {}".format(np.std(delta_m[delta_m != 0])))
    # print("Median change (not 0): {}".format(np.median(delta_m[delta_m != 0])))
    #
    # gene_diff = pd.DataFrame({"n_changed": np.sum(delta_m != 0, axis=1), "sum_change": np.sum(delta_m, axis=1)}, index=overlapping_genes)
    # gene_diff = gene_diff.loc[gene_diff["n_changed"] != 0, :]
    # gene_diff.sort_values(by="sum_change", inplace=True)
    # print(gene_diff)
    # del gene_diff
    #
    # print("Sample delta:")
    # sample_diff = pd.DataFrame({"n_changed": np.sum(delta_m != 0, axis=0), "sum_change": np.sum(delta_m, axis=0)}, index=overlapping_samples)
    # sample_diff = sample_diff.loc[sample_diff["n_changed"] != 0, :]
    # sample_diff.sort_values(by="sum_change", inplace=True)
    # print(sample_diff)
    # del delta_m, sample_diff

    # Set the mean to 0.
    if args.zero_mean:
        print("\tSetting row average to zero.")
        expr_m1 = expr_m1 - expr_m1.mean(axis=1, keepdims=True)
        expr_m2 = expr_m2 - expr_m2.mean(axis=1, keepdims=True)

    # print("Plotting some genes")
    # # Plotting.
    # for gene_index, gene in enumerate(overlapping_genes):
    #     if gene not in ["AC063944.3", "ANGPTL5", "ARMT1", "CD109", "MIR4500HG", "MYOZ2", "SCN3A", "SCRN2", "TMEM14A", "WDPCP", "WFDC2", "ZNF518A", "ZNF579"]:
    #         continue
    #     plot_regplot(
    #         df=pd.DataFrame({"expr1": expr_m1[gene_index, :], "expr2": expr_m2[gene_index, :], "hue": datasets_s.to_numpy()}),
    #         x="expr1",
    #         y="expr2",
    #         hue="hue",
    #         palette=palette,
    #         xlabel=args.label1,
    #         ylabel=args.label2,
    #         title=gene,
    #         filename=gene + "_expression"
    #     )
    #     # if gene_index > 10:
    #     #     break
    # # exit()


    print("Correlating per sample")

    # Correlate the samples.
    sample_corr_data = []
    n_correlations = (n_samples * n_samples) if args.corr_all_samples else n_samples
    corr_index = 0
    for sample_index1, sample1 in enumerate(overlapping_samples):
        for sample_index2, sample2 in enumerate(overlapping_samples):
            if not args.corr_all_samples and sample1 != sample2:
                continue
            if corr_index == 0 or corr_index % 1000 == 0:
                print("\tCalculated {:,} / {:,} correlations".format(corr_index, n_correlations), end='\r')

            coef = np.nan
            if np.std(expr_m1[:, sample_index1]) > 1e-7 and np.std(expr_m2[:, sample_index2]) > 1e-7:
                coef, _ = correlate(x=expr_m1[:, sample_index1], y=expr_m2[:, sample_index2])

            sample_corr_data.append([sample_to_dataset_dict[sample1], sample1, sample_to_dataset_dict[sample2], sample2, n_genes, coef])
            corr_index += 1
    print("\tCalculated {:,} / {:,} correlations".format(corr_index, n_correlations))

    print("Correlating per gene")

    # Correlate the genes per dataset.
    gene_corr_data = []
    datasets_a = datasets_s.to_numpy()
    n_dataset = len(np.unique(datasets_a))
    n_correlations = (n_dataset + (1 if n_dataset > 1 else 0)) * n_genes
    corr_index = 0
    for dataset in list(np.unique(datasets_a)) + ["all"]:
        # if dataset != "MS":
        #     continue
        if n_dataset == 1 and dataset == "all":
            continue

        # Create a mask for the samples.
        sample_mask = np.ones(len(overlapping_samples), dtype=bool)
        if dataset != "all":
            sample_mask = datasets_a == dataset
        n_samples = np.sum(sample_mask)

        # Correlate the genes.
        for gene_index, gene in enumerate(overlapping_genes):
            if corr_index == 0 or corr_index % 1000 == 0:
                print("\tCalculated {:,} / {:,} correlations".format(corr_index, n_correlations), end='\r')

            coef = np.nan
            expr_a1 = expr_m1[gene_index, :][sample_mask]
            expr_a2 = expr_m2[gene_index, :][sample_mask]
            if np.std(expr_a1) > 1e-7 and np.std(expr_a2) > 1e-7:
                coef, _ = correlate(x=expr_a1, y=expr_a2)

            gene_corr_data.append([dataset, gene, n_samples, coef])
            corr_index += 1
    print("\tCalculated {:,} / {:,} correlations".format(corr_index, n_correlations))

    sample_corr_df = pd.DataFrame(sample_corr_data, columns=["dataset1", "sample1", "dataset2", "sample2", "n_genes", args.method + "r"])
    sample_corr_df.sort_values(by=args.method + "r", inplace=True)
    print(sample_corr_df)
    sample_corr_df.to_csv(sample_corr_fpath, sep="\t", header=True, index=False)

    gene_corr_df = pd.DataFrame(gene_corr_data, columns=["dataset", "gene", "n_samples", args.method + "r"])
    gene_corr_df.dropna(inplace=True)
    gene_corr_df.sort_values(by=args.method + "r", inplace=True)
    print(gene_corr_df)
    gene_corr_df.to_csv(gene_corr_fpath, sep="\t", header=True, index=False)

sample_corr_df = pd.read_csv(sample_corr_fpath, sep="\t", header=0, index_col=None)
print(sample_corr_df)
gene_corr_df = pd.read_csv(gene_corr_fpath, sep="\t", header=0, index_col=None)
print(gene_corr_df)

################################################

same_sample_df = sample_corr_df.loc[sample_corr_df["sample1"] == sample_corr_df["sample2"], :].copy()
same_sample_df.sort_values(by=args.method + "r", inplace=True)
print(same_sample_df)

plot_histplot(
    df=same_sample_df,
    x=args.method + "r",
    hue="dataset1",
    palette=palette,
    xlabel=args.method + " r",
    ylabel="count",
    title=args.label1 + " vs " + args.label2,
    filename="sample_" + args.method + "r"
)
# for dataset in sample_corr_df["dataset1"].unique():
#     if not "Columbia" in dataset:
#         continue
#     plot_heatmap(
#         df=sample_corr_df.loc[(sample_corr_df["dataset1"] == dataset) & (sample_corr_df["dataset2"] == dataset), :],
#         x=args.method + "r",
#         xlabel=args.method + " r",
#         ylabel="count",
#         title=args.label1 + " vs " + args.label2,
#         filename="sample_" + dataset + "_" + args.method + "r"
#     )

plot_histplot(
    df=gene_corr_df,
    x=args.method + "r",
    hue="dataset",
    palette=palette,
    xlabel=args.method + " r",
    ylabel="count",
    title=args.label1 + " vs " + args.label2,
    filename="gene_" + args.method + "r"
)

print("END")