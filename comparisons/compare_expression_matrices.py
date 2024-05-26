#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse


"""
Syntax: 
./compare_expression_matrices.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--expr1", type=str, required=True, help="")
parser.add_argument("--label1", type=str, required=True, help="")
parser.add_argument("--expr2", type=str, required=True, help="")
parser.add_argument("--label2", type=str, required=True, help="")
parser.add_argument("--genes", nargs="*", type=str, required=False, help="")
parser.add_argument("--gte", type=str, required=False, help="")
parser.add_argument("--exclude", type=str, required=False, help="")
parser.add_argument("--palette", type=str, required=False, help="")
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

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

print("Loading data")
expr_df1 = pd.read_csv(args.expr1, sep="\t", header=0, index_col=0)
print("\tExpression matrix 1: {}".format(expr_df1.shape))

expr_df2 = pd.read_csv(args.expr2, sep="\t", header=0, index_col=0)
print("\tExpression matrix 2: {}".format(expr_df2.shape))

gte_df = None
dataset_dict = {}
if args.gte is not None:
    gte_df = pd.read_csv(args.gte, sep="\t", header=0, index_col=None)
    dataset_dict = dict(zip(gte_df["individual_id"], gte_df["dataset"]))

overlapping_genes = list(set(expr_df2.index).intersection(set(expr_df1.index)))
n_genes = len(overlapping_genes)

overlapping_samples = list(set(expr_df2.columns).intersection(set(expr_df1.columns)))
if args.exclude is not None:
    overlapping_samples = [sample for sample in overlapping_samples if dataset_dict[sample] != args.exclude]

n_samples = len(overlapping_samples)
print("\tOverlapping genes: {}".format(n_genes))
print("\tOverlapping samples: {}".format(n_samples))

expr_m1 = expr_df1.loc[overlapping_genes, overlapping_samples].to_numpy()
expr_m2 = expr_df2.loc[overlapping_genes, overlapping_samples].to_numpy()
del expr_df1, expr_df2

print("Correlating data")
correlations = []
i = 0
for i, gene in enumerate(overlapping_genes):
    if i ==0 or i % 1000 == 0:
        print("\tParsed {:,} / {:,} genes".format(i, n_genes), end='\r')

    # Pearson r.
    pearson_coef, _ = stats.pearsonr(x=expr_m1[i, :], y=expr_m2[i, :])

    # Spearman r.
    spearman_coef, _ = stats.spearmanr(a=expr_m1[i, :], b=expr_m2[i, :])

    # Residual sum of squares.
    res = expr_m1[i, :] - expr_m2[i, :]
    res_squared = res * res
    rss = np.sum(res_squared)

    correlations.append([i, gene, pearson_coef, abs(pearson_coef), spearman_coef, abs(spearman_coef), rss])
print("\tParsed {:,} / {:,} genes".format(i, n_genes))

corr_df = pd.DataFrame(correlations, columns=["index", "gene", "pearson_coef", "abs_pearson_coef", "spearman_coef", "abs_spearman_coef", "rss"])
corr_df.sort_values(by="abs_pearson_coef", inplace=True)
print(corr_df)
for col in corr_df.columns:
    if col in ["index", "gene"]:
        continue
    print("{}:\tmin:{:.4f}\tmax:{:.4f}\tmean:{:.4f}\tstd:{:.4f}\tmedian:{:.4f}".format(col, corr_df[col].min(), corr_df[col].max(), corr_df[col].mean(), corr_df[col].std(), corr_df[col].median()))


def plot_histplot(df, x="x", xlabel="", ylabel="", title="", filename="plot"):
    sns.set_style("ticks")

    sns.set(rc={'figure.figsize': (9, 12)})
    sns.set_style("ticks")
    fig, ax = plt.subplots()
    sns.despine(fig=fig, ax=ax)

    g = sns.histplot(x=x,
                    data=df,
                     color="#808080",
                    ax=ax)

    ax.set_title(title,
                 fontsize=22,
                 fontweight='bold')
    ax.set_ylabel(ylabel,
                  fontsize=14,
                  fontweight='bold')
    ax.set_xlabel(xlabel,
                  fontsize=14,
                  fontweight='bold')

    ax.annotate(
        "N = {:,}".format(df.shape[0]),
        xy=(0.03, 0.95),
        xycoords=ax.transAxes,
        color="#808080",
        fontsize=14,
        fontweight='bold'
    )
    ax.annotate(
        "mean = {:.2f}".format(df[x].mean()),
        xy=(0.03, 0.90),
        xycoords=ax.transAxes,
        color="#808080",
        fontsize=14,
        fontweight='bold'
    )

    plt.tight_layout()
    fpath = os.path.join(args.outdir, "plot", "{}_{}_{}.png".format(args.label1, args.label2, filename))
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

    facecolors = "#808080"
    if hue is not None:
        facecolors = df["hue"].map(palette)

    sns.regplot(x=x, y=y, data=df, ci=None,
                scatter_kws={'facecolors': facecolors,
                             # 's': 10,
                             # 'alpha': 0.2,
                             'linewidths': 0},
                line_kws={"color": "#000000"},
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
        for key, color in palette.items():
            if key in df[hue].values.tolist():
                handles.append(mpatches.Patch(color=color, label=key))
        ax2.legend(handles=handles, loc="center")
    ax2.set_axis_off()

    plt.tight_layout()
    fpath = os.path.join(args.outdir, "plot", "{}_{}_{}.png".format(args.label1, args.label2, filename))
    print("Saving {}".format(os.path.basename(fpath)))
    fig.savefig(fpath)
    plt.close()

for column in corr_df.columns:
    if column in ["index", "gene", "rss"]:
        continue

    plot_histplot(
        df=corr_df,
        x=column,
        xlabel=column,
        ylabel="count",
        title=args.label1 + " vs " + args.label2,
        filename=column
    )

if args.genes is not None:
    for gene in args.genes:
        index = overlapping_genes.index(gene)
        plot_df = pd.DataFrame({"expr1": expr_m1[index, :], "expr2": expr_m2[index, :]}, index=overlapping_samples)
        if dataset_dict:
            plot_df["hue"] = plot_df.index.map(dataset_dict)
        plot_regplot(
            df=plot_df,
            x="expr1",
            y="expr2",
            hue="hue" if dataset_dict else None,
            palette=palette if dataset_dict else None,
            xlabel=args.label1,
            ylabel=args.label2,
            title=gene,
            filename=gene + "_expression"
        )

print("END")