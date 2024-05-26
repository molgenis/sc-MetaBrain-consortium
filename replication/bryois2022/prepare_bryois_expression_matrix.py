#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse

"""
Syntax: 
./prepare_bryois_expression_matrix.py -h
"""

parser = argparse.ArgumentParser(description="")
parser.add_argument("--indir", type=str, required=True, help="")
parser.add_argument("--gte", type=str, required=True, help="")
parser.add_argument("--palette", type=str, required=False, default=None, help="A color palette file.")
parser.add_argument("--outdir", type=str, required=True, help="")
parser.add_argument("--extensions", nargs="+", type=str, choices=["png", "pdf", "eps"], default=["png"], help="The figure file extension.. Default: 'png'.")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import numpy as np
import pandas as pd
import json
import os
from statsmodels.regression.linear_model import OLS
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plot_outdir = os.path.join(args.outdir, 'plot')
if not os.path.exists(plot_outdir):
    os.makedirs(plot_outdir)

# Set the right pdf font for exporting.
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Load  the color palette.
palette = None
if args.palette is not None:
    with open(args.palette) as f:
        palette = json.load(f)
    f.close()


def create_indicator_df(df, column):
    value_counts = list(zip(*np.unique(df[column], return_counts=True)))
    value_counts.sort(key=lambda x: -x[1])
    values = [csc[0] for csc in value_counts]

    indicator_df = pd.DataFrame(0, index=df.index, columns=values)
    for value in values:
        indicator_df.loc[(df[column] == value).values, value] = 1
    indicator_df.index.name = "-"

    return indicator_df.iloc[:, 1:]

def remove_covariates(Y, X):
    # print(Y)
    # print(X.loc[Y.columns, :])
    Y_m = Y.to_numpy()
    X_m = np.hstack((np.ones((Y_m.shape[1], 1)), X.loc[Y.columns, :].to_numpy()))

    Y_m_resid = np.empty(Y_m.shape, dtype=np.float64)
    n_rows = Y_m.shape[0]
    i = 0
    for i in range(n_rows):
        if i % 5000 == 0:
            print("\t\t{:,}/{:,} rows processed [{:.2f}%]".format(i, (n_rows - 1), (100 / (n_rows - 1)) * i))

        # Mask the Nan values.
        sample_mask = ~np.isnan(Y_m[i, :])

        # Calculate residuals using OLS.
        Y_m_resid[i, ~sample_mask] = np.nan
        Y_m_resid[i, sample_mask] = OLS(Y_m[i, sample_mask], X_m[sample_mask, :]).fit().resid
    print("\t\t{:,}/{:,} rows processed [{:.2f}%]".format(i, (n_rows - 1), (100 / (n_rows - 1)) * i))
    return pd.DataFrame(Y_m_resid, index=Y.index, columns=Y.columns)

def calculate_pca(df, n_components=2):
    # df = (df - df.mean(axis=0)) / df.std(axis=0)

    pca = PCA(n_components=n_components)
    pca.fit(df)

    pca_names = ["PC{}".format(i) for i in range(1, n_components + 1)]

    pca_df = pd.DataFrame(pca.components_, index=pca_names, columns=df.columns).T
    var_df = pd.Series(pca.explained_variance_, index=pca_names)
    return pca_df, var_df

def plot_scatterplot(df, x="x", y="y", legend=False, hue=None, palette=None,
                     xlabel="", ylabel="", title="", filename="scatterplot"):
        fig, (ax, legend_ax) = plt.subplots(nrows=1, ncols=2, figsize=(8, 6), gridspec_kw={"width_ratios": [0.9, 0.1]})
        sns.set(color_codes=True)
        sns.set_style("ticks")

        sns.despine(fig=fig, ax=ax)
        sns.scatterplot(x=x,
                        y=y,
                        hue=hue,
                        data=df,
                        palette=palette,
                        linewidth=0,
                        legend=legend,
                        ax=ax)

        ax.set_title(title,
                     fontsize=20,
                     fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=12,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=12,
                      fontweight='bold')

        legend_ax.set_axis_off()
        if palette is not None and hue is not None:
            handles = []
            for key, value in palette.items():
                if key in df[hue].unique():
                    handles.append(mpatches.Patch(color=value, label="{} [N={:,}]".format(key, df.loc[df[hue] == key].shape[0])))
            legend_ax.legend(handles=handles, fontsize=10)

        plt.tight_layout()
        for extension in args.extensions:
            fig.savefig(os.path.join(plot_outdir, "{}.{}".format(filename, extension)))
        plt.close()


print("Loading GTE")
gte_df = pd.read_csv(args.gte, sep="\t", header=None, index_col=None)
gte_df.columns = ["genotype_id", "expression_id", "dataset"]
print(gte_df)

print("Loading sample keys")
ms_key_df = pd.read_csv(os.path.join(args.indir, "meta_ms_public_with_key.txt"), sep="\t", header=0, index_col=None)
ad_key_df = pd.read_csv(os.path.join(args.indir, "meta_ad_public_with_key.txt"), sep="\t", header=0, index_col=None)
key_df = pd.concat([ms_key_df, ad_key_df], axis=0)
print(key_df)
del ms_key_df, ad_key_df

# Merge the key info to the gte frame.
key_df = key_df[["individual_id", "individual_id_anon"]].drop_duplicates().dropna()
gte_df = gte_df.merge(key_df, left_on="expression_id", right_on="individual_id_anon", how="left")
sample_trans = dict(zip(gte_df["individual_id"], gte_df["expression_id"]))
print(gte_df)

print("Loading bed")
bed_df = pd.read_csv(os.path.join(args.indir, "Excitatory.neurons.bed.gz"), sep="\t", header=0, index_col=None)
bed_df.index = [gene.split("_")[0] for gene in bed_df["ID"]]
loc_df = bed_df.iloc[:, :4].copy()
expr_df = bed_df.iloc[:, 4:].copy()
if len(set(expr_df.index)) != expr_df.shape[0]:
    print("Warning, indices are not unique.")

print(expr_df)
print(bed_df)

print("Loading and plotting covariate matrix")
cov_df = pd.read_csv(os.path.join(args.indir, "Excitatory.neurons.cov.txt.gz"), sep="\t", header=0, index_col=0).T
for column in cov_df.columns:
    if column in ["study", "diagnosis"]:
        continue
    cov_df[column] = cov_df[column].astype(float)

plot_scatterplot(
    df=cov_df,
    x="PC1",
    y="PC2",
    hue="study",
    palette=palette,
    xlabel="PC1",
    ylabel="PC2",
    title="Bryois Excitatory.neurons Geno",
    filename="Excitatory.neurons.cov.pca.geno"
)
plot_scatterplot(
    df=cov_df,
    x="PC1_exp",
    y="PC2_exp",
    hue="study",
    palette=palette,
    xlabel="PC1_exp",
    ylabel="PC2_exp",
    title="Bryois Excitatory.neurons",
    filename="Excitatory.neurons.cov.pca"
)

print("Correcting expression matrix")
correct_df = cov_df.loc[:, [column for column in cov_df.columns if column not in ["study", "diagnosis"]]].copy()
study_indicator_df = create_indicator_df(df=cov_df, column="study")
diagnosis_indicator_df = create_indicator_df(df=cov_df, column="diagnosis")
correct_df = correct_df.merge(study_indicator_df, left_index=True, right_index=True).merge(diagnosis_indicator_df, left_index=True, right_index=True)
print(correct_df)
del cov_df

res_df = remove_covariates(Y=expr_df, X=correct_df)

print("Reformatting columns")
expr_df.columns = [col if col in gte_df["expression_id"].values else sample_trans[col] for col in expr_df.columns]
res_df.columns = [col if col in gte_df["expression_id"].values else sample_trans[col] for col in expr_df.columns]

for df, label in [(expr_df, "Raw"), (res_df, "CovariatesRemovedOLS")]:
    overlap = set(gte_df.iloc[:, 0]).intersection(set(df.columns))
    print("{} overlap: {:,} / {:,}".format(label, len(overlap), gte_df.shape[0]))
    missing = set(gte_df.iloc[:, 0]).difference(set(df.columns))
    print("{} missing: {}".format(label, ", ".join(missing)))

print("Saving files")
gte_df = gte_df.loc[gte_df["expression_id"].isin(overlap), ["genotype_id", "expression_id", "dataset"]]
print(gte_df)
gte_df.to_csv(os.path.join(args.outdir, "Excitatory.neurons.gte.txt"), sep="\t", header=False, index=False)
gte_df["dataset"] = "dataset"
gte_df.to_csv(os.path.join(args.outdir, "Excitatory.neurons.gte.nodataset.txt"), sep="\t", header=False, index=False)

bed_df.to_csv(os.path.join(args.outdir, "Excitatory.neurons.bed.txt.gz"), sep="\t", header=True, index=True, compression="gzip")
expr_df.loc[:, gte_df["expression_id"]].to_csv(os.path.join(args.outdir, "Excitatory.neurons.Exp.txt.gz"), sep="\t", header=True, index=True, compression="gzip")
res_df.loc[:, gte_df["expression_id"]].to_csv(os.path.join(args.outdir, "Excitatory.neurons.Exp.CovariatesRemovedOLS.txt.gz"), sep="\t", header=True, index=True, compression="gzip")

##########

print("Plotting PCA")
for df, label in [(expr_df, "Raw"), (res_df, "CovariatesRemovedOLS")]:
    pca_df, var_df = calculate_pca(df=df)
    plot_df = pca_df.merge(gte_df, left_index=True, right_on="expression_id")

    plot_scatterplot(
        df=plot_df,
        x="PC1",
        y="PC2",
        hue="dataset",
        palette=palette,
        xlabel="PC1 [{:.0f}%]".format(var_df["PC1"]),
        ylabel="PC2 [{:.0f}%]".format(var_df["PC2"]),
        title="Bryois Excitatory.neurons",
        filename="Excitatory.neurons.expr.{}.pca".format(label)
    )