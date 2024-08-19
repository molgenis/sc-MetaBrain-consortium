#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse


"""
Syntax: 
./pseudobulk_pca.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--work_dir", type=str, required=True, help="")
parser.add_argument("--in_dir", type=str, required=True, help="")
parser.add_argument("--ancestry", type=str, required=False, default="EUR", help="")
parser.add_argument("--cell_level", type=str, required=False, default="L1", help="")
parser.add_argument("--cell_types", nargs="*", type=str, required=False, default=["AST", "END", "EX", "IN", "MIC", "OLI", "OPC", "PER"], help="")
parser.add_argument("--individual_aggregate", type=str, required=False, default="IID", help="")
parser.add_argument("--gte", type=str, required=True, help="")
parser.add_argument("--correct_datasets", action='store_true', help="")
parser.add_argument("--palette", type=str, required=False, default=None, help="A color palette file.")
parser.add_argument("--extensions", nargs="+", type=str, choices=["png", "pdf", "eps"], default=["png"], help="The figure file extension.. Default: 'png'.")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import numpy as np
import pandas as pd
import os
import json
from statsmodels.regression.linear_model import OLS
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

outdir = os.path.join(args.work_dir, 'plots', 'pseudobulk_pca')
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Set the right pdf font for exporting.
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Load  the color palette.
palette = None
if args.palette is not None:
    with open(args.palette) as f:
        palette = json.load(f)
    f.close()

def construct_dataset_df(gte_df):
    dataset_sample_counts = list(zip(*np.unique(gte_df.iloc[:, 2], return_counts=True)))
    dataset_sample_counts.sort(key=lambda x: -x[1])
    datasets = [csc[0] for csc in dataset_sample_counts]

    dataset_df = pd.DataFrame(0, index=gte_df.iloc[:, 1], columns=datasets)
    for dataset in datasets:
        dataset_df.loc[(gte_df.iloc[:, 2] == dataset).values, dataset] = 1
    dataset_df.index.name = "-"

    return dataset_df

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
    x = df.T.to_numpy()
    x = (x - np.mean(x, axis=0)) / np.sqrt(np.sum(x ** 2, axis=0) / max(1, (x.shape[0] - 1)))

    pca = PCA(n_components=n_components)
    pca_names = ["PC{}".format(i) for i in range(1, n_components + 1)]
    pca_df = pd.DataFrame(pca.fit_transform(x), columns=pca_names, index=df.index)
    var_df = pd.Series(pca.explained_variance_ratio_, index=pca_names)
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
            fig.savefig(os.path.join(outdir, "{}.{}".format(filename, extension)))
        plt.close()

print("Loading GTE data...")
gte_df = pd.read_csv(args.gte, sep="\t", header=None)
gte_df.columns = ["genotype_id", "expression_id", "dataset"]

print("Loading metadata and applying QC filters")
metadata = pd.read_csv(os.path.join(args.in_dir, "..", "..", "metadata.tsv.gz"), sep="\t", header=0, index_col=None, low_memory=False)
# THIS HAS TO BE IDENTICAL TO create_expression_matrices.R filter except for the cell type filter!!!!!!!!!!!
# ps. the fact that individual-sample aggregate is missing here is fine since we only use that to completely remove
# samples and those would not be in this metadata file anymore.
metadata = metadata.loc[(metadata["DropletType"] == "singlet") &
                        (metadata["tag"] == "NotOutlier") &
                        (metadata["Provided_Ancestry"] == args.ancestry) &
                        (metadata["cell_treatment"] == "UT"), :]

# Cell type counts per sample.
ct_count_per_ind = metadata.loc[:, [args.individual_aggregate, args.cell_level, "Barcode"]].groupby([args.individual_aggregate, args.cell_level]).count()
ct_count_per_ind.reset_index(drop=False, inplace=True)
ct_count_per_ind.columns = [args.individual_aggregate, "cell_type", "cell_count"]
ct_count_per_ind["expression_id"] = [individual_aggregate.split(";;")[0] for individual_aggregate in ct_count_per_ind[args.individual_aggregate]]
ct_count_per_ind.index = ct_count_per_ind["expression_id"] + "_" + ct_count_per_ind["cell_type"]
print(ct_count_per_ind)

dataset_df = None
# full_dataset_df = None
if args.correct_datasets:
    print("Creating dataset df...")
    dataset_df = construct_dataset_df(gte_df=gte_df)

    # full_dataset_df_list = []
    # for cell_type in args.cell_types:
    #     tmp_dataset_df = dataset_df.copy()
    #     tmp_dataset_df.index = ["{}_{}".format(index, cell_type) for index in tmp_dataset_df.index]
    #     full_dataset_df_list.append(tmp_dataset_df)
    # full_dataset_df = pd.concat(full_dataset_df_list, axis=0)
    # del full_dataset_df_list, tmp_dataset_df

print("Loading expression data...")
df_list = []
for cell_type in args.cell_types:
    fpath = os.path.join(args.in_dir, cell_type, "PreQC", cell_type + ".qtlInput.txt")
    if not os.path.exists(fpath):
        print("Error, could not find '{}'".format(fpath))

    df = pd.read_csv(fpath, sep="\t", header=0, index_col=0)
    print("\tLoaded '{}' with shape {}".format(os.path.basename(fpath), df.shape))

    if args.correct_datasets:
        df = remove_covariates(Y=df, X=dataset_df.iloc[:, 1:])
    pca_df, var_df = calculate_pca(df=df)
    plot_df = pca_df.merge(gte_df, left_index=True, right_on="expression_id").merge(ct_count_per_ind.loc[ct_count_per_ind["cell_type"] == cell_type], on="expression_id", how="left")
    print(plot_df)

    for hue, pal, legend in [("dataset", palette, False), ("cell_count", None, True)]:
        plot_scatterplot(
            df=plot_df,
            x="PC1",
            y="PC2",
            legend=legend,
            hue=hue,
            palette=pal,
            xlabel="PC1 [{:.0f}%]".format(var_df["PC1"]),
            ylabel="PC2 [{:.0f}%]".format(var_df["PC2"]),
            title="Merged Pseudobulk",
            filename="pseudobulk_expression_pca_{}_{}{}".format(cell_type, hue, "_dataset_corrected" if args.correct_datasets else "")
        )

    df.columns = [column + "_" + cell_type for column in df.columns]
    df_list.append(df)
df = pd.concat(df_list, axis=1).dropna()
print(df)

# if args.correct_datasets:
#     df = remove_covariates(Y=df, X=full_dataset_df.iloc[:, 1:])
pca_df, var_df = calculate_pca(df=df)
expression_ids = []
cell_types = []
for index in pca_df.index:
    expression_id, cell_type = index.split("_")
    expression_ids.append(expression_id)
    cell_types.append(cell_type)
pca_df["expression_id"] = expression_ids
pca_df["cell_type"] = cell_types
print(pca_df)
plot_df = pca_df.merge(ct_count_per_ind[["cell_count"]], left_index=True, right_index=True, how="left").merge(gte_df, on="expression_id")
print(plot_df)

print(plot_df)
plot_scatterplot(
    df=plot_df,
    x="PC1",
    y="PC2",
    hue="dataset",
    palette=palette,
    xlabel="PC1 [{:.0f}%]".format(var_df["PC1"]),
    ylabel="PC2 [{:.0f}%]".format(var_df["PC2"]),
    title="Merged Pseudobulk",
    filename="pseudobulk_expression_pca_dataset{}".format("_dataset_corrected" if args.correct_datasets else "")
)
plot_scatterplot(
    df=plot_df,
    x="PC1",
    y="PC2",
    hue="cell_type",
    palette=palette,
    xlabel="PC1 [{:.0f}%]".format(var_df["PC1"]),
    ylabel="PC2 [{:.0f}%]".format(var_df["PC2"]),
    title="Merged Pseudobulk",
    filename="pseudobulk_expression_pca_cell_type{}".format("_dataset_corrected" if args.correct_datasets else "")
)
plot_scatterplot(
    df=plot_df,
    x="PC1",
    y="PC2",
    hue="cell_count",
    legend=True,
    xlabel="PC1 [{:.0f}%]".format(var_df["PC1"]),
    ylabel="PC2 [{:.0f}%]".format(var_df["PC2"]),
    title="Merged Pseudobulk",
    filename="pseudobulk_expression_pca_cellcount{}".format("_dataset_corrected" if args.correct_datasets else "")
)


