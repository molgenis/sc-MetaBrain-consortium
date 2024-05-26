#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse


"""
Syntax: 
./check_replication_i2.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--work_dir", type=str, required=True, help="")
parser.add_argument("--merged", type=str, required=True, help="")
parser.add_argument("--replication", type=str, required=True, help="")
parser.add_argument("--cell_type", type=str, required=True, help="")
parser.add_argument("--palette", type=str, required=False, default=None, help="A color palette file.")
parser.add_argument("--extensions", nargs="+", type=str, choices=["png", "pdf", "eps"], default=["png"], help="The figure file extension.. Default: 'png'.")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import json
import gzip
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

outdir = os.path.join(args.work_dir, 'plots', 'check_replication_i2')
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

def plot(df, type, x="variable", y="value", hue=None, palette=None, xlabel="", ylabel="", title="", filename=""):
    sns.set(rc={'figure.figsize': (12, 9)})
    sns.set_style("ticks")
    fig, (ax1, ax2) = plt.subplots(nrows=1,
                                   ncols=2,
                                   gridspec_kw={"width_ratios": [0.99, 0.01]})
    sns.despine(fig=fig, ax=ax1)

    if type == "box":
        sns.violinplot(x=x,
                       y=y,
                       hue=hue,
                       data=df,
                       palette=palette,
                       cut=0,
                       dodge=True,
                       ax=ax1)

        plt.setp(ax1.collections, alpha=.75)

        sns.boxplot(x=x,
                    y=y,
                    hue=hue,
                    data=df,
                    palette=palette,
                    dodge=True,
                    ax=ax1)

        for patch in ax1.artists:
            r, g, b, a = patch.get_facecolor()
            patch.set_facecolor((r, g, b, .75))

        if ax1.get_legend() is not None:
            ax1.get_legend().remove()

        ax1.set_xticklabels(ax1.get_xticklabels(), rotation=0)
    elif type == "kde":
        sns.kdeplot(
            data=df,
            x=x,
            hue=hue,
            fill=True,
            palette=palette,
            ax=ax1
        )
    else:
        raise NotImplementedError("")

    ax2.set_axis_off()
    if palette is not None and hue is not None:
        handles = []
        for key, value in palette.items():
            if key in df[hue].unique():
                handles.append(mpatches.Patch(color=value, label=key))
        ax2.legend(handles=handles, fontsize=10)

    plt.tight_layout()
    for extension in args.extensions:
        fig.savefig(os.path.join(outdir, "{}.{}".format(filename, extension)))
    plt.close()

print("Loading merged data")
merged_df = pd.read_csv(args.merged, sep="\t", header=0, index_col=0)
print("\tLoaded '{}' with shape {}".format(os.path.basename(args.merged), merged_df.shape))
selection = set(merged_df.index.tolist())

print("Loading replication data")
column_dict = None
dataset_r_column = None
dataset_z_column = None
dataset_n_column = None
columns = ["Gene", "SNP", "MetaI2"]
data = []
nrows = 0
with gzip.open(args.replication, "rt") as f:
    for i, line in enumerate(f):
        if i % 1e6 == 0:
            print("  Parsed {:,} rows, saved {:,}".format(i, nrows))
        values = line.rstrip("\n").split("\t")
        if i == 0:
            column_dict = {}
            for i, value in enumerate(values):
                column_dict[value] = i
                if value.startswith("DatasetCorrelationCoefficients"):
                    dataset_r_column = value
                    columns.extend(["CorrCoef-" + value for value in value.replace("DatasetCorrelationCoefficients(", "").replace(")", "").split(";")])
                elif value.startswith("DatasetZScores"):
                    dataset_z_column = value
                    columns.extend(["ZScore-" + value for value in value.replace("DatasetZScores(", "").replace(")", "").split(";")])
                elif value.startswith("DatasetSampleSizes"):
                    dataset_n_column = value
                    columns.extend(["SampleSize-" + value for value in value.replace("DatasetSampleSizes(", "").replace(")", "").split(";")])
            continue

        # if len(data) == 100:
        #     break

        gene_snp = values[column_dict["Gene"]] + "_" + values[column_dict["SNP"]]
        if gene_snp not in selection:
            continue

        extended_values = []

        data.append([values[column_dict["Gene"]], values[column_dict["SNP"]], values[column_dict["MetaI2"]]] + values[column_dict[dataset_r_column]].split(";") + values[column_dict[dataset_z_column]].split(";") + values[column_dict[dataset_n_column]].split(";"))
        nrows += 1
    print("  Parsed {:,} rows, saved {:,}".format(i, nrows))
f.close()


replication_df = pd.DataFrame(data, columns=columns).replace('-', np.nan)
replication_df.index = replication_df["Gene"] + "_" + replication_df["SNP"]
for column in replication_df.columns[2:]:
    replication_df[column] = replication_df[column].astype(float)
    if column.startswith("ZScore-"):
        replication_df["ABS-" + column] = replication_df[column].abs()
print("\tLoaded '{}' with shape {}".format(os.path.basename(args.replication), replication_df.shape))
print(replication_df)

df = merged_df.merge(replication_df[["MetaI2"] + [column for column in replication_df.columns if column.startswith("ABS-ZScore-")]], left_index=True, right_index=True, how="left")
df = df.loc[~df["MetaI2"].isna(), ]
df["signif"] = "Neither"
df.loc[(df["Bryois FDR"] < 0.05) & (df["mbQTL FDR"] >= 0.05), "signif"] = "Discovery"
df.loc[(df["Bryois FDR"] >= 0.05) & (df["mbQTL FDR"] < 0.05), "signif"] = "Replication"
df.loc[(df["Bryois FDR"] < 0.05) & (df["mbQTL FDR"] < 0.05), "signif"] = "Both"
print(df)

dfm = df.melt(id_vars=["signif"], value_vars=[col for col in df.columns if col.startswith("ABS-ZScore-")])
dfm["variable"] = [value.replace("ABS-ZScore-", "") for value in dfm["variable"]]
plot(
    df=dfm,
    type="box",
    x="signif",
    hue="variable",
    palette=palette,
    xlabel="dataset",
    ylabel="Z-score",
    title="eQTL z-score - replicating or not",
    filename=args.cell_type + "_eqtl_zscore_replication_subset"
)

plot(
    df=df,
    type="kde",
    x="MetaI2",
    hue="signif",
    xlabel="MetaI2",
    ylabel="Density",
    title="eQTL I2 - replicating or not",
    filename=args.cell_type + "_eqtl_i2_replication"
)

df = df.loc[df["signif"] != "Neither", :]
plot(
    df=df,
    type="kde",
    x="MetaI2",
    hue="signif",
    xlabel="MetaI2",
    ylabel="Density",
    title="eQTL I2 - replicating or not",
    filename=args.cell_type + "_eqtl_i2_replication_subset"
)


