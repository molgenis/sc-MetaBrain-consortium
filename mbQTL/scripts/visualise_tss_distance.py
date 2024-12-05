#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--data", required=True, type=str, help="")
parser.add_argument("--nbins", required=False, type=int, default=5, help="")
parser.add_argument("--signif_col", required=False, type=str, default=None, help="")
parser.add_argument("--alpha", required=False, type=float, default=0.05, help="")
parser.add_argument("--extensions", nargs="+", type=str, choices=["png", "pdf", "eps"], default=["png"], help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

os.makedirs(os.path.dirname(args.out), exist_ok=True)

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def plot(data, x="x", hue=None, hline=None, vline=None, palette=None, title="", xlabel="", filename="plot"):
    plt.rcParams["figure.figsize"] = (12, 9)
    fig, ax = plt.subplots()

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    if hue is not None and palette is not None and not isinstance(palette, dict):
        hue_values = data[hue].unique()
        cmap = matplotlib.colormaps.get_cmap(palette)
        palette = {}
        for i, hue_value in enumerate(hue_values):
            palette[hue_value] = matplotlib.colors.rgb2hex(cmap(1 / (i + 1)))

    if hue is not None:
        for hue_value in data[hue].unique():
            hue_color = palette[hue_value]
            s = data.loc[data[hue]==hue_value, x]
            if s.shape[0] <= 1:
                continue
            s.plot.kde(color=hue_color)
    else:
        data[x].plot.kde(color="#000000")

    # Clip the distribution to the nearest round number.
    min_x = data[x].min()
    max_x = data[x].max()
    n_numbers = max(len(str(abs(min_x))), len(str(abs(min_x))))
    plt.xlim(round(min_x, -1 * n_numbers), round(max_x, -1 * n_numbers))

    if vline is not None:
        for value in vline:
            ax.axvline(value, ls='--', color="#808080", alpha=0.5, zorder=-1, linewidth=2)
    if hline is not None:
        for value in hline:
            ax.axhline(value, ls='--', color="#808080", alpha=0.5, zorder=-1, linewidth=2)

    y_pos = 0.95
    y_pos_step = -0.04
    ax.annotate(
        'N = {:,}'.format(data.shape[0]),
        xy=(0.03, y_pos),
        xycoords=ax.transAxes,
        color="#000000",
        fontsize=12,
        fontweight='bold')
    y_pos += y_pos_step

    if hue is not None:
        for hue_value in data[hue].unique():
            color = "#000000"
            if palette is not None and hue_value in palette:
                color = palette[hue_value]
            n = (data[hue] == hue_value).sum()
            ax.annotate(
                'N {} = {:,} [{:.2f}%]'.format(hue_value, n, (100 / data.shape[0]) * n),
                xy=(0.03, y_pos),
                xycoords=ax.transAxes,
                color=color,
                fontsize=12,
                fontweight='bold')
            y_pos += y_pos_step

    ax.set_title(title,
                 fontsize=16,
                 color="#000000",
                 weight='bold')
    ax.set_xlabel(xlabel,
                  fontsize=12,
                  color="#000000",
                  weight='bold')
    ax.set_ylabel("density",
                  fontsize=12,
                  color="#000000",
                  weight='bold')

    fig.tight_layout()
    for extension in args.extensions:
        plt.savefig(args.out + filename + '.' + extension, bbox_inches="tight")


######################################################

print("Loading data ...")
df = pd.read_csv(args.data, sep="\t", header=0, index_col=None)
print("\tLoaded dataframe: {} with shape: {}".format(os.path.basename(args.data), df.shape))

if args.signif_col is not None and args.signif_col in df.columns:
    df = df.loc[df[args.signif_col] < args.alpha, :]
    print("\tFiltering on {}<{} leaves {:,} effects".format(args.signif_col, args.alpha, df.shape[0]))

if df.shape[0] == 0:
    for extension in args.extensions:
        plt.savefig(args.out + '-TSSDistance-VariantType.' + extension)
        plt.savefig(args.out + '-TSSDistance-PValBins.' + extension)
    exit()

# Adding variant type column.
df["VariantType"] = df["SNPAlleles"].str.len().map({0: "NA", 1: "NA", 2: "NA", 3: "SNP"}, na_action=None)
df["VariantType"] = df["VariantType"].fillna("INDEL")
n_snp = df["VariantType"] == "SNP"
n_indel = df["VariantType"] == "INDEL"
print("\tSNP: {:,} / {:,} [{:.2f}%] of variants".format(n_snp.sum(), df.shape[0], (100 / df.shape[0]) * n_snp.sum()))
print("\tINDEL: {:,} / {:,} [{:.2f}%] of effects".format(n_indel.sum(), df.shape[0], (100 / df.shape[0]) * n_indel.sum()))

# Adding significance bins.
p_value_col = None
if "MetaP" in df.columns:
    p_value_col = "MetaP"
elif "MetaP-Random" in df.columns:
    p_value_col = "MetaP-Random"
else:
    print("Error, could not identify p-value column.")
df["Pbin"] = pd.qcut(df[p_value_col], q=args.nbins)

print("Plotting the distance to the TSS ...")
# only works for cis-eQTLs
cis_mask = (df["QTLType"] == "CIS") & (df["GeneChr"] == df["SNPChr"])
print("\t{:,} / {:,} [{:.2f}%] cis-eQTL effects".format(cis_mask.sum(), df.shape[0], (100 / df.shape[0]) * cis_mask.sum()))
tss_df = df.loc[cis_mask, ["VariantType", "Pbin", "GenePos", "SNPPos"]].copy()
tss_df["distance"] = tss_df["GenePos"] - tss_df["SNPPos"]
del df

plot(
    data=tss_df,
    x="distance",
    hue="VariantType",
    xlabel="distance",
    title="Distance to TSS (variant type)",
    palette={"NA": "#808080", "SNP": "#56B4E9", "INDEL": "#D55E00"},
    vline=[0],
    filename="-TSSDistance-VariantType"
)

plot(
    data=tss_df,
    x="distance",
    hue="Pbin",
    xlabel="distance",
    title="Distance to TSS ({} bins)".format(p_value_col),
    palette="viridis",
    vline=[0],
    filename="-TSSDistance-PValBins"
)

print("\nEND")