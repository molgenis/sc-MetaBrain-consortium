#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import pandas as pd
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--data", required=True, type=str, help="")
parser.add_argument("--transpose", action="store_true", default=False, help="")
parser.add_argument("--gte", required=False, type=str, default=None, help="")
parser.add_argument("--center", action="store_true", default=False, help="")
parser.add_argument("--scale", action="store_true", default=False, help="")
parser.add_argument("--eval_n_pcs", required=False, type=int, default=3, help="")
parser.add_argument("--zscore", required=False, type=int, default=None, help="")
parser.add_argument("--data_out", required=True, type=str, help="")
parser.add_argument("--plot_out", required=False, type=str, default=None, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

if args.plot_out is None:
    args.plot_out = args.data_out

os.makedirs(os.path.dirname(args.data_out), exist_ok=True)
os.makedirs(os.path.dirname(args.plot_out), exist_ok=True)

COLORMAP = ["#000000", "#E69F00",  "#56B4E9",  "#CC79A7",
            "#F0E442",  "#0072B2",  "#D55E00",  "#009E73",
            "#E8E8E8",  "#D3D3D3",  "#A9A9A9", "#808080"]

def prcomp(x, retx=True, center=True, scale=False):
    """
    https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/prcomp
    """
    n_samples, n_features = x.shape

    mean_ = False
    if center:
        mean_ = np.mean(x, axis=0)
        x -= mean_

    scale_ = False
    if scale:
        scale_ = np.sqrt(np.sum(x ** 2, axis=0) / max(1, (n_samples - 1)))
        x /= scale_

    # Perform singular value decomposition. FUll matrix is false since the it is not square.
    U, S, Vt = linalg.svd(x, full_matrices=False)

    # flip eigenvectors' sign to enforce deterministic output
    # columns of u, rows of v
    max_abs_cols = np.argmax(np.abs(U), axis=0)
    signs = np.sign(U[max_abs_cols, range(U.shape[1])])
    U *= signs
    Vt *= signs[:, np.newaxis]

    out = {
        "sdev": S / np.sqrt(n_samples - 1),
        "rotation": np.transpose(Vt),
        "center": mean_,
        "scale": scale_
    }
    if retx:
        out["x"] = U * S

    return out

def plot_scree(data, x, y1, y2, lines=None, xlabel='', ylabel1='', ylabel2='', title='', filename='PCA'):
    plt.rcParams["figure.figsize"] = (12, 9)
    fig, ax = plt.subplots()

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.scatter(data[x], data[y1], c="#000000", alpha=0.5, linewidths=0)
    plt.plot(data[x], data[y1], c="#000000", alpha=0.5)

    for (axis, value, color) in lines:
        if axis == "x":
            ax.axvline(value, ls='--', color=color, alpha=0.5, zorder=-1, linewidth=2)
        elif axis == "y1":
            ax.axhline(value, ls='--', color=color, alpha=0.5, zorder=-1, linewidth=2)

    ax.set_title(title,
                 fontsize=16,
                 color="#000000",
                 weight='bold')
    ax.set_xlabel(xlabel,
                  fontsize=12,
                  color="#000000",
                  weight='bold')
    ax.set_ylabel(ylabel1,
                  fontsize=12,
                  color="#000000",
                  weight='bold')

    if y2 is not None:
        ax2 = ax.twinx()
        plt.scatter(data[x], data[y2], c="#000000", alpha=0.5, linewidths=0)
        plt.plot(data[x], data[y2], c="#000000", alpha=0.5)

        for (axis, value, color) in lines:
            if axis == "y2":
                ax2.axhline(value, ls='--', color=color, alpha=0.5, zorder=-1, linewidth=2)

        ax2.set_ylabel(ylabel2,
                       fontsize=12,
                       color="#000000",
                       weight='bold')

    fig.tight_layout()
    plt.savefig(args.plot_out + filename + '.png', bbox_inches="tight")


def plot_embedding(data, z=None, label_mask=None, annot=None, title='', filename='PCA'):
    if annot is None:
        annot = {}

    columns = [column for column in data.columns if column not in [z, label_mask]]
    ncolumns = len(columns)

    c = None
    cmap = None
    if z is not None:
        values = data[z].unique()
        if len(values) <= len(COLORMAP):
            cmap = dict(zip(values, COLORMAP[:len(values)]))
            c = data[z].map(cmap)

    plt.rcParams["figure.figsize"] = (4 * ncolumns, 4 * ncolumns)
    fig, axs = plt.subplots(nrows=ncolumns, ncols=ncolumns, sharex='col', sharey='row')

    for i, column1 in enumerate(columns):
        for j, column2 in enumerate(columns):
            ax = axs[i, j]
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)

            if i == j:
                ax.set_axis_off()
                ax.annotate(
                    annot[column1] if column1 in annot else column1,
                    xy=(0.5, 0.5),
                    xycoords=ax.transAxes,
                    color="#000000",
                    fontsize=16,
                    ha = 'center',
                    fontweight='bold')
                continue
            elif i == 0 and j == len(columns) - 1:
                ax.set_axis_off()
                if cmap is not None:
                    handles = []
                    for key, color in cmap.items():
                        handles.append(mpatches.Patch(color=color, label="{} [n={:,}]".format(key, sum(data[z] == key))))
                    ax.legend(handles=handles, loc="center")
            elif j > i:
                ax.set_axis_off()
                continue
            else:
                ax.scatter(data[column2], data[column1], c=c, alpha=0.5, linewidths=0)
                ax.axvline(0, ls='--', color="#000000", alpha=0.5, zorder=-1, linewidth=2)
                ax.axhline(0, ls='--', color="#000000", alpha=0.5, zorder=-1, linewidth=2)

                if label_mask is not None:
                    for index, point in data.loc[data[label_mask], :].iterrows():
                        color = "#b22222"
                        if cmap is not None and cmap[point[z]] == "#D55E00":
                            color = "#000000"

                        ax.annotate(
                            str(index),
                            xy=(point[column2], point[column1]),
                            color=color,
                            fontsize=8,
                            ha='center',
                            va='bottom',
                            fontweight='bold')

    fig.suptitle(title,
                 fontsize=16,
                 color="#000000",
                 weight='bold')
    fig.tight_layout()
    plt.savefig(args.plot_out + filename + '.png', bbox_inches="tight")

print("Loading data...")
df = pd.read_csv(args.data, sep="\t", header=0, index_col=0)

if df.shape[0] == 0 or df.isnull().values.any() or df.isin([np.inf, -np.inf]).values.any():
    print("Error, input must not be empty or contain infs / NaNs")
    exit()

# Transform to numpy for speed.
samples = df.index.to_numpy()
features = df.columns.to_numpy()
m = df.to_numpy()
del df

if args.transpose:
    print("\nTransposing data...")
    m = np.transpose(m)
    tmp_samples = samples
    samples = features
    features = tmp_samples
    del tmp_samples

gte_df = None
if args.gte is not None:
    print("\nFiltering data...")
    gte_df = pd.read_csv(args.gte, sep="\t", header=None, index_col=None)
    keep = set(gte_df.iloc[:, 1].values)
    print("  Loaded {:,} samples".format(len(keep)))
    samples_mask = np.array([sample in keep for sample in samples], dtype=bool)
    features_mask = np.array([feature in keep for feature in features], dtype=bool)
    if np.sum(samples_mask) != 0 and np.sum(features_mask) == 0:
        print("  Kept {:,} / {:,} samples".format(np.sum(samples_mask), len(samples)))
        m = m[samples_mask, :]
        samples = samples[samples_mask]
    elif np.sum(samples_mask) == 0 and np.sum(features_mask) != 0:
        print("  Kept {:,} / {:,} features".format(np.sum(features_mask), len(features)))
        m = m[:, features_mask]
        features = features[features_mask]
    else:
        print("Error, GTE file does not match indices nor colnames.")
        exit()

print("\nRemoving zero variance columns...")
features_mask = np.std(m, axis=0) != 0
m = m[:, features_mask]
features = features[features_mask]
print("  Kept {:,} / {:,} features".format(np.sum(features_mask), np.size(features_mask)))

print("\nCalculating Principal Components...")
pca = prcomp(x=m, center=args.center, scale=args.scale)
pca_indices = ["PC{}_exp".format(i) for i in range(1, min(len(samples), len(features)) + 1)]
projection_df = pd.DataFrame(pca["x"], index=samples, columns=pca_indices).astype(float).T
rotation_df = pd.DataFrame(pca["rotation"], columns=pca_indices)
expl_var_df = pd.DataFrame({"expl var %": ((pca["sdev"] ** 2) / np.sum(pca["sdev"] ** 2)) * 100}, index=pca_indices)
expl_var_df["cum expl var %"] = expl_var_df["expl var %"].cumsum()

# print(projection_df)
# print(rotation_df)
# print(expl_var_df)

print("\nIdentifying outlier samples (abs z-score >= {})...".format(args.zscore))
proj_eval_df = projection_df.iloc[:args.eval_n_pcs, :].T
proj_eval_df = (proj_eval_df - proj_eval_df.mean(axis=0)) / proj_eval_df.std(axis=0)
proj_eval_df["include"] = True
proj_eval_df.loc[proj_eval_df.abs().max(axis=1) >= args.zscore, "include"] = False
for i in range(args.eval_n_pcs):
    index = proj_eval_df.columns[i]
    n_outliers = (proj_eval_df.iloc[:, i].abs() >= args.zscore).sum()
    print("\t{}: {:,} outlier samples".format(index, n_outliers))
print("\ttotal: {:,} outlier samples".format(proj_eval_df.shape[0] - proj_eval_df["include"].sum()))

include_samples = proj_eval_df.index[proj_eval_df["include"]]
exclude_samples = proj_eval_df.index[~proj_eval_df["include"]]
if gte_df is None:
    include_df = pd.DataFrame({0: include_samples})
    exclude_df = pd.DataFrame({0: exclude_samples})
else:
    include_df = gte_df.loc[gte_df.iloc[:, 1].isin(include_samples), :]
    exclude_df = gte_df.loc[gte_df.iloc[:, 1].isin(exclude_samples), :]
del proj_eval_df, include_samples

print("\nSaving results...")
projection_df.to_csv(args.data_out + "Pcs.txt.gz", sep="\t", header=True, index=True, compression="gzip")
rotation_df.to_csv(args.data_out + "Pcs_rot.txt.gz", sep="\t", header=True, index=True, compression="gzip")
expl_var_df.to_csv(args.data_out + "Pcs_var.txt.gz", sep="\t", header=True, index=True, compression="gzip")
include_df.to_csv(args.data_out + "samples_include.txt", sep="\t", header=False, index=False)
exclude_df.to_csv(args.data_out + "samples_exclude.txt", sep="\t", header=False, index=False)

print("\nPlotting embedding...")
expl_var_df["index"] = range(len(expl_var_df))
plot_scree(
    data=expl_var_df,
    x="index",
    y1="expl var %",
    y2="cum expl var %",
    xlabel="Component Number",
    ylabel1="% of variance explained",
    ylabel2="cummulative explained variance",
    lines=[("y1", 0, "#000000"), ("y2", 100, "#000000")],
    title='PCA Scree - {}{}Counts'.format('Centered ' if args.scale else '', 'Scaled ' if args.scale else ''),
    filename="Scree"
)

expl_var = dict(zip(expl_var_df.index, expl_var_df["expl var %"]))
annot = {}
for label in projection_df.index:
    annot[label] = "{}\n{:.2f}%".format(label, expl_var[label])

projection_df = projection_df.iloc[:args.eval_n_pcs, :].T
projection_df["outlier"] = False
projection_df.loc[exclude_samples, "outlier"] = True
z = None
if gte_df is not None:
    projection_df = projection_df.merge(gte_df[[1, 2]].rename(columns={1: "sample", 2: "dataset"}), left_index=True, right_on="sample", how="left").set_index("sample")
    z = "dataset"

plot_embedding(
    data=projection_df,
    z=z,
    label_mask="outlier",
    annot=annot,
    title='PCA - {}{}Counts'.format('Centered ' if args.scale else '', 'Scaled ' if args.scale else ''),
    filename="Pcs"
)

print("Done")
