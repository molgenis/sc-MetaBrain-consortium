#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--data", required=True, type=str, help="")
parser.add_argument("--transpose", action="store_true", default=False, help="")
parser.add_argument("--gte", required=False, type=str, default=None, help="")
parser.add_argument("--scale", action="store_true", default=False, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import pandas as pd
import numpy as np
from scipy import linalg


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

print("Loading data...")
df = pd.read_csv(args.data, sep="\t", header=0, index_col=0)

# Transform to numpy for speed.
samples = df.index.to_numpy()
features = df.columns.to_numpy()
m = df.to_numpy()
del df

if args.transpose:
    print("Transposing data...")
    m = np.transpose(m)
    tmp_samples = samples
    samples = features
    features = tmp_samples
    del tmp_samples

if args.gte is not None:
    print("Filtering data...")
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

print("Removing zero variance columns...")
features_mask = np.std(m, axis=0) != 0
m = m[:, features_mask]
features = features[features_mask]

print("Calculating Principal Components...")
pca = prcomp(x=m, scale=args.scale)
pca_indices = ["PC{}_exp".format(i) for i in range(1, len(samples) + 1)]
projection_df = pd.DataFrame(pca["x"], index=samples, columns=pca_indices).astype(float).T
rotation_df = pd.DataFrame(pca["rotation"], columns=pca_indices)
expl_var_df = pd.Series((pca["sdev"] ** 2) / np.sum(pca["sdev"] ** 2), index=pca_indices, name="Explained Variance")

# print(projection_df)
# print(rotation_df)
# print(expl_var_df)

print("Saving results...")
projection_df.to_csv(args.out + "Pcs.txt.gz", sep="\t", header=True, index=True, compression="gzip")
rotation_df.to_csv(args.out + "Pcs_rot.txt.gz", sep="\t", header=True, index=True, compression="gzip")
expl_var_df.to_csv(args.out + "Pcs_var.txt.gz", sep="\t", header=True, index=True, compression="gzip")

print("Done")
