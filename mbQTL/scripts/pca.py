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

# def calculate_pca(X):
#     pca = PCA()
#     projection = pca.fit_transform(X)
#     return projection, pca.components_, pca.explained_variance_

def calculate_pca(X):
    n_samples, n_features = X.shape

    # Center data
    X -= np.mean(X, axis=0)

    U, S, Vt = linalg.svd(X, full_matrices=False)

    # flip eigenvectors' sign to enforce deterministic output
    # columns of u, rows of v
    max_abs_cols = np.argmax(np.abs(U), axis=0)
    signs = np.sign(U[max_abs_cols, range(U.shape[1])])
    U *= signs
    Vt *= signs[:, np.newaxis]

    # Get variance explained by singular values
    explained_variance_ = (S ** 2) / (n_samples - 1)
    total_var = explained_variance_.sum()
    explained_variance_ratio_ = explained_variance_ / total_var

    U *= S

    return U, Vt, explained_variance_ratio_

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

if args.scale:
    print("Scaling data...")
    m = (m - np.mean(m, axis=0)) / np.std(m, axis=0)

print("Calculating Principal Components...")
x, rotation, expl_var = calculate_pca(X=m)
index = ["PC{}_exp".format(i) for i in range(1, len(samples) + 1)]
projection_df = pd.DataFrame(x, index=samples, columns=index).astype(float).T
rotation_df = pd.DataFrame(rotation, index=index, columns=features)
expl_var_df = pd.Series(expl_var, index=index, name="Explained Variance")

# print(projection_df)
# print(rotation_df)
# print(expl_var_df)

print("Saving results...")
projection_df.to_csv(args.out + "Pcs.txt.gz", sep="\t", header=True, index=True, compression="gzip")
rotation_df.to_csv(args.out + "Pcs_rot.txt.gz", sep="\t", header=True, index=True, compression="gzip")
expl_var_df.to_csv(args.out + "Pcs_var.txt.gz", sep="\t", header=True, index=True, compression="gzip")

print("Done")
