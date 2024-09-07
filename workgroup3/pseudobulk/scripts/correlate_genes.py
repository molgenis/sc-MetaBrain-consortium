#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import numpy as np
import scanpy
import gzip
import os

import cProfile
import pstats
profiler = cProfile.Profile()
profiler.enable()

parser = argparse.ArgumentParser(description="")
parser.add_argument("--counts", required=True, type=str, help="")
parser.add_argument("--feature_name", required=False, type=str, default="HGNC", choices=["HGNC", "ENSG", "HGNC_ENSG"], help="")
parser.add_argument("--min_obs", required=False, type=int, default=10, help="")
parser.add_argument("--outfile", required=True, type=str, help="")
args = parser.parse_args()

dir = os.path.dirname(args.outfile)
if dir != "":
    os.makedirs(dir, exist_ok=True)

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")


def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)

def load_counts(counts_fpath):
    adata = scanpy.read_10x_h5(counts_fpath)
    count_matrix = adata.X
    n_barcodes, n_genes = count_matrix.shape
    genes = adata.var_names.to_numpy()
    barcodes = adata.obs_names.to_numpy()
    if np.size(genes) != n_genes:
        print("Error, matrix size does not match gene annotation.")
        exit()
    if np.size(barcodes) != n_barcodes:
        print("Error, matrix size does not match barcodes annotation.")
        exit()
    print("\tLoaded raw matrix with {:,} barcodes and {:,} genes.".format(n_barcodes, n_genes))

    # Check for unique barcodes.
    if len(np.unique(barcodes)) != np.size(barcodes):
        print("Error, not all barcodes are unique.")
        exit()

    # Parse the gene info
    gene_ids = adata.var['gene_ids']
    gene_ids = gene_ids.reset_index(drop=False)
    gene_ids = gene_ids.to_numpy()
    if not (genes == gene_ids[:, 0]).all():
        print("Error, 'count_data.var_names' are expected to have the same order as 'count_data.var['gene_ids''.")
        exit()

    # Set gene names.
    features = get_features(m=gene_ids, indices={"HGNC": 0, "ENSG": 1})

    # Remove duplicates.
    u, c = np.unique(features, return_counts=True)
    dup = u[c > 1]
    features_mask = np.in1d(features, dup)
    if np.sum(features_mask) > 0:
        print("\tRemoving {:,} duplicate genes.".format(np.sum(features_mask)))
        adata = adata[:, ~features_mask]
        adata.var_names = features[~features_mask]

    # Parse the barcode info.
    # barcode_ids = count_data.var['feature_types'].to_numpy()

    print("\tDuplicate filtered matrix has shape: {}".format(count_matrix.shape))
    return adata

def get_features(m, indices):
    if args.feature_name == "HGNC":
        return m[:, indices["HGNC"]]
    elif args.feature_name == "ENSG":
        return m[:, indices["ENSG"]]
    elif args.feature_name == "HGNC_ENSG":
        return np.char.add(np.char.add(m[:, indices["HGNC"]], "_"), m[:, indices["ENSG"]])
    else:
        print("Unexpected feature name '{}'".format(args.feature_name))
        exit()

def weight_matrix(m, weight):
    """
    https://stackoverflow.com/questions/38641691/weighted-correlation-coefficient-with-pandas
    """
    # wtd_cor
    weight = weight / np.mean(weight)

    # stdz
    sw = np.sum(weight)
    wtd_mean = np.sum(m * weight[:, np.newaxis], axis=0) / sw
    m = m - wtd_mean
    return m

def pre_weighted_cov(x, y, w, sw):
    """
    https://stackoverflow.com/questions/38641691/weighted-correlation-coefficient-with-pandas
    """
    return np.sum(w * x * y) / sw

def pre_weighted_corr(cov_x, cov_y, cov_xy):
    """
    https://stackoverflow.com/questions/38641691/weighted-correlation-coefficient-with-pandas
    """
    return cov_xy / np.sqrt(cov_x * cov_y)

#############################################

print("\nLoading counts ...")
adata = load_counts(counts_fpath=args.counts)

# TODO: check for NaN or INF

print("\nFiltering genes ...")
feature_mask = np.sum(adata.X > 0, axis=0) >= args.min_obs
print("\tExcluding {:,} features due to min_obs >={} filter".format(np.size(feature_mask) - np.sum(feature_mask), args.min_obs))
adata = adata[:, feature_mask]
del feature_mask

# Converting to dense matrix.
m = adata.X.A
features = adata.var_names
del adata

print("\nWeighing matrix ...")
w = m.sum(axis=1)
sw = np.sum(w)
weighted_m = weight_matrix(m=m, weight=w)

print("\nCalculating correlation ...")
n_features = weighted_m.shape[1]
n_correlations = int(((n_features * n_features) / 2) - n_features)
print("\tCalculating {:,} correlations for {:,} features".format(n_correlations, n_features))

fh = gzopen(args.outfile, "wt")
fh.write("Feature1\tFeature2\tCorrelation\n")

index = 0
for i in range(n_features):
    featurei = features[i]

    x = weighted_m[:, i]
    cov_x = pre_weighted_cov(x, x, w, sw)
    for j in range(n_features):
        if j >= i:
            continue
        if index == 0 or index % 1e5 == 0:
            print("\tCalculated {:,} / {:,} correlations".format(index, n_correlations), end='\r')

        featurej = features[j]
        y = weighted_m[:, j]
        cov_y = pre_weighted_cov(y, y, w, sw)
        cov_xy = pre_weighted_cov(x, y, w, sw)
        beta = pre_weighted_corr(cov_x=cov_x, cov_y=cov_y, cov_xy=cov_xy)

        fh.write(f"{featurei}\t{featurej}\t{beta}\n")

        index += 1
fh.close()
print("\tCalculated {:,} / {:,} correlations".format(index, n_correlations))

print("Done")

profiler.disable()
stats = pstats.Stats(profiler).sort_stats('cumtime')
stats.print_stats(40)