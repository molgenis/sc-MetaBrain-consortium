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
parser.add_argument("--chunk_size", required=False, type=int, default=10000, help="")
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

# Normalise the weights.
norm_w = w / np.mean(w)

# weight the matrix: m - (np.sum(m * weight[:, np.newaxis], axis=0) / np.sum(weight))
weighted_m = m - (np.einsum('ij,i->j', m, norm_w) / np.sum(norm_w))

# calculate the weighted cov matrix: np.sum(w[:, np.newaxis] * weighted_m * weighted_m) / np.sum(w)
cov_m = np.einsum('ij,i->j', (weighted_m * weighted_m), w) / np.sum(w)

print("\nCalculating correlation ...")
n_features = weighted_m.shape[1]
n_correlations = int(((n_features * n_features) / 2) - n_features)
print("\tCalculating {:,} correlations for {:,} features".format(n_correlations, n_features))

fh = gzopen(args.outfile, "wt")
fh.write("Feature1\tFeature2\tCorrelation\n")

# Mimicks: R, C = np.triu_indices(n_features,1) to get the upper triangle but
# creating that mask for 178M correlations takes too long and too much memory.
# Instead, we create it ourselves in chunks and calculate the correlations for
# one chunk all at once. This way we can balance memory and speed by adjusting
# the chunk size.
R_chunk = np.empty(args.chunk_size, dtype=int)
C_chunk = np.empty(args.chunk_size, dtype=int)
chunk_index = 0
total_index = 0
for i in range(n_features):
    for j in range(n_features):
        if j >= i:
            continue

        R_chunk[chunk_index] = i
        C_chunk[chunk_index] = j

        chunk_index += 1
        total_index += 1
        if chunk_index < args.chunk_size:
            continue

        # The chunk is full. We can now calculate the pairwise cov
        # between each x, y combination in the upper triangle indices
        # of the current chunk. We can extract the cov's for x * x and y * y
        # from the matrix uses the triangle indices.
        cov_x = cov_m[R_chunk]
        cov_y = cov_m[C_chunk]
        cov_xy = np.einsum('ij,i->j', (weighted_m[:, R_chunk] * weighted_m[:, C_chunk]), w) / np.sum(w)

        # Now calculate the weighted pearson correlations for the whole chunk at once.
        betas = cov_xy / np.sqrt(cov_x * cov_y)

        # Look up with genes we were processing and write the results to file.
        # Is there a faster way to do this? Should I open and close the file
        # for each chunk?
        featuresi = features[R_chunk]
        featuresj = features[C_chunk]
        for featurei, featurej, beta in zip(featuresi, featuresj, betas):
            fh.write(f"{featurei}\t{featurej}\t{beta}\n")

        print("\tCalculated {:,} / {:,} correlations".format(total_index, n_correlations), end='\r')
        chunk_index = 0

# Do not forget the final chunk.
if chunk_index > 0:
    # Cut off the end part of the chunk that we did not fill up.
    R_chunk = R_chunk[:chunk_index]
    C_chunk = C_chunk[:chunk_index]

    # Process similar as above.
    cov_x = cov_m[R_chunk]
    cov_y = cov_m[C_chunk]
    cov_xy = np.einsum('ij,i->j', (weighted_m[:, R_chunk] * weighted_m[:, C_chunk]), w) / np.sum(w)
    betas = cov_xy / np.sqrt(cov_x * cov_y)

    featuresi = features[R_chunk]
    featuresj = features[C_chunk]
    for featurei, featurej, beta in zip(featuresi, featuresj, betas):
        fh.write(f"{featurei}\t{featurej}\t{beta}\n")
fh.close()
print("\tCalculated {:,} / {:,} correlations".format(total_index, n_correlations))

print("Done")

profiler.disable()
stats = pstats.Stats(profiler).sort_stats('cumtime')
stats.print_stats(40)