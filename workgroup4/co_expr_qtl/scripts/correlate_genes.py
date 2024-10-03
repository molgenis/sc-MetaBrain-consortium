#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import pandas as pd
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
parser.add_argument("--weights", required=True, type=str, help="")
parser.add_argument("--feature_name", required=False, type=str, default="HGNC", choices=["HGNC", "ENSG", "HGNC_ENSG"], help="")
parser.add_argument("--min_obs", required=False, type=int, default=10, help="")
parser.add_argument("--chunk_size", required=False, type=int, default=10000, help="")
parser.add_argument("--outfile", required=True, type=str, help="")
parser.add_argument("--feature_list", required=True, type=str, help="")
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
weights = pd.read_csv(args.weights, sep="\t", header=0, index_col=0)

if not weights.index.equals(adata.obs_names):
    print("Error, weight cell barcodes do not match matrix cell barcodes")
    exit()

feature_mask = np.sum(adata.X > 0, axis=0) >= args.min_obs
print("\tExcluding {:,} features due to min_obs >={} filter".format(np.size(feature_mask) - np.sum(feature_mask), args.min_obs))
adata = adata[:, feature_mask]
del feature_mask

print("\nFiltering features using provided feature list ...")
gene_list = pd.read_csv(args.feature_list,sep="\t").iloc[:,0].tolist()
# # Sort var_names
adata = adata[:, adata.var_names.isin(gene_list)]
print(f"Remaining features: {len(adata.var_names)}")
del gene_list

# Converting to dense matrix.
m = adata.X.A
if np.isnan(m).any():
    print("Error, count matrix contains NaN values")
    exit()

w = weights.iloc[:,0].to_numpy(dtype=int)
if np.isnan(w).any():
    print("Error, weights array contains NaN values")
    exit()

features = adata.var_names

# # Test example.
# x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
# y = np.array([1, 2, 1, 3, 3, 2, 2, 3, 4, 1, 3, 3, 1, 1, 3, 3, 3, 3, 0, 1])
# w = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0])
# m = np.vstack((x, y)).T
# features = np.array(["GeneA", "GeneB"])
# #   correlation   std.err   t.value   p.value
# # Y  -0.2694926 0.2269819 -1.187287 0.2505473

print("\nWeighing matrix ...")
# Normalise the weights.
norm_w = w / np.mean(w)

# weight the matrix: m - (np.sum(m * weight[:, np.newaxis], axis=0) / np.sum(weight))
weighted_m = m - (np.einsum('ij,i->j', m, norm_w) / np.sum(norm_w))
del m

# calculate the weighted cov matrix: np.sum(w[:, np.newaxis] * weighted_m * weighted_m) / np.sum(w)
sw = np.sum(w)
weighted_m_sq = np.einsum('ij,ij->ij', weighted_m, weighted_m)
cov_a = np.einsum('ij,i->j', weighted_m_sq, w) / sw
del weighted_m_sq

print("\nCalculating correlation ...")
n_features = weighted_m.shape[1]
n_correlations = int(((n_features * n_features) - n_features) / 2)
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

        # Save the upper triangle indices.
        R_chunk[chunk_index] = i
        C_chunk[chunk_index] = j

        # Increment counters.
        chunk_index += 1
        total_index += 1

        # Check if the chunk is full, if not
        # grap the next upper triangle indices.
        if chunk_index < args.chunk_size:
            continue

        # The chunk is full. We can now calculate the pairwise cov
        # between each x, y combination in the upper triangle indices
        # of the current chunk. We can extract the cov's for x * x and y * y
        # from the matrix uses the triangle indices.
        cov_xx = cov_a[R_chunk]
        cov_yy = cov_a[C_chunk]
        cov_xx_yy = np.sqrt(np.einsum('i,i->i', cov_xx, cov_yy))
        weighted_m_sq_pairs = np.einsum('ij,ij->ij', weighted_m[:, R_chunk], weighted_m[:, C_chunk])
        cov_xy = np.einsum('ij,i->j', weighted_m_sq_pairs, w) / sw
        del weighted_m_sq_pairs

        # Now calculate the weighted pearson correlations for the whole chunk at once.
        betas = cov_xy / cov_xx_yy

        # Look up with genes we were processing and write the results to a
        # file line by line.
        featuresi = features[R_chunk]
        featuresj = features[C_chunk]
        for chunk_i in range(chunk_index):
            fh.write(f"{featuresi[chunk_i]}\t{featuresj[chunk_i]}\t{betas[chunk_i]}\n")

        print("\tCalculated {:,} / {:,} correlations".format(total_index, n_correlations), end='\r')
        chunk_index = 0

# Do not forget the final chunk.
if chunk_index > 0:
    # Cut off the end part of the chunk that we did not fill up.
    R_chunk = R_chunk[:chunk_index]
    C_chunk = C_chunk[:chunk_index]

    # Process similar as above.
    cov_xx = cov_a[R_chunk]
    cov_yy = cov_a[C_chunk]
    cov_xx_yy = np.sqrt(np.einsum('i,i->i', cov_xx, cov_yy))
    weighted_m_sq_pairs = np.einsum('ij,ij->ij', weighted_m[:, R_chunk], weighted_m[:, C_chunk])
    cov_xy = np.einsum('ij,i->j', weighted_m_sq_pairs, w) / sw
    del weighted_m_sq_pairs
    betas = cov_xy / cov_xx_yy

    featuresi = features[R_chunk]
    featuresj = features[C_chunk]
    for chunk_i in range(chunk_index):
        fh.write(f"{featuresi[chunk_i]}\t{featuresj[chunk_i]}\t{betas[chunk_i]}\n")
fh.close()
print("\tCalculated {:,} / {:,} correlations".format(total_index, n_correlations))

print("Done")

profiler.disable()
stats = pstats.Stats(profiler).sort_stats('cumtime')
stats.print_stats(40)