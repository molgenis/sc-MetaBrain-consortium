#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import pandas as pd
import numpy as np
import scanpy
import gzip
import os
import struct

parser = argparse.ArgumentParser(description="")
parser.add_argument("--counts", required=True, type=str, help="File with counts")
parser.add_argument("--weights", required=True, type=str, help="File with weights")
parser.add_argument("--feature_name", required=False, type=str, default="HGNC", choices=["HGNC", "ENSG", "HGNC_ENSG"], help="What type of feature name to use in the output.")
parser.add_argument("--geneannotation", required=False, type=str, help="Gene chromosome annotation")
parser.add_argument("--chr", required=False, type=str, help="Limit correlation calculation to gene pairs where the first gene is on specified chromosome")
parser.add_argument("--egenelist", required=False, type=str, help="List of all valid eQTL genes")
parser.add_argument("--coegenelist", required=False, type=str, help="List of all valid co-eQTL genes")
parser.add_argument("--min_obs", required=False, type=int, default=10, help="Minimum non-zero observations per gene")
parser.add_argument("--chunk_size", required=False, type=int, default=10000, help="Chunk size")
parser.add_argument("--binary", action="store_true", default=False, help="Output the results in binary format")
parser.add_argument("--out", required=True, type=str, help="Output filename")
args = parser.parse_args()

if os.path.dirname(args.out) != "":
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

if args.geneannotation is not None and args.chr is None:
    print("Warning, --geneannotation is ignored if --chr is not given.")
if args.chr is not None and args.geneannotation is None:
    print("Error, --geneannotation must be given in if --chr is given.")
    exit()

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        if mode == "w":
            return gzip.open(file, mode + 't', 4)
        else:
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
        print("Error, 'count_data.var_names' are expected to have the same order as 'count_data.var['gene_ids']'.")
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

def load_annotation(annotation, chr=None):
    print("  Loading gene annotation...")
    annot_genes = set()
    pos = {}
    with gzopen(annotation, mode='r') as f:
        for i, line in enumerate(f):
            values = line.rstrip().split("\t")
            if i == 0:
                pos = {value: index for index, value in enumerate(values)}
                continue
            if chr is not None and values[pos["Chr"]] != chr:
                continue

            if args.feature_name == "HGNC":
                annot_genes.add(values[pos["GeneSymbol"]])
            elif args.feature_name == "ENSG":
                annot_genes.add(values[pos["Gene"]])
            elif args.feature_name == "HGNC_ENSG":
                annot_genes.add(values[pos["Gene"] + "_" + values[pos["GeneSymbol"]]])
            else:
                print("Unexpected feature name '{}'".format(args.feature_name))
                exit()
    f.close()
    print("\tLoaded {:,} genes{}".format(len(annot_genes), " for chr '{}'".format(chr) if chr is not None else ""))

    return annot_genes

def read_as_set(fpath):
    fh = gzopen(fpath,'r')
    data = set()
    for line in fh:
        data.add(line.strip())
    return data


#############################################

print("\nLoading counts ...")
adata = load_counts(counts_fpath=args.counts)
weights = pd.read_csv(args.weights, sep="\t", header=0, index_col=0)
geneset = set(adata.var_names)

print("\nFiltering genes to correlate ...")
print("\tInput data has {:,} genes".format(len(geneset)))

# read a set of genes that we want to limit to (if any)
if args.egenelist is None:
    egeneset = set(geneset)
else:
    egeneset = read_as_set(fpath=args.egenelist)
    print("\tInput eGene set has {:,} genes".format(len(egeneset)))

if args.coegenelist is None:
    coegeneset = set(geneset)
else:
    coegeneset = read_as_set(fpath=args.coegenelist)
    print("\tInput co-eGene set has {:,} genes".format(len(coegeneset)))

if args.geneannotation is not None and args.chr is not None:
    # limit the eGene set to only genes on the chromosome of interest.
    # Load the gene to chromosome annotation
    chr_genes = load_annotation(annotation=args.geneannotation, chr=args.chr)
    egeneset = egeneset.intersection(chr_genes)
    print("\tLimit eGenes to {:,} genes on chromosome '{}'".format(len(egeneset), args.chr))
    if len(egeneset) == 0:
        print("Specified to run on chromosome {}, but no genes in the data match to this string.".format(args.chr))
        exit()
    del chr_genes

    args.out = f"{args.out}.chr{args.chr}"

# Filter the input data.
feature_mask = np.array([gene in egeneset or gene in coegeneset for gene in adata.var_names])
adata = adata[:, feature_mask]
print("\tExcluding {:,} features due to eGene / co-eGene / chromosome filter".format(np.size(feature_mask) - np.sum(feature_mask)))
if np.sum(feature_mask) == 0:
    print("Error, no data")
    exit()
del feature_mask

# TODO: check for NaN or INF in input count matrix?

print("\nFiltering genes where counts sum to zero...")
feature_mask = np.sum(adata.X > 0, axis=0) >= args.min_obs
print("\tExcluding {:,} additional features due to min_obs >={} filter".format(np.size(feature_mask) - np.sum(feature_mask), args.min_obs))
adata = adata[:, feature_mask]
geneset = set(adata.var_names)
egeneset = egeneset.intersection(geneset)
coegeneset = coegeneset.intersection(geneset)
if np.sum(feature_mask) == 0:
    print("Error, no data")
    exit()
del feature_mask

print("\nPreprocessing...")
# Converting to dense matrix.
m = adata.X.A
w = weights.loc[adata.obs_names, ].iloc[:, 0].to_numpy()
features = adata.var_names
del adata

# # Test example.
# x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
# y = np.array([1, 2, 1, 3, 3, 2, 2, 3, 4, 1, 3, 3, 1, 1, 3, 3, 3, 3, 0, 1])
# w = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0])
# m = np.vstack((x, y)).T
# features = np.array(["GeneA", "GeneB"])
# #   correlation   std.err   t.value   p.value
# # Y  -0.2694926 0.2269819 -1.187287 0.2505473

# Determine how many values we have and will calculate.
n_features = m.shape[1]
n_egenes = len(egeneset)
n_coegenes = len(coegeneset)
n_correlations = int((n_egenes * n_coegenes) - n_egenes - (((n_egenes * n_egenes) - n_egenes) / 2))
print("\tCalculating {:,} correlations for {:,} features".format(n_correlations, n_features))

for feature in features:
    if "," in feature:
        print("Error, features contain comma.")
        exit()

features_str = ",".join(features)
n_feature_chars = len(features_str)
is_egene = [feature in egeneset for feature in features]

# Write the header to the output file.
if args.binary:
    if (n_features * n_features) > 2147483647:
        print("Error, number of values exceeds the max integer within 4 bytes.")
        exit()
    if n_feature_chars > 2147483647:
        print("Error, the feature string exceeds the max integer within 4 bytes.")
        exit()

    fh = gzopen(args.out + ".corr.dat", "wb")
    fh.write(struct.pack('>4i', n_egenes, n_coegenes, n_correlations, n_feature_chars))
    fh.write(struct.pack(f'>{n_feature_chars}s', features_str.encode('ascii')))
    fh.write(struct.pack(f'>{n_features}?', *is_egene))
else:
    fh = gzopen(args.out + ".corr.txt.gz", "w")
    fh.write(f"{n_egenes}\n{n_coegenes}\n{n_correlations}\n{n_feature_chars}\n")
    fh.write(f"{features_str}\n")
    fh.write(",".join([str(value) for value in is_egene]) + "\n")

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

# Mimicks: R, C = np.triu_indices(n_features,1) to get the upper triangle but
# creating that mask for 178M correlations takes too long and too much memory.
# Instead, we create it ourselves in chunks and calculate the correlations for
# one chunk all at once. This way we can balance memory and speed by adjusting
# the chunk size.
R_chunk = np.empty(args.chunk_size, dtype=int) # row indices
C_chunk = np.empty(args.chunk_size, dtype=int) # column indices
V_chunk = np.empty(args.chunk_size, dtype=int) # value indices
if args.binary:
    O_chunk = bytearray(12 * args.chunk_size) # output array; initialize byte buffer once, 12 bytes = 1 int and 1 double
else:
    O_chunk = ["NaN"] * args.chunk_size # output array
chunk_index = 0
total_index = 0
for i, feature1 in enumerate(features):
    if feature1 not in egeneset:
        continue
    for j, feature2 in enumerate(features):
        if feature2 not in coegeneset:
            continue
        if i == j:
            continue
        if i > j and feature2 in egeneset:
            continue

        # Save the upper triangle indices.
        R_chunk[chunk_index] = i
        C_chunk[chunk_index] = j
        V_chunk[chunk_index] = (i * n_features) + j

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
            if args.binary:
                byteidx = chunk_i * 12
                byteidxend = byteidx + 12
                O_chunk[byteidx:byteidxend] = struct.pack('>id', V_chunk[chunk_i], betas[chunk_i])
            else:
                O_chunk[chunk_i] = f"{V_chunk[chunk_i]}\t{betas[chunk_i]}"

        if args.binary:
            fh.write(O_chunk)
        else:
            fh.write("\n".join(O_chunk) + "\n")

        print("\tCalculated {:,} / {:,} correlations".format(total_index, n_correlations), end='\r')
        chunk_index = 0

# Do not forget the final chunk.
if chunk_index > 0:
    # Cut off the end part of the chunk that we did not fill up.
    R_chunk = R_chunk[:chunk_index]
    C_chunk = C_chunk[:chunk_index]
    V_chunk = V_chunk[:chunk_index]
    if args.binary:
        O_chunk = O_chunk[:(12 * chunk_index)]
    else:
        O_chunk = O_chunk[:chunk_index]

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
        if args.binary:
            byteidx = chunk_i * 12
            byteidxend = byteidx + 12
            O_chunk[byteidx:byteidxend] = struct.pack('>id', V_chunk[chunk_i], betas[chunk_i])
        else:
            O_chunk[chunk_i] = f"{V_chunk[chunk_i]}\t{betas[chunk_i]}"

    if args.binary:
        fh.write(O_chunk)
    else:
        fh.write("\n".join(O_chunk) + "\n")
fh.close()
print("\tCalculated {:,} / {:,} correlations".format(total_index, n_correlations))

print("Done")