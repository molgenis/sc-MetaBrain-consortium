#!/usr/bin/env python
# Author: A. Kooijmans

import argparse
import gzip
import os
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description="")
parser.add_argument("--input", required=True,  nargs="*", help="Individual correlation files per sample")
parser.add_argument("--gene_pairs", required=True, type=str, help="List of gene pairs for a chromosome")
parser.add_argument("--out", required=True, type=str, help="Output file")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

def get_sample_id(filepath):
    file = filepath.split('/')[-1]
    pool_sample = file.split(".")[0:2]
    sample_id = ".".join(pool_sample)
    return sample_id                

print(f"Reading in gene pairs...")
with gzip.open(args.gene_pairs,"rt") as f:
    for line in f:
        gene_pairs = line.strip().split("\t")
        break
f.close()

n_gene_pairs = len(gene_pairs)

print("Creating gene indexing dictionary...")
gene_pair_idx = {gene_pair: i for i, gene_pair in enumerate(gene_pairs)}

print(f"Writing aggregated file...")
fileOut = gzip.open(args.out,"at")
fileOut.write("\t")
fileOut.write("\t".join(gene_pairs))
fileOut.write("\n")

header = None
counter = 0
for file in args.input:
    corr_array = np.full(n_gene_pairs, np.nan)
    sample_id = get_sample_id(file)
    print(f"Processing {sample_id}...")
    with gzip.open(file,"rt") as f:
        for line in f:
            counter += 1
            if header == None:
                header = line
                continue
            values = line.strip().split("\t")
            gene1 = values[0]
            gene2 = values[1]
            corr = values[2]

            if f"{gene1}_{gene2}" in gene_pair_idx:
                row_idx = gene_pair_idx[f"{gene1}_{gene2}"]
                corr_array[row_idx] = corr
                continue
            elif f"{gene2}_{gene1}" in gene_pair_idx:
                row_idx = gene_pair_idx[f"{gene2}_{gene1}"]
                corr_array[row_idx] = corr
                continue

    fileOut.write(f"{sample_id}\t")
    fileOut.write("\t".join(map(str, corr_array))) 
    fileOut.write("\n")
fileOut.close()
print("Done.")