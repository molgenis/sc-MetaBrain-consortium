#!/usr/bin/env python
# Author: A. Kooijmans

import time
import sys
import argparse
import gzip
import os
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description="")
parser.add_argument("--annot", required=False, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="")
parser.add_argument("--feature_path", required=True, type=str, help="")
parser.add_argument("--chr", required=True, type=str, help="")
parser.add_argument("--autosomes_only", action="store_true", default=False, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

os.makedirs(os.path.dirname(args.out), exist_ok=True)

print("Loading feature list...")
features = pd.read_csv(args.feature_path,sep="\t")
feature_list = features.iloc[:, 0].sort_values().tolist()

print("Loading chromosome annotation...")
annot = pd.read_csv(f"{args.annot}", sep="\t")
if args.autosomes_only:
    annot = annot[~annot["chromosome"].isin(["X", "Y"])]
annot = annot[annot["feature_id"].isin(feature_list)] # Filter using feature list
all_features = np.array(annot.feature_id) # All genes

annot = annot[annot["chromosome"] == args.chr] # Select chromosome
annot = annot.sort_values(by='feature_id') # Sort in alphabetical order
chr_features = np.array(annot.feature_id) # Extract genes 
print(len(chr_features))

#all_features = set(all_features) - set(chr_features) 
all_features = np.array(list(all_features))
print(len(all_features))
del annot

n_all = len(all_features)
n_chr = len(chr_features)
n_pairs = n_all * n_chr

print(f"Writing {n_pairs} gene pairs for chromosome {args.chr}...")
fileOut = gzip.open(f"{args.out}gene.pairs.chr.{args.chr}.txt.gz","wt")
fileOut.write("\t")
idx = 0
for i in range(n_chr):
    gene_i = chr_features[i]
    for j in range(n_all):
        gene_j = all_features[j]
        if gene_i == gene_j:
            continue
        gene_pair = f"{gene_i}_{gene_j}\t"
        fileOut.write(gene_pair)
        idx += 1
        if idx % 1000000 == 1:
            print(f"Lines written: {idx}",end="\r")

fileOut.write("\n")
fileOut.close()
print("Done")