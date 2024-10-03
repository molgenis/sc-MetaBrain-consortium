#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import pandas as pd
import numpy as np
import gzip
import os
import cProfile
import pstats

profiler = cProfile.Profile()
profiler.enable()

parser = argparse.ArgumentParser(description="")
parser.add_argument("--in_gtf", required=True, type=str, help="")
parser.add_argument("--gene_pairs", required=True, type=str, help="")
parser.add_argument("--chr", required=False, type=str, default="gene_name", help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

os.makedirs(os.path.dirname(args.out), exist_ok=True)

# Read in correlation gene pairs
print(f"Reading feature pairs corresponding to chromosome {args.chr} ...")
with gzip.open(args.gene_pairs,"rt") as f:
    for line in f:
        feature_pairs = line.strip().split("\t")

print(f"Processing {len(feature_pairs)} feature pairs ...")
feature_dict = {}
for x in feature_pairs:
    feature = x.split("_")[0]
    feature_dict[feature] = np.empty(4,dtype="U50")

# Extract info from gtf file 
print("Reading gtf annotation file ...")
id_dict = {}
with open(args.in_gtf, mode="r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        seqname, source, feature, start, end, score, strand, frame, attribute = line.rstrip("\n").split("\t")
        if feature != "gene":
            continue
        attributes = {key: val.replace('"', '') for key, val in (attr.split(' ') for attr in attribute.split("; "))}
        if attributes["gene_name"] in feature_dict.keys():
            feature_dict[attributes["gene_name"]][0] = source
            feature_dict[attributes["gene_name"]][1] = start
            feature_dict[attributes["gene_name"]][2] = end
            feature_dict[attributes["gene_name"]][3] = strand
        id_dict[ attributes["gene_name"]] = attributes["gene_id"]

# Write annotation for each gene pair
print(f"Writing annotation file ...")
fileOut = gzip.open(f"{args.out}chr.{args.chr}.annotation.txt.gz","wt")
fileOut.write("Platform\tArrayAddress\tSymbol\tChr\tChrStart\tChrEnd\tProbe\tStrand\n")
for pair in feature_pairs:
    array_adress = pair
    ft0,ft1 = x.split("_")
    platform = feature_dict[ft0][0]
    chr_start = feature_dict[ft0][1]
    chr_end = feature_dict[ft0][2]
    strand = feature_dict[ft0][3]
    symbol = f"{id_dict[ft0]}_{id_dict[ft1]}" # Get ensg symbol
    fileOut.write(f"{platform}\t{array_adress}\t{symbol}\t{args.chr}\t{chr_start}\t{chr_end}\t{symbol}\t{strand}\n")
fileOut.close()

print("Done.")

profiler.disable()
stats = pstats.Stats(profiler).sort_stats('cumtime')
stats.print_stats(40)