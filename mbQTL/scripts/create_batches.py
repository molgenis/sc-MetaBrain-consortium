#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import os
import gzip

parser = argparse.ArgumentParser(description="")
parser.add_argument("--exp", required=True, type=str, help="")
parser.add_argument("--n_genes", required=False, type=int, default=200, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)

fhin = gzopen(args.exp, mode='r')

print("Filtering genes...")
pos = {}
genes = set()
chunk_id = 0
chunk = []
for i, line in enumerate(fhin):
    values = line.rstrip("\n").split("\t")
    if i == 0:
        pos = {label: index for index, label in enumerate(values)}
        continue

    # Save the gene.
    gene = values[0]
    if gene in genes:
        print("Warning, duplicate genes.")
    chunk.append(gene)
    genes.add(gene)

    # Check if the gene is in range.
    if len(chunk) == args.n_genes:
        fhout = gzopen(os.path.join(args.out + "-chunk" + str(chunk_id) + "-genes.txt"), mode='w')
        for gene in chunk:
            fhout.write(gene + "\n")
        fhout.close()
        chunk = []
        chunk_id += 1

fhin.close()

if len(chunk) > 0:
    fhout = gzopen(os.path.join(args.out + "-chunk" + str(chunk_id) + "-genes.txt"), mode='w')
    for gene in chunk:
        fhout.write(gene + "\n")
    fhout.close()
    chunk = []
    chunk_id += 1

print("  Saved {:,} genes into {:,} chunks".format(len(genes), chunk_id))
print("Done")