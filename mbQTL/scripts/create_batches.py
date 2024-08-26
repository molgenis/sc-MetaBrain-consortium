#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import os
import gzip

parser = argparse.ArgumentParser(description="")
parser.add_argument("--genes", required=True, type=str, help="")
parser.add_argument("--annotation", required=False, type=str, default=None, help="")
parser.add_argument("--header_index", required=False, type=int, default=None, help="")
parser.add_argument("--gene_index", required=False, type=int, default=0, help="")
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

annot_genes = None
if args.annotation is not None:
    print("Loading annotation genes...")
    annot_genes = set()
    pos = {}
    with gzopen(args.annotation, mode='r') as f:
        for i, line in enumerate(f):
            values = line.rstrip().split("\t")
            if i == 0:
                pos = {value: index for index, value in enumerate(values)}
                continue
            annot_genes.add(values[pos["ArrayAddress"]])
    f.close()
    print("  Loaded {:,} genes".format(len(annot_genes)))

fhin = gzopen(args.genes, mode='r')

print("Splitting genes in chunks ...")
genes = set()
n_skipped = 0
n_duplicated = 0
chunk_id = 0
chunk = []
for i, line in enumerate(fhin):
    values = line.rstrip("\n").split("\t")
    if args.header_index is not None and i == args.header_index:
        continue

    # Save the gene.
    gene = values[args.gene_index]
    if annot_genes is not None and gene not in annot_genes:
        print("Warning, gene '{}' is not in annotation. Skipping gene.".format(gene))
        n_skipped += 1
        continue
    if gene in genes:
        print("Warning, gene '{}' is duplicated. Skipping gene.".format(gene))
        n_duplicated += 1
        continue

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
if args.annotation is not None:
    print("  {:,} genes are skipped because of --annotation".format(n_skipped))
print("  {:,} genes are duplicated".format(n_duplicated))

if chunk_id == 0:
    print("Error, no chunks created.")
    exit(1)

print("Done")