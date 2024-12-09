#!/usr/bin/env python
# Author: A. Kooijmans, H.J. Westra, and M. Vochteloo

import argparse
import os
import sys
import gzip

parser = argparse.ArgumentParser(description="")
parser.add_argument("--chr", required=True, type=str,  help="")
parser.add_argument("--corr", required=True, type=str,  help="")
parser.add_argument("--top_effects", required=True, type=str,  help="")
parser.add_argument("--chr_col", type=str, default="GeneChr", help="")
parser.add_argument("--significance_col", type=str, default="qval", help="")
parser.add_argument("--gene_col", type=str, default="Gene", help="")
parser.add_argument("--snp_col", type=str, default="SNP", help="")
parser.add_argument("--alpha", type=float, default=0.05,  help="")
parser.add_argument("--out", required=True, type=str,  help="")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

sys.stdout.flush()

os.makedirs(args.out, exist_ok=True) 

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        if mode == "w":
            return gzip.open(file, mode + 't', 4)
        else:
            return gzip.open(file, mode + 't')
    else:
        return open(file, mode)

print("Getting gene pairs IDs from correlation row names")
fin = gzopen(args.corr,"r")
gene_pairs = set() 
counter = 0
fin.readline()
for line in fin:
    counter += 1
    gp = line.strip().split("\t")[0] 
    if "_" not in gp:
        print(f"Error, invalid delimiter found in gene pair '{gp}', gene pair delimeter should be '_'")
        exit()
    gene_pairs.add(gp)    
    if counter % 100000 == 0:
        print(f"{counter} lines processed",end="\r")

print(f"{len(gene_pairs)} gene pairs found")

fin.close()

sys.stdout.flush()

print(f"\nGetting significant top eQTL - SNP combinations using a threshold of: {args.significance_col} < {args.alpha}")
gene_snp_map = {}
fin = gzopen(args.top_effects,"r")
header = fin.readline().strip().split("\t")
for line in fin:
    values = line.strip().split("\t")
    chr = values[header.index(args.chr_col)]
    if chr == args.chr:    
        qval = float(values[header.index(args.significance_col)])
        if qval < args.alpha:
            gene = values[header.index(args.gene_col)]
            snp = values[header.index(args.snp_col)]
            gene_snp_map[gene] = snp
fin.close()
print(f"{len(gene_snp_map)} significant eQTL - SNP combinations found")

sys.stdout.flush()

print("\nGetting SNP gene pair triplets to test ")
fout_genepairs = gzopen(f"{args.out}genepairs.txt","w")
fout_triplest = gzopen(f"{args.out}snp_genepair_triplets_to_test.txt","w")
fout_groups = gzopen(f"{args.out}genepair_egene_groups.txt","w")

eqtl_genes = gene_snp_map.keys()

counter = 0
for gp in sorted(gene_pairs):
    gene1,gene2 = gp.split("_")
    if gene1 in eqtl_genes:
        counter += 1
        snp = gene_snp_map[gene1]
        fout_genepairs.write(f"{gp}\n")
        fout_triplest.write(f"{snp}\t{gp}\n")
        fout_groups.write(f"{gp}\t{gene1}\n")
fout_genepairs.close()
fout_triplest.close()
fout_groups.close()
print(f"{counter} SNP genepair triplets found")

sys.stdout.flush()
print("\nDone.\n")