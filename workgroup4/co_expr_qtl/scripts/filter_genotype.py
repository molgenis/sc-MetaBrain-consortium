#!/usr/bin/env python
# Author: A. Kooijmans, H.J. Westra, and M. Vochteloo

import gzip
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("--snp_genepair", required=True, type=str,  help="eQTL top effects file")
# parser.add_argument("--linkfile", required=True, type=str,  help="eQTL top effects file")
parser.add_argument("--genotype", required=True, type=str,  help="genotype file")
parser.add_argument("--output", required=True, type=str,  help="output file")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        if mode == "w":
            return gzip.open(file, mode + 't', 4)
        else:
            return gzip.open(file, mode + 't')
    else:
        return open(file, mode)

print("Loading SNPs")
fin = gzopen(args.snp_genepair,"r")
top_snps = set()
for line in fin:
    snp = line.strip().split("\t")[0]
    top_snps.add(snp)
fin.close()
print(f"{len(top_snps)} SNPs loaded")

# print("\nLoading samples")
# fin = gzopen(args.linkfile,"r")
# samples = set()
# for line in fin:
#     id = line.strip().split("\t")[1]
#     samples.add(id)
# fin.close()
# print(f"{len(samples)} samples loaded")

# ADD: Filter on sample ids

# Load genotype file
fin = gzopen(args.genotype,"r")
fout = gzopen(args.output,"w")

print("\nFiltering genotype file")
counter = 0
found_snps = set()
for line in fin:
    counter += 1 
    if line.startswith("#"):
        fout.write(line)
    else:
        snp = line.split("\t")[2]
        if snp in top_snps:
            found_snps.add(snp)
            fout.write(line)
            print(f"{len(found_snps)} / {len(top_snps)} lines written",end="\r")
            if len(found_snps) == len(top_snps):
                break
fin.close()
fout.close()
print(f"{len(found_snps)} / {len(top_snps)} lines written")

if len(found_snps) < len(top_snps):
    print(f"Warning, {len(top_snps - found_snps)} not found in genotype file")

print("\nDone.")
