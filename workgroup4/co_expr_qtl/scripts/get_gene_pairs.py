#!/usr/bin/env python

import gzip
import argparse
import sys
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--input", required=True, type=str,  help="")
parser.add_argument("--out", required=True, type=str,  help="")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

os.makedirs(args.out, exist_ok=True) 

fin = gzip.open(args.input,"rt")
fout = open(f"{args.out}genepairs.txt","wt")

print("Extracting gene pairs from correlation file ...")
counter = 0
for line in fin:
    gene_pair = line.strip().split("\t")[0]
    fout.write(f"{gene_pair}\n")
    counter += 1
    if counter % 100000 == 0:
        print(f"{counter} lines processed",end="\r")
fin.close()
fout.close()
print("Done.")
