#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import os
import gzip
import re

parser = argparse.ArgumentParser(description="")
parser.add_argument("--annotation", required=True, type=str, help="")
parser.add_argument("--genomic_range", required=False, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

chromosome = None
start = None
end = None
if args.genomic_range is not None:
    match = re.match("([0-9]{1,2}|X|Y|MT):([0-9]+)-([0-9]+)", args.genomic_range)
    chromosome = str(match.group(1))
    start = int(match.group(2))
    end = int(match.group(3))

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)

fhin = gzopen(args.annotation, mode='r')
fhout = gzopen(os.path.join(args.out + "-genes.txt"), mode='w')

print("Filtering genes...")
pos = {}
genes = set()
for i, line in enumerate(fhin):
    values = line.rstrip("\n").split("\t")
    if i == 0:
        pos = {label: index for index, label in enumerate(values)}
        continue

    # Check if the gene is in range.
    if chromosome is not None and values[pos["Chr"]] != chromosome:
        continue
    mean = (float(values[pos["ChrStart"]]) + float(values[pos["ChrEnd"]])) / 2
    if start is not None and mean < start:
        continue
    if end is not None and mean > end:
        continue

    # Save the gene.
    gene = values[pos["ArrayAddress"]]
    if gene in genes:
        print("Warning, duplicate genes.")
    fhout.write(gene + "\n")
    genes.add(gene)

fhout.close()

print("  Saved {:,} genes".format(len(genes)))
print("Done")