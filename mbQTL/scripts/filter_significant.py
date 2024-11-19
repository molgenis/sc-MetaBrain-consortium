#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import gzip
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--data", required=True, type=str, help="")
parser.add_argument("--signif_col", required=False, type=str, default="MetaP", help="")
parser.add_argument("--alpha", required=False, type=float, default=0.05, help="")
parser.add_argument("--outfile", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

dir = os.path.dirname(args.outfile)
if dir != "":
    os.makedirs(dir, exist_ok=True)

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)

print("Filtering results file...")
data = []
header = None
with gzopen(args.data, 'r') as f:
    for i, line in enumerate(f):
        elems = line.rstrip("\n").split("\t")
        if i == 0:
            header = line
            pos = {column: index for index, column in enumerate(elems)}
            if args.signif_col not in pos:
                print("Error, could not find column '{}' in results file header.".format(args.signif_col))
                exit()
            continue

        try:
            value = float(elems[pos[args.signif_col]])
        except ValueError:
            continue

        if value < args.alpha:
            data.append((value, line))
f.close()
data.sort(key=lambda x: x[0])
print("  Loaded {:,} results of which {:,} were kept".format(i, len(data)))
if header is None:
    print("  Error, no data in input file.")
    exit()

print("Writing significant file...")
fhout = gzopen(args.outfile, 'w')
fhout.write(header)
for _, line in data:
    fhout.write(line)
fhout.close()

print("\nDone")
