#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import os
import gzip

parser = argparse.ArgumentParser(description="")
parser.add_argument("--gte", required=True, type=str, help="")
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

print("Loading GTE...")
lines = [[""], ["Dataset"]]
with gzopen(args.gte, mode='r') as fhin:
    for i, line in enumerate(fhin):
        values = line.rstrip("\n").split("\t")
        lines[0].append(values[1])
        lines[1].append(values[2])
fhin.close()


print("\nSaving dataset as covariates...")
fhout = gzopen(args.outfile, mode='w')
for line in lines:
    fhout.write("\t".join(line) + "\n")
fhout.close()

print("Done")