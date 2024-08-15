#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import os
import gzip

parser = argparse.ArgumentParser(description="")
parser.add_argument("--gte", required=True, type=str, help="")
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

print("Loading GTE...")
datasets = {}
with gzopen(args.gte, mode='r') as fhin:
    for i, line in enumerate(fhin):
        values = line.rstrip("\n").split("\t")
        if values[2] not in datasets:
            datasets[values[2]] = []
        datasets[values[2]].append(line)
fhin.close()


print("\nSplitting on dataset column...")
for dataset, lines in datasets.items():
    print("\tDataset {} has {:,} samples".format(dataset, len(lines)))
    fhout = gzopen(os.path.join(args.out + "gte_" + dataset + ".txt"), mode='w')
    for line in lines:
        fhout.write(line)
    fhout.close()

print("Done")