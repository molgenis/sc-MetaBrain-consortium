#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import glob
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--input", required=True, type=str, help="")
parser.add_argument("--cell_type", required=True, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

if not args.out.endswith(os.sep):
    args.out += os.sep
os.makedirs(os.path.dirname(args.out), exist_ok=True)


print("Checking for sample with multiple files ...")
sample_fpaths = {}
for fpath in glob.glob(args.input):
    # expression/{pool}/data/{pool}.{sample}.{cell_type}.raw.counts.h5"
    pool, sample, cell_type = os.path.basename(fpath).split(".")[:3]
    if cell_type != args.cell_type:
        continue

    if sample in sample_fpaths:
        sample_fpaths[sample].append(fpath)
    else:
        sample_fpaths[sample] = [fpath]

print("Calculating average read counts ...")
for sample, fpaths in sample_fpaths.items():
    counts_outpath = os.path.join(args.out, "{sample}.raw.counts.h5".format(sample=sample))
    weights_outpath = os.path.join(args.out, "{sample}.raw.weights.txt.gz".format(sample=sample))
    if len(fpaths) == 0:
        # Copy file.
        fh = open(counts_outpath, "w")
        fh.close()
        fh = open(weights_outpath, "w")
        fh.close()
    else:
        print("\t Combining: {}".format(", ".join([os.path.basename(fpath) for fpath in fpaths])))
        # Combine matrices.
        fh = open(counts_outpath, "w")
        fh.close()
        fh = open(weights_outpath, "w")
        fh.close()

print("Done.")
