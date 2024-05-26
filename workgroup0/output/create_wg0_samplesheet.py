#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse


"""
Syntax:
./create_wg0_samplesheet.py -h 
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--work_dir", type=str, required=True, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import pandas as pd
import glob
import os
import re

data = []

print("Creating WGO samplesheet file")
for fpath in glob.glob(os.path.join(args.work_dir, "fastq", "*")):
    match = re.match("(.+)_S[0-9]_L00[0-4]_R[12]_001\.fastq\.gz", os.path.basename(fpath))
    if match is None:
        print("\tSkipped '{}'".format(os.path.basename(fpath)))
        continue
    sample = match.group(1)
    data.append([sample, sample])

df = pd.DataFrame(data, columns=["Fastq", "Sample"])
df.sort_values(by="Fastq", inplace=True)
df.drop_duplicates(inplace=True)

df.to_csv(os.path.join(args.work_dir, "samplesheet_wg0.tsv"), sep="\t", header=True, index=False)
print("Saved file of shape: {}".format(df.shape))