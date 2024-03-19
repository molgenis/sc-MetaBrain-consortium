#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse


"""
Syntax: 
./rename_fastq_files.py -h
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

print("Renaming fastq files.")
data = []
for fpath in glob.glob(os.path.join(args.work_dir, "*")):
    old_name = os.path.basename(fpath)
    match = re.match("(.+)_(L00[0-9])_001\.(R[12])\.fastq\.gz", old_name)
    if match is None:
        print("\tSkipped '{}'".format(old_name))
        continue
    sample = match.group(1)
    lane = match.group(2)
    read_type = match.group(3)
    new_name = "{sample}_S1_{lane}_{read_type}_001.fastq.gz".format(sample=sample, lane=lane, read_type=read_type)
    os.rename(fpath, os.path.join(args.work_dir, new_name))
    print("\tRenamed '{}' into '{}'".format(old_name, new_name))
    data.append([old_name, new_name])

df = pd.DataFrame(data, columns=["OldName", "NewName"])
df.sort_values(by="OldName", inplace=True)
# for _, row in df.iterrows():
#     print(row["OldName"], row["NewName"])
df.to_csv(os.path.join(args.work_dir, "rename_table.tsv"), sep="\t", header=True, index=False)
print("Saved file of shape: {}".format(df.shape))