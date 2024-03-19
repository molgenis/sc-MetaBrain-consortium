#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse


"""
Syntax: 
./create_wg0_flle_directories.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--link_table", type=str, required=True, help="")
parser.add_argument("--cellranger_dir", type=str, required=True, help="")
parser.add_argument("--cellbender_dir", type=str, required=True, help="")
parser.add_argument("--out", type=str, required=True, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import pandas as pd
import os

link_table_df = pd.read_csv(args.link_table, sep=",", header=0, index_col=None)
data = []
for _, row in link_table_df.iterrows():
    pool = row[1]
    counts = args.cellbender_dir + str(pool) + "/cellbender_remove_background_output_filtered.h5"
    barcodes = args.cellbender_dir + str(pool) + "/cellbender_remove_background_output_cell_barcodes.csv"
    bam = args.cellranger_dir + str(pool) + "/outs/possorted_genome_bam.bam"
    for filepath in [counts, barcodes, bam]:
        if not os.path.exists(filepath):
            print("Error, {} does not exist.".format(filepath))
            exit()
    data.append([pool, counts, barcodes, bam])
df = pd.DataFrame(data, columns=["Pool", "Counts", "Barcodes", "Bam"])
df.to_csv(args.out, sep="\t", header=True, index=False)