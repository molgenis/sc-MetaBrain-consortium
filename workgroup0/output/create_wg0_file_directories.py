#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse
import glob
import os
import pandas as pd

"""
Syntax: 
./create_wg0_file_directories.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--indirs", nargs="+", type=str, required=True, help="")
parser.add_argument("--cellranger", action="store_true", default=False, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

def load_cellbender_runs(indir, filenames):
    pool_fpaths = {}
    for filename in filenames:
        fpaths = glob.glob(os.path.join(indir, "CellBender", "*", filename))
        pools = [fpath.split(os.sep)[-2].split("Run")[0] for fpath in fpaths]
        if len(set(pools)) != len(pools):
            print("Error, not all pools are unique.")
            exit()
        if len(set(pools).intersection(pool_fpaths)) > 0:
            print("Error, not all pools are unique.")
            exit()
        pool_fpaths.update(dict(zip(pools, fpaths)))
    return pool_fpaths

def load_cellranger(indir, filename):
    fpaths = glob.glob(os.path.join(indir, "CellRanger", "*", "outs", filename))
    pools = [fpath.split(os.sep)[-2] for fpath in fpaths]
    return dict(zip(pools, fpaths))

########################################################################

for indir in args.indirs:
    print("\nProcessing {}".format(indir))
    indir = os.path.join(str(os.path.dirname(os.path.abspath(__file__))), indir)

    out = os.path.join(indir, "Combine_Results", "wg0_file_directories.tsv")
    outdir = os.path.dirname(out)
    if outdir != "":
        os.makedirs(outdir, exist_ok=True)

    count_files = {}
    barcode_files = {}
    if args.cellranger:
        print("Loading CellRanger count matrices")
        count_files = load_cellranger(indir=indir, filename="filtered_feature_bc_matrix.h5")
        n_count_files = len(count_files)
        print(f"  Found {n_count_files} CellRanger count matrices.")

        print("Loading CellRanger barcodes")
        barcode_files = load_cellranger(indir=indir, filename=os.sep.join(["filtered_feature_bc_matrix", "barcodes.tsv.gz"]))
        n_barcode_files = len(barcode_files)
        print(f"  Found {n_barcode_files} CellRanger barcodes.")

        if n_count_files == 0 or n_barcode_files == 0:
            exit()
    else:
        print("Loading CellBender count matrices")
        count_files = load_cellbender_runs(
            indir=indir,
            filenames=[
                "cellbender_remove_background_output_filtered.h5",
                "cellbender_feature_bc_matrix_filtered.h5"
            ])
        n_count_files = len(count_files)
        print(f"  Found {n_count_files} CellBender count matrices.")

        print("Loading CellBender barcodes")
        barcode_files = load_cellbender_runs(
            indir=indir,
            filenames=[
                "cellbender_remove_background_output_cell_barcodes.csv",
                "cellbender_feature_bc_matrix_cell_barcodes.csv"
            ])
        n_barcode_files = len(barcode_files)
        print(f"  Found {n_barcode_files} CellBender barcodes files.")

        if n_count_files == 0 or n_barcode_files == 0:
            exit()

    print("Loading CellRanger BAM files")
    bam_files = load_cellranger(indir=indir, filename="possorted_genome_bam.bam")
    n_bam_files = len(bam_files)
    print(f"  Found {n_bam_files} CellRanger BAM files.")

    print("Constructing poolsheet")
    df = pd.concat([
        pd.DataFrame(count_files, index=["Counts"]),
        pd.DataFrame(barcode_files, index=["Barcodes"]),
        pd.DataFrame(bam_files, index=["Bam"])
    ]).T.reset_index().rename(columns={"index": "Pool"})
    df = df.sort_values(by="Pool").reset_index(drop=True)
    print(df)

    print("Saving output")
    df.to_csv(out, sep="\t", header=True, index=False)

    print("Done.")
