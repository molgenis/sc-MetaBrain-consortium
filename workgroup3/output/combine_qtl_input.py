#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse
import os

"""
Syntax: 
./combine_qtl_input.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--work_dir", type=str, required=True, help="")
parser.add_argument("--datasets", nargs="+", type=str, required=False, default=["2024-06-27-Cain2023", "2024-06-27-Mathys2019", "2024-06-27-RocheAD2022", "2024-06-27-RocheColumbia2022", "2024-06-27-RocheMS2022", "2024-06-27-Zhou2020"], help="")
parser.add_argument("--ancestry", type=str, required=False, default="EUR", help="")
parser.add_argument("--cell_level", type=str, required=False, default="L1", help="")
parser.add_argument("--cell_types", nargs="+", type=str, required=False, default=["AST", "END", "EX", "IN", "MIC", "OLI", "OPC", "PER"], help="")
parser.add_argument("--out_folder", type=str, required=True, help="")
parser.add_argument("--dryrun", action='store_true', help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import pandas as pd

def path_exists(fpath):
    parts = fpath.split(os.sep)
    for i in range(2, len(parts) + 1):
        sub_fpath = os.sep.join(parts[:i])
        if not os.path.exists(sub_fpath):
            print("\tError, '{}' does not exist.".format(sub_fpath))
            return False

    return True

for cell_type in args.cell_types:
    smf_df_list = []
    expr_df_list = []
    cov_df_list = []

    print("Processing '{}'".format(cell_type))

    for dataset in args.datasets:
        # Loading sample mapping file.
        smf_fpath = os.path.join(args.work_dir, dataset, "expression_input", args.ancestry, args.cell_level, cell_type, "PostQC", cell_type + ".smf.txt")
        if not path_exists(fpath=smf_fpath):
            continue
        print("\tLoading '{}'".format(smf_fpath))
        smf_df = pd.read_csv(smf_fpath, sep="\t", header=None, index_col=None)
        smf_df["dataset"] = dataset
        smf_df_list.append(smf_df)
        del smf_fpath, smf_df

        # Loading expression matrix.
        expr_fpath = os.path.join(args.work_dir, dataset, "expression_input", args.ancestry, args.cell_level, cell_type, "PostQC", cell_type + ".qtlInput.txt")
        if not path_exists(fpath=expr_fpath):
            continue
        print("\tLoading '{}'".format(expr_fpath))
        expr_df = pd.read_csv(expr_fpath, sep="\t", header=0, index_col=0)
        expr_df_list.append(expr_df)
        del expr_fpath, expr_df

        # Loading expression matrix.
        cov_fpath = os.path.join(args.work_dir, dataset, "expression_input", args.ancestry, args.cell_level, cell_type, "PostQC", cell_type + ".covariates.txt")
        if not path_exists(fpath=cov_fpath):
            continue
        print("\tLoading '{}'".format(cov_fpath))
        cov_df = pd.read_csv(cov_fpath, sep="\t", header=0, index_col=None)
        cov_df_list.append(cov_df)
        del cov_fpath, cov_df

    print("Combining files")
    smf_df = pd.concat(smf_df_list, axis=0)
    del smf_df_list
    print(smf_df)

    out_smf_path = os.path.join(args.work_dir, args.out_folder, "expression_input", args.ancestry, args.cell_level, cell_type, "PostQC", cell_type + ".gte.txt")
    os.makedirs(os.path.dirname(out_smf_path), exist_ok=True)
    if not args.dryrun:
        print("\tSaving '{}'".format(out_smf_path))
        smf_df.to_csv(out_smf_path, sep="\t", header=False, index=False)
    del smf_df, out_smf_path

    expr_df = pd.concat(expr_df_list, axis=1).dropna(axis=0)
    del expr_df_list
    print(expr_df)

    out_expr_path = os.path.join(args.work_dir, args.out_folder, "expression_input", args.ancestry, args.cell_level, cell_type, "PostQC", cell_type + ".qtlInput.txt")
    os.makedirs(os.path.dirname(out_expr_path), exist_ok=True)
    if not args.dryrun:
        print("\tSaving '{}'".format(out_expr_path))
        expr_df.to_csv(out_expr_path, sep="\t", header=True, index=True)
    del expr_df, out_expr_path

    cov_df = pd.concat(cov_df_list, axis=0)
    del cov_df_list
    print(cov_df)
    print("N cells: {:,}".format(cov_df["CellCount"].sum()))

    out_cov_path = os.path.join(args.work_dir, args.out_folder, "expression_input", args.ancestry, args.cell_level, cell_type, "PostQC", cell_type + ".covariates.txt")
    os.makedirs(os.path.dirname(out_cov_path), exist_ok=True)
    if not args.dryrun:
        print("\tSaving '{}'".format(out_cov_path))
        cov_df.to_csv(out_cov_path, sep="\t", header=True, index=True)
    del cov_df, out_cov_path

print("Done")