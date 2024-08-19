#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse
import os

"""
Syntax: 
./exclude_duplicate_samples.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--work_dir", type=str, required=True, help="")
parser.add_argument("--datasets", nargs="+", type=str, required=False, default=["2024-06-27-Cain2023", "2024-06-27-Mathys2019", "2024-06-27-RocheAD2022", "2024-06-27-RocheColumbia2022", "2024-06-27-RocheMS2022", "2024-06-27-Zhou2020"], help="")
parser.add_argument("--ancestry", type=str, required=False, default="EUR", help="")
parser.add_argument("--cell_level", type=str, required=False, default="L1", help="")
parser.add_argument("--cell_types", nargs="+", type=str, required=False, default=["AST", "END", "EX", "IN", "MIC", "OLI", "OPC", "PER"], help="")
parser.add_argument("--individual_aggregate", type=str, required=False, default="Assignment", help="")
parser.add_argument("--sample_aggregate", type=str, required=False, default="Assignment_Run_Lane", help="")
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
    ind_sample_cellcount_info = {}

    print("Processing '{}'".format(cell_type))

    print("\nLoading covariate info from datasets.")
    for dataset in args.datasets:
        sample_aggregate = args.sample_aggregate
        if dataset in ["2024-06-27-RocheAD2022", "2024-06-27-RocheMS2022"]:
            sample_aggregate = args.individual_aggregate

        # Loading covariates
        cov_fpath = os.path.join(args.work_dir, dataset, "expression_input", args.ancestry, args.cell_level, cell_type, "PreQC", cell_type + ".covariates.txt")
        if not path_exists(fpath=cov_fpath):
            continue
        print("\tLoading '{}'".format(cov_fpath))
        cov_df = pd.read_csv(cov_fpath, sep="\t", header=0, index_col=None, dtype=str)
        cov_df["id"] = cov_df[args.individual_aggregate] + "+" + cov_df[sample_aggregate]
        print("\t\tLoaded '{}' samples".format(cov_df.shape[0]))
        del cov_fpath

        # Loading exclusions.
        excl_fpath = os.path.join(args.work_dir, dataset, "expression_input", args.ancestry, args.cell_level, cell_type, "manual_selection", cell_type + "_PreQC_exclude_smf.txt")
        if not path_exists(fpath=excl_fpath):
            continue

        print("\tLoading '{}'".format(excl_fpath))
        try:
            excl_df = pd.read_csv(excl_fpath, sep="\t", header=None, index_col=None, dtype=str)
            excl_df["id"] = excl_df.iloc[:, 0] + "+" + excl_df.iloc[:, 1]
            print("\t\tExcluding '{}' samples".format(excl_df.shape[0]))
            cov_df = cov_df.loc[~cov_df["id"].isin(excl_df["id"]), :]
            del excl_df
        except pd.errors.EmptyDataError:
            pass

        # Saving the individuals.
        for _, row in cov_df.iterrows():
            ind = row[args.individual_aggregate]
            sample = row[sample_aggregate]
            cell_count = float(row["CellCount"])
            if ind not in ind_sample_cellcount_info:
                ind_sample_cellcount_info[ind] = {}
            ind_sample_cellcount_info[ind][dataset] = (ind, sample, cell_count)

        del sample_aggregate, cov_df, excl_fpath

    print("\nFinding best sample per individual.")
    exclude_samples = {}
    n_issues = 0
    for ind, ind_samples in ind_sample_cellcount_info.items():
        if len(ind_samples) > 1:
            n_issues += 1

            # Means multiple datasets.
            max_cells = 0
            max_dataset = None
            info = []
            for dataset, (ind, sample, cell_count) in ind_samples.items():
                info.append(dataset + ": " + str(cell_count))
                if cell_count > max_cells:
                    max_cells = cell_count
                    max_dataset = dataset

            print("\tIndividual '{}' occurs in dataset '{}'. Keeping '{}'.".format(ind, ", ".join(info), max_dataset))
            for dataset, (ind, sample, cell_count) in ind_samples.items():
                if dataset == max_dataset:
                    continue
                if dataset not in exclude_samples:
                    exclude_samples[dataset] = []
                exclude_samples[dataset].append([ind, sample])
    print("\tFound {} issues.".format(n_issues))
    del ind_sample_cellcount_info

    print("\nWriting exclusion files.")
    for dataset, exclude_ids in exclude_samples.items():
        print("\tProcessing '{}'".format(dataset))

        excl_fpath = os.path.join(args.work_dir, dataset, "expression_input", args.ancestry, args.cell_level, cell_type, "manual_selection", cell_type + "_PreQC_exclude_smf.txt")

        exclude_rows = []
        try:
            excl_df = pd.read_csv(excl_fpath, sep="\t", header=None, index_col=None, dtype=str)
            if not args.dryrun:
                ori_excl_fpath = excl_fpath.replace(".txt", "_original.txt")
                print("\tSaving '{}'".format(ori_excl_fpath))
                excl_df.to_csv(ori_excl_fpath, sep="\t", header=False, index=False)
            for _, row in excl_df.iterrows():
                exclude_rows.append([row[0], row[1]])
        except pd.errors.EmptyDataError:
            pass

        print("\t\t{} loaded".format(len(exclude_rows)))

        for exclude_id in exclude_ids:
            exclude_rows.append(exclude_id)

        if not args.dryrun:
            print("\tSaving '{}'".format(excl_fpath))
            pd.DataFrame(exclude_rows).to_csv(excl_fpath, sep="\t", header=False, index=False)
        print("\t\t{} written".format(len(exclude_rows)))

        del exclude_rows