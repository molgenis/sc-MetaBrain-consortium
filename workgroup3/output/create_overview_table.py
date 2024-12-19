#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse
import os

"""
Syntax: 
./create_overview_table.py -h
"""

parser = argparse.ArgumentParser(description="")
parser.add_argument("--workdir", type=str, required=True, help="")
parser.add_argument("--cell_qc_setting", type=str, required=False, default="500nCountRNA_0nFeatureRNA_100Complex_100PcntRB_5PcntMT_0MALAT1_DynamicCapBarcodes_FalseCRBarcodes", help="")
parser.add_argument("--merge_setting", type=str, required=False, default="5Cells", help="")
parser.add_argument("--norm_setting", type=str, required=False, default="10Obs_1CPM_BryoisNorm", help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import glob
import numpy as np
import pandas as pd

cell_types = [os.path.basename(fpath) for fpath in glob.glob(os.path.join(args.workdir, "mbQTL", "*"))]

exclude_cols = ["Sample", "Cell type", "ncells_pool"]
with pd.ExcelWriter(os.path.join(args.workdir, 'overview.xlsx')) as writer:
    # Load the eQTL summary stats.
    overview = pd.read_csv(os.path.join(args.workdir, "mbQTL", "output", "mbQTL-results.txt"), sep="\t", header=0, index_col=0)
    overview.index = [value.split(".")[0] for value in overview["Prefix"]]
    if args.norm_setting.endswith("BryoisNorm"):
        overview.insert(loc=2, column="NGenesAligned", value=np.nan)
        overview.insert(loc=3, column="NVariableGenes", value=np.nan)
        overview.insert(loc=4, column="NObservedGenes", value=np.nan)
        overview.insert(loc=5, column="NExpressedGenes", value=np.nan)
    else:
        print("Error, not implemented.")
        exit()
    overview.insert(loc=9, column="ExpressionN", value=np.nan)
    overview.insert(loc=10, column="NPCAOutliers", value=np.nan)
    overview.insert(loc=11, column="eQTLMaxN", value=np.nan)
    overview.insert(loc=14, column="nCells", value=np.nan)

    # Save.
    overview.to_excel(writer,
                      sheet_name="overview",
                      na_rep="NA",
                      index=True,
                      header=True)

    for ct in cell_types:
        if ct == "output":
            continue
        print(f'Processing cell type: {ct}')

        # Load the number of genes.
        gene_stats = pd.read_csv(os.path.join(args.workdir, "expression", args.cell_qc_setting, args.merge_setting, args.norm_setting, ct + ".pseudobulk.TMM.stats.tsv"), sep="\t", header=0, index_col=1)

        # Add the summary stats to the overview sheet.
        overview.loc[ct, "NGenesAligned"] = gene_stats.loc["Input", "ngenes"]
        overview.loc[ct, "NVariableGenes"] = gene_stats.loc["Varfilter", "ngenes"]
        overview.loc[ct, "NObservedGenes"] = gene_stats.loc["MinObsFilter", "ngenes"]
        overview.loc[ct, "NExpressedGenes"] = gene_stats.loc["CPMFilter", "ngenes"]

        # Load the number of cells.
        pb_stats = pd.read_csv(os.path.join(args.workdir, "expression", args.cell_qc_setting, args.merge_setting, ct + ".pseudobulk.stats.tsv"), sep="\t", header=0, index_col=0)
        dataset_count = pb_stats["Dataset"].value_counts().to_frame().rename(columns={"count": "NSamples"})
        df = dataset_count.merge(pb_stats[[col for col in pb_stats.columns if col not in exclude_cols]].groupby("Dataset").sum(), left_index=True, right_index=True)
        del dataset_count

        # Add the PCA outlier removal info.
        for suffix, colname, colname_suffix in [("exclude", "PCAOutlier", "Removed"), ("include", "eQTL", "Kept")]:
            # Add the number of samples that were exclude by PCA outlier removal.
            gte_subset = pd.read_csv(os.path.join(args.workdir, "mbQTL", ct, "pca", "all", "AllSamples",f"{ct}.pseudobulk.TMM.samples_{suffix}.txt"), sep="\t", header=None, index_col=None)
            dataset_count = gte_subset[2].value_counts().to_frame().rename(columns={"count": f"N{colname}Samples"})
            df = df.merge(dataset_count, left_index=True, right_index=True, how="left")
            df[f"N{colname}Samples"] = df[f"N{colname}Samples"].fillna(0)
            del dataset_count

            # Also add how many cells were removed and kept.
            pb_sats_removed = pb_stats.loc[pb_stats["Sample"].isin(gte_subset[0]), ["Dataset", "ncells"]].groupby("Dataset").sum().rename(columns={"ncells": f"{colname}NCells{colname_suffix}"})
            df = df.merge(pb_sats_removed, left_index=True, right_index=True, how="left")
            df[f"{colname}NCells{colname_suffix}"] = df[f"{colname}NCells{colname_suffix}"].fillna(0)

        # Add the total row.
        df.loc["total", :] = df.sum(axis=0)
        print(df)

        # Add the summary stats to the overview sheet.
        overview.loc[ct, "ExpressionN"] = df.loc["total", "NSamples"]
        overview.loc[ct, "NPCAOutliers"] = df.loc["total", "NPCAOutlierSamples"]
        overview.loc[ct, "eQTLMaxN"] = df.loc["total", "NeQTLSamples"]
        overview.loc[ct, "nCells"] = df.loc["total", "eQTLNCellsKept"]

        # Save.
        df.to_excel(writer,
                    sheet_name=ct,
                    na_rep="NA",
                    index=True,
                    header=True)

    # Update overview.
    overview.to_excel(writer,
                      sheet_name="overview",
                      na_rep="NA",
                      index=True,
                      header=True)

print("Done")