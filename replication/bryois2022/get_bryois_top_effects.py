#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse

"""
Syntax: 
./get_bryois_top_effects.py -h
"""

parser = argparse.ArgumentParser(description="")
parser.add_argument("--work_dir", type=str, required=True, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import pandas as pd
import os
import gzip

print("Loading excel file")
df = pd.read_excel(os.path.join(args.work_dir, "41593_2022_1128_MOESM3_ESM.xlsx"), sheet_name="Table S2", skiprows=3, dtype=str)
print(df)

# This is basically the same order as the original excel matrix but I add the 'pval' column, fix the 'symbol' column (excel transfered some symbols to dates), and combine symbol + ensembl into symbol_ensembl since the all files also use that format.
col_order = ["symbol_ensembl", "SNP", "effect_allele", "other_allele", "dist_TSS", "beta", "pval", "bpval", "adj_p","beta_metabrain", "p_metabrain", "Replication"]
for ct in df["cell_type"].unique():
    print("Processing '{}'".format(ct))

    # Select the cell type we are intersted in.
    ct_top_df = df.loc[df["cell_type"] == ct, :].copy()

    # Remove columns that are wrong are no longer needed.
    ct_top_df.drop(["cell_type", "symbol"], axis=1, inplace=True)

    # Round  the beta to 6 decimals to ensure a match.
    ct_top_df["beta"] = ct_top_df["beta"].astype(float).round(6).astype(str)

    # Save the order to preserve it in the end.
    row_order = list((ct_top_df["ensembl"] + "_" + ct_top_df["SNP"]).values)
    print(ct_top_df)

    # Parse the excel cell type to the file cell type.
    cell_type = ct.replace(" / ", "...").replace(" ", ".")

    # Create a dict of the eQTL effects we want to keep.
    effects = dict(zip(ct_top_df["ensembl"], ct_top_df["SNP"]))

    # Loop through each chromosome and save the effects that were in the top matrix.
    lines = []
    for chr in range(1, 23):
        print("\tLoading chromosome {}".format(chr))
        ct_chr_path = os.path.join(args.work_dir, "{}.{}.gz".format(cell_type, chr))
        if not os.path.exists(ct_chr_path):
            print("Error, could not find file '{}'".format(ct_chr_path))
            exit()

        # header is missing but equals: ["HGNC_ENSEMBL", "SNP", "Distance to TSS", "Nominal p-value", "Beta"]
        with gzip.open(ct_chr_path, "rt") as f:
            for line in f:
                symbol_ensembl, snp, dist_tss, pval, beta = line.rstrip("\n").split(" ")
                symbol, ensembl = symbol_ensembl.split("_")
                if ensembl in effects and effects[ensembl] == snp:
                    lines.append([symbol_ensembl, symbol, ensembl, snp, dist_tss, pval, beta])
        f.close()
    ct_all_df = pd.DataFrame(lines, columns=["symbol_ensembl", "symbol", "ensembl", "SNP", "dist_TSS", "pval", "beta"])

    # Round  the beta to 6 decimals to ensure a match.
    ct_all_df["beta"] = ct_all_df["beta"].astype(float).round(6).astype(str)

    # Merge the matrices and check if no rows were removed.
    ct_df = ct_top_df.merge(ct_all_df, how="inner")
    if ct_df.shape[0] != ct_top_df.shape[0]:
        print("Error, could not find all effects.")
        exit()

    # Set the index to reorder the rows.
    ct_df.index = ct_df["ensembl"] + "_" + ct_df["SNP"]
    ct_df = ct_df.loc[row_order, col_order]

    # Store the output.
    ct_df.to_csv(os.path.join(args.work_dir, "{}.top.gz".format(cell_type)), sep="\t", header=True, index=False, compression="gzip")
    print("\tSaved {}.top.gz with shape: {}".format(cell_type, ct_df.shape))

    del ct_top_df, row_order, effects, lines, ct_all_df, ct_df
