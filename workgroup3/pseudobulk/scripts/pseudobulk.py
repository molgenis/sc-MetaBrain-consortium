#!/usr/bin/env python
# Author: M. Vochteloo

import warnings
warnings.simplefilter("ignore", UserWarning)

import numpy as np
import pandas as pd
import scanpy
import argparse
import gzip
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--poolsheet", required=True, type=str, help="")
parser.add_argument("--pool", required=True, type=str, help="")
parser.add_argument("--psam", required=True, type=str, help="")
parser.add_argument("--cell_annotation", required=True, type=str, help="")
parser.add_argument("--droplet_type_annotation", required=True, type=str, help="")
parser.add_argument("--cell_type_annotation", required=True, type=str, help="")
parser.add_argument("--ct_pairing", required=False, type=str, default=None, help="")
parser.add_argument("--rb_genes", required=True, type=str, help="")
parser.add_argument("--mt_genes", required=True, type=str, help="")
parser.add_argument("--sample_aggregate", required=False, type=str, default="Assignment", help="")
parser.add_argument("--ancestry", required=False, type=str, default="EUR", help="")
parser.add_argument("--cell_level", required=False, type=str, default="L1", help="")
parser.add_argument("--ncount_rna", required=False, type=int, default=500, help="")
parser.add_argument("--nfeature_rna", required=False, type=int, default=0, help="")
parser.add_argument("--percent_rb", required=False, type=int, default=0, help="")
parser.add_argument("--percent_mt", required=False, type=int, default=5, help="")
parser.add_argument("--malat1", required=False, type=int, default=0, help="")
parser.add_argument("--feature_name", required=False, type=str, default="HGNC", choices=["HGNC", "ENSG", "HGNC_ENSG"], help="")
parser.add_argument("--aggregate_method", required=False, type=str, default="sum", choices=["sum", "mean"], help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

if not os.path.isdir(args.out):
    os.makedirs(args.out, exist_ok=True)

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


def load_file(fpath, sep="\t", header=0, index_col=None, must_contain=None):
    data = []
    columns = None
    with gzopen(fpath, mode="r") as fh:
        for index, line in enumerate(fh):
            if index == 0 or (index + 1) % 1e5 == 0:
                print("\tParsed {:,} lines".format(index + 1), end='\r')

            if index != header and must_contain is not None and must_contain not in line:
                continue

            values = line.rstrip("\n").split(sep)
            if index == header:
                columns = values
                continue
            data.append(values)
    print("\tParsed {:,} lines".format(index + 1), end='\n')

    df = pd.DataFrame(data, columns=columns)
    if index_col is not None:
        if len(set(df.columns)) != df.shape[1]:
            print("Error, not all columns are unique.")
            exit()

        column = df.columns[index_col]
        df.index = df[column]
        df.drop(column, axis=1, inplace=True)

    print("\tLoaded {} with shape: {}".format(os.path.basename(fpath), df.shape))
    return df

def load_metadata(barcodes_fpath):
    print("  Loading Barcodes")
    barcodes = load_file(barcodes_fpath, header=None).rename(columns={0: "Barcode"})
    barcodes.reset_index(drop=False, inplace=True)
    barcodes["Barcode"] = barcodes["Barcode"].str.split("-", n=1, expand=True)[0] + "_" + str(args.pool)
    barcodes["Pool"] = args.pool

    print("  Loading cell annotations")
    cell_annot = load_file(args.cell_annotation, must_contain=args.pool)

    print("  Loading droplet type annotations")
    droplet_annot = load_file(args.droplet_type_annotation, must_contain=args.pool)

    print("  Loading cell type annotations")
    type_annot = load_file(args.cell_type_annotation, must_contain=args.pool)

    print("  Loading PSAM")
    psam = load_file(args.psam).rename(columns={"#FID": "FID"})

    ct_pairs = None
    if args.ct_pairing is not None:
        print("  Loading cell type pairing")
        ct_pairs = load_file(args.ct_pairing, sep=";")

    print("\nMerging metadata ...")
    metadata = barcodes.merge(cell_annot, on="Barcode", how="inner").merge(droplet_annot, on=["Pool", "Barcode"], how="inner").merge(type_annot, on=["Pool", "Barcode"], how="left").merge(psam, left_on="Assignment", right_on="IID", how="left")
    if ct_pairs is not None:
        metadata = metadata.merge(ct_pairs, on=ct_pairs.columns[1], how="left")

    # Rorder back into original order.
    metadata = metadata.sort_values(by="index", ascending=True).drop(["index"], axis=1)

    # Check if all barcodes are unique.
    if len(metadata["Barcode"].uniuqe()) != metadata.shape[0]:
        print("Error, not all barcodes are unique in the metadata.")

    print("\tLoaded metadata with shape: {}".format(metadata.shape))
    return metadata

def load_counts(counts_fpath):
    count_data = scanpy.read_10x_h5(counts_fpath)
    count_matrix = count_data.X
    n_barcodes, n_genes = count_matrix.shape
    genes = count_data.var_names.to_numpy()
    barcodes = count_data.obs_names.to_numpy()
    if np.size(genes) != n_genes:
        print("Error, matrix size does not match gene annotation.")
        exit()
    if np.size(barcodes) != n_barcodes:
        print("Error, matrix size does not match barcodes annotation.")
        exit()
    print("\tLoaded matrix with {:,} barcodes and {:,} genes.".format(n_barcodes, n_genes))

    # Check for unique barcodes.
    if len(np.unique(barcodes)) != np.size(barcodes):
        print("Error, not all barcodes are unique.")
        exit()

    # Parse the gene info
    gene_ids = count_data.var['gene_ids']
    gene_ids = gene_ids.reset_index(drop=False)
    gene_ids = gene_ids.to_numpy()
    if not (genes == gene_ids[:, 0]).all():
        print("Error, 'count_data.var_names' are expected to have the same order as 'count_data.var['gene_ids''.")
        exit()

    # Set gene names.
    features = None
    if args.feature_name == "HGNC":
        features = gene_ids[:, 0]
    elif args.feature_name == "ENSG":
        features = gene_ids[:, 1]
    elif args.feature_name == "HGNC_ENSG":
        features = np.char.add(np.char.add(gene_ids[:, 0], "_"), gene_ids[:, 1])
    else:
        print("Unexpected feature name '{}'".format(args.feature_name))
        exit()

    # Remove duplicates.
    u, c = np.unique(features, return_counts=True)
    dup = u[c > 1]
    features_mask = np.in1d(features, dup)
    if np.sum(features_mask) > 0:
        print("\tRemoving {:,} duplicate genes.".format(np.sum(features_mask)))
        count_matrix = count_matrix[:, ~features_mask]
        features = features[~features_mask]


    # Parse the barcode info.
    # barcode_ids = count_data.var['feature_types'].to_numpy()

    print("\tLoaded expression with shape: {}".format(count_matrix.shape))
    return count_matrix, barcodes, features

def load_features(fpath):
    gene_ids = load_file(fpath=fpath).to_numpy()

    if args.feature_name == "HGNC":
        features = gene_ids[:, 1]
    elif args.feature_name == "ENSG":
        features = gene_ids[:, 0]
    elif args.feature_name == "HGNC_ENSG":
        features = np.char.add(np.char.add(gene_ids[:, 1], "_"), gene_ids[:, 0])
    else:
        print("Unexpected feature name '{}'".format(args.feature_name))
        exit()

    return features

def get_malat1_feature():
    if args.feature_name == "HGNC":
        return "MALAT1"
    elif args.feature_name == "ENSG":
        return "ENSG00000251562"
    elif args.feature_name == "HGNC_ENSG":
        return "MALAT1_ENSG00000251562"
    else:
        print("Unexpected feature name '{}'".format(args.feature_name))
        exit()

def calc_barcode_qc(counts, barcodes, features):
    # First calculate the total number of counts and features per barcode.
    info = pd.DataFrame(np.hstack((np.sum(counts, axis=1), np.sum(counts != 0, axis=1))), columns=["nCount_RNA", "nFeature_RNA"], index=barcodes)

    print("  Adding complexity")
    info["complexity"] = np.nan
    mask = info["nCount_RNA"] > 0
    info.loc[mask, "complexity"] = info.loc[mask, "nFeature_RNA"] / info.loc[mask, "nCount_RNA"]
    del mask

    print("  Adding ribosomal %")
    rb_features = load_features(args.rb_genes)
    rb_mask = np.in1d(features, rb_features)
    info["rb"] = np.sum(counts[:, rb_mask], axis=1)

    info["percent.rb"] = np.nan
    mask = info["nCount_RNA"] > 0
    info.loc[mask, "percent.rb"] = (info.loc[mask, "rb"] / info.loc[mask, "nCount_RNA"]) * 100
    del rb_features, rb_mask, mask

    print("  Adding mitochondiral %")
    mt_features = load_features(args.mt_genes)
    mt_mask = np.in1d(features, mt_features)
    info["mt"] = np.sum(counts[:, mt_mask], axis=1)

    info["percent.mt"] = np.nan
    mask = info["nCount_RNA"] > 0
    info.loc[mask, "percent.mt"] = (info.loc[mask, "mt"] / info.loc[mask, "nCount_RNA"]) * 100
    del mt_features, mt_mask, mask

    print("  Adding MALAT1")
    info["MALAT1"] = 0
    malat1_feature = get_malat1_feature()
    malat1_mask = features == malat1_feature
    info["MALAT1"] = counts[:, malat1_mask].toarray()

    return info

def apply_barcode_qc(barcode_qc):
    n_barcodes = barcode_qc.shape[0]
    print("  Input: {:,} barcodes".format(n_barcodes))

    ncount_rna_mask = barcode_qc["nCount_RNA"] >= args.ncount_rna
    print("\tnCount_RNA >= {} yields {:,} [{:.2f}%] barcodes".format(args.ncount_rna, np.sum(ncount_rna_mask), (100 / n_barcodes) * np.sum(ncount_rna_mask)))

    nfeature_rna_mask = barcode_qc["nFeature_RNA"] >= args.nfeature_rna
    print("\tnFeature_RNA >= {} yields {:,} [{:.2f}%] barcodes".format(args.nfeature_rna, np.sum(nfeature_rna_mask), (100 / n_barcodes) * np.sum(nfeature_rna_mask)))

    percent_rb_mask = barcode_qc["percent.rb"] <= args.percent_rb
    print("\tpercent_rb <= {} yields {:,} [{:.2f}%] barcodes".format(args.percent_rb, np.sum(percent_rb_mask), (100 / n_barcodes) * np.sum(percent_rb_mask)))

    percent_mt_mask = barcode_qc["percent.mt"] <= args.percent_mt
    print("\tpercent_rb <= {} yields {:,} [{:.2f}%] barcodes".format(args.percent_mt, np.sum(percent_mt_mask), (100 / n_barcodes) * np.sum(percent_mt_mask)))

    malat1_mask = barcode_qc["MALAT1"] >= args.malat1
    print("\tMALAT1 >= {} yields {:,} [{:.2f}%] barcodes".format(args.malat1, np.sum(malat1_mask), (100 / n_barcodes) * np.sum(malat1_mask)))

    mask = ncount_rna_mask & nfeature_rna_mask & percent_rb_mask & percent_mt_mask & malat1_mask
    print("  Pass-QC: {:,} [{:.2f}%] barcodes".format(np.sum(mask), (100 / n_barcodes) * np.sum(mask)))

    barcode_qc["tag"] = "Outlier"
    barcode_qc.loc[mask, "tag"] = "NotOutlier"

    return barcode_qc

def pseudobulk_per_ct(counts, features, metadata):
    ncells = counts.shape[0]

    expr_data = []
    expr_columns = []
    cell_data = []
    samples = metadata[args.sample_aggregate].unique()
    for sample in samples:
        if sample == "doublet" or sample == "unassigned":
            continue

        sample_mask = metadata[args.sample_aggregate] == sample
        n_sample_cells = np.sum(sample_mask)
        print("  {} has {:,} input cells [{:.2f}%]".format(sample, n_sample_cells, (100 / ncells) * n_sample_cells))

        if n_sample_cells == 0:
            continue

        cell_types = metadata.loc[sample_mask, args.cell_level].unique()
        for cell_type in cell_types:
            if pd.isnull(cell_type):
                continue

            sample_ct_mask = (metadata[args.sample_aggregate] == sample) & (metadata[args.cell_level] == cell_type)
            n_sample_ct_cells = np.sum(sample_ct_mask)
            print("\t{} has {:,} input cells".format(cell_type, n_sample_ct_cells))

            mask1 = sample_ct_mask & (metadata["DropletType"] == "singlet")
            mask2 = sample_ct_mask & (metadata["tag"] == "NotOutlier")
            mask3 = sample_ct_mask & (metadata["Provided_Ancestry"] == args.ancestry)
            mask4 = sample_ct_mask & (metadata["cell_treatment"] == "UT")
            mask = sample_ct_mask & mask1 & mask2 & mask3 & mask4

            print("\t  DropletType - Singlet: N = {:,} [{:.2f}%]".format(np.sum(mask1), (100 / n_sample_ct_cells) * np.sum(mask1)))
            print("\t  tag - NotOutlier: N = {:,} [{:.2f}%]".format(np.sum(mask2), (100 / n_sample_ct_cells) * np.sum(mask2)))
            print("\t  Provided_Ancestry - {}: N = {:,} [{:.2f}%]".format(args.ancestry, np.sum(mask3), (100 / n_sample_ct_cells) * np.sum(mask3)))
            print("\t  cell_treatment - UT: N = {:,} [{:.2f}%]".format(np.sum(mask4), (100 / n_sample_ct_cells) * np.sum(mask4)))
            print("\t{} has {:,} cells pass QC [{:.2f}%]\n".format(cell_type, np.sum(mask), (100 / n_sample_ct_cells) * np.sum(mask)))

            n_sample_ct_pass_cells = np.sum(mask)
            cell_data.append([sample, cell_type, n_sample_ct_pass_cells])
            if n_sample_ct_pass_cells == 0:
                continue

            if args.aggregate_method == "sum":
                expr_data.append(np.sum(counts[mask, :], axis=0))
            elif args.aggregate_method == "mean":
                expr_data.append(np.mean(counts[mask, :], axis=0))
            else:
                print("Error, unexpected aggregate_method type {}".format(args.aggregate_method))
                exit()

            expr_columns.append(sample + "_" + cell_type)

    expr_df = None
    if len(expr_data) > 0:
        expr_df = pd.DataFrame(np.vstack(expr_data), columns=features, index=expr_columns)

    cell_df = None
    if len(cell_data) > 0:
        cell_df = pd.DataFrame(cell_data, columns=["sample", "cell type", "cells"])
        print("\tCombined {:,} cells over {:,} samples.".format(cell_df["cells"].sum(), expr_df.shape[1] if expr_df is not None else 0))
        print(cell_df)

        for cell_type in cell_df["cell type"].unique():
            ncells = cell_df.loc[cell_df["cell type"] == cell_type, "cells"].sum()
            print("\t  Cell type: {} has {:,} cells over all samples".format(cell_type, ncells))

    return expr_df, cell_df

#############################################

print("\nExtracting input fpaths from poolsheet ...")
fpaths = load_file(args.poolsheet, must_contain=args.pool).to_dict("index")
if len(fpaths) != 1:
    print("Error, pool is not unique.")
    exit()
if fpaths[0]["Pool"] != args.pool:
    print("Error, unexpected Pool value.")
    exit()

counts_fpath = fpaths[0]["Counts"]
barcodes_fpath = fpaths[0]["Barcodes"]
print("  Counts fpath: " + counts_fpath)
print("  Barcodes fpath: " + barcodes_fpath)

print("\nCreating metadata ...")
metadata = load_metadata(barcodes_fpath=barcodes_fpath)

print("\tSaving file")
metadata.to_csv(os.path.join(args.out, str(args.pool) + ".full.metadata.tsv.gz"), sep="\t", header=True, index=False, compression="gzip")

print("\nLoading counts matrix ...")
counts, barcodes, features = load_counts(counts_fpath=counts_fpath)

print("\nCalculating barcode QC stats ...")
barcode_qc = calc_barcode_qc(counts=counts, barcodes=barcodes, features=features)

print("\tSaving file")
barcode_qc.to_csv(os.path.join(args.out, str(args.pool) + ".qc_metrics.tsv.gz"), sep="\t", header=True, index=True, compression="gzip")

print("\nApplying barcode QC ...")
barcode_qc = apply_barcode_qc(barcode_qc=barcode_qc)

print("\nMerging metadata and barcode QC ...")
# Add the Barcode_Pool column.
barcode_qc.reset_index(drop=False, inplace=True)
barcode_qc = barcode_qc.rename(columns={"index": "barcode"})
barcode_qc["Barcode"] = barcode_qc["barcode"].str.split("-", n=1, expand=True)[0] + "_" + str(args.pool)

# Merge on Barcode and reorder by the count matrix (barcode column).
metadata = metadata.merge(barcode_qc, on="Barcode", how="inner").set_index("barcode").loc[barcodes, :]

print("\tSaving file")
metadata.to_csv(os.path.join(args.out, str(args.pool) + ".full.metadata.tsv.gz"), sep="\t", header=True, index=False, compression="gzip")

print("\nPseudobulking per cell type ...")
expr, cell = pseudobulk_per_ct(counts=counts, features=features, metadata=metadata)
print(expr)
print(cell)

print("\tSaving file")
if expr is not None:
    expr.to_csv(os.path.join(args.out, str(args.pool) + ".pseudobulk.tsv.gz"), sep="\t", header=True, index=True, compression="gzip")
if cell is not None:
    cell.to_csv(os.path.join(args.out, str(args.pool) + ".pseudobulk.cells.tsv.gz"), sep="\t", header=True, index=True, compression="gzip")

print("Done")