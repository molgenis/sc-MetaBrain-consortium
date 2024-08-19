#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse

"""
Syntax: 
./create_pseudobulk_v2.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--metadata", type=str, required=True, help="")
parser.add_argument("--wg0_dir", type=str, required=True, help="")
parser.add_argument("--bryois_dir", type=str, required=True, help="")
parser.add_argument("--metadata_dir", type=str, required=True, help="")
parser.add_argument("--datasets", nargs="+", type=str, required=False, default=["2023-09-10-Cain2023", "2023-09-10-Mathys2019", "2023-09-10-RocheAD2022", "2023-09-10-RocheColumbia2022", "2023-09-10-RocheMS2022", "2023-09-07-Zhou2020"], help="")
parser.add_argument("--cell_types", nargs="+", type=str, required=False, default=["Astrocytes", "Endothelial.cells", "Excitatory.neurons", "Inhibitory.neurons", "Microglia", "Oligodendrocytes", "OPCs...COPs", "Pericytes"], help="")
parser.add_argument("--out_dir", type=str, required=True, help="")
parser.add_argument("--cellbender", action='store_true', help="")
parser.add_argument("--unfiltered", action='store_true', help="")
parser.add_argument("--save_annotation", action='store_true', help="")
parser.add_argument("--remove_dupl", action='store_true', help="")
parser.add_argument("--force", action='store_true', help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import warnings
warnings.simplefilter("ignore", UserWarning)

import os
import glob
import numpy as np
import pandas as pd
import scanpy
import gzip

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + "t")
    else:
        return open(file, mode)

def load_bryois_barcodes():
    ad_barcodes_df = pd.read_csv(os.path.join(args.bryois_dir, "ad.cell_type.labels.txt.gz"), sep="\t", header=0, index_col=None)
    ms_barodes_df = pd.read_csv(os.path.join(args.bryois_dir, "ms.cell_type.labels.txt.gz"), sep="\t", header=0, index_col=None)
    barcodes_df = pd.concat([ad_barcodes_df, ms_barodes_df], axis=0)
    del ad_barcodes_df, ms_barodes_df

    barcodes_df["label"] = barcodes_df["label"].str.replace(" / ", "...")
    barcodes_df["label"] = barcodes_df["label"].str.replace(" ", ".")

    barcodes_df = barcodes_df.rename(columns={"individual": "individual_id", "sample": "sample_id"})

    return barcodes_df

def aggregate_by_individual_id(in_fpath, out_fpath, key_columns=None, value_columns=None, sep="\t", header=0, index_col=None):
    # Load the sample_id data.
    df = pd.read_csv(in_fpath, sep=sep, header=header, index_col=index_col)
    print("\tLoaded dataframe: {} with shape: {}".format(os.path.basename(in_fpath), df.shape))
    print(df)

    columns = key_columns + value_columns
    df = df.loc[:, columns]
    order = df.columns
    df = df.groupby(key_columns).sum()
    df.reset_index(drop=False, inplace=True)
    df = df.loc[:, order]
    if "perc_expressed" in order:
        df["perc_expressed"] = ((df["n_expressed"] / df["n_cells"]) * 100).round(2)
    print(df)

    # Save the output.
    df.to_csv(out_fpath, sep=sep, header=False if header is None else True, index=False if index_col is None else True, compression="gzip")
    print("\tSaved dataframe: {} with shape:  {}".format(os.path.basename(out_fpath), df.shape))


def save_df(df, fpath, sep="\t", header=True, index=True):
    compression = "infer"
    if fpath.endswith(".gz"):
        compression = "gzip"

    os.makedirs(os.path.dirname(fpath), exist_ok=True)
    df.to_csv(fpath, sep=sep, header=header, index=index, compression=compression)

##########################################

print("Loading metadata")
metadata_df = pd.read_csv(args.metadata, sep="\t", header=0, index_col=None)
for column in ["sample_id", "individual_id", "genotype_id"]:
    na_sample_id_mask = metadata_df[column].isna()
    metadata_df.loc[na_sample_id_mask, column] = "NA" + metadata_df.loc[na_sample_id_mask, :].index.astype(str)
metadata_sample_ids = set(metadata_df["sample_id"].unique())

print("Loading Bryois et al. barcodes")
bryois_barcodes_df = load_bryois_barcodes()
print(bryois_barcodes_df)
barcode_sample_ids = set(bryois_barcodes_df["sample_id"].unique())

print("Check for missing sample IDs.")
missing_barcode_sample_ids = metadata_sample_ids.difference(barcode_sample_ids)
missing = False
if len(missing_barcode_sample_ids) > 0:
    print("\tWarning, {:,} samples have no barcodes".format(len(missing_barcode_sample_ids)))
    missing = True
missing_metadata_sample_ids = barcode_sample_ids.difference(metadata_sample_ids)
if len(missing_metadata_sample_ids) > 0:
    print("\tWarning, {:,} samples have no metadata".format(len(missing_metadata_sample_ids)))
    missing = True
if not missing:
    print("\tall sample IDs are present.")

# Create dictionaries from the metadata.
sample_accession_id_to_sample_id = dict(zip(metadata_df["sample_accession_id"], metadata_df["sample_id"]))
sample_accession_id_to_individual_id = dict(zip(metadata_df["sample_accession_id"], metadata_df["individual_id"]))
sample_accession_id_to_genotype_id = dict(zip(metadata_df["sample_accession_id"], metadata_df["genotype_id"]))
sample_id_to_individual_id = dict(zip(metadata_df["sample_id"], metadata_df["individual_id"]))

################################################

print("Generating summed expression file")
for dataset in args.datasets:
    print("Processing '{}'".format(dataset))

    # Create the output file.
    out_dir = os.path.join(args.out_dir, "CellBender" if args.cellbender else "CellRanger", "unfiltered" if args.unfiltered else "filtered", str(dataset))
    os.makedirs(out_dir, exist_ok=True)

    sample_stats_fpath = os.path.join(out_dir, "filtering_stats.sample_id.tsv.gz")
    sample_expr_fpath = os.path.join(out_dir, "sum_expression.sample_id.tsv.gz")
    cell_annot_fpath = os.path.join(out_dir, "cell_annotations.tsv.gz")
    wg1_metadata_fpath = os.path.join(out_dir, "Final_Assignments_demultiplexing_doublets.tsv.gz")
    wg2_metadata_fpath = os.path.join(out_dir, "azimuth_all.metadata.tsv.gz")
    if (not os.path.exists(sample_stats_fpath) or
            not os.path.exists(sample_expr_fpath) or
            (args.save_annotation and not os.path.exists(cell_annot_fpath)) or
            (args.save_annotation and not os.path.exists(wg1_metadata_fpath)) or
            (args.save_annotation and not os.path.exists(wg2_metadata_fpath)) or
            args.force):

        # Find the input files.
        count_inpath = os.path.join(args.wg0_dir, dataset, "CellRanger", "*", "outs", ("raw" if args.unfiltered else "filtered") + "_feature_bc_matrix.h5")
        if args.cellbender:
            count_inpath = os.path.join(args.wg0_dir, dataset, "CellBender", "*", "cellbender_remove_background_output" + ("" if args.unfiltered else "_filtered") + ".h5")
        print("\tLoading {}".format(count_inpath))

        count_fpaths = glob.glob(count_inpath)
        n_count_fpaths = len(count_fpaths)
        print("\t\tFound {:,} count matrices".format(n_count_fpaths))
        if n_count_fpaths == 0:
            continue

        # Open the files and write the header.
        stats_fho = gzopen(sample_stats_fpath, mode="w")
        stats_columns = ["dataset", "sample_accession_id", "individual_id", "sample_id", "n_barcodes", "n_genes", "cell_type", "n_cells_query", "n_cells_kept"]
        stats_fho.write("\t".join(stats_columns) + "\n")

        expr_fho = gzopen(sample_expr_fpath, mode="w")
        expr_columns = ["cell_type", "sample_id", "individual_id", "ensembl", "symbol", "counts", "n_expressed", "n_cells", "perc_expressed"]
        expr_fho.write("\t".join(expr_columns) + "\n")

        annot_fho = gzopen(cell_annot_fpath, mode="w")
        wg1_meta_fho = gzopen(wg1_metadata_fpath, mode="w")
        wg2_meta_fho = gzopen(wg2_metadata_fpath, mode="w")
        annot_columns = ["Barcode", "sequencing_platform", "sequencing_run", "sequencing_lane", "scrna_platform", "plate_based", "umi_based", "biomaterial", "sorting", "cell_treatment", "sample_condition"]
        wg1_meta_columns = ["Pool", "Barcode", "Assignment", "DropletType"]
        wg2_meta_columns = ["Barcode", "Pool", "L1"]

        if args.save_annotation:
            annot_fho.write("\t".join(annot_columns) + "\n")
            wg1_meta_fho.write("\t".join(wg1_meta_columns) + "\n")
            wg2_meta_fho.write("\t".join(wg2_meta_columns) + "\n")

        ################################################

        n_rows_stats = 0
        n_rows_annot = 0
        n_rows_expr = 0
        for count_index, count_fpath in enumerate(count_fpaths):
            # Extract the sample accession id from the input path (i.e. pool).
            sample_accession_id = str(count_fpath.split(os.sep)[-3])
            if args.cellbender:
                sample_accession_id = str(count_fpath.split(os.sep)[-2])
                if sample_accession_id.endswith("Run1"):
                    sample_accession_id = str(sample_accession_id[:-4])

            # Find the matching sample_id and individual_id.
            if sample_accession_id not in sample_accession_id_to_sample_id or sample_accession_id not in sample_accession_id_to_individual_id:
                print("\tError, sample_accession_id '{}' not found in metadata.".format(sample_accession_id))
                exit()
            sample_id = str(sample_accession_id_to_sample_id[sample_accession_id])
            individual_id = str(sample_accession_id_to_individual_id[sample_accession_id])
            genotype_id = sample_accession_id_to_genotype_id[sample_accession_id]

            print("\tProcessing: {}/{}\t sample_accession_id '{}'\tsample_id '{}'\tindividual_id '{}'\tgenotype_id '{}'.".format(count_index, n_count_fpaths - 1, sample_accession_id, sample_id, individual_id, genotype_id))

            # Read in the data.
            count_data = scanpy.read_10x_h5(count_fpath)
            # count_matrix = count_data.X.toarray() # Warning this makes a dense matrix and could use much RAM.
            count_matrix = count_data.X
            n_barcodes, n_genes = count_matrix.shape
            genes = count_data.var_names.to_numpy()
            barcodes = count_data.obs_names.to_numpy()
            print("\t\tLoaded matrix with {:,} barcodes and {:,} genes.".format(n_barcodes, n_genes))

            # HGNC to ensembl dict.
            gene_ids = count_data.var['gene_ids']
            gene_ids = gene_ids.reset_index(drop=False)
            gene_ids = gene_ids.to_numpy()
            # barcode_ids = count_data.var['feature_types'].to_numpy()
            if not (genes == gene_ids[:, 0]).all():
                print("Error, 'count_data.var_names' are expected to have the same order as 'count_data.var['gene_ids''.")
                exit()

            # Save the annotation info.
            if args.save_annotation:
                print("\t\tSaving annotation files.")
                sample_df = bryois_barcodes_df.loc[bryois_barcodes_df["sample_id"] == sample_id, :]
                n_cells_query = sample_df.shape[0]
                sample_ct_counts_query = sample_df["label"].value_counts()
                sample_barcode_to_ct = dict(zip(sample_df["barcode"], sample_df["label"]))

                for barcode in barcodes:
                    barcode_sample_acc_id = barcode.split("-1")[0] + "_" + sample_accession_id
                    droplet_type = "singlet" if barcode in sample_barcode_to_ct else "exclude"
                    cell_type = sample_barcode_to_ct[barcode] if barcode in sample_barcode_to_ct else "exclude"

                    annot_fho.write("\t".join([barcode_sample_acc_id, "NONE", sample_accession_id, "1", "NONE", "N", "Y", "NONE", "NONE", "UT", "NONE"]) + "\n")
                    wg1_meta_fho.write("\t".join([sample_accession_id, barcode_sample_acc_id, genotype_id, droplet_type]) + "\n")
                    wg2_meta_fho.write("\t".join([barcode_sample_acc_id, sample_accession_id, cell_type]) + "\n")

                    n_rows_annot += 1

                del sample_df, n_cells_query, sample_ct_counts_query, sample_barcode_to_ct

            # Filter the data.
            for cell_type in args.cell_types:
                # Select the barcodes that are from the current sample and the current cell type.
                keep_barcodes = set(bryois_barcodes_df.loc[(bryois_barcodes_df["sample_id"] == sample_id) & (bryois_barcodes_df["label"] == cell_type), "barcode"].values)
                barcode_mask = np.array([barcode in keep_barcodes for barcode in barcodes], dtype=bool)
                n_cells_query = len(keep_barcodes)

                # Calculate the info.
                str_cell_type = cell_type.replace("...", " / ").replace(".", " ")
                counts = np.sum(count_matrix[barcode_mask, :], axis=0)
                n_expressed = np.sum(count_matrix[barcode_mask, :] > 0, axis=0)
                n_cells = np.sum(barcode_mask)

                # Save the filter stats data.
                stats_data = [dataset, sample_accession_id, individual_id, sample_id, str(n_barcodes), str(n_genes), cell_type, str(n_cells_query), str(n_cells)]
                stats_fho.write("\t".join(stats_data) + "\n")
                n_rows_stats += 1

                # Print the info.
                print("\t\tFiltered expression matrix for {:,} {} barcodes and found {:,} barcodes.".format(n_cells_query, cell_type, n_cells))

                # No need to continue if we are not keeping any cells.
                if n_cells == 0:
                    continue

                # Calculate the other info.
                perc_expressed = np.round((n_expressed / n_cells) * 100, 2)

                # Save the pseudobulk data.
                seen = set()
                for i in range(n_genes):
                    if args.remove_dupl and genes[i] in seen:
                        continue
                    if genes[i] != gene_ids[i, 0]:
                        print("\t\tError, genes[{i}] did not match gene_ids[{i}, 0].".format(i=i))
                        exit()
                    # expr_fho.write("\t".join([str_cell_type, sample_id, individual_id, gene_ids[i, 1], genes[i], str(counts[i]), str(n_expressed[i]), str(n_cells), str(perc_expressed[i])]) + "\n") # Use in case of dense matrix.
                    expr_fho.write("\t".join([str_cell_type, sample_id, individual_id, gene_ids[i, 1], genes[i], str(counts[0, i]), str(n_expressed[0, i]), str(n_cells), str(perc_expressed[0, i])]) + "\n")
                    n_rows_expr += 1
                    seen.add(genes[i])

        stats_fho.close()
        expr_fho.close()
        annot_fho.close()
        wg1_meta_fho.close()
        wg2_meta_fho.close()
        print("\n\tSaved dataframe: {} with shape: ({}, {})".format(os.path.basename(sample_stats_fpath), n_rows_stats, len(stats_columns)))
        print("\tSaved dataframe: {} with shape: ({}, {})".format(os.path.basename(sample_expr_fpath), n_rows_expr, len(expr_columns)))
        print("\tSaved dataframe: {} with shape: ({}, {})".format(os.path.basename(cell_annot_fpath), n_rows_annot, len(annot_columns)))
        print("\tSaved dataframe: {} with shape: ({}, {})".format(os.path.basename(wg2_metadata_fpath), n_rows_annot, len(wg1_meta_columns)))
        print("\tSaved dataframe: {} with shape: ({}, {})".format(os.path.basename(wg1_metadata_fpath), n_rows_annot, len(wg2_meta_columns)))

    ################################################

    ind_stats_fpath = os.path.join(out_dir, "filtering_stats.individual_id.tsv.gz")
    ind_expr_fpath = os.path.join(out_dir, "sum_expression.individual_id.tsv.gz")
    if (not os.path.exists(ind_stats_fpath) or not os.path.exists(ind_expr_fpath)) or args.force:
        print("\nAggregating summed expression on individual_id")
        aggregate_by_individual_id(
            in_fpath=sample_stats_fpath,
            out_fpath=ind_stats_fpath,
            key_columns=["dataset", "individual_id", "cell_type"],
            value_columns=["n_barcodes", "n_cells_query", "n_cells_kept"],
        )
        aggregate_by_individual_id(
            in_fpath=sample_expr_fpath,
            out_fpath=ind_expr_fpath,
            key_columns=["cell_type", "individual_id", "ensembl", "symbol"],
            value_columns=["counts",  "n_expressed", "n_cells", "perc_expressed"],
        )

print("\nDone")