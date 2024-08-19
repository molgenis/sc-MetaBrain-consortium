#!/usr/bin/env python
# Author: M. Vochteloo
import argparse

"""
Syntax: 
./create_cell_annotation_file.py -h
"""

parser = argparse.ArgumentParser(description="")
parser.add_argument("--wg0_indir", required=True, type=str, help="")
parser.add_argument("--wg1_indir", required=True, type=str, help="")
parser.add_argument("--metadata_indir", required=True, type=str, help="")
parser.add_argument("--metadata_outdir", required=True, type=str, help="")
parser.add_argument("--bryois_indir", required=False, type=str, help="")
parser.add_argument("--dataset", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")


import numpy as np
import pandas as pd
import glob
import os
import gzip


def gzopen(file):
    if file.endswith(".gz"):
        return gzip.open(file, 'rt')
    else:
        return open(file, 'r')

def load_file(fpath, header=None, sep="\t"):
    data = {}
    columns = []
    with gzopen(fpath) as f:
        for line_index, line in enumerate(f):
            values = line.rstrip("\n").split(sep)
            if header is not None and line_index == header:
                columns = values
                continue

            split_columns = []
            split_values = []
            for value_index, value in enumerate(values):
                if ";" not in value:
                    split_columns.append(value_index)
                    split_values.append(value)
                    continue

                for subvalue in value.split(";"):
                    if subvalue == "":
                        continue

                    split_subvalue = subvalue.split("=")
                    split_columns.append(split_subvalue[0])
                    split_values.append(split_subvalue[1])

            if not len(set(split_columns)) == len(split_columns):
                print("Error 1")
                exit()

            if not len(split_columns) == len(split_values):
                print("Error 2")
                exit()

            for column, value in zip(split_columns, split_values):
                if column not in data:
                    data[column] = []
                data[column].append(value)
    f.close()

    return pd.DataFrame(data, columns=columns if columns else None)

def load_mathys_metadata():
    # df1 = pd.read_csv(os.path.join(args.metadata_indir, args.dataset, "metadata", "snRNAseqPFC_BA10_Sample_key.csv"), sep=",", header=0, index_col=None, dtype=str)
    # print(df1)

    df2 = pd.read_csv(os.path.join(args.metadata_indir, args.dataset, "metadata", "snRNAseqPFC_BA10_assay_scRNAseq_metadata.csv"), sep=",", header=0, index_col=None, dtype=str)
    df2.drop(["readLength", "fileName"], axis=1, inplace=True)
    df2 = df2.drop_duplicates()

    df3 = pd.read_csv(os.path.join(args.metadata_indir, args.dataset, "metadata", "snRNAseqPFC_BA10_biospecimen_metadata.csv"), sep=",", header=0, index_col=None, dtype=str)

    df4 = pd.read_csv(os.path.join(args.metadata_indir, args.dataset, "metadata", "snRNAseqPFC_BA10_id_mapping.csv"), sep=",", header=0, index_col=None, dtype=str)
    df4.drop(["fastq"], axis=1, inplace=True)
    df4 = df4.drop_duplicates()

    df = df3.merge(df2, how="left").merge(df4, how="left")
    df["scrna_platform"] = "10x Chromium Single Cell 3′ Reagent Kits v.2"

    return df

def load_zhou_metadata():
    df2 = pd.read_csv(os.path.join(args.metadata_indir, args.dataset, "metadata", "snRNAseqAD_TREM2_assay_scRNAseq_metadata.csv"), sep=",", header=0, index_col=None, dtype=str)

    df3 = pd.read_csv(os.path.join(args.metadata_indir, args.dataset, "metadata", "snRNAseqAD_TREM2_biospecimen_metadata.csv"), sep=",", header=0, index_col=None, dtype=str)

    df = df3.merge(df2, how="left")
    df["scrna_platform"] = "10x Chromium Single Cell 5’"

    return df

def load_cain_metadata():
    df1 = pd.read_csv(os.path.join(args.metadata_indir, "AMP-AD", "metadata", "ROSMAP_biospecimen_metadata.csv"), sep=",", header=0, index_col=None)
    df1 = df1.loc[(df1["assay"] == "scrnaSeq"), :].dropna(axis=1, how='all')
    df1.drop(['assay'], axis=1, inplace=True)
    # print(df1)
    # print(len(df1["individualID"].unique()))
    # print(len(df1["specimenID"].unique()))
    # print(df1.loc[df1["individualID"] == "R1721075", :])

    df2 = pd.read_csv(os.path.join(args.metadata_indir, "AMP-AD", "metadata", "ROSMAP_assay_scrnaSeq_metadata.csv"), sep=",", header=0, index_col=None)
    df2.drop(['assay'], axis=1, inplace=True)
    # print(df2)
    # print(len(df2["specimenID"].unique()))

    df = df1.merge(df2, how="left")
    # print(df)
    # print(df.loc[df["individualID"] == "R1721075", :])
    # exit()

    return df

def load_fujita_metadata():
    df = load_cain_metadata()
    df["id"] = df["rnaBatch"] + "_" + df["platformLocation"]
    return df

def load_rosmap_metadata():
    if args.dataset == "Mathys2019":
        df = load_mathys_metadata()
    elif args.dataset == "Zhou2020":
        df = load_zhou_metadata()
    elif args.dataset == "Cain2023":
        df = load_cain_metadata()
    elif args.dataset == "Fujita2022":
        df = load_fujita_metadata()
    else:
        print("Error, unepxected --dataset in load_rosmap_metadata()")
        exit()

    df3 = pd.read_csv(os.path.join(args.metadata_indir, "AMP-AD", "metadata", "ROSMAP_clinical.csv"), sep=",", header=0, index_col=None)
    df3["projid"] = df3["projid"].astype(str)

    df = df.merge(df3, how="left")
    print(df)

    df["sample_condition"] = "healthy"
    df.loc[~df["age_first_ad_dx"].isna(), "sample_condition"] = "Alzheimer's disease"
    # TODO: this is not resulting in 24;24 for Mathys2019

    return df

def load_bryois_metadata():
    df1_list = []
    # df2_list = []
    df3_list = []
    df4_list = []
    for study in ["Bryois_RocheMS_EGAD00001009169", "Bryois_RocheAD_EGAD00001009166", "Bryois_Columbia_EGAD00001009168"]:
        if study == "Bryois_Columbia_EGAD00001009168":
            df1 = pd.read_csv(os.path.join(args.bryois_indir, "key_columbia_rosmap.txt"), sep="\t", header=0, index_col=None)
            if df1.shape[0] > 0:
                df1["case_or_control"] = "control"
                df1.loc[df1["diagnosis"] != "Ctrl", "case_or_control"] = "case"
                df1["phenotype"] = df1["diagnosis"].map({"Ctrl": "healthy", "AD": "Alzheimer's disease"})
                df1["organism_part"] = "Brain"
                df1 = df1.rename(columns={"tissue": "organism_region"})
                print(df1)
                df1_list.append(df1)

                # Missing: SAMEA, gender, region, ENA-CHECKLIST
        else:
            df1 = load_file(os.path.join(args.metadata_indir, "Roche", study, "metadata", "delimited_maps", "Analysis_Sample_meta_info.map"))
            if df1.shape[0] > 0:
                df1 = df1.rename(columns={0: "SAMEA", 1: "sample_id_anon", "subject_id": "individual_id_anon"})
                df1[["organism_part", "organism_region"]] = df1["organism_part"].str.split(" - ", n=1, expand=True)
                print(df1)
                df1_list.append(df1)

        # df2 = load_file(os.path.join(args.metadata_indir, "Roche", study, "metadata", "delimited_maps", "Run_Sample_meta_info.map")).drop_duplicates()
        # if df2.shape[0] > 0:
        #     df2_list.append(df2)

        df3 = pd.read_csv(os.path.join(args.metadata_indir, "Roche", study, "metadata", "delimited_maps", "Sample_File.map"), sep="\t", header=None, index_col=None)
        if df3.shape[0] > 0:
            df3 = df3.rename(columns={0: "sample_alias", 1: "sample_accession_id", 2: "filename", 3: "file_unique_accession_id"})
            df3 = df3[["sample_alias", "sample_accession_id"]].drop_duplicates()
            df3["individual_id_anon"] = "Ind" + df3["sample_alias"].str.split("_", n=1, expand=True)[0]
            df3["sample_id_anon"] = df3["individual_id_anon"] + "-Sample" + df3["sample_alias"].str.split("_", n=1, expand=True)[1]
            df3_list.append(df3)

        df4 = pd.read_csv(os.path.join(args.metadata_indir, "Roche", study, "metadata", "delimited_maps", "Study_Experiment_Run_sample.map"), sep="\t", header=None, index_col=None)
        if df4.shape[0] > 0:
            df4 = df4.rename(columns={0: "study_ref_secondary_id", 1: "study_title", 3: "platform", 4: "instrument_model", 5: "library_layout", 6: "library_name", 7: "library_strategy", 8: "library_source", 9: "library_selection", 10: "experiment_id", 11: "run_accession_id", 12: "submitter_id", 14: "sample_accession_id"})
            df4.drop([2, 13], axis=1, inplace=True)
            df4_list.append(df4)

        del df1, df3, df4

    df1 = pd.concat(df1_list, axis=0)
    # df2 = pd.concat(df2_list, axis=0)
    df3 = pd.concat(df3_list, axis=0).drop_duplicates()
    df4 = pd.concat(df4_list, axis=0)
    del df1_list, df3_list, df4_list

    return df4.merge(df3, how="left").merge(df1, how="left")

def load_bryois_keys():
    ad_keys_df = pd.read_csv(os.path.join(args.bryois_indir, "excitatory_neurons_eqtl", "meta_ad_public_with_key.txt"), sep="\t", header=0, index_col=None)
    ms_keys_df = pd.read_csv(os.path.join(args.bryois_indir, "excitatory_neurons_eqtl", "meta_ms_public_with_key.txt"), sep="\t", header=0, index_col=None)
    keys_df = pd.concat([ad_keys_df, ms_keys_df], axis=0)
    del ad_keys_df, ms_keys_df

    keys_df = keys_df.rename(columns={"individual_id_anon": "subject_id"})

    return keys_df

def load_barcodes_file(method="CellBender", pool_name="Pool"):
    dataset = glob.glob(os.path.join(args.wg0_indir, "*" + args.dataset))
    if len(dataset) == 1:
        dataset = dataset[0]
    else:
        print("Error in load_barcodes_file()")
        exit()

    inpath = None
    if method == "CellBender":
        inpath = os.path.join(args.wg0_indir, dataset, "CellBender", "*", "cellbender_feature_bc_matrix_cell_barcodes.csv")
    elif method == "CellRanger":
        inpath = os.path.join(args.wg0_indir, dataset, "CellRanger", "*", "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
    else:
        print("Error, unexpected method input '{}'".format(method))
        exit()

    barcodes_dfs = []
    for fpath in glob.glob(inpath):
        pool = fpath.replace(os.path.join(args.wg0_indir, dataset), "").split(os.sep)[2].replace("Run1", "")
        barcodes_fpath = os.path.join(fpath)
        barcodes_df = pd.read_csv(barcodes_fpath, sep="\t", header=None, index_col=None)
        barcodes_df.columns = ["barcode"]
        barcodes_df[pool_name] = pool
        barcodes_df["method"] = method
        barcodes_dfs.append(barcodes_df)

    return pd.concat(barcodes_dfs, axis=0)

def load_barcodes(pool_name="Pool"):
    cellranger_barcodes_df = load_barcodes_file(method="CellRanger", pool_name=pool_name)
    cellranger_barcodes_df["CellRanger"] = True
    cellbender_barcodes_df = load_barcodes_file(method="CellBender", pool_name=pool_name)
    cellbender_barcodes_df["CellBender"] = True
    barcodes_df = cellbender_barcodes_df.merge(cellranger_barcodes_df, how="outer", on=["barcode", pool_name])
    barcodes_df.drop(["method_x", "method_y"], axis=1, inplace=True)
    barcodes_df["CellRanger"] = barcodes_df["CellRanger"].fillna(False)
    barcodes_df["CellRanger"] = barcodes_df["CellRanger"].fillna(False)
    barcodes_df["Barcode"] = barcodes_df["barcode"].str.split("-", n=1, expand=True)[0] + "_" + barcodes_df[pool_name].astype(str)
    return barcodes_df

def load_demultiplex():
    dataset = glob.glob(os.path.join(args.wg1_indir, "*" + args.dataset + "*"))
    if len(dataset) == 1:
        dataset = dataset[0]
    else:
        print("Error in load_demultiplex()")
        exit()

    df = pd.read_csv(os.path.join(args.wg1_indir, dataset, "CombinedResults", "Final_Assignments_demultiplexing_doublets.tsv.gz"), sep="\t", header=None, index_col=None)
    print(df)
    exit()

def load_bryois_barcodes():
    ad_barcodes_df = pd.read_csv(os.path.join(args.bryois_indir, "ad.cell_type.labels.txt.gz"), sep="\t", header=0, index_col=None)
    ms_barodes_df = pd.read_csv(os.path.join(args.bryois_indir, "ms.cell_type.labels.txt.gz"), sep="\t", header=0, index_col=None)
    barcodes_df = pd.concat([ad_barcodes_df, ms_barodes_df], axis=0)
    del ad_barcodes_df, ms_barodes_df

    barcodes_df["label"] = barcodes_df["label"].str.replace(" / ", "_")
    barcodes_df["label"] = barcodes_df["label"].str.replace(" ", "_")
    barcodes_df["Bryois2022"] = True

    return barcodes_df

def load_mathys_barcodes():
    mathys_df = pd.read_csv(os.path.join(args.metadata_indir, "Mathys2019", "processed", "filtered_column_metadata.txt"), sep="\t", header=0, index_col=None, dtype=str)
    mathys_df["barcode"] = mathys_df["TAG"].str.split(".", n=1, expand=True)[0] + "-1"
    mathys_df["Mathys2019"] = True
    mathys_df = mathys_df[["barcode", "projid", "broad.cell.type", "Subcluster", "Mathys2019"]]
    mathys_df = mathys_df.rename({"broad.cell.type": "Mathys_ct_label", "Subcluster": "Mathys_ct_subcluster"})
    return mathys_df

print("Loading full link table")
full_link_df = pd.read_csv(os.path.join(args.metadata_outdir, args.dataset, args.dataset + "_full_link_table.tsv"), sep="\t", header=0, index_col=None, dtype=str)
if args.dataset == "RocheMS2022" or args.dataset == "RocheAD2022":
    full_link_df = full_link_df.rename(columns={"Ind": "subject_id", "sampleID": "sample_accession_id"})
    full_link_df["sample_id_anon"] = full_link_df["id"].str.split("_", n=1, expand=True)[0]
elif args.dataset == "RocheColumbia2022":
    full_link_df = full_link_df.rename(columns={"individual_id_anon": "subject_id", "SampleID": "sample_accession_id"})
elif args.dataset == "Cain2023":
    full_link_df = full_link_df.rename(columns={"scrnaSeqID": "specimenID"})
else:
    pass
print(full_link_df)

print("\nMerging metadata")
if args.dataset in ["RocheMS2022", "RocheAD2022", "RocheColumbia2022"]:
    print("\nLoading Bryois metadata")
    bryois_keys_df = load_bryois_keys()
    bryois_metadata_df = load_bryois_metadata()
    metadata_df = full_link_df.merge(bryois_metadata_df, how='left').merge(bryois_keys_df, how='left')
elif args.dataset in ["Mathys2019", "Zhou2020", "Cain2023", "Fujita2022"]:
    print("\nLoading ROSMAP metadata")
    rosmap_metadata_df = load_rosmap_metadata()
    metadata_df = full_link_df.merge(rosmap_metadata_df, how='left')
else:
    print("Unexpected --dataset 1")
    exit()
print(metadata_df)

print("\nSaving metadata file")
metadata_df.to_csv(os.path.join(args.metadata_outdir, args.dataset, args.dataset + "_full_metadata_table.tsv"), sep="\t", index=False, header=True)
print("\tSaved data frame with shape: {}".format(metadata_df.shape))

#################################################################

pool_name = None
if args.dataset in ["RocheMS2022", "RocheAD2022", "RocheColumbia2022"]:
    pool_name = "sample_accession_id"
elif args.dataset in ["Mathys2019"]:
    pool_name = "projid"
elif args.dataset in ["Zhou2020", "Cain2023"]:
    pool_name = "individualID"
elif args.dataset in ["Fujita2022"]:
    pool_name = "id"
else:
    print("Unexpected --dataset 2")
    exit()

print("\nLoading barcodes data")
barcodes_df = load_barcodes(pool_name=pool_name)
print(barcodes_df)
# print(barcodes_df["id"].value_counts())
# print(barcodes_df["CellRanger"].value_counts())
# print(barcodes_df["CellBender"].value_counts())

# print("\nLoading demultiplexing output")
# demultiplex_df = load_demultiplex()

# print("\nLoading Bryois barcodes data")
# bryois_barcodes_df = load_bryois_barcodes()
# bryois_barcodes_df = bryois_barcodes_df.rename(columns={"individual": "individual_id", "sample": "sample_id", "label": "Bryois_ct_label"})
# print(bryois_barcodes_df)
#
# print("\nLoading Mathys 2019 barcodes data")
# mathys_df = load_mathys_barcodes()
# mathys_df = mathys_df.rename({"broad.cell.type": "Mathys_ct_label", "Subcluster": "Mathys_ct_subcluster"})
# print(mathys_df)

print("\nBuilding cell annotation file")
if args.dataset in ["RocheMS2022", "RocheAD2022", "RocheColumbia2022"]:
    # df = barcodes_df.merge(metadata_df, how="left").merge(bryois_barcodes_df, how="left")
    df = barcodes_df.merge(metadata_df, how="left")
    df = df.rename(columns={
        "instrument_model": "sequencing_platform",
        "run_accession_id": "sequencing_run",
        "experiment_id": "sequencing_lane",
        "library_name": "scrna_platform",
        "organism_part": "biomaterial",
        "phenotype": "sample_condition"
    })
    df["plate_based"] = "N"
    df["umi_based"] = "Y"
    df["cell_treatment"] = "UT"
elif args.dataset in ["Mathys2019", "Zhou2020", "Cain2023", "Fujita2022"]:
    # df = barcodes_df.merge(metadata_df, how="left").merge(mathys_df, how="left")
    df = barcodes_df.merge(metadata_df, how="left")
    df = df.rename(columns={
        "platform": "sequencing_platform",
        "tissue": "biomaterial",
        "phenotype": "sample_condition"
    })
    df["sequencing_run"] = df["projid"].astype(str) + "_run" + df["sequencingBatch"].astype(str)
    df["sequencing_lane"] = df["projid"].astype(str) + "_lane1"
    df["plate_based"] = "N"
    df["umi_based"] = "Y"
    df["cell_treatment"] = "UT"
else:
    print("Unexpected --dataset 3")
    exit()

default_columns = [
    ("Barcode", None),
    ("sequencing_platform", "NONE"),
    ("sequencing_run", "NONE"),
    ("sequencing_lane", "NONE"),
    ("scrna_platform", "NONE"),
    ("plate_based", "NONE"),
    ("umi_based", "NONE"),
    ("biomaterial", "NONE"),
    ("sorting", "NONE"),
    ("cell_treatment", "NONE"),
    ("sample_condition", "NONE")
]
order = []
for default_column, default_value in default_columns:
    order.append(default_column)
    if default_column not in df.columns:
        if default_value is None:
            print("\tError, missing {}".format(default_column))
            exit()
        else:
            print("\tAdding: {} = {}".format(default_column, default_value))
            df[default_column] = default_value
df = df[order + ["CellBender"]]
print(df)

print("Saving cell annotation file")
df.loc[df["CellBender"], :].to_csv(os.path.join(args.metadata_outdir, args.dataset, args.dataset + "_cell_annotations_CellBender.tsv.gz"), sep="\t", index=False, header=True, compression="gzip")
print("Saved data frame with shape: {}".format(df.shape))

df.loc[df["CellRanger"], :].to_csv(os.path.join(args.metadata_outdir, args.dataset, args.dataset + "_cell_annotations_CellRanger.tsv.gz"), sep="\t", index=False, header=True, compression="gzip")
print("Saved data frame with shape: {}".format(df.shape))

#################################################################