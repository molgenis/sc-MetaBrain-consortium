#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse
import os

"""
Syntax: 
./parse_bryois_ampad_metadata.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--wg0_dir", type=str, required=True, help="")
parser.add_argument("--wg3_dir", type=str, required=True, help="")
parser.add_argument("--bryois_dir", type=str, required=True, help="")
parser.add_argument("--metadata_dir", type=str, required=True, help="")
parser.add_argument("--out_folder", type=str, required=True, help="")
parser.add_argument("--dryrun", action='store_true', help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import warnings
warnings.simplefilter("ignore", UserWarning)

import pandas as pd
import glob
import gzip


def load_scmetabrain_pools():
    data = []
    for dataset in ["2023-09-10-Cain2023", "2023-09-10-Mathys2019", "2023-09-10-RocheAD2022", "2023-09-10-RocheColumbia2022","2023-09-10-RocheMS2022", "2023-09-07-Zhou2020"]:
        study = dataset.split("-")[-1]
        for count_fpath in glob.glob(os.path.join(args.wg0_dir, dataset, "CellRanger", "*", "outs", "filtered_feature_bc_matrix.h5")):
            data.append([count_fpath.split(os.sep)[-3], study])
    return pd.DataFrame(data, columns=["sample_accession_id", "Study"])

def load_bryois_links():
    df = pd.read_csv(os.path.join(args.bryois_dir, "cell_type_labels_MV", "cell_type.labels.sample_id.individual_id.txt"), sep="\t", header=0, index_col=None)
    return df

def load_bryois_keys():
    ad_keys_df = pd.read_csv(os.path.join(args.bryois_dir, "excitatory_neurons_eqtl", "meta_ad_public_with_key.txt"), sep="\t", header=0, index_col=None)
    ms_keys_df = pd.read_csv(os.path.join(args.bryois_dir, "excitatory_neurons_eqtl", "meta_ms_public_with_key.txt"), sep="\t", header=0, index_col=None)
    return pd.concat([ad_keys_df, ms_keys_df], axis=0)


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
                print("Error, not all values are unique.")
                exit()

            if not len(split_columns) == len(split_values):
                print("Error, keys are values do not have the same length.")
                exit()

            for column, value in zip(split_columns, split_values):
                if column not in data:
                    data[column] = []
                data[column].append(value)
    f.close()

    return pd.DataFrame(data, columns=columns if columns else None)

def load_bryois_delimited_maps():
    datasets = ["Bryois_RocheMS_EGAD00001009169", "Bryois_RocheAD_EGAD00001009166", "Bryois_Columbia_EGAD00001009168"]

    study_exp_run_sample_list = []
    for fpath in glob.glob(os.path.join(args.metadata_dir, "Roche", "*", "metadata", "delimited_maps", "Study_Experiment_Run_sample.map")):
        if fpath.split(os.sep)[-4] not in datasets:
            continue
        study_exp_run_sample_df = pd.read_csv(fpath, sep="\t", header=None, index_col=None)

        if study_exp_run_sample_df.shape[0] > 0:
            study_exp_run_sample_df = study_exp_run_sample_df.rename(columns={0: "study_ref_secondary_id", 1: "study_title", 3: "platform", 4: "instrument_model", 5: "library_layout", 6: "library_name", 7: "library_strategy", 8: "library_source", 9: "library_selection", 10: "experiment_id", 11: "run_accession_id", 12: "submitter_id", 14: "sample_accession_id"})
            study_exp_run_sample_df.drop([2, 13], axis=1, inplace=True)
            study_exp_run_sample_list.append(study_exp_run_sample_df)
    study_exp_run_sample_df = pd.concat(study_exp_run_sample_list, axis=0).drop_duplicates()
    del study_exp_run_sample_list

    sample_file_list = []
    for fpath in glob.glob(os.path.join(args.metadata_dir, "Roche", "*", "metadata", "delimited_maps", "Sample_File.map")):
        if fpath.split(os.sep)[-4] not in datasets:
            continue
        sample_file_df = pd.read_csv(fpath, sep="\t", header=None, index_col=None)

        if sample_file_df.shape[0] > 0:
            sample_file_df = sample_file_df.rename(columns={0: "sample_alias", 1: "sample_accession_id", 2: "filename", 3: "file_unique_accession_id"})
            sample_file_df = sample_file_df[["sample_alias", "sample_accession_id"]].drop_duplicates()
            sample_file_df["individual_id_anon"] = "Ind" + sample_file_df["sample_alias"].str.split("_", n=1, expand=True)[0]
            sample_file_df["sample_id_anon"] = sample_file_df["individual_id_anon"] + "-Sample" + sample_file_df["sample_alias"].str.split("_", n=1, expand=True)[1]
            sample_file_list.append(sample_file_df)
    sample_file_df = pd.concat(sample_file_list, axis=0).drop_duplicates()
    del sample_file_list

    return study_exp_run_sample_df.merge(sample_file_df, how="left")

def load_roche_metadata():
    key_df = load_bryois_keys()
    print(key_df)
    delimited_maps_df = load_bryois_delimited_maps()
    print(delimited_maps_df)
    print([column for column in key_df.columns if column in delimited_maps_df.columns])
    df = delimited_maps_df.merge(key_df, how="left")
    del key_df, delimited_maps_df
    return df[["sample_accession_id", "sample_id", "individual_id"]]

def load_ampad_metadata():
    datasets = ["Mathys2019", "Zhou2020"]

    bioscp_metadata_list = []
    for fpath in glob.glob(os.path.join(args.metadata_dir, "*", "metadata", "*biospecimen_metadata.csv")):
        dataset = fpath.split(os.sep)[-3]
        if dataset not in datasets:
            continue
        bioscp_metadata_df = pd.read_csv(fpath, sep=",", header=0, index_col=None, dtype=str)
        if "assay" in bioscp_metadata_df.columns:
            bioscp_metadata_df = bioscp_metadata_df.loc[(bioscp_metadata_df["assay"] == "scrnaSeq"), :].dropna(axis=1, how='all')
        bioscp_metadata_df = bioscp_metadata_df.dropna(axis=1, how='all').drop_duplicates()
        if dataset == "Mathys2019":
            bioscp_metadata_df["sample_accession_id"] = bioscp_metadata_df["projid"]
        elif dataset == "Zhou2020":
            bioscp_metadata_df["sample_accession_id"] = bioscp_metadata_df["individualID"]
        bioscp_metadata_list.append(bioscp_metadata_df)
    bioscp_metadata_df = pd.concat(bioscp_metadata_list, axis=0).drop_duplicates()
    del bioscp_metadata_list

    sc_metadata_list = []
    for fpath in glob.glob(os.path.join(args.metadata_dir, "*", "metadata", "*_metadata.csv")):
        if fpath.split("_")[-2].lower() != "scrnaseq":
            continue
        if fpath.split(os.sep)[-3] not in datasets:
            continue

        sc_metadata_df = pd.read_csv(fpath, sep=",", header=0, index_col=None, dtype=str)
        if fpath.split(os.sep)[-3] == "AMP-AD":
            sc_metadata_df = sc_metadata_df.loc[sc_metadata_df["dataContributionBatch"].isin(["DLPFC Experiment 1"]), :]
            # cain_df = load_cain_metadata()
            # print(cain_df)
            # print([col for col in cain_df.columns if col in sc_metadata_df.columns])
            # sc_metadata_df = sc_metadata_df.merge(cain_df, how="left")
            # print(sc_metadata_df)
            # exit()

        remove_cols = [col for col in ["readLength", "fileName"] if col in sc_metadata_df]
        if len(remove_cols) > 0:
            sc_metadata_df.drop(remove_cols, axis=1, inplace=True)

        sc_metadata_df = sc_metadata_df.dropna(axis=1, how='all').drop_duplicates()
        sc_metadata_list.append(sc_metadata_df)
    sc_metadata_df = pd.concat(sc_metadata_list, axis=0)
    del sc_metadata_list

    id_mapping_list = []
    for fpath in glob.glob(os.path.join(args.metadata_dir, "*", "metadata", "*id_mapping.csv")):
        if fpath.split(os.sep)[-3] not in datasets:
            continue
        id_mapping_df = pd.read_csv(fpath, sep=",", header=0, index_col=None, dtype=str)
        id_mapping_df = id_mapping_df.dropna(axis=1, how='all').drop_duplicates()
        id_mapping_list.append(id_mapping_df)
    id_mapping_df = pd.concat(id_mapping_list, axis=0).drop_duplicates()
    id_mapping_df.drop(["fastq"], axis=1, inplace=True)
    del id_mapping_list

    # TODO use raw input.
    link_table_list = []
    for fpath in glob.glob(os.path.join(args.metadata_dir, "..", "..", "processeddata", "single-cell", "metadata", "*", "*_full_link_table.tsv")):
        if fpath.split(os.sep)[-2] not in datasets:
            continue
        link_table_df = pd.read_csv(fpath, sep="\t", header=0, index_col=None, dtype=str)
        print(link_table_df)
        link_table_list.append(link_table_df)
    link_table_df = pd.concat(link_table_list, axis=0).drop_duplicates()
    del link_table_list

    df = bioscp_metadata_df.merge(link_table_df[["specimenID", "individualID", "wholeGenomeSeqID"]], how="left").merge(sc_metadata_df, how="left").merge(id_mapping_df, how="left")
    df["sample_id"] = df["specimenID"].str.split("_", n=1, expand=True)[0]
    df = df[["sample_accession_id", "sample_id", "wholeGenomeSeqID"]].drop_duplicates()
    df = df.rename(columns={"wholeGenomeSeqID": "individual_id"})
    return df

##########################################

print("Loading scMetaBrain Pools")
scmb_pool_df = load_scmetabrain_pools()
print(scmb_pool_df)

print("Loading Bryois links")
bryois_link_df = load_bryois_links()
print(bryois_link_df)

print("Loading Bryois et al. metadata")
roche_covariate_df = load_roche_metadata()
roche_metadata_df = load_roche_metadata()
print(roche_metadata_df)

print("Loading AMP-AD metadata")
ampad_metadata_df = load_ampad_metadata()
print(ampad_metadata_df)

# for sample_id in bryois_link_df["sample_id"]:
#     if sample_id not in roche_metadata_df["sample_id"].values and sample_id not in ampad_metadata_df["sample_id"].values:
#         print(sample_id)
# print("-----")
# for _, row in scmb_pool_df.iterrows():
#     if row["Pool"] not in roche_metadata_df["sample_accession_id"].values and row["Pool"] not in ampad_metadata_df["sample_accession_id"].values:
#         print(row["Pool"], row["Study"])

print("Merging data")
metadata_df = pd.concat([roche_metadata_df, ampad_metadata_df], axis=0)
metadata_df.drop(["individual_id"], axis=1, inplace=True)
df = scmb_pool_df.merge(metadata_df, how="left").merge(bryois_link_df, how="left")
print(df)

if not args.dryrun:
    print("Saving data")
    df.to_csv(os.path.join(args.wg3_dir, args.out_folder, "Bryois2022", "broiys_ampad_metadata_v2.txt.gz"), sep="\t", header=True, index=False, compression="gzip")
