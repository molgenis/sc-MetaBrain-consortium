#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse

"""
Syntax: 
./compare_cell_annotation.py -h
"""

parser = argparse.ArgumentParser(description="")
parser.add_argument("--scmbd_indir", type=str, required=True, help="")
parser.add_argument("--scmd_metadata", type=str, required=True, help="")
parser.add_argument("--cell_mapping", type=str, required=True, help="")
parser.add_argument("--bryois_indir", type=str, required=True, help="")
parser.add_argument("--gte", type=str, required=True, help="")
parser.add_argument("--outdir", type=str, required=True, help="")
parser.add_argument("--extensions", nargs="+", type=str, choices=["png", "pdf", "eps"], default=["png"], help="The figure file extension.. Default: 'png'.")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#############################################

print("Loading Bryois data ...")
bryois_ad_df = pd.read_csv(os.path.join(args.bryois_indir, "ad.cell_type.labels.txt.gz"), sep="\t", header=0, index_col=None)
bryois_ms_df = pd.read_csv(os.path.join(args.bryois_indir, "ms.cell_type.labels.txt.gz"), sep="\t", header=0, index_col=None)
bryois_df = pd.concat([bryois_ad_df, bryois_ms_df], axis=0)
bryois_df.rename(columns={"individual": "individual_id", "sample": "sample_id"}, inplace=True)
del bryois_ad_df, bryois_ms_df

print(bryois_df)
print(bryois_df["label"].value_counts())
# Oligodendrocytes     507.786
# Excitatory neurons   442.186
# Astrocytes           165.416
# Inhibitory neurons   118.709
# Microglia             92.062
# OPCs / COPs           59.597
# Endothelial cells     27.120
# Pericytes             13.978

print(bryois_df[["individual_id", "sample_id"]].drop_duplicates())

#############################################

# Final_Assignments_demultiplexing_doublets.tsv.gz
# /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2024-01-10-scMetaBrain-Workgroup3DownstreamAnalyses/2024-04-04-Freeze1Bryois/expression_input/metadata.tagged.tsv.gz
# azimuth_all.metadata.tsv.gz

# Final_Assignments_demultiplexing_doublets_no_cellbender.tsv.gz
# /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2024-01-10-scMetaBrain-Workgroup3DownstreamAnalyses/2024-05-06-Freeze1Bryois-NoCellBender/expression_input/metadata.tagged.tsv.gz
# azimuth_all.metadata.no_cellbender.tsv.gz

print("Loading scMetaBrain data ...")
scmb_wg1_df = pd.read_csv(os.path.join(args.scmbd_indir, "Final_Assignments_demultiplexing_doublets_no_cellbender.tsv.gz"), sep="\t", header=0, index_col=None)
scmb_wg1_tag_df = pd.read_csv("/groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2024-01-10-scMetaBrain-Workgroup3DownstreamAnalyses/2024-05-06-Freeze1Bryois-NoCellBender/expression_input/metadata.tagged.tsv.gz", sep="\t", header=0, index_col=None)
scmb_wg2_df = pd.read_csv(os.path.join(args.scmbd_indir, "azimuth_all.metadata.no_cellbender.tsv.gz"), sep="\t", header=0, index_col=None)
l1_df = pd.read_csv(args.cell_mapping, sep=";", header=0, index_col=None)
scmb_df = scmb_wg1_df.merge(scmb_wg1_tag_df, on=["Pool", "Barcode"], how="left").merge(scmb_wg2_df, on=["Pool", "Barcode"], how="inner").merge(l1_df, on="predicted.subclass", how="left")
del scmb_wg1_df, scmb_wg1_tag_df, scmb_wg2_df, l1_df

print(scmb_df)
print(scmb_df["L1"].value_counts())
# EX   1.696.773
# OLI    786.618
# AST    211.294
# IN     188.120
# MIC    134.512
# OPC     62.573
# PER     36.400
# END     36.028

print(scmb_df.loc[scmb_df["DropletType"] == "singlet", "L1"].value_counts())
print(scmb_df.loc[scmb_df["tag"] == "NotOutlier", "L1"].value_counts())
print(scmb_df.loc[(scmb_df["DropletType"] == "singlet") & (scmb_df["tag"] == "NotOutlier"), "L1"].value_counts())

print(scmb_df[["Pool", "Assignment"]].drop_duplicates())

#############################################

# print("Loading Bryois metadata ...")
# bryois_ct_labels_smf_df = bryois_df[["individual_id", "sample_id"]].drop_duplicates()
# # bryois_ct_labels_smf_df.to_csv(os.path.join(args.bryois_indir, "cell_type.labels.sample_id.individual_id.txt"), sep="\t", header=True, index=False)
# # bryois_ct_labels_smf_df = pd.read_csv(os.path.join(args.bryois_indir, "cell_type.labels.sample_id.individual_id.txt"), sep="\t", header=0, index_col=None)
# print(bryois_ct_labels_smf_df)
#
# ms_key_df = pd.read_csv(os.path.join(args.bryois_indir, "excitatory_neurons_eqtl", "meta_ms_public_with_key.txt"), sep="\t", header=0, index_col=None)
# ad_key_df = pd.read_csv(os.path.join(args.bryois_indir, "excitatory_neurons_eqtl", "meta_ad_public_with_key.txt"), sep="\t", header=0, index_col=None)
# key_df = pd.concat([ms_key_df, ad_key_df], axis=0)
# del ms_key_df, ad_key_df
# print(key_df)
#
# smf_df = bryois_ct_labels_smf_df.merge(key_df, on=["individual_id", "sample_id"], how="left")
# del bryois_ct_labels_smf_df, key_df
# print(smf_df)
#
# trans_sample_anon = dict(zip(smf_df["sample_id_anon"], smf_df["sample_id"]))
#
# print("Loading scMetaBrain metadata ...")
# bryois_to_scmbd_sample_id = {}
# bryois_to_scmbd_individual_id = {}
# for dataset in ["Mathys2019", "RocheAD2022", "RocheColumbia2022", "RocheMS2022", "Zhou2020", "Cain2023"]:
#     full_link_df = pd.read_csv(os.path.join(args.scmd_metadata, dataset, dataset + "_full_link_table.tsv"), sep="\t", header=0, index_col=None)
#     keep = []
#     if dataset == "Mathys2019":
#         bryois_to_scmbd_sample_id.update(dict(zip(
#             full_link_df["specimenID"].str.split("_", n=1, expand=True)[0],
#             full_link_df["projid"].astype(str)
#         )))
#         bryois_to_scmbd_individual_id.update(dict(zip(
#             full_link_df["specimenID"].str.split("_", n=1, expand=True)[0],
#             full_link_df["wholeGenomeSeqID"]
#         )))
#     elif dataset in ["RocheAD2022", "RocheMS2022"]:
#         bryois_to_scmbd_sample_id.update(dict(zip(
#             [trans_sample_anon[sample_id_anon] if sample_id_anon in trans_sample_anon else sample_id_anon for sample_id_anon in full_link_df["id"].str.split("_", n=1, expand=True)[0]],
#             full_link_df["sampleID"]
#         )))
#         bryois_to_scmbd_individual_id.update(dict(zip(
#             [trans_sample_anon[sample_id_anon] if sample_id_anon in trans_sample_anon else sample_id_anon for sample_id_anon in full_link_df["id"].str.split("_", n=1, expand=True)[0]],
#             full_link_df["Ind"]
#         )))
#     elif dataset == "RocheColumbia2022":
#         bryois_to_scmbd_sample_id.update(dict(zip(
#             [trans_sample_anon[sample_id_anon] if sample_id_anon in trans_sample_anon else sample_id_anon for sample_id_anon in full_link_df["sample_id_anon"]],
#             full_link_df["SampleID"]
#         )))
#         bryois_to_scmbd_individual_id.update(dict(zip(
#             [trans_sample_anon[sample_id_anon] if sample_id_anon in trans_sample_anon else sample_id_anon for sample_id_anon in full_link_df["sample_id_anon"]],
#             full_link_df["individual_id"]
#         )))
#     elif dataset == "Zhou2020":
#         full_link_df.loc[full_link_df["specimenID"] == "P7", "wholeGenomeSeqID"] = "MAP50104134"
#
#         bryois_to_scmbd_sample_id.update(dict(zip(
#             full_link_df["specimenID"],
#             full_link_df["individualID"]
#         )))
#         bryois_to_scmbd_individual_id.update(dict(zip(
#             full_link_df["specimenID"],
#             full_link_df["wholeGenomeSeqID"]
#         )))
#     elif dataset == "Cain2023":
#         cain_sample_id_dict = {}
#         for sample_id in smf_df["sample_id"]:
#             if not sample_id.startswith("MFC"):
#                 continue
#             splitted_sample_id = sample_id.split("-")
#             cain_sample_id_dict["-".join(splitted_sample_id[:-1])] = splitted_sample_id[-1]
#
#         full_link_df["sample_id"] = [scrna_seq_id + "-" + cain_sample_id_dict[scrna_seq_id] if scrna_seq_id in cain_sample_id_dict else scrna_seq_id for scrna_seq_id in full_link_df["scrnaSeqID"]]
#
#         bryois_to_scmbd_sample_id.update(dict(zip(
#             full_link_df["sample_id"],
#             full_link_df["individualID"]
#         )))
#         bryois_to_scmbd_individual_id.update(dict(zip(
#             full_link_df["sample_id"],
#             full_link_df["wholeGenomeSeqID"]
#         )))
#
# smf_df["scmbd_sample_id"] = smf_df["sample_id"].map(bryois_to_scmbd_sample_id)
# smf_df["scmbd_individual_id"] = smf_df["sample_id"].map(bryois_to_scmbd_individual_id)
# print(smf_df)
#
# smf_df.to_csv(os.path.join(args.outdir, "smf.txt"), sep="\t", header=True, index=False)
smf_df = pd.read_csv(os.path.join(args.outdir, "smf.txt"), sep="\t", header=0, index_col=None)
print(smf_df)

###############################################################

# scmb_ct_labels_smf_df = scmb_df[["Pool", "Assignment"]].drop_duplicates().dropna()
# scmb_ct_labels_smf_df.columns = ["scmbd_sample_id", "scmbd_individual_id"]
# scmb_ct_labels_smf_df.to_csv(os.path.join(args.outdir, "cell_type.labels.pool.assignment.txt"), sep="\t", header=True, index=False)
# # scmb_ct_labels_smf_df = pd.read_csv(os.path.join(args.outdir, "cell_type.labels.pool.assignment.txt"), sep="\t", header=0, index_col=None)
# print(scmb_ct_labels_smf_df)

# smf_index = smf_df["scmbd_sample_id"] + "_" + smf_df["scmbd_individual_id"]
# scmb_index = scmb_ct_labels_smf_df["scmbd_sample_id"] + "_" + scmb_ct_labels_smf_df["scmbd_individual_id"]
# for index in scmb_index:
#     if index in smf_index.values:
#         # print(index, " FOUND")
#         pass
#     else:
#         print(index, " NOT FOUND")

###############################################################


def create_confusion_matrix(df, x, y):
    row_counts = list(zip(*np.unique(df[y], return_counts=True)))
    row_counts.sort(key=lambda x: x[0])
    row_labels = ["{}\n[n={:,.0f}]".format(label, size) for label, size in row_counts]

    col_counts = list(zip(*np.unique(df[x], return_counts=True)))
    col_counts.sort(key=lambda x: x[0])
    col_labels = ["{}\n[n={:,.0f}]".format(label, size) for label, size in col_counts]

    ct_trans = {
        "Oligodendrocytes": "OLI",
        "Excitatory neurons": "EX",
        "Astrocytes": "AST",
        "Inhibitory neurons": "IN",
        "Microglia": "MIC",
        "OPCs / COPs": "OPC",
        "Endothelial cells": "END",
        "Pericytes": "PER"
    }

    confusion_df = pd.DataFrame(np.nan, index=row_labels, columns=col_labels)
    annotation_df = pd.DataFrame("", index=row_labels, columns=col_labels)
    tp_count = 0
    tp_total = 0
    for row_label, (row_value, row_count) in zip(row_labels, row_counts):
        for col_label, (col_value, _) in zip(col_labels, col_counts):
            n_overlap = df.loc[(df[y] == row_value) & (df[x] == col_value), :].shape[0]
            frac_overlap = np.nan
            if n_overlap > 0:
                frac_overlap = n_overlap / row_count
                annotation_df.loc[row_label, col_label] = "{:.2f}\nn={:,.0f}".format(frac_overlap, n_overlap)

            if row_value == ct_trans[col_value]:
                tp_count += n_overlap
                tp_total += row_count
            confusion_df.loc[row_label, col_label] = frac_overlap

    tpr = np.nan
    if tp_total > 0:
        tpr = tp_count / tp_total

    return confusion_df, annotation_df, tpr


def plot_heatmap(df, annot_df, xlabel="", ylabel="", title="", outfile="plot"):
    cmap = sns.diverging_palette(246, 24, as_cmap=True)

    fig, axes = plt.subplots(nrows=2,
                             ncols=2,
                             figsize=(1 * df.shape[1] + 5, 1 * df.shape[0] + 5),
                             gridspec_kw={"width_ratios": [0.2, 0.8],
                                          "height_ratios": [0.8, 0.2]})
    sns.set(color_codes=True)

    row_index = 0
    col_index = 0
    for _ in range(4):
        ax = axes[row_index, col_index]
        if row_index == 0 and col_index == 1:

            sns.heatmap(df, cmap=cmap, vmin=-1, vmax=1, center=0,
                        square=True, annot=annot_df, fmt='',
                        cbar=False, annot_kws={"size": 12},
                        ax=ax)

            plt.setp(ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=20, rotation=0))
            plt.setp(ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=20, rotation=90))

            ax.set_xlabel(xlabel, fontsize=14)
            ax.xaxis.set_label_position('top')

            ax.set_ylabel(ylabel, fontsize=14)
            ax.yaxis.set_label_position('right')

            ax.set_title(title, fontsize=40)
        else:
            ax.set_axis_off()

        col_index += 1
        if col_index > 1:
            col_index = 0
            row_index += 1

    for extension in args.extensions:
        fig.savefig(os.path.join(args.outdir, "{}.{}".format(outfile, extension)))
    plt.close()

###############################################################

print("Merging data data ...")

bryois_df = bryois_df.merge(smf_df[["sample_id", "scmbd_sample_id", "scmbd_individual_id"]], on="sample_id", how="left")
# print(bryois_df[["sample_id", "individual_id"]].drop_duplicates())
# print(bryois_df[["scmbd_sample_id", "scmbd_individual_id"]].drop_duplicates())

scmb_df["barcode"] = scmb_df["Barcode"].str.split("_", n=1, expand=True)[0] + "-1"
print(scmb_df)

df = bryois_df.merge(scmb_df[["Pool", "Assignment", "barcode", "DropletType", "tag", "L1"]], left_on=["scmbd_sample_id", "scmbd_individual_id", "barcode"], right_on=["Pool", "Assignment", "barcode"], how="inner")
print(df)
print(df.loc[df["label"] == "Excitatory neurons", ["DropletType"]].value_counts())
print(df.loc[df["label"] == "Excitatory neurons", ["tag"]].value_counts())
print(df.loc[df["label"] == "Excitatory neurons", ["L1"]].value_counts())
# print(df[["sample_id", "scmbd_sample_id"]].drop_duplicates())

df.to_csv(os.path.join(args.outdir, "merged_cell_assignments.txt.gz"), sep="\t", header=True, index=False, compression="gzip")
# df = pd.read_csv(os.path.join(args.outdir, "merged_cell_assignments.txt.gz"), sep="\t", header=0, index_col=None)
print(df)

labels = df["label"].unique()
labels.sort()

confusion_df, annotation_df, tpr = create_confusion_matrix(df=df, x="label", y="L1")
print(confusion_df)

plot_heatmap(
    df=confusion_df,
    annot_df=annotation_df,
    title="TPR: {:.3f}".format(tpr),
    outfile="Bryois_label_vs_scMetaBrain_L1_confusion_matrix")
