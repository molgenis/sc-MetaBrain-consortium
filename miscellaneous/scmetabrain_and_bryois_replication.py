#!/usr/bin/env python3

"""
File:         scmetabrain_and_bryois_replication.py
Created:      2023/04/20
Last Changed: 2023/04/21
Author:       M.Vochteloo

Copyright (C) 2020 M.Vochteloo
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""

# Standard imports.
from __future__ import print_function
from datetime import datetime
import argparse
import glob
import re
import os

# Third party imports.
import numpy as np
import pandas as pd
import h5py
from statsmodels.stats import multitest
import rpy2.robjects as robjects
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adjustText import adjust_text

# Local application imports.

"""
Syntax:
./scmetabrain_and_bryois_replication.py \
    --work_dir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-20-ReplicateInBryois \
    --dataset_outdir 2023-04-21-Mathys2019 \
    --wg3_folder /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-11-WorkGroup3eQTLAndDEA/2023-04-12-Mathys2019 \
    --exclude_ct EX IN END AST OLI OPC MIC \
    --bryois_folder /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-20-ReplicateInBryois/Bryois2022
    
"""

# Metadata
__program__ = "scMetaBrain and Bryois et al. 2022 Replication"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@rug.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.work_dir = getattr(arguments, 'work_dir')
        dataset_outdir = getattr(arguments, 'dataset_outdir')
        self.wg3_folder = getattr(arguments, 'wg3_folder')
        self.annotation_level = getattr(arguments, 'annotation_level')
        self.exclude_ct = getattr(arguments, 'exclude_ct')
        self.bryois_folder = getattr(arguments, 'bryois_folder')
        self.bryois_n = 196
        self.extensions = getattr(arguments, 'extensions')
        self.force = getattr(arguments, 'force')
        self.verbose = getattr(arguments, 'verbose')
        self.qvalues_script = getattr(arguments, 'qvalues')
        self.rb_script = getattr(arguments, 'rb')

        # Pre-process the dataset output directory.
        date_str = datetime.now().strftime("%Y-%m-%d")
        if not re.match("^[0-9]{4}-(0[1-9]|1[0-2])-(0[1-9]|[1-2][0-9]|3[0-1])", dataset_outdir):
            dataset_outdir = "{}-{}".format(date_str, dataset_outdir)
        self.dataset_outdir = dataset_outdir
        self.dataset_plot_outdir = os.path.join(self.dataset_outdir, "plot")

        for dir in [self.dataset_outdir, self.dataset_plot_outdir]:
            if not os.path.exists(dir):
                os.makedirs(dir)

        ########################################################################

        self.bryois_ct_dict = {
            "Astrocytes": "AST",
            "Endothelial.cells": "END",
            "Excitatory.neurons": "EX",
            "Inhibitory.neurons": "IN",
            "Microglia": "MIC",
            "Oligodendrocytes": "OLI",
            "OPCs...COPs": "OPC",
            "Pericytes": "PER"
        }

        self.palette = {
            "AST": "#D55E00",
            "END": "#CC79A7",
            "EX": "#0072B2",
            "IN": "#56B4E9",
            "MIC": "#E69F00",
            "OLI": "#009E73",
            "OPC": "#F0E442",
            "PER": "#808080"
        }

        self.shared_xlim = None
        self.shared_ylim = None

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit")
        parser.add_argument("--work_dir",
                            type=str,
                            required=True,
                            default=None,
                            help="The working directory.")
        parser.add_argument("--dataset_outdir",
                            type=str,
                            required=True,
                            default=None,
                            help="The name of the output directory where you "
                                 "would like all outputs/results saved.")
        parser.add_argument("--wg3_folder",
                            type=str,
                            default=None,
                            help="The path to the outputs from WG3")
        parser.add_argument("--annotation_level",
                            type=str,
                            choices=["L1", "L2"],
                            default="L1",
                            help="The annotation level to use. "
                                 "Default: 'L1'.")
        parser.add_argument("--exclude_ct",
                            nargs="+",
                            type=str,
                            default=[],
                            help="The cell types results to exclude.")
        parser.add_argument("--bryois_folder",
                            type=str,
                            default=None,
                            help="The path to the Bryois et al. 2022 data")
        parser.add_argument("-e",
                            "--extensions",
                            nargs="+",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")
        parser.add_argument("--force",
                            action='store_true',
                            help="Force a rerun of all files.")
        parser.add_argument("--verbose",
                            action='store_true',
                            help="Print all info.")

        # Required external scripts.
        parser.add_argument("--qvalues",
                            type=str,
                            default="qvalue_truncp.R",
                            help="The path to the qvalues script")
        parser.add_argument("--rb",
                            type=str,
                            default="Rb.R",
                            help="The path to the Rb script")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data")
        plot_data = {}
        for filename, cell_type in self.bryois_ct_dict.items():
            if cell_type in self.exclude_ct:
                continue

            print("  Working on '{}'".format(cell_type))

            ct_workdir = os.path.join(self.dataset_outdir, cell_type)
            if not os.path.exists(ct_workdir):
                os.makedirs(ct_workdir)

            bryois_df_path = os.path.join(ct_workdir, "bryois2022.txt.gz")
            scmetabrain_df_path = os.path.join(ct_workdir, "scmetabrain.txt.gz")
            merged_df_path = os.path.join(ct_workdir, "merged.txt.gz")

            if os.path.exists(merged_df_path) and not self.force:
                df = self.load_file(merged_df_path, header=0, index_col=0)
            else:
                if os.path.exists(bryois_df_path):
                # if os.path.exists(bryois_df_path) and not self.force:
                    bryois_df = self.load_file(bryois_df_path, header=0, index_col=0)
                else:
                    bryois_df = self.load_bryois_data(outpath=bryois_df_path, filename=filename)

                bryois_df = bryois_df.loc[~bryois_df["pos_hg38"].isna(), :]
                bryois_df["pos_hg38"] = bryois_df["pos_hg38"].astype(int)
                bryois_df.index = bryois_df["ENSG"] + "_" + bryois_df["pos_hg38"].astype(str) + "_" + bryois_df["alleles"]

                if os.path.exists(scmetabrain_df_path):
                # if os.path.exists(scmetabrain_df_path) and not self.force:
                    scmetabrain_df = self.load_file(scmetabrain_df_path, header=0, index_col=0)
                else:
                    scmetabrain_df = self.load_wg3_data(cell_type=cell_type,
                                                        outpath=scmetabrain_df_path,
                                                        ensg_hits=set(bryois_df["ENSG"].values),
                                                        index_hits=set(bryois_df.index.values))
                    if scmetabrain_df is None:
                        print("  Warning, no overlapping eQTLs found.")
                        continue

                scmetabrain_df.index = scmetabrain_df["ENSG"] + "_" + scmetabrain_df["snp_position"].astype(int).astype(str) + "_" + scmetabrain_df["alleles"]

                df = self.merge_data(bryois_df=bryois_df, scmetabrain_df=scmetabrain_df, outpath=merged_df_path)

            plot_data[cell_type] = df

        print("\nVisualizing comparison")
        replication_stats_df = self.visualise_data(plot_data=plot_data)

        print("\nReplication stats:")
        for label in replication_stats_df["label"].unique():
            print("\t{}".format(label))
            stats_df = replication_stats_df.loc[replication_stats_df["label"] == label, :]
            stats_df_mean = stats_df[["variable", "value"]].groupby("variable").mean()
            for index, row in stats_df_mean.iterrows():
                print("\t  {}: {:.2f}".format(index, row["value"]))

            stats_df_sum = stats_df[["variable", "value"]].groupby("variable").sum()
            print("\t  Overall concordance: {:,}/{:,} [{:.2f}%]".format(
                stats_df_sum.loc["N concordant", "value"],
                stats_df_sum.loc["N", "value"],
                (100 / stats_df_sum.loc["N", "value"]) * stats_df_sum.loc["N concordant", "value"]))
            print("")

        self.save_file(df=replication_stats_df,
                       outpath=os.path.join(self.dataset_outdir,
                                            "replication_stats.txt.gz"))

    def load_bryois_data(self, filename, outpath):
        print("  Loading Bryois et al. 2022 data")

        dfs = []
        for chr in range(1, 23):
            print("   Loading chromosome {}".format(chr))
            ct_chr_path = os.path.join(self.bryois_folder, "{}.{}.gz".format(filename, chr))
            if not os.path.exists(ct_chr_path):
                continue

            ct_chr_df = self.load_file(ct_chr_path, sep=" ", header=None, index_col=None)
            ct_chr_df.columns = ["Gene_id", "SNP_id", "Distance to TSS", "Nominal p-value", "Beta"]
            ct_chr_df = ct_chr_df.loc[ct_chr_df.groupby('Gene_id')["Nominal p-value"].idxmin(), :]
            dfs.append(ct_chr_df)

            del ct_chr_df

        if len(dfs) == 0:
            print("  Warning, no data loaded for this cell type")
            return None

        df = pd.concat(dfs, axis=0)

        if df.shape[0] == 0:
            print("  Warning, no data loaded for this cell type")
            return None

        gene_id_df = df["Gene_id"].str.split("_", n=None, expand=True)
        gene_id_df.columns = ["HGNC", "ENSG"]
        df = pd.concat([gene_id_df, df], axis=1)
        del gene_id_df

        print("\tLoading SNP position file")
        snp_pos_df = self.load_file(inpath=os.path.join(self.bryois_folder, "snp_pos.txt"), header=0, index_col=None)
        snp_pos_df['allele1'] = np.minimum(snp_pos_df['effect_allele'], snp_pos_df['other_allele'])
        snp_pos_df['allele2'] = np.maximum(snp_pos_df['effect_allele'], snp_pos_df['other_allele'])
        snp_pos_df['alleles'] = snp_pos_df['allele1'] + snp_pos_df['allele2']
        snp_pos_df.drop(['allele1', 'allele2'], axis=1, inplace=True)
        snp_pos_df.columns = ["SNP_id" if col == "SNP" else col for col in snp_pos_df.columns]

        print("  Merging with SNP position file")
        df = df.merge(snp_pos_df, on="SNP_id", how="left")

        print("  Saving file")
        self.save_file(df=df, outpath=outpath)
        # self.save_file(df=df, outpath=outpath.replace(".txt.gz", ".xlsx"))

        return df

    def load_wg3_data(self, cell_type, outpath, ensg_hits, index_hits):
        print("  Loading scMetaBrain data")

        results_file_prefix = "qtl_results_"

        h5_files = glob.glob(os.path.join(self.wg3_folder, "output", self.annotation_level, cell_type, "qtl", "{}*.h5".format(results_file_prefix)))
        h5_files = [(int(os.path.basename(h5_file).split("_")[2]), h5_file) for h5_file in h5_files]
        h5_files.sort(key=lambda x: x[0])
        if len(h5_files) == 0:
            print("  Warning, no h5 files found for this cell type")
            return None

        dfs = []
        prev_chr = None
        for chr, h5_filepath in h5_files:
            if prev_chr is None or chr != prev_chr:
                print("   Loading chromosome {}".format(os.path.basename(h5_filepath).split("_")[2]))
            h5_df = self.load_h5_file(filepath=h5_filepath,
                                      cell_type=cell_type,
                                      results_file_prefix=results_file_prefix,
                                      ensg_hits=ensg_hits,
                                      index_hits=index_hits)
            if h5_df is None:
                continue

            dfs.append(h5_df)
            prev_chr = chr
            del h5_df

        if len(dfs) == 0:
            print("  Warning, no data loaded for this cell type")
            return None

        df = pd.concat(dfs, axis=0)

        if df.shape[0] == 0:
            print("  Warning, no data loaded for this cell type")
            return None

        print("  Saving file")
        self.save_file(df=df, outpath=outpath)
        # self.save_file(df=df, outpath=outpath.replace(".txt.gz", ".xlsx"))

        return df

    def load_h5_file(self, filepath, cell_type, results_file_prefix, ensg_hits, index_hits):
        analysis_subset = os.path.basename(filepath).replace(results_file_prefix, "").replace(".h5", "")
        snp_metadata_file = os.path.join(self.wg3_folder, "output", self.annotation_level, cell_type, "qtl", "snp_metadata_{}.txt".format(analysis_subset))
        feature_metadata_file = os.path.join(self.wg3_folder, "output", self.annotation_level, cell_type, "qtl", "feature_metadata_{}.txt".format(analysis_subset))

        if os.path.exists(feature_metadata_file):
            ffea_df = self.load_file(feature_metadata_file, header=0, index_col=None)
        elif os.path.exists(feature_metadata_file + ".gz"):
            ffea_df = self.load_file(feature_metadata_file + ".gz", header=0, index_col=None)
        else:
            print("  Warning, skipping: '{}' missing feature metadata file.".format(analysis_subset))
            return None

        if os.path.exists(snp_metadata_file):
            fsnp_df = self.load_file(snp_metadata_file, header=0, index_col=None)
        elif os.path.exists(snp_metadata_file + ".gz"):
            fsnp_df = self.load_file(snp_metadata_file + ".gz", header=0, index_col=None)
        else:
            print("  Warning, skipping: '{}' missing SNP metadata file.".format(analysis_subset))
            return None

        ffea_df = ffea_df.rename(index=str,
                                 columns={"chromosome": "feature_chromosome",
                                          "start": "feature_start",
                                          "end": "feature_end"})
        fsnp_df = fsnp_df.rename(index=str,
                                 columns={"chromosome": "snp_chromosome",
                                          "position": "snp_position"})

        ########################################################################

        snp_id_df = fsnp_df["snp_id"].str.split(":", n=None, expand=True)
        snp_id_df.columns = ["snp_chromosome", "snp_position", "alleleA", "alleleB"]
        snp_id_df.drop(["snp_chromosome", "snp_position"], axis=1, inplace=True)
        fsnp_df = pd.concat([fsnp_df, snp_id_df], axis=1)
        fsnp_df['allele1'] = np.minimum(fsnp_df['alleleA'], fsnp_df['alleleB'])
        fsnp_df['allele2'] = np.maximum(fsnp_df['alleleA'], fsnp_df['alleleB'])
        fsnp_df['alleles'] = fsnp_df['allele1'] + fsnp_df['allele2']
        fsnp_df.drop(['alleleA', 'alleleB', 'allele1', 'allele2'], axis=1, inplace=True)

        ########################################################################

        frez = h5py.File(filepath, 'r')
        frez_keys = [k.replace('_i_', '') for k in list(frez.keys())]

        # Filter on Bryois ENSG hits.
        hgnc_to_ensg_dict = dict(zip(ffea_df["feature_id"], ffea_df["ENSG"]))
        frez_keys = [frez_key for frez_key in frez_keys if frez_key in hgnc_to_ensg_dict and hgnc_to_ensg_dict[frez_key] in ensg_hits]
        if len(frez_keys) == 0:
            return None

        dfs = []
        for frez_key in frez_keys:
            frez_df = pd.DataFrame(np.array(frez[frez_key]))
            frez_df['feature_id'] = frez_key
            dfs.append(frez_df)
        df = pd.concat(dfs, axis=0)
        df['snp_id'] = df['snp_id'].astype(str)

        if self.verbose:
            print("\tLoaded h5 file: {} with shape: {}".format(os.path.basename(filepath), df.shape))

        df = pd.merge(df, ffea_df, on='feature_id', how='left')

        ########################################################################

        if(len(glob.glob(os.path.join(self.wg3_folder, "output", self.annotation_level, cell_type, "qtl", "snp_qc_metrics_naContaining_feature_*.txt"))) > 0):
            print("  Error, code not implemented yet")
            exit()
            # tmp_df = pd.DataFrame(columns=df.columns)
            # for key in frez_keys:
            #     qc_metrics_inpath = os.path.join(self.wg3_folder, "output", self.annotation_level, cell_type, "qtl", "snp_qc_metrics_naContaining_feature_{}.txt".format(key))
            #     if os.path.isfile(qc_metrics_inpath):
            #         fsnp_rel = self.load_file(qc_metrics_inpath, header=0, index_col=None)
            #         tmp_t = df.loc[df["feature_id"] == key]
            #         fsnp_t = fsnp_df.loc[:,["snp_id", "snp_chromosome", "snp_position", "assessed_allele"]]
            #         fsnp_t = pd.merge(fsnp_t, fsnp_rel, on='snp_id', how='right')
            #         tmp_t = pd.merge(tmp_t, fsnp_t, on='snp_id', how='left')
            #         tmp_df = tmp_df.append(tmp_t, sort=False)
            #     else:
            #         tmp_t = df.loc[df["feature_id"] == key]
            #         tmp_t = pd.merge(tmp_t, fsnp_df, on='snp_id', how='left')
            #         tmp_df = tmp_df.append(tmp_t,sort=False)
            #     data[key]=np.zeros(len(np.unique(list(frez_keys))),dtype='object')+np.nan
            # df = tmp_df
            # del tmp_df
        else:
            df = pd.merge(df, fsnp_df, on='snp_id', how='left')

        df.index = df["ENSG"] + "_" + df["snp_position"].astype(str) + "_" + df["alleles"]
        df = df.loc[[index for index in df.index if index in index_hits], :].copy()

        df['empirical_feature_p_value'] = df['empirical_feature_p_value'].astype(float)
        df['p_value'] = df['p_value'].astype(float)

        del ffea_df, fsnp_df

        return df

    def merge_data(self, bryois_df, scmetabrain_df, outpath):
        # Make columns unique.
        bryois_df.columns = ["bryois_{}".format(col) for col in bryois_df.columns]
        scmetabrain_df.columns = ["scmetabrain_{}".format(col) for col in scmetabrain_df.columns]

        print("  Merging discovery and replication eQTLs")
        df = bryois_df.merge(scmetabrain_df, left_index=True, right_index=True)
        print("\t{} overlapping entries".format(df.shape[0]))

        df["flip"] = df["scmetabrain_assessed_allele"] != df["bryois_effect_allele"]
        df["scmetabrain_beta"] = df["scmetabrain_beta"] * df["flip"].map({True: -1, False: 1})

        print("\tSaving file")
        self.save_file(df=df, outpath=outpath)
        # self.save_file(df=df, outpath=outpath.replace(".txt.gz", ".xlsx"))

        return df

    def load_file(self, inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)

        if self.verbose:
            print("\tLoaded dataframe: {} "
                  "with shape: {}".format(os.path.basename(inpath),
                                          df.shape))
        return df

    def save_file(self, df, outpath, header=True, index=True, sep="\t", na_rep="NA",
                  sheet_name="Sheet1"):
        if outpath.endswith('xlsx'):
            df.to_excel(outpath,
                        sheet_name=sheet_name,
                        na_rep=na_rep,
                        header=header,
                        index=index)
        else:
            compression = 'infer'
            if outpath.endswith('.gz'):
                compression = 'gzip'

            df.to_csv(outpath, sep=sep, index=index, header=header,
                      compression=compression)

        if self.verbose:
            print("\tSaved dataframe: {} "
                  "with shape: {}".format(os.path.basename(outpath),
                                          df.shape))

    def visualise_data(self, plot_data):
        cell_types = list(plot_data.keys())
        cell_types.sort()

        nrows = 3
        ncols = len(cell_types)

        self.shared_ylim = {i: (0, 1) for i in range(nrows)}
        self.shared_xlim = {i: (0, 1) for i in range(ncols)}

        replication_stats = []

        sns.set(rc={'figure.figsize': (ncols * 8, nrows * 6)})
        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharex='col',
                                 sharey='row')

        for col_index, cell_type in enumerate(cell_types):
            print("\tWorking on '{}'".format(cell_type))

            # Select the required columns.
            df = plot_data[cell_type]
            df["bryois_n_samples"] = self.bryois_n
            plot_df = df.loc[:, ["bryois_HGNC",
                                 "scmetabrain_n_samples",
                                 "scmetabrain_maf",
                                 "scmetabrain_empirical_feature_p_value",
                                 "scmetabrain_beta",
                                 "scmetabrain_beta_se",
                                 "bryois_n_samples",
                                 "bryois_Nominal p-value",
                                 "bryois_Beta"
                                 ]].copy()
            plot_df.columns = ["Gene symbol",
                               "scMetaBrain N",
                               "scMetaBrain MAF",
                               "scMetaBrain pvalue",
                               "scMetaBrain beta",
                               "scMetaBrain beta se",
                               "Bryois N",
                               "Bryois pvalue",
                               "Bryois beta"]
            plot_df = plot_df.loc[~plot_df["Bryois pvalue"].isna(), :]
            plot_df.sort_values(by="scMetaBrain pvalue", inplace=True)
            print(plot_df)

            # Calculate the replication standard error.
            self.pvalue_to_zscore(df=plot_df,
                                  beta_col="Bryois beta",
                                  p_col="Bryois pvalue",
                                  prefix="Bryois ")
            self.zscore_to_beta(df=plot_df,
                                z_col="Bryois z-score",
                                maf_col="scMetaBrain MAF",
                                n_col="Bryois N",
                                prefix="Bryois zscore-to-")

            # Convert the interaction beta to log scale.
            plot_df["scMetaBrain log beta"] = [np.log(abs(beta) + 1) * np.sign(beta) for beta in plot_df["scMetaBrain beta"]]
            plot_df["Bryois log beta"] = [np.log(abs(beta) + 1) * np.sign(beta) for beta in plot_df["Bryois beta"]]
            print(plot_df)

            include_ylabel = False
            if col_index == 0:
                include_ylabel = True

            if col_index == 0:
                for row_index, panel in enumerate(["A", "B", "C"]):
                    axes[row_index, col_index].annotate(
                        panel,
                        xy=(-0.3, 0.9),
                        xycoords=axes[row_index, col_index].transAxes,
                        color="#000000",
                        fontsize=40
                    )

            print("\tPlotting row 1.")
            xlim, ylim, stats1 = self.scatterplot(
                df=plot_df,
                fig=fig,
                ax=axes[0, col_index],
                x="Bryois log beta",
                y="scMetaBrain log beta",
                xlabel="",
                ylabel="scMetaBrain log beta",
                title=cell_type,
                color=self.palette[cell_type],
                include_ylabel=include_ylabel
            )
            self.update_limits(xlim, ylim, 0, col_index)

            print("\tPlotting row 2.")
            xlim, ylim, stats2 = self.scatterplot(
                df=plot_df.loc[plot_df["Bryois pvalue"] <= 0.05, :],
                fig=fig,
                ax=axes[1, col_index],
                x="Bryois log beta",
                y="scMetaBrain log beta",
                xlabel="",
                ylabel="scMetaBrain log beta",
                title="",
                color=self.palette[cell_type],
                include_ylabel=include_ylabel,
                pi1_column="scMetaBrain pvalue",
                rb_columns=[("Bryois zscore-to-beta", "Bryois zscore-to-se"), ("scMetaBrain beta", "scMetaBrain beta se")]
            )
            self.update_limits(xlim, ylim, 1, col_index)

            print("\tPlotting row 3.")
            xlim, ylim, stats3 = self.scatterplot(
                df=plot_df.loc[(plot_df["Bryois pvalue"] <= 0.05) & (plot_df["scMetaBrain pvalue"] <= 0.05), :],
                fig=fig,
                ax=axes[2, col_index],
                x="Bryois log beta",
                y="scMetaBrain log beta",
                label="Gene symbol",
                xlabel="Bryois log beta",
                ylabel="scMetaBrain log beta",
                title="",
                color=self.palette[cell_type],
                include_ylabel=include_ylabel
            )
            self.update_limits(xlim, ylim, 2, col_index)
            print("")

            for stats, label in zip([stats1, stats2, stats3], ["all", "discovery significant", "both significant"]):
                stats_m = stats.melt()
                stats_m["label"] = label
                stats_m["cell type"] = cell_type
                replication_stats.append(stats_m)

        for (m, n), ax in np.ndenumerate(axes):
            (xmin, xmax) = self.shared_xlim[n]
            (ymin, ymax) = self.shared_ylim[m]

            xmargin = (xmax - xmin) * 0.05
            ymargin = (ymax - ymin) * 0.05

            ax.set_xlim(xmin - xmargin - 1, xmax + xmargin)
            ax.set_ylim(ymin - ymargin, ymax + ymargin)

        # Add the main title.
        fig.suptitle("scMetaBrain replication in Bryois et al. 2022",
                     fontsize=40,
                     color="#000000",
                     weight='bold')

        for extension in self.extensions:
            fig.savefig(os.path.join(self.dataset_plot_outdir, "replication_plot.{}".format(extension)))
        plt.close()

        # Construct the replication stats data frame.
        replication_stats_df = pd.concat(replication_stats, axis=0)
        replication_stats_df.dropna(inplace=True)

        return replication_stats_df

    @staticmethod
    def pvalue_to_zscore(df, beta_col, p_col, prefix=""):
        p_values = df[p_col].to_numpy()
        zscores = stats.norm.ppf(p_values / 2)
        mask = np.ones_like(p_values)
        mask[df[beta_col] > 0] = -1
        df["{}z-score".format(prefix)] = zscores * mask
        df.loc[df[p_col] == 1, "{}z-score".format(prefix)] = 0
        df.loc[df[p_col] == 0, "{}z-score".format(prefix)] = -40.

    @staticmethod
    def zscore_to_beta(df, z_col, maf_col, n_col, prefix=""):
        chi = df[z_col] * df[z_col]
        a = 2 * df[maf_col] * (1 - df[maf_col]) * (df[n_col] + chi)
        df["{}beta".format(prefix)] = df[z_col] / a ** (1/2)
        df["{}se".format(prefix)] = 1 / a ** (1/2)

    def scatterplot(self, df, fig, ax, x="x", y="y", facecolors=None,
                    label=None, max_labels=15, xlabel="", ylabel="", title="",
                    color="#000000", ci=95, include_ylabel=True,
                    pi1_column=None, rb_columns=None):
        sns.despine(fig=fig, ax=ax)

        if not include_ylabel:
            ylabel = ""

        if facecolors is None:
            facecolors = "#808080"
        else:
            facecolors = df[facecolors]

        n = df.shape[0]
        concordance = np.nan
        n_concordant = np.nan
        coef = np.nan
        pi1 = np.nan
        rb = np.nan

        if n > 0:
            lower_quadrant = df.loc[(df[x] < 0) & (df[y] < 0), :]
            upper_quadrant = df.loc[(df[x] > 0) & (df[y] > 0), :]
            n_concordant = lower_quadrant.shape[0] + upper_quadrant.shape[0]
            concordance = (100 / n) * n_concordant

            if n > 1:
                coef, p = stats.pearsonr(df[x], df[y])

                if pi1_column is not None:
                    pi1 = self.calculate_p1(p=df[pi1_column])

                if rb_columns is not None:
                    rb_est = self.calculate_rb(
                        b1=df[rb_columns[0][0]],
                        se1=df[rb_columns[0][1]],
                        b2=df[rb_columns[1][0]],
                        se2=df[rb_columns[1][1]],
                        )
                    rb = rb_est[0]

            sns.regplot(x=x, y=y, data=df, ci=ci,
                        scatter_kws={'facecolors': facecolors,
                                     'edgecolors': "#808080"},
                        line_kws={"color": color},
                        ax=ax
                        )

            if label is not None:
                texts = []
                for i, (_, point) in enumerate(df.iterrows()):
                    if i > max_labels:
                        continue
                    texts.append(ax.text(point[x],
                                         point[y],
                                         str(point[label]),
                                         color=color))

                adjust_text(texts,
                            ax=ax,
                            only_move={'points': 'x',
                                       'text': 'xy',
                                       'objects': 'x'},
                            autoalign='x',
                            expand_text=(1., 1.),
                            expand_points=(1., 1.),
                            arrowprops=dict(arrowstyle='-', color='#808080'))

        ax.axhline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)
        ax.axvline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)

        y_pos = 0.9
        if n > 0:
            ax.annotate(
                'N = {:,}'.format(n),
                xy=(0.03, 0.9),
                xycoords=ax.transAxes,
                color=color,
                fontsize=14,
                fontweight='bold'
            )
            y_pos -= 0.05

        if not np.isnan(coef):
            ax.annotate(
                'r = {:.2f}'.format(coef),
                xy=(0.03, y_pos),
                xycoords=ax.transAxes,
                color=color,
                fontsize=14,
                fontweight='bold'
            )
            y_pos -= 0.05

        if not np.isnan(concordance):
            ax.annotate(
                'concordance = {:.0f}%'.format(concordance),
                xy=(0.03, y_pos),
                xycoords=ax.transAxes,
                color=color,
                fontsize=14,
                fontweight='bold'
            )
            y_pos -= 0.05

        if not np.isnan(pi1):
            ax.annotate(
                '\u03C01 = {:.2f}'.format(pi1),
                xy=(0.03, y_pos),
                xycoords=ax.transAxes,
                color=color,
                fontsize=14,
                fontweight='bold'
            )
            y_pos -= 0.05

        if not np.isnan(rb):
            ax.annotate(
                'Rb = {:.2f}'.format(rb),
                xy=(0.03, y_pos),
                xycoords=ax.transAxes,
                color=color,
                fontsize=14,
                fontweight='bold'
            )

        ax.set_title(title,
                     fontsize=22,
                     color=color,
                     weight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        stats_df = pd.DataFrame([[n, n_concordant, concordance, coef, pi1, rb]],
                                columns=["N", "N concordant", "concordance", "pearsonr", "pi1", "Rb"],
                                index=[0])

        return (df[x].min(), df[x].max()), (df[y].min(), df[y].max()), stats_df

    def update_limits(self, xlim, ylim, row, col):
        row_ylim = self.shared_ylim[row]
        if ylim[0] < row_ylim[0]:
            row_ylim = (ylim[0], row_ylim[1])
        if ylim[1] > row_ylim[1]:
            row_ylim = (row_ylim[0], ylim[1])
        self.shared_ylim[row] = row_ylim

        col_xlim = self.shared_xlim[col]
        if xlim[0] < col_xlim[0]:
            col_xlim = (xlim[0], col_xlim[1])
        if xlim[1] > col_xlim[1]:
            col_xlim = (col_xlim[0], xlim[1])
        self.shared_xlim[col] = col_xlim

    @staticmethod
    def calculate_p1(p):
        robjects.r("source('qvalue_truncp.R')")
        p = robjects.FloatVector(p)
        qvalue_truncp = robjects.globalenv['qvalue_truncp']
        pi0 = qvalue_truncp(p)[0]
        return 1 - np.array(pi0)

    @staticmethod
    def calculate_rb(b1, se1, b2, se2, theta=0):
        robjects.r("source('Rb.R')")
        b1 = robjects.FloatVector(b1)
        se1 = robjects.FloatVector(se1)
        b2 = robjects.FloatVector(b2)
        se2 = robjects.FloatVector(se2)
        calcu_cor_true = robjects.globalenv['calcu_cor_true']
        rb = calcu_cor_true(b1, se1, b2, se2, theta)
        return np.array(rb)[0]

    def print_arguments(self):
        print("Arguments:")
        print("  > Working directory:           {}".format(self.work_dir))
        print("  > Dataset output directory:    {}".format(self.dataset_outdir))
        print("  > Workgroup 3 folder:          {}".format(self.wg3_folder))
        print("  > Annotation level:            {}".format(self.annotation_level))
        print("  > Exclude cell type:           {}".format(", ".join(self.exclude_ct)))
        print("  > Bryois data folder:          {}".format(self.bryois_folder))
        print("  > Plot extensions:             {}".format(", ".join(self.extensions)))
        print("  > Force:                       {}".format(self.force))
        print("  > Verbose:                     {}".format(self.verbose))
        print("  > Qvalues script:              {}".format(self.qvalues_script))
        print("  > Rb script:                   {}".format(self.rb_script))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
