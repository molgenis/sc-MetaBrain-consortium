#!/usr/bin/env python3

"""
File:         bulk_metabrain_replication.py
Created:      2022/02/10
Last Changed: 2022/04/24
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
import os
import re

# Third party imports.
import numpy as np
import pandas as pd
import h5py
from statsmodels.stats import multitest
import rpy2.robjects as robjects
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from adjustText import adjust_text

# Local application imports.

"""
Syntax:
./bulk_metabrain_replication.py \
    --work_dir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-20-ReplicateInBryois \
    --wg3_folder /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-11-WorkGroup3eQTLAndDEA/2023-04-12-Mathys2019 \
    --dataset_outdir 2023-04-24-Mathys2019-vs-BulkMetaBrain
    
"""

# Metadata
__program__ = "Bulk MetaBrain Replication"
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
        arguments = self.create_argument_parser()
        self.work_dir = getattr(arguments, 'work_dir')
        dataset_outdir = getattr(arguments, 'dataset_outdir')
        self.wg3_folder = getattr(arguments, 'wg3_folder')
        self.annotation_level = getattr(arguments, 'annotation_level')
        self.exclude_ct = getattr(arguments, 'exclude_ct')
        self.extensions = getattr(arguments, 'extensions')
        self.force = getattr(arguments, 'force')
        self.verbose = getattr(arguments, 'verbose')
        self.qvalues_script = getattr(arguments, 'qvalues')
        self.rb_script = getattr(arguments, 'rb')

        # Define the discovery data.
        self.discovery_path = os.path.join(self.work_dir, "deKlein2023", "merged_decon_results.txt.gz")

        # Pre-process the dataset output directory.
        date_str = datetime.now().strftime("%Y-%m-%d")
        if not re.match("^[0-9]{4}-(0[1-9]|1[0-2])-(0[1-9]|[1-2][0-9]|3[0-1])", dataset_outdir):
            dataset_outdir = "{}-{}".format(date_str, dataset_outdir)

        self.dataset_outdir = dataset_outdir
        self.dataset_plot_outdir = os.path.join(self.dataset_outdir, "plot")

        for dir in [dataset_outdir, self.dataset_outdir, self.dataset_plot_outdir]:
            if not os.path.exists(dir):
                os.makedirs(dir)

        self.metabrain_ct_dict = {
            "Astrocyte": "AST",
            "EndothelialCell": "END",
            "Excitatory": "EX",
            "Inhibitory": "IN",
            "Microglia": "MIC",
            "Oligodendrocyte": "OLI"
        }

        self.palette = {
            "AST": "#D55E00",
            "END": "#CC79A7",
            "EX": "#0072B2",
            "IN": "#56B4E9",
            "MIC": "#E69F00",
            "OLI": "#009E73"
        }

        self.shared_xlim = None
        self.shared_ylim = None

        matplotlib.rcParams['pdf.fonttype'] = 42

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
        parser.add_argument("-e",
                            "--extensions",
                            nargs="+",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")
        parser.add_argument("--verbose",
                            action='store_true',
                            help="Print all info.")
        parser.add_argument("--force",
                            action='store_true',
                            help="Force a rerun of all files.")


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

        print("Loading bulk MetaBrain data")
        full_bulk_metabrain_df = self.load_file(self.discovery_path, header=0, index_col=None)
        full_bulk_metabrain_df["ENSG"] = full_bulk_metabrain_df["Gene"].str.split(".", n=None, expand=True)[0]
        snp_id_df = full_bulk_metabrain_df["SNP"].str.split(":", n=None, expand=True)
        snp_id_df.columns = ["CHR", "Position", "RS_id", "alleles"]
        alleles_df = snp_id_df["alleles"].str.split("_", n=None, expand=True)
        alleles_df.columns = ["alleleA", "alleleB"]
        alleles_df['allele1'] = np.minimum(alleles_df['alleleA'], alleles_df['alleleB'])
        alleles_df['allele2'] = np.maximum(alleles_df['alleleA'], alleles_df['alleleB'])
        alleles_df['alleles'] = alleles_df['allele1'] + alleles_df['allele2']
        full_bulk_metabrain_df = pd.concat([full_bulk_metabrain_df, snp_id_df[["Position"]], alleles_df[["alleles"]]], axis=1)
        full_bulk_metabrain_df.index = full_bulk_metabrain_df["ENSG"] + "_" + full_bulk_metabrain_df["Position"].astype(str) + "_" + full_bulk_metabrain_df["alleles"]
        full_bulk_metabrain_df.drop(["alleles"], axis=1, inplace=True)
        del snp_id_df, alleles_df

        print("Loading single-nucleus MetaBrain data")
        plot_data = {}
        for bulk_cell_type, sn_cell_type in self.metabrain_ct_dict.items():
            if sn_cell_type in self.exclude_ct:
                continue

            print("  Working on '{}'".format(sn_cell_type))

            ct_workdir = os.path.join(self.dataset_outdir, sn_cell_type)
            if not os.path.exists(ct_workdir):
                os.makedirs(ct_workdir)

            sn_metabrain_df_path = os.path.join(ct_workdir, "scmetabrain.txt.gz")
            merged_df_path = os.path.join(ct_workdir, "merged.txt.gz")

            if os.path.exists(merged_df_path) and not self.force:
                df = self.load_file(merged_df_path, header=0, index_col=0)
            else:
                # Subset bulk.
                bulk_metabrain_df = full_bulk_metabrain_df.loc[:, ["Gene", "ENSG", "Gene symbol", "SNP", "Alleles", "Allele assessed", "N", "HW pval", "Minor allele", "MAF", "Overall z-score", "{} pvalue".format(bulk_cell_type), "{} beta".format(bulk_cell_type), "{} interaction beta".format(bulk_cell_type), "{} BH-FDR".format(bulk_cell_type)]].copy()
                bulk_metabrain_df.columns = ["Gene", "ENSG", "Gene symbol", "SNP", "Alleles", "Allele assessed", "N", "HW pval", "Minor allele", "MAF", "Overall z-score", "Pvalue", "Beta", "Interaction beta", "BH-FDR"]
                print(bulk_metabrain_df)

                if os.path.exists(sn_metabrain_df_path) and not self.force:
                    sn_metabrain_df = self.load_file(sn_metabrain_df_path,
                                                     header=0,
                                                     index_col=0)
                else:
                    sn_metabrain_df = self.load_wg3_data(cell_type=sn_cell_type,
                                                         outpath=sn_metabrain_df_path,
                                                         top_effect=False,
                                                         ensg_hits=set(bulk_metabrain_df["ENSG"].values),
                                                         index_hits=set(bulk_metabrain_df.index.values))
                    if sn_metabrain_df is None:
                        print("  Warning, no overlapping eQTLs found.")
                        return None

                print(sn_metabrain_df)

                df = self.merge_data(bulk_metabrain_df=bulk_metabrain_df,
                                     sn_metabrain_df=sn_metabrain_df,
                                     outpath=merged_df_path)


            plot_data[sn_cell_type] = df

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


    def load_wg3_data(self, cell_type, outpath, remove_indels=True,
                      top_effect=False, ensg_hits=None, index_hits=None):
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
            prev_chr = chr

            h5_df = self.load_h5_file(filepath=h5_filepath,
                                      cell_type=cell_type,
                                      results_file_prefix=results_file_prefix,
                                      ensg_hits=ensg_hits,
                                      index_hits=index_hits)
            if h5_df is None:
                continue

            if remove_indels:
                h5_df = h5_df.loc[h5_df["alleles"].str.len() == 2, :]

            if top_effect:
                h5_df = h5_df.loc[h5_df.groupby('feature_id')["empirical_feature_p_value"].idxmin(), :]

            dfs.append(h5_df)
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

    def load_h5_file(self, filepath, cell_type, results_file_prefix,
                     ensg_hits=None, index_hits=None):
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

        if ensg_hits is not None:
            # Filter on Bryois ENSG hits.
            hgnc_to_ensg_dict = dict(zip(ffea_df["feature_id"], ffea_df["ENSG"]))
            frez_keys = [frez_key for frez_key in frez_keys if frez_key in hgnc_to_ensg_dict and hgnc_to_ensg_dict[frez_key] in ensg_hits]

        if len(frez_keys) == 0:
            print("  Warning, skipping: '{}', no overlap in genes.".format(analysis_subset))
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
            df = pd.merge(df, fsnp_df, on='snp_id', how='inner')

        df.index = df["ENSG"] + "_" + df["snp_position"].astype(str) + "_" + df["alleles"]

        if index_hits is not None:
            df = df.loc[[index for index in df.index if index in index_hits], :]

        df['empirical_feature_p_value'] = df['empirical_feature_p_value'].astype(float)
        df['p_value'] = df['p_value'].astype(float)

        del ffea_df, fsnp_df

        return df

    def merge_data(self, bulk_metabrain_df, sn_metabrain_df, outpath):
        # Make columns unique.
        bulk_metabrain_df.columns = ["bulk_{}".format(col) for col in bulk_metabrain_df.columns]
        sn_metabrain_df.columns = ["sn_{}".format(col) for col in sn_metabrain_df.columns]

        print("  Merging discovery and replication eQTLs")
        df = bulk_metabrain_df.merge(sn_metabrain_df, left_index=True, right_index=True)
        print("\t{} overlapping entries".format(df.shape[0]))

        df["flip"] = df["bulk_Allele assessed"] != df["sn_assessed_allele"]
        df["sn_beta"] = df["sn_beta"] * df["flip"].map({True: -1, False: 1})

        df["sn_BH-FDR"] = np.nan
        repl_mask = np.logical_and((df["bulk_BH-FDR"] <= 0.05).to_numpy(), (~df["sn_empirical_feature_p_value"].isna()).to_numpy())
        n_overlap = np.sum(repl_mask)
        if n_overlap > 1:
            df.loc[repl_mask, "sn_BH-FDR"] = multitest.multipletests(df.loc[repl_mask, "sn_empirical_feature_p_value"], method='fdr_bh')[1]

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

        if ncols == 1:
            axes = axes.reshape(nrows, 1)

        for col_index, cell_type in enumerate(cell_types):
            print("\tWorking on '{}'".format(cell_type))

            # Select the required columns.
            df = plot_data[cell_type]
            plot_df = df.loc[:, ["bulk_Gene symbol",
                                 "bulk_N",
                                 "bulk_MAF",
                                 "bulk_Pvalue",
                                 "bulk_BH-FDR",
                                 "bulk_Interaction beta",
                                 "sn_n_samples",
                                 "sn_empirical_feature_p_value",
                                 "sn_BH-FDR",
                                 "sn_beta",
                                 "sn_beta_se"
                                 ]].copy()
            plot_df.columns = ["Gene symbol",
                               "bulkMetaBrain N",
                               "bulkMetaBrain MAF",
                               "bulkMetaBrain pvalue",
                               "bulkMetaBrain FDR",
                               "bulkMetaBrain beta",
                               "snMetaBrain N",
                               "snMetaBrain pvalue",
                               "snMetaBrain FDR",
                               "snMetaBrain beta",
                               "snMetaBrain beta se"]
            plot_df = plot_df.loc[(~plot_df["bulkMetaBrain pvalue"].isna()) & (~plot_df["snMetaBrain pvalue"].isna()), :]
            plot_df.sort_values(by="bulkMetaBrain pvalue", inplace=True)

            # Calculate the replication standard error.
            self.pvalue_to_zscore(df=plot_df,
                                  beta_col="bulkMetaBrain beta",
                                  p_col="bulkMetaBrain pvalue",
                                  prefix="bulkMetaBrain ")
            self.zscore_to_beta(df=plot_df,
                                z_col="bulkMetaBrain z-score",
                                maf_col="bulkMetaBrain MAF",
                                n_col="bulkMetaBrain N",
                                prefix="bulkMetaBrain zscore-to-")

            # Convert the interaction beta to log scale.
            plot_df["bulkMetaBrain log beta"] = [np.log(abs(beta) + 1) * np.sign(beta) for beta in plot_df["bulkMetaBrain beta"]]
            plot_df["snMetaBrain log beta"] = [np.log(abs(beta) + 1) * np.sign(beta) for beta in plot_df["snMetaBrain beta"]]

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
                x="bulkMetaBrain log beta",
                y="snMetaBrain log beta",
                xlabel="",
                ylabel="snMetaBrain log beta",
                title=cell_type,
                color=self.palette[cell_type],
                include_ylabel=include_ylabel
            )
            self.update_limits(xlim, ylim, 0, col_index)

            print("\tPlotting row 2.")
            xlim, ylim, stats2 = self.scatterplot(
                df=plot_df.loc[plot_df["bulkMetaBrain FDR"] <= 0.05, :],
                fig=fig,
                ax=axes[1, col_index],
                x="bulkMetaBrain log beta",
                y="snMetaBrain log beta",
                xlabel="",
                ylabel="snMetaBrain log beta",
                title="",
                color=self.palette[cell_type],
                include_ylabel=include_ylabel,
                pi1_column="snMetaBrain pvalue",
                rb_columns=[("bulkMetaBrain zscore-to-beta", "bulkMetaBrain zscore-to-se"), ("snMetaBrain beta", "snMetaBrain beta se")]
            )
            self.update_limits(xlim, ylim, 1, col_index)

            print("\tPlotting row 3.")
            xlim, ylim, stats3 = self.scatterplot(
                df=plot_df.loc[(plot_df["bulkMetaBrain FDR"] <= 0.05) & (plot_df["snMetaBrain FDR"] <= 0.05), :],
                fig=fig,
                ax=axes[2, col_index],
                x="bulkMetaBrain log beta",
                y="snMetaBrain log beta",
                label="Gene symbol",
                xlabel="bulkMetaBrain log beta",
                ylabel="snMetaBrain log beta",
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
        fig.suptitle("de Klein 2023 replication in scMetaBrain",
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

            if n > 2:
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

    def calculate_p1(self, p):
        robjects.r("source('{}')".format(self.qvalues_script))
        p = robjects.FloatVector(p)
        qvalue_truncp = robjects.globalenv['qvalue_truncp']
        pi0 = qvalue_truncp(p)[0]
        return 1 - np.array(pi0)

    def calculate_rb(self, b1, se1, b2, se2, theta=0):
        robjects.r("source('Rb.R')".format(self.rb_script))
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
        print("  > Plot extensions:             {}".format(", ".join(self.extensions)))
        print("  > Force:                       {}".format(self.force))
        print("  > Verbose:                     {}".format(self.verbose))
        print("  > Qvalues script:              {}".format(self.qvalues_script))
        print("  > Rb script:                   {}".format(self.rb_script))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
