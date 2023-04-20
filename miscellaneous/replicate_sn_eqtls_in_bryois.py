#!/usr/bin/env python3

"""
File:         replicate_sn_eqtls_in_bryois.py
Created:      2023/04/20
Last Changed:
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
import gzip
import glob
import re
import os

# Third party imports.
import numpy as np
import pandas as pd
from statsmodels.stats import multitest
# import rpy2.robjects as robjects
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adjustText import adjust_text

# Local application imports.

"""
Syntax:
./replicate_sn_eqtls_in_bryois.py \
    --work_dir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-20-ReplicateInBryois \
    --dataset_outdir 2023-04-12-Mathys2019 \
    --wg3_folder /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-11-WorkGroup3eQTLAndDEA/2023-04-12-Mathys2019 \
    --exclude_ct END \
    --bryois_folder /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-20-ReplicateInBryois/Bryois2022
    
"""

# Metadata
__program__ = "Replicate sneQTLs in Bryois et al. 2022"
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

        self.dist_to_repl_ct_dict = {
            "AST": "Astrocytes",
            "END": "Endothelial.cells",
            "EX": "Excitatory.neurons",
            "IN": "Inhibitory.neurons",
            "MIC": "Microglia",
            "OLI": "Oligodendrocytes",
            "OPC": "OPCs...COPs",
            "PER": "Pericytes"
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
                            default=None,
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

        disc_df_path = os.path.join(self.dataset_outdir, "discovery.txt.gz")
        repl_df_path = os.path.join(self.dataset_outdir, "replication.txt.gz")
        df_path = os.path.join(self.dataset_outdir, "merged.txt.gz")

        print("Loading data")
        if os.path.exists(df_path) and not self.force:
            df = self.load_file(df_path, header=0, index_col=0)
        else:
            if os.path.exists(disc_df_path) and os.path.exists(repl_df_path) and not self.force:
                disc_df = self.load_file(disc_df_path, header=0, index_col=0)
                repl_df = self.load_file(repl_df_path, header=0, index_col=0)
            else:
                disc_df, ids = self.load_disc_data(disc_df_path=disc_df_path)
                repl_df = self.load_repl_data(repl_df_path=repl_df_path, ids=ids)

            print("Merging discovery and replication")
            df = self.merge_data(disc_df=disc_df, repl_df=repl_df, df_path=df_path)

        print("\nVisualizing comparison")
        replication_stats_df = self.visualise_data(df=df)

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

    def load_disc_data(self, disc_df_path):
        disc_key_columns = ["snp_id",
                            "feature_chromosome",
                            "feature_start",
                            "feature_end",
                            "ENSG",
                            "biotype",
                            "snp_chromosome",
                            "snp_position",
                            "Alleles",
                            "feature_id"]

        disc_df = None
        for cell_type_path in glob.glob(os.path.join(self.wg3_folder, "output", self.annotation_level, "*")):
            cell_type = os.path.basename(cell_type_path)
            cell_type_qtl_path = os.path.join(cell_type_path, "top_qtl_results_all.txt.gz")
            if not os.path.isdir(cell_type_path) or not os.path.exists(cell_type_qtl_path) or (self.exclude_ct is not None and cell_type in self.exclude_ct):
                continue

            print("\tProcessing '{}'".format(cell_type))

            ####################################################################

            print("\t  Loading discovery eQTLs")
            ct_disc_df = self.load_file(cell_type_qtl_path, header=0, index_col=None)
            effect_allele_data = []
            other_allele_data = []
            alleles_data = []
            mask = []
            for _, row in ct_disc_df.iterrows():
                splitted_snp_id = row["snp_id"].split(":")
                if len(splitted_snp_id) != 4:
                    effect_allele_data.append(np.nan)
                    other_allele_data.append(np.nan)
                    alleles_data.append(np.nan)
                    mask.append(False)
                    continue

                (_, _, allele1, allele2) = splitted_snp_id

                mask_value = True
                if row["assessed_allele"] == allele1:
                    effect_allele_data.append(allele1)
                    other_allele_data.append(allele2)
                elif row["assessed_allele"] == allele2:
                    effect_allele_data.append(allele2)
                    other_allele_data.append(allele1)
                else:
                    effect_allele_data.append(np.nan)
                    other_allele_data.append(np.nan)
                    mask_value = False

                alleles = [allele1, allele2]
                alleles.sort()
                alleles_data.append("/".join(alleles))

                if len(allele1) != 1 or len(allele2) != 1:
                    mask_value = False

                mask.append(mask_value)

            ct_disc_df["EA"] = effect_allele_data
            ct_disc_df["OA"] = other_allele_data
            ct_disc_df["Alleles"] = alleles_data
            ct_disc_df = ct_disc_df.loc[mask, :]
            ct_disc_df.drop(["assessed_allele"], axis=1, inplace=True)
            ct_disc_df.columns = [col if col in disc_key_columns else "{} {}".format(cell_type, col) for col in ct_disc_df.columns]

            if disc_df is None:
                disc_df = ct_disc_df
            else:
                disc_df = disc_df.merge(ct_disc_df,
                                        on=disc_key_columns,
                                        how="outer")

            del ct_disc_df

        disc_df.index = disc_df["ENSG"] + "_" + disc_df["snp_position"].astype(int).astype(str) + "_" + disc_df["Alleles"]

        col_order = disc_key_columns + [col for col in disc_df.columns if col not in disc_key_columns]
        disc_df = disc_df.loc[:, col_order]
        disc_df.columns = ["scMetaBrain " + col for col in disc_df.columns]
        print(disc_df)

        print("\tSaving data")
        self.save_file(df=disc_df, outpath=disc_df_path)

        return disc_df, set((disc_df["scMetaBrain ENSG"] + "_" + disc_df["scMetaBrain snp_position"].astype(int).astype(str)).values)


    def load_repl_data(self, repl_df_path, ids):
        repl_key_columns = ["HGNC",
                            "Ensembl",
                            "SNP",
                            "Distance to TSS"]

        print("\tLoading replication SNP position file")
        repl_snp_pos_df = self.load_file(inpath=os.path.join(self.bryois_folder, "snp_pos.txt"), header=0, index_col=None)
        repl_snp_pos_df.columns = ["SNP", "chr", "pos hg19", "pos hg38", "EA", "OA"]
        snp_info_dict = dict(zip(repl_snp_pos_df["SNP"], repl_snp_pos_df["pos hg38"]))

        repl_df = None
        for cell_type, repl_cell_type in self.dist_to_repl_ct_dict.items():
            print("\t  Loading replication eQTLs")
            repl_data = []
            for chr in range(1, 23):
                repl_file_path = os.path.join(self.bryois_folder, "{}.{}.gz".format(repl_cell_type, chr))
                if not os.path.exists(repl_file_path):
                    continue
                print("\t\t{}".format(os.path.basename(repl_file_path)))
                repl_data.extend(self.read_repl_file(repl_file_path, ids, snp_info_dict))

            if len(repl_data) == 0:
                print("\t  Warning, no overlap found in replication datasets.")
                continue

            print("\t  {:,.0f} / {:,.0f} of eQTLs found in replication datasets.".format(len(repl_data), len(ids)))

            ct_repl_df = pd.DataFrame(repl_data,
                                      columns=["HGNC",
                                               "Ensembl",
                                               "SNP",
                                               "Distance to TSS",
                                               "{} p-value".format(cell_type),
                                               "{} beta".format(cell_type)])
            if repl_df is None:
                repl_df = ct_repl_df
            else:
                repl_df = repl_df.merge(ct_repl_df,
                                        on=repl_key_columns,
                                        how="outer")

            del ct_repl_df

        ####################################################################

        if repl_df.shape[0] == 0:
            print("Error, replication data is empty.")
            exit()

        print("\tMerging with SNP position file")
        repl_df = repl_df.merge(repl_snp_pos_df, on="SNP", how="left")

        alleles_data = []
        for index, row in repl_df.iterrows():
            alleles = [row["EA"], row["OA"]]
            alleles.sort()
            alleles_data.append("/".join(alleles))
        repl_df["Alleles"] = alleles_data

        repl_df.index = repl_df["Ensembl"] + "_" + repl_df["pos hg38"].astype(int).astype(str) + "_" + repl_df["Alleles"]
        repl_df["N"] = self.bryois_n

        repl_key_columns = repl_key_columns + ["chr", "pos hg19", "pos hg38", "EA", "OA", "Alleles", "N"]
        col_order = repl_key_columns + [col for col in repl_df.columns if col not in repl_key_columns]
        repl_df = repl_df.loc[:, col_order]
        repl_df.columns = ["Bryois " + col for col in repl_df.columns]
        print(repl_df)

        self.save_file(df=repl_df, outpath=repl_df_path)

        return repl_df

    def merge_data(self, disc_df, repl_df, df_path):
        # Check what cell types we have.
        discovery_cell_types = [col.replace("scMetaBrain ", "").replace(" beta", "") for col in disc_df.columns if col.startswith("scMetaBrain ") and col.endswith(" beta")]
        replication_cell_types = [col.replace("Bryois ", "").replace(" beta", "") for col in repl_df.columns if col.startswith("Bryois ") and col.endswith(" beta")]

        print("\tMerging discovery and replication eQTLs")
        df = disc_df.merge(repl_df, left_index=True, right_index=True)

        print("\tMatching the direction of effect")
        drop_columns = []
        for cell_type in discovery_cell_types:
            df["{} flip".format(cell_type)] = df["scMetaBrain {} EA".format(cell_type)] != df["Bryois EA"]
            df.loc[:, "scMetaBrain {} beta".format(cell_type)] = df["scMetaBrain {} beta".format(cell_type)] * df["{} flip".format(cell_type)].map({True: -1, False: 1})
            drop_columns.extend(["{} flip".format(cell_type), "scMetaBrain {} EA".format(cell_type), "scMetaBrain {} OA".format(cell_type)])
        df.drop(drop_columns, axis=1, inplace=True)

        # print("\tAdding BH-FDR for the replication.")
        # overlapping_cell_types = list(set(discovery_cell_types).intersection(set(replication_cell_types)))
        # for cell_type in overlapping_cell_types:
        #     print("\t  {}".format(cell_type))
        #     df["Bryois {} BH-FDR".format(cell_type)] = np.nan
        #     discovery_mask = (df["scMetaBrain {} BH-FDR".format(cell_type)] <= 0.05).to_numpy()
        #     print("\t\tDiscovery N-ieqtls: {:,}".format(np.sum(discovery_mask)))
        #     replication_mask = (~df["Bryois {} p-value".format(cell_type)].isna()).to_numpy()
        #     mask = np.logical_and(discovery_mask, replication_mask)
        #     n_overlap = np.sum(mask)
        #     if n_overlap > 1:
        #         df.loc[mask, "Bryois {} BH-FDR".format(cell_type)] = \
        #         multitest.multipletests(df.loc[mask, "Bryois {} p-value".format(cell_type)], method='fdr_bh')[1]
        #     n_replicating = df.loc[df["Bryois {} BH-FDR".format(cell_type)] <= 0.05, :].shape[0]
        #     print("\t\tReplication N-ieqtls: {:,} / {:,} [{:.2f}%]".format(n_replicating, n_overlap, (100 / n_overlap) * n_replicating))

        # print("\tReordering columns")
        # TODO

        print("\tSaving output")
        self.save_file(df=df,
                       outpath=df_path,
                       index=False)
        # self.save_file(df=df,
        #                outpath=df_path.replace(".txt.gz", "xlsx"),
        #                index=False,
        #                sheet_name="Bryois et al. 2022 replication")

        return df

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def read_repl_file(filepath, ids, snp_info_dict):
        lines = []
        with gzip.open(filepath, 'rt') as f:
            for line in f:
                (gene_id, snp_id, dist_to_tss, nom_pval, beta) = line.strip("\n").split(" ")
                (hgnc_name, ensembl_name) = gene_id.split("_")
                if snp_id not in snp_info_dict.keys():
                    continue
                line_id = ensembl_name + "_" + str(int(snp_info_dict[snp_id]))
                if line_id not in ids:
                    continue

                lines.append([hgnc_name, ensembl_name, snp_id, dist_to_tss, nom_pval, beta])
        f.close()

        return lines

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t", na_rep="NA",
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
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    def visualise_data(self, df):
        print(df)
        print(df.columns.tolist())

        discovery_cell_types = [col.replace("scMetaBrain ", "").replace(" beta", "") for col in df.columns if col.startswith("scMetaBrain ") and col.endswith(" beta")]
        replication_cell_types = [col.replace("Bryois ", "").replace(" beta", "") for col in df.columns if col.startswith("Bryois ") and col.endswith(" beta")]
        cell_types = list(set(discovery_cell_types).intersection(set(replication_cell_types)))
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

        for col_index, ct in enumerate(cell_types):
            print("\tWorking on '{}'".format(ct))

            # Select the required columns.
            plot_df = df.loc[:, ["Bryois HGNC",
                                 "scMetaBrain {} n_samples".format(ct),
                                 "scMetaBrain {} maf".format(ct),
                                 "scMetaBrain {} empirical_feature_p_value".format(ct),
                                 "scMetaBrain {} beta".format(ct),
                                 "scMetaBrain {} beta_se".format(ct),
                                 "Bryois N",
                                 "Bryois {} p-value".format(ct),
                                 "Bryois {} beta".format(ct)
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
            plot_df = plot_df.loc[(~plot_df["scMetaBrain pvalue"].isna()) & (~plot_df["Bryois pvalue"].isna()), :]
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
                x="scMetaBrain log beta",
                y="Bryois log beta",
                xlabel="",
                ylabel="Bryois log eQTL beta",
                title=ct,
                color=self.palette[ct],
                include_ylabel=include_ylabel
            )
            self.update_limits(xlim, ylim, 0, col_index)

            print("\tPlotting row 2.")
            xlim, ylim, stats2 = self.scatterplot(
                df=plot_df.loc[plot_df["scMetaBrain pvalue"] <= 0.05, :],
                fig=fig,
                ax=axes[1, col_index],
                x="scMetaBrain log beta",
                y="Bryois log beta",
                xlabel="",
                ylabel="Bryois log beta",
                title="",
                color=self.palette[ct],
                include_ylabel=include_ylabel,
                #pi1_column="Bryois pvalue",
                #rb_columns=[("scMetaBrain beta", "scMetaBrain beta se"), ("Bryois zscore-to-beta", "Bryois zscore-to-se")]
            )
            self.update_limits(xlim, ylim, 1, col_index)

            print("\tPlotting row 3.")
            xlim, ylim, stats3 = self.scatterplot(
                df=plot_df.loc[plot_df["Bryois pvalue"] <= 0.05, :],
                fig=fig,
                ax=axes[2, col_index],
                x="scMetaBrain log beta",
                y="Bryois log beta",
                label="Gene symbol",
                xlabel="scMetaBrain log beta",
                ylabel="Bryois log beta",
                title="",
                color=self.palette[ct],
                include_ylabel=include_ylabel
            )
            self.update_limits(xlim, ylim, 2, col_index)
            print("")

            for stats, label in zip([stats1, stats2, stats3], ["all", "discovery significant", "both significant"]):
                stats_m = stats.melt()
                stats_m["label"] = label
                stats_m["cell type"] = ct
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
        print("  > Qvalues script:              {}".format(self.qvalues_script))
        print("  > Rb script:                   {}".format(self.rb_script))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
