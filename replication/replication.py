#!/usr/bin/env python3

"""
File:         replication.py
Created:      2024/02/28
Last Changed: 2024/05/01
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
import argparse
import gzip
import glob
import os
import re
import json

# Third party imports.
# import tabix
import h5py
import numpy as np
import pandas as pd
from natsort import natsort_keygen
from statsmodels.stats import multitest
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.rinterface_lib.embedded import RRuntimeError
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adjustText import adjust_text

# Local application imports.

"""
Syntax:
./replication.py -h
"""

# Metadata
__program__ = "Replication Plot"
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

CHROMOSOMES = [str(chr) for chr in range(1, 23)]
METHODS = ["LIMIX", "LIMIX_REDUCED", "mbQTL", "mbQTL_MetaBrain", "eQTLMappingPipeline", "eQTLgenPhase2", "Bryois"]

class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.discovery_method = getattr(arguments, 'discovery_method')
        self.discovery_path = getattr(arguments, 'discovery_path')
        self.discovery_name = getattr(arguments, 'discovery_name')
        self.discovery_cell_type = getattr(arguments, 'discovery_cell_type')
        self.replication_method = getattr(arguments, 'replication_method')
        self.replication_path = getattr(arguments, 'replication_path')
        self.replication_name = getattr(arguments, 'replication_name')
        self.replication_cell_type = getattr(arguments, 'replication_cell_type')
        self.gene = getattr(arguments, 'gene')
        self.snp = getattr(arguments, 'snp')
        self.effect_size = getattr(arguments, 'effect_size')
        self.pvalue = getattr(arguments, 'pvalue')
        self.alpha = getattr(arguments, 'alpha')
        self.fdr_calc_method = getattr(arguments, 'fdr_calc_method')
        self.log_modulus = getattr(arguments, 'log_modulus')
        self.palette_path = getattr(arguments, 'palette')
        self.extensions = getattr(arguments, 'extensions')
        outdir = getattr(arguments, 'outdir')
        self.force = getattr(arguments, 'force')
        self.save = getattr(arguments, 'save')
        self.qvalue_truncp_script = getattr(arguments, 'qvalue_truncp')
        self.rb_script = getattr(arguments, 'rb')

        if self.log_modulus:
            print("Warning, applying log modulus changes how to plot looks but not how the"
                  "replication statistics are calculated!")

        # Check if the names are identical.
        if self.discovery_name == self.replication_name:
            print("Warning, discovery and replication name are identical")
            self.discovery_name = self.discovery_name + "1"
            self.replication_name = self.replication_name + "2"

        # Set variables.
        if outdir is None:
            outdir = str(os.path.dirname(os.path.abspath(__file__)))
        self.data_outdir = os.path.join(outdir, "replication_data", "{}discovery_{}replication".format(self.discovery_name, self.replication_name))
        self.plot_outdir = os.path.join(outdir, "replication_plot", "{}discovery_{}replication".format(self.discovery_name, self.replication_name))
        for outdir in [self.data_outdir, self.plot_outdir]:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

        # Loading palette.
        self.palette = None
        if self.palette_path is not None:
            with open(self.palette_path) as f:
                self.palette = json.load(f)
            f.close()

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

        # Matching cell types. key = standard name, value = dataset specific name
        self.cell_type_translate = {
            "Astrocytes": "AST",
            "Endothelial.cells": "END",
            "Excitatory.neurons": "EX",
            "Inhibitory.neurons": "IN",
            "Microglia": "MIC",
            "Oligodendrocytes": "OLI",
            "OPCs...COPs": "OPC",
            "Pericytes": "PER"
        }

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
        parser.add_argument("--discovery_method",
                            type=str,
                            required=True,
                            choices=METHODS,
                            help="")
        parser.add_argument("--discovery_path",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--discovery_name",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--discovery_cell_type",
                            type=str,
                            required=False,
                            default=None,
                            help="")
        parser.add_argument("--replication_method",
                            type=str,
                            required=True,
                            choices=METHODS,
                            help="")
        parser.add_argument("--replication_path",
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--replication_name",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--replication_cell_type",
                            type=str,
                            required=False,
                            default=None,
                            help="")
        parser.add_argument("--gene",
                            type=str,
                            choices=["hgnc", "ensembl"],
                            default="hgnc",
                            help="")
        parser.add_argument("--snp",
                            type=str,
                            choices=["rsid", "chr:pos"],
                            default="rsid",
                            help="")
        parser.add_argument("--effect_size",
                            type=str,
                            choices=["beta", "z-score"],
                            default="beta",
                            help="")
        parser.add_argument("--pvalue",
                            type=str,
                            choices=["permuted", "bonferroni", "nominal"],
                            default="permuted",
                            help="")
        parser.add_argument("--alpha",
                            type=float,
                            required=False,
                            default=0.05,
                            help="")
        parser.add_argument("--fdr_calc_method",
                            type=str,
                            choices=["qvalues", "bh_fdr", "none"],
                            default="qvalues",
                            help="The multiple testing correction method "
                                 "to use. Default: 'qvalues'.")
        parser.add_argument("--log_modulus",
                            action='store_true',
                            help="")
        parser.add_argument("--palette",
                            type=str,
                            required=False,
                            default=None,
                            help="A color palette file.")
        parser.add_argument("--extensions",
                            nargs="+",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")
        parser.add_argument("--outdir",
                            type=str,
                            default=None,
                            help="")
        parser.add_argument("--force",
                            action='store_true',
                            help="")
        parser.add_argument("--save",
                            action='store_true',
                            help="")

        # Required external scripts.
        parser.add_argument("--qvalue_truncp",
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

        disc = self.get_method_class(method=self.discovery_method)(
            type="discovery",
            name=self.discovery_name,
            path=self.discovery_path,
            cell_type=self.discovery_cell_type,
            gene=self.gene,
            snp=self.snp,
            effect_size=self.effect_size,
            pvalue=self.pvalue
        )
        repl = self.get_method_class(method=self.replication_method)(
            type="replication",
            name=self.replication_name,
            path=self.replication_path,
            cell_type=self.replication_cell_type,
            gene=self.gene,
            snp=self.snp,
            effect_size=self.effect_size,
            pvalue=self.pvalue
        )

        # First get the top effect per gene in the discovery dataset.
        print("### Loading discovery data... ###")
        discovery_top_filepath = os.path.join(self.data_outdir, disc.get_id() + "_TopEffects.txt.gz")
        if os.path.exists(discovery_top_filepath) and not self.force:
            discovery_top_df = self.load_file(discovery_top_filepath)
        else:
            discovery_top_df = disc.get_top_effects()
            if self.save:
                self.save_file(df=discovery_top_df, outpath=discovery_top_filepath)
        print(discovery_top_df)
        print("\n")

        if discovery_top_df.shape[0] == 0:
            print("Error, discovery dataframe is empty.")
            exit()

        # Create a dict of top effects with keys being HGNC names and values being RS IDs.
        discovery_eqtls = disc.get_eqtl_effects(df=discovery_top_df)

        # Select the specific top effects from the discovery datasets in the replication dataset.
        print("### Loading replication data... ###")
        replication_filepath = os.path.join(self.data_outdir, disc.get_id() + "_TopEffects_LookupIn_" + repl.get_id() + ".txt.gz")
        if os.path.exists(replication_filepath) and not self.force:
            replication_df = self.load_file(replication_filepath)
        else:
            replication_df = repl.get_specific_effects(effects=discovery_eqtls)
            if self.save:
                self.save_file(df=replication_df, outpath=replication_filepath)
        print(replication_df)
        print("\n")

        if replication_df.shape[0] == 0:
            print("Error, replication dataframe is empty.")
            exit()

        # Overlap the data frames and harmonise direction of effect.
        print("### Overlapping summary statistics... ###")
        overlapped_df = self.overlap_summary_stats(
            disc_df=disc.standardize_format(df=discovery_top_df),
            repl_df=repl.standardize_format(df=replication_df),
            disc_name=disc.get_name(),
            repl_name=repl.get_name()
        )
        self.save_file(df=overlapped_df, outpath=os.path.join(self.data_outdir, disc.get_id() + "_Disc_" + repl.get_id() + "_Repl_MergedEffects.txt.gz"))
        print(overlapped_df)
        print("\n")

        # Remove missingness.
        overlapped_df = overlapped_df.loc[~overlapped_df[repl.get_name() + " effect_size"].isna(), :]

        # Calculate the replication statistics.
        print("### Calculating replication statistics... ###")
        replication_stats_df = self.calculate_replication_stats(
            df=overlapped_df,
            disc_name=disc.get_name(),
            repl_name=repl.get_name()
        )
        self.save_file(df=replication_stats_df, outpath=os.path.join(self.data_outdir, disc.get_id() + "_Disc_" + repl.get_id() + "_Repl_ReplicationStats.txt.gz"))
        print(replication_stats_df)
        print("\n")

        # Plot the comparison.
        print("### Plotting... ###")
        self.plot_replication(
            df=overlapped_df,
            replication_stats_df=replication_stats_df,
            disc_name=disc.get_name(),
            repl_name=repl.get_name(),
            disc_cell_type=disc.get_cell_type(),
            repl_cell_type=repl.get_cell_type(),
            log_modulus=self.log_modulus
        )
        print("### End ###")

    def print_arguments(self):
        print("Arguments:")
        print("  > Discovery method: {}".format(self.discovery_method))
        print("  > Discovery path: {}".format(self.discovery_path))
        print("  > Discovery name: {}".format(self.discovery_name))
        print("  > Discovery cell type: {}".format(self.discovery_cell_type))
        print("  > Replication method: {}".format(self.replication_method))
        print("  > Replication path: {}".format(self.replication_path))
        print("  > Replication name: {}".format(self.replication_name))
        print("  > Replication cell type: {}".format(self.replication_cell_type))
        print("  > Palette path: {}".format(self.palette_path))
        print("  > Gene: {}".format(self.gene))
        print("  > SNP: {}".format(self.snp))
        print("  > Effect size: {}".format(self.effect_size))
        print("  > P-value: {}".format(self.pvalue))
        print("  > Alpha: {}".format(self.alpha))
        print("  > FDR calc. method: {}".format(self.fdr_calc_method))
        print("  > Log modulus transform: {}".format(self.log_modulus))
        print("  > Data directory: {}".format(self.data_outdir))
        print("  > Plot directory: {}".format(self.plot_outdir))
        print("  > Extensions: {}".format(self.extensions))
        print("  > Force: {}".format(self.force))
        print("  > Save: {}".format(self.save))
        print("  > qvalues truncp script: {}".format(self.qvalue_truncp_script))
        print("  > Rb script: {}".format(self.rb_script))
        print("")

    @staticmethod
    def get_method_class(method):
        if method == "LIMIX":
            return LIMIX
        elif method == "LIMIX_REDUCED":
            return LIMIX_REDUCED
        elif method == "mbQTL":
            return mbQTL
        elif method == "mbQTL_MetaBrain":
            return mbQTL_MetaBrain
        elif method == "eQTLMappingPipeline":
            return eQTLMappingPipeline
        elif method == "eQTLgenPhase2":
            return eQTLgenPhase2
        elif method == "Bryois":
            return Bryois
        else:
            print("Error, unexpected method '{}'".format(method))
            exit()

    @staticmethod
    def load_file(inpath, header=0, index_col=0, sep="\t", low_memory=True,
                  nrows=None, skiprows=None, usecols=None):
        print("Loading file...")
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows,
                         usecols=usecols)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t", na_rep="NA",
                  sheet_name="Sheet1"):
        print("Saving file...")
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

    def overlap_summary_stats(self, disc_df, repl_df, disc_name, repl_name):
        print("Discovery dataframe '{}' has shape: {}".format(disc_name, disc_df.shape))
        print("Replication dataframe '{}' has shape: {}".format(repl_name, repl_df.shape))

        df = disc_df.merge(repl_df, left_index=True, right_index=True, how='left')

        na_mask = (df[repl_name + " pvalue"].isna()).to_numpy()
        print("Number of effects overlapping: {:,}".format(df.shape[0] - np.sum(na_mask)))
        if np.sum(na_mask) == df.shape[0]:
            return None

        allele_columns = [disc_name + " EA", disc_name + " OA", repl_name + " EA", repl_name + " OA"]
        if len(set(allele_columns).difference(df.columns)) == 0:
            alleles_match_mask = (((df[disc_name + " EA"] == df[repl_name + " EA"]) & (df[disc_name + " OA"] == df[repl_name + " OA"])) |
                                  ((df[disc_name + " EA"] == df[repl_name + " OA"]) & (df[disc_name + " OA"] == df[repl_name + " EA"])))
            mismatches_mask = ~na_mask & ~alleles_match_mask
            print("Number of variants with mismatched alleles: {:,}".format(np.sum(mismatches_mask)))
            if np.sum(mismatches_mask) > 0:
                # print(df.loc[mismatches_mask, [disc_name + " SNP", repl_name + " SNP"] + allele_columns])
                df.loc[mismatches_mask, [column for column in df.columns if column.startswith(repl_name)]] = np.nan
                na_mask = (df[repl_name + " pvalue"].isna()).to_numpy()
                print("Number of effects overlapping after mismatch removal: {:,}".format(df.shape[0] - np.sum(na_mask)))
                if np.sum(na_mask) == df.shape[0]:
                    return None
        else:
            print("Warning, could not verify that the alleles match. Assuming it is fine.")

        df[repl_name + " flip"] = False
        df.loc[~na_mask, repl_name + " flip"] = df.loc[~na_mask, repl_name + " EA"] != df.loc[~na_mask, disc_name + " EA"]
        df[repl_name + " effect_size"] = df[repl_name + " effect_size"] * df[repl_name + " flip"].map({True: -1, False: 1})
        df[repl_name + " EA"] = df[disc_name + " EA"]
        if repl_name + " OA" in df.columns:
            df[repl_name + " OA"] = df[disc_name + " OA"]
        print("Number of effects flipped: {:,}".format(df[repl_name + " flip"].sum()))

        fdr_calc_function = None
        if self.fdr_calc_method == "qvalues":
            fdr_calc_function = self.qvalues
        elif self.fdr_calc_method == "bh_fdr":
            fdr_calc_function = self.bh_fdr
        elif self.fdr_calc_method == "none":
            fdr_calc_function = lambda p: p
        else:
            print("Error, FDR calculation method '{}' is not implemented".format(self.fdr_calc_method))
            exit()

        if disc_name + " FDR" not in df.columns:
            print("Warning, discovery summary statistics did not include a dataset level multiple testing corrected p-value. "
                  "Calculating '{}' over '{}' p-values instead.".format(self.fdr_calc_method, self.pvalue))
            df[disc_name + " FDR"] = fdr_calc_function(p=df[disc_name + " pvalue"])
        discovery_mask = (df[disc_name + " FDR"] <= self.alpha).to_numpy()
        print("\nNumber of discovery significant effects: {:,}".format(np.sum(discovery_mask)))

        na_mask = (df[repl_name + " pvalue"].isna()).to_numpy()
        mask = np.logical_and(discovery_mask, ~na_mask)
        df[repl_name + " FDR"] = np.nan
        if np.sum(mask) > 1:
            df.loc[mask, repl_name + " FDR"] = fdr_calc_function(p=df.loc[mask, repl_name + " pvalue"])

        replication_mask = (df[repl_name + " FDR"] <= self.alpha).to_numpy()
        print("Number of replication significant effects: {:,}".format(np.sum(replication_mask)))
        
        return df
    
    def calculate_replication_stats(self, df, disc_name, repl_name):
        repl_stats_df = df.loc[~df[repl_name + " pvalue"].isna(), :].copy()
        repl_stats_df.sort_values(by=disc_name + " pvalue", inplace=True)
        
        disc_rb_beta, disc_rb_se = self.prepare_rb_info(df=repl_stats_df, name=disc_name)
        repl_rb_beta, repl_rb_se = self.prepare_rb_info(df=repl_stats_df, name=repl_name)

        replication_stats = []
        for disc_signif, repl_signif in [(False, False), (True, False), (True, True)]:
            mask = pd.Series(True, index=repl_stats_df.index)
            if disc_signif:
                mask = mask & (repl_stats_df[disc_name + " FDR"] <= self.alpha)
            if repl_signif:
                mask = mask & (repl_stats_df[repl_name + " FDR"] <= self.alpha)
            
            n = mask.sum()
            ac = np.nan
            coef = np.nan
            coefp = np.nan
            pi1 = np.nan
            rb = np.nan

            if n > 0:
                ac = self.ac(df=repl_stats_df.loc[mask, :], x=disc_name + " effect_size", y=repl_name + " effect_size")
            if n > 1:
                coef, coefp = self.pearson_r(df=repl_stats_df.loc[mask, :], x=disc_name + " effect_size", y=repl_name + " effect_size")

                if disc_signif and not repl_signif:
                    pi1 = min(self.pi1(df=repl_stats_df.loc[mask, :], x=repl_name + " pvalue"), 1)
                    rb = self.rb(
                        df=repl_stats_df.loc[mask, :],
                        b1=disc_rb_beta,
                        se1=disc_rb_se,
                        b2=repl_rb_beta,
                        se2=repl_rb_se,
                    )
                    if not np.isnan(rb):
                        rb = max(-1., min(rb, 1.))

            replication_stats.append([disc_signif, repl_signif, n, coef, coefp, ac, pi1, rb])
        replication_stats_df = pd.DataFrame(replication_stats,
                                            columns=["Disc significant",
                                                     "Repl significant",
                                                     "N",
                                                     "Coef",
                                                     "CoefP",
                                                     "AC",
                                                     "pi1",
                                                     "Rb"])
        replication_stats_df.index = replication_stats_df["Disc significant"].astype(str) + "_" + replication_stats_df["Repl significant"].astype(str)
        return replication_stats_df

    @staticmethod
    def bh_fdr(p):
        return multitest.multipletests(p, method='fdr_bh')[1]

    @staticmethod
    def qvalues(p):
        qvalue = importr("qvalue")
        pvals = robjects.FloatVector(p)
        qobj = robjects.r['qvalue'](pvals)
        return np.array(qobj.rx2('qvalues'))

    def prepare_rb_info(self, df, name):
        zscore_col = name + " effect_size"
        if self.effect_size == "beta":
            if name + " se" in df.columns:
                return name + " effect_size", name + " se"

            # Convert beta + p-value to z-score.
            zscore_col = name + " derived z-score"
            self.pvalue_to_zscore(df=df,
                                  beta_col=name + " effect_size",
                                  p_col=name + " pvalue",
                                  zscore_col=zscore_col)

        if not name + " MAF" in df or not name + " N" in df:
            return None, None

        self.zscore_to_beta(df=df,
                            zscore_col=zscore_col,
                            maf_col=name + " MAF",
                            n_col=name + " N",
                            beta_col=name + " derived beta",
                            se_col=name + " derived se")

        return name + " derived beta", name + " derived se"

    @staticmethod
    def pvalue_to_zscore(df, beta_col, p_col, zscore_col):
        p_values = df[p_col].to_numpy()
        zscores = stats.norm.ppf(p_values / 2)
        mask = np.ones_like(p_values)
        mask[df[beta_col] > 0] = -1
        df[zscore_col] = zscores * mask
        df.loc[df[p_col] == 1, zscore_col] = 0
        df.loc[df[p_col] == 0, zscore_col] = -40.

    @staticmethod
    def zscore_to_beta(df, zscore_col, maf_col, n_col, beta_col, se_col):
        chi = df[zscore_col] * df[zscore_col]
        a = 2 * df[maf_col] * (1 - df[maf_col]) * (df[n_col] + chi)
        df[beta_col] = df[zscore_col] / a ** (1/2)
        df[se_col] = 1 / a ** (1/2)

    @staticmethod
    def pearson_r(df, x, y):
        return stats.pearsonr(df[x], df[y])

    @staticmethod
    def ac(df, x, y):
        lower_quadrant = df.loc[(df[x] < 0) & (df[y] < 0), :]
        upper_quadrant = df.loc[(df[x] > 0) & (df[y] > 0), :]
        n_concordant = lower_quadrant.shape[0] + upper_quadrant.shape[0]
        concordance = (100 / df.shape[0]) * n_concordant
        return concordance

    def pi1(self, df, x):
        robjects.r("source('{}')".format(self.qvalue_truncp_script))
        p = robjects.FloatVector(df[x])
        qvalue_truncp = robjects.globalenv['qvalue_truncp']
        try:
            pi0 = qvalue_truncp(p)[0]
        except RRuntimeError as e:
            print("Warning, could not calculate pi1.")
            print(e)
            return np.nan

        return 1 - np.array(pi0)

    def rb(self, df, b1, se1, b2, se2, theta=0):
        if None in [b1, se1, b2, se2]:
            return np.nan

        robjects.r("source('{}')".format(self.rb_script))
        b1 = robjects.FloatVector(df[b1])
        se1 = robjects.FloatVector(df[se1])
        b2 = robjects.FloatVector(df[b2])
        se2 = robjects.FloatVector(df[se2])
        calcu_cor_true = robjects.globalenv['calcu_cor_true']
        rb = calcu_cor_true(b1, se1, b2, se2, theta)
        return np.array(rb)[0][0]

    def translate_cell_type(self, cell_type):
        if cell_type in self.cell_type_translate:
            return self.cell_type_translate[cell_type]
        return cell_type

    def plot_replication(self, df, replication_stats_df, disc_name, repl_name, disc_cell_type, repl_cell_type, log_modulus=False, title=""):
        plot_df = df[[
            disc_name + " label",
            disc_name + " effect_size",
            disc_name + " FDR",
            repl_name + " effect_size",
            repl_name + " FDR"
        ]].copy()

        color = "#000000"
        disc_standard_cell_type = self.translate_cell_type(cell_type=disc_cell_type)
        repl_standard_cell_type = self.translate_cell_type(cell_type=repl_cell_type)
        if (self.palette is not None) and (disc_standard_cell_type == repl_standard_cell_type) and (disc_standard_cell_type in self.palette.keys()):
            color = self.palette[disc_standard_cell_type]

        nrows = 1
        ncols = 3

        shared_xlim = {i: (0, 1) for i in range(ncols)}
        shared_ylim = {i: (0, 1) for i in range(nrows)}

        sns.set(rc={'figure.figsize': (ncols * 8, nrows * 6)})
        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharex='none',
                                 sharey='none')
        if nrows == 1:
            axes = axes[np.newaxis, ...]

        prefix = ""
        if log_modulus:
            prefix = " log"
            plot_df["{}{} effect_size".format(disc_name, prefix)] = self.calc_log_modulus(df=plot_df, column=disc_name + " effect_size")
            plot_df["{}{} effect_size".format(repl_name, prefix)] = self.calc_log_modulus(df=plot_df, column=repl_name + " effect_size")

        plot_df["facecolors"] = "#808080"
        plot_df.loc[(plot_df[disc_name + " FDR"] <= self.alpha) & (plot_df[repl_name + " FDR"] <= self.alpha), "facecolors"] = color

        print("Plotting column 1.")
        xlim, ylim = self.scatterplot(
            df=plot_df,
            fig=fig,
            ax=axes[0, 0],
            x="{}{} effect_size".format(disc_name, prefix),
            y="{}{} effect_size".format(repl_name, prefix),
            xlabel="{} {}{} eQTL {}".format(disc_name, disc_standard_cell_type, prefix, self.effect_size),
            ylabel="{} {}{} eQTL {}".format(repl_name, repl_standard_cell_type, prefix, self.effect_size),
            title="All overlapping",
            color=color,
            include_ylabel=True,
            annotations=["N = {:,}".format(replication_stats_df.loc["False_False", "N"]),
                         "r = {:.2f}".format(replication_stats_df.loc["False_False", "Coef"]),
                         "AC = {:.0f}%".format(replication_stats_df.loc["False_False", "AC"])]
        )
        shared_xlim, shared_ylim = self.update_limits(
            shared_xlim=shared_xlim,
            shared_ylim=shared_ylim,
            xlim=xlim,
            ylim=ylim,
            row=0,
            col=2
        )

        print("Plotting column 2.")
        xlim, ylim = self.scatterplot(
            df=plot_df.loc[plot_df[disc_name + " FDR"] <= self.alpha, :],
            fig=fig,
            ax=axes[0, 1],
            x="{}{} effect_size".format(disc_name, prefix),
            y="{}{} effect_size".format(repl_name, prefix),
            facecolors="facecolors",
            xlabel="{} {}{} eQTL {}".format(disc_name, disc_standard_cell_type, prefix, self.effect_size),
            ylabel="",
            title=disc_name + " signif.",
            color=color,
            include_ylabel=False,
            annotations=["N = {:,}".format(replication_stats_df.loc["True_False", "N"]),
                         "r = {:.2f}".format(replication_stats_df.loc["True_False", "Coef"]),
                         "AC = {:.0f}%".format(replication_stats_df.loc["True_False", "AC"]),
                         "\u03C01 = {:.2f}".format(replication_stats_df.loc["True_False", "pi1"]),
                         "Rb = {:.2f}".format(replication_stats_df.loc["True_False", "Rb"])]
        )
        shared_xlim, shared_ylim = self.update_limits(
            shared_xlim=shared_xlim,
            shared_ylim=shared_ylim,
            xlim=xlim,
            ylim=ylim,
            row=0,
            col=2
        )

        print("Plotting column 3.")
        xlim, ylim = self.scatterplot(
            df=plot_df.loc[(plot_df[disc_name + " FDR"] <= self.alpha) & (plot_df[repl_name + " FDR"] <= self.alpha), :],
            fig=fig,
            ax=axes[0, 2],
            x="{}{} effect_size".format(disc_name, prefix),
            y="{}{} effect_size".format(repl_name, prefix),
            label=disc_name + " label",
            xlabel="{} {}{} eQTL {}".format(disc_name, disc_standard_cell_type, prefix, self.effect_size),
            ylabel="",
            title="Both signif.",
            color=color,
            include_ylabel=False,
            annotations=["N = {:,}".format(replication_stats_df.loc["True_True", "N"]),
                         "r = {:.2f}".format(replication_stats_df.loc["True_True", "Coef"]),
                         "AC = {:.0f}%".format(replication_stats_df.loc["True_True", "AC"])]
        )
        shared_xlim, shared_ylim = self.update_limits(
            shared_xlim=shared_xlim,
            shared_ylim=shared_ylim,
            xlim=xlim,
            ylim=ylim,
            row=0,
            col=2
        )

        for (m, n), ax in np.ndenumerate(axes):
            (xmin, xmax) = shared_xlim[n]
            (ymin, ymax) = shared_ylim[m]

            xmargin = (xmax - xmin) * 0.05
            ymargin = (ymax - ymin) * 0.05

            ax.set_xlim(xmin - xmargin - 1, xmax + xmargin)
            ax.set_ylim(ymin - ymargin, ymax + ymargin)

        # Add the main title.
        fig.suptitle(title,
                     fontsize=25,
                     color="#000000",
                     weight='bold')

        plt.tight_layout()
        for extension in self.extensions:
            filename = "{}_{}_Disc_{}_{}_Repl.{}".format(disc_name, disc_cell_type, repl_name, repl_cell_type, extension)
            fig.savefig(os.path.join(self.plot_outdir, filename))
            print("Saved {}".format(filename))
        plt.close()

    @staticmethod
    def calc_log_modulus(df, column):
        return df[column].apply(lambda value: np.log(abs(value) + 1) * np.sign(value))

    def scatterplot(self, df, fig, ax, x="x", y="y", facecolors=None,
                    label=None, max_labels=15, xlabel="", ylabel="", title="",
                    color="#000000", ci=95, include_ylabel=True,
                    annotations=None):
        sns.despine(fig=fig, ax=ax)

        if not include_ylabel:
            ylabel = ""

        scatter_kws_facecolors = "#808080"
        if facecolors is not None:
            scatter_kws_facecolors = df[facecolors]

        n = df.shape[0]
        if n > 1:
            sns.regplot(x=x, y=y, data=df, ci=ci,
                        scatter_kws={'facecolors': scatter_kws_facecolors,
                                     'edgecolors': scatter_kws_facecolors},
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
                            arrowprops=dict(arrowstyle='-', color='#808080'))

        ax.axhline(0, ls='--', color="#808080", alpha=0.3, zorder=-1)
        ax.axvline(0, ls='--', color="#808080", alpha=0.3, zorder=-1)

        y_pos = 0.95
        for annotation in annotations:
            if annotation.startswith("N = "):
                if "N = {:,}".format(n) != annotation:
                    print("Warning, are you sure these stats match the data you are plotting?")

            ax.annotate(
                annotation,
                xy=(0.03, y_pos),
                xycoords=ax.transAxes,
                color=color,
                fontsize=14,
                fontweight='bold'
            )
            y_pos -= 0.05

        ax.set_title(title,
                     fontsize=18,
                     color=color,
                     weight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        return (df[x].min(), df[x].max()), (df[y].min(), df[y].max())

    @staticmethod
    def update_limits(shared_xlim, shared_ylim, xlim, ylim, row, col):
        row_ylim = shared_ylim[row]
        if ylim[0] < row_ylim[0]:
            row_ylim = (ylim[0], row_ylim[1])
        if ylim[1] > row_ylim[1]:
            row_ylim = (row_ylim[0], ylim[1])
        shared_ylim[row] = row_ylim

        col_xlim = shared_xlim[col]
        if xlim[0] < col_xlim[0]:
            col_xlim = (xlim[0], col_xlim[1])
        if xlim[1] > col_xlim[1]:
            col_xlim = (col_xlim[0], xlim[1])
        shared_xlim[col] = col_xlim
        return shared_xlim, shared_ylim

##############################################################################################################

class Dataset:
    def __init__(self, type, name, path, cell_type, gene, snp, effect_size, pvalue):
        self.class_name = None
        self.type = type
        self.name = name
        self.path = path
        self.cell_type = cell_type
        self.gene = gene
        self.snp = snp
        self.effect_size = effect_size
        self.pvalue = pvalue

        # Set the class variables.
        self.hgnc_gene = None
        self.ensembl_gene = None

        self.rsid_snp = None
        self.chr_pos_snp = None

        self.beta_effect_size = None
        self.zscore_effect_size = None

        self.nominal_pvalue = None
        self.bonferroni_pvalue = None
        self.permuted_pvalue = None

        # Option 0: data does not exist [(None, None, sep)]
        # Option 1: data is a full singular column [("A", None, sep)]
        # Option 2: data is part of a singular column [("A", "(a-zA-Z]+)_", sep)], make sure that group 1 of the regex captures the info of interest
        # Option 3: data is two or more full columns [("A", None, sep), ("B", None, sep)]
        # Option 4: data is a combination of option 2 and option 3 [("A", "(a-zA-Z]+)_", sep), ("B", None, sep)]
        # Option 5: data needs info from other file [("A", {"A": "a"}, sep)], prepare info as a translate dictionary
        self.columns = {
            "gene": [(None, None, None)],
            "SNP": [(None, None, None)],
            "EA": [(None, None, None)],
            "OA": [(None, None, None)],
            "effect_size": [(None, None, None)],
            "se": [(None, None, None)],
            "pvalue": [(None, None, None)],
            "FDR": [(None, None, None)],
            "N": [(None, None, None)],
            "MAF": [(None, None, None)],
            "label": [(None, None, None)]
        }

        # Default input files.
        self.top_effects_path = None
        self.all_effects_path = None

    def get_class_name(self):
        return self.class_name

    def get_type(self):
        return self.type

    def get_name(self):
        return self.name

    def get_path(self):
        return self.path
    
    def get_cell_type(self):
        return self.cell_type
    
    def get_id(self):
        if self.cell_type is None:
            return self.get_name()
        return self.name + "_" + self.cell_type

    def get_gene(self):
        return self.gene

    def get_gene_label(self):
        if self.gene == "hgnc":
            if self.hgnc_gene is None:
                print("Error, class '{}' does not have HGNC gene info.".format(self.class_name))
                raise ValueError()
            return self.hgnc_gene
        elif self.gene == "ensembl":
            if self.ensembl_gene is None:
                print("Error, class '{}' does not have ENSEMBL gene info.".format(self.class_name))
                raise ValueError()
            return self.ensembl_gene

        print("Error, unexpected problem in get_gene_label() for class '{}' with gene = {}, hgnc_gene = {}, and ensembl_gene = {}".format(self.name, self.gene, self.hgnc_gene, self.ensembl_gene))
        exit()

    def get_snp(self):
        return self.snp

    def get_snp_label(self):
        if self.snp == "rsid":
            if self.rsid_snp is None:
                print("Error, class '{}' does not have RSID SNP info.".format(self.class_name))
                raise ValueError()
            return self.rsid_snp
        elif self.snp == "chr:pos":
            if self.chr_pos_snp is None:
                print("Error, class '{}' does not have CHR:POS SNP info.".format(self.class_name))
                raise ValueError()
            return self.chr_pos_snp

        print("Error, unexpected problem in get_gene_label() for class '{}' with snp = {}, rsid_snp = {}, and chr_pos_snp = {}".format(self.name, self.snp, self.rsid_snp, self.chr_pos_snp))
        exit()

    def get_effect_size(self):
        return self.effect_size

    def get_effect_size_label(self):
        if self.effect_size == "beta":
            if self.beta_effect_size is None:
                print("Error, class '{}' does not have beta info.".format(self.class_name))
                raise ValueError()
            return self.beta_effect_size
        elif self.effect_size == "z-score":
            if self.zscore_effect_size is None:
                print("Error, class '{}' does not have z-score info.".format(self.class_name))
                raise ValueError()
            return self.zscore_effect_size

        print("Error, unexpected problem in get_effect_size_label() for class '{}' with effect_size = {}, beta = {}, and z-score = {}".format(self.name, self.effect_size, self.beta_effect_size, self.zscore_effect_size))
        exit()

    def get_pvalue(self):
        return self.pvalue

    def get_pvalue_label(self):
        if self.pvalue == "nominal" or self.type == "replication":
            if self.nominal_pvalue is None:
                print("Error, class '{}' does not have a nominal p-value info.".format(self.class_name))
                raise ValueError()
            return self.nominal_pvalue

        if self.type == "discovery" and self.pvalue == "bonferroni":
            if self.bonferroni_pvalue is None:
                print("Error, class '{}' does not have a bonferroni p-value info.".format(self.class_name))
                raise ValueError()
            return self.bonferroni_pvalue

        if self.type == "discovery" and self.pvalue == "permuted":
            if self.permuted_pvalue is None:
                print("Error, class '{}' does not have a permuted p-value info.".format(self.class_name))
                raise ValueError()
            return self.permuted_pvalue

        print("Error, unexpected problem in get_pvalue_label() for class '{}' with pvalue = {}, nominal_pvalue = {}, bonferroni_pvalue = {}, and permuted_pvalue = {}".format(self.name, self.pvalue, self.nominal_pvalue, self.bonferroni_pvalue, self.permuted_pvalue))
        exit()

    def get_eqtl_effects(self, df):
        return self.create_dict_from_df(df=df, key=self.columns["gene"], value=self.columns["SNP"])

    def create_dict_from_df(self, df, key, value):
        keys = df.apply(lambda row: self.extract_info(data=row, query=key), axis=1)
        values = df.apply(lambda row: self.extract_info(data=row, query=value), axis=1)
        return dict(zip(keys, values))

    @staticmethod
    def extract_info(data, query):
        info = ""
        for (column, pattern, sep) in query:
            if column is not None:
                if pattern is None:
                    info += str(data[column])
                elif isinstance(pattern, str):
                    try:
                        info += re.match(pattern, str(data[column])).group(1)
                    except AttributeError:
                        print("Error, pattern did not match in extract_info():\tre.match({}, {}).group(1)".format(pattern, str(data[column])))
                        exit()
                elif isinstance(pattern, dict):
                    if not data[column] in pattern:
                        return None
                    info += pattern[data[column]]
                else:
                    print("Error, unexpected input in extract_info()")
                    exit()

            if sep is not None:
                info += sep

        if info == "":
            info = None

        return info

    def trans_label_to_index(self, inpath, label_dict):
        colname_to_index_dict = self.get_file_colname_to_index(inpath)
        return self.update_label_dict(label_dict=label_dict, colname_to_index_dict=colname_to_index_dict)

    def get_file_colname_to_index(self, inpath, header=0, sep="\t"):
        colnames = self.get_file_header(inpath=inpath, header=header, sep=sep)
        return dict(zip(colnames, range(len(colnames))))

    @staticmethod
    def get_file_header(inpath, header=0, sep="\t"):
        if inpath.endswith(".gz"):
            fhin = gzip.open(inpath, 'rt')
        else:
            fhin = open(inpath, 'r')

        colnames = None
        for i, line in enumerate(fhin):
            if i != header:
                continue

            colnames = line.rstrip("\n").split(sep)
            break
        fhin.close()

        return colnames

    @staticmethod
    def update_label_dict(label_dict, colname_to_index_dict):
        indices_dict = {}
        for argument, query in label_dict.items():
            index_query = []
            for (column, pattern, sep) in query:
                index_query.append((colname_to_index_dict[column], pattern, sep))
            indices_dict[argument] = index_query
            del index_query

        return indices_dict

    def get_top_effects(self):
        df = self.get_effects(mode="top")
        return self.add_missing_info(df=df)

    def get_all_effects(self):
        df = self.get_effects(mode="all")
        return self.add_missing_info(df=df)

    def get_specific_effects(self, effects):
        df = self.get_effects(effects=effects, mode="specific")
        return self.add_missing_info(df=df)

    def get_effects(self, effects=None, mode="top"):
        all_effects_path = self.all_effects_path
        if all_effects_path and not os.path.exists(all_effects_path) and os.path.exists(all_effects_path + ".gz"):
            all_effects_path = all_effects_path + ".gz"

        top_effects_path = self.top_effects_path
        if top_effects_path and not os.path.exists(top_effects_path) and os.path.exists(top_effects_path + ".gz"):
            top_effects_path = top_effects_path + ".gz"

        df = None
        if mode == "all":
            df = self.load_file(all_effects_path)
        elif mode == "top":
            df = self.load_file(top_effects_path)
        elif mode == "specific":
            specific_entries_cols = self.trans_label_to_index(
                inpath=all_effects_path,
                label_dict={"key": self.columns["gene"], "value": self.columns["SNP"]}
            )
            df = self.load_partial_file(
                all_effects_path,
                effects=effects,
                specific_entries_cols=specific_entries_cols
            )
        else:
            print("Error, mode '{}' not implemented in get_effects() - {}".format(mode, self.class_name))
            exit()

        return df

    @staticmethod
    def load_file(inpath, header=0, index_col=None, sep="\t", low_memory=True,
                  nrows=None, skiprows=None, usecols=None):
        if not os.path.exists(inpath):
            print("Error, '{}' file not found".format(os.path.basename(inpath)))
            exit()

        print("Loading file...")
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows,
                         usecols=usecols)
        print("Loaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def load_excel(inpath, sheet_name, skiprows=None, dtype=None):
        if not os.path.exists(inpath):
            print("Error, '{}' file not found".format(os.path.basename(inpath)))
            exit()

        print("Loading file...")
        df = pd.read_excel(inpath, sheet_name=sheet_name, skiprows=skiprows, dtype=dtype)
        print("Loaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    def load_partial_file(self, inpath, header=0, index_col=None, sep="\t", nrows=None,
                          skiprows=None, usecols=None, top_entries_cols=None,
                          specific_entries_cols=None, effects=None, ignore_tabix=False):
        if not os.path.exists(inpath):
            print("Error, '{}' file not found".format(os.path.basename(inpath)))
            raise FileNotFoundError()

        if self.snp == "chr:pos" and os.path.exists(inpath + ".tbi") and not ignore_tabix:
            df = self.load_partial_file_w_tabix(
                inpath=inpath,
                header=header,
                index_col=index_col,
                sep=sep,
                usecols=usecols,
                top_entries_cols=top_entries_cols,
                specific_entries_cols=specific_entries_cols,
                effects=effects
            )
        else:
            df = self.load_partial_file_per_line(
                inpath=inpath,
                header=header,
                index_col=index_col,
                sep=sep,
                nrows=nrows,
                skiprows=skiprows,
                usecols=usecols,
                top_entries_cols=top_entries_cols,
                specific_entries_cols=specific_entries_cols,
                effects=effects
            )

        print("Loaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    def load_partial_file_w_tabix(self, inpath, header=0, index_col=None, sep="\t", usecols=None,
                                  top_entries_cols=None, specific_entries_cols=None, effects=None):
        if not os.path.exists(inpath):
            print("Error, '{}' file not found".format(os.path.basename(inpath)))
            raise FileNotFoundError()

        print("Loading partial file with tabix...")
        tb = tabix.open(inpath)
        columns = self.get_file_header(inpath=inpath, header=header, sep=sep)
        lines = []
        indices = None if index_col is None else []
        top_entries = {}
        n_effects = len(effects)
        i = 0
        for i, snp in enumerate(effects.values()):
            if i % 1e6 == 0:
                print("  Parsed {:,} / {:,} effects".format(i, n_effects), end='\r')
            chr, pos = snp.split(":")
            try:
                records = tb.query(chr, int(pos) - 1, int(pos) + 1)
            except tabix.TabixError:
                continue

            for values in records:
                if index_col is not None:
                    indices.append(values[index_col])
                    del values[index_col]

                if usecols is not None:
                    values = [values[i] for i in usecols]

                if specific_entries_cols is not None and effects is not None:
                    key = self.extract_info(data=values, query=specific_entries_cols["key"])
                    if key not in effects:
                        continue

                    if "value" in specific_entries_cols.keys():
                        value = self.extract_info(data=values, query=specific_entries_cols["value"])
                        if effects[key] != value:
                            continue

                if top_entries_cols is None:
                    lines.append(values)
                else:
                    key = values[top_entries_cols["key"]]
                    value = values[top_entries_cols["value"]]
                    if key not in top_entries:
                        top_entries[key] = (value, values)
                    elif key in top_entries and value < top_entries[key][0]:
                        top_entries[key] = (value, values)
                    else:
                        pass
        print("  Parsed {:,} / {:,} effects".format(i, n_effects))

        if top_entries_cols is not None:
            for _, (_, values) in top_entries.items():
                lines.append(values)

        return pd.DataFrame(lines, columns=columns, index=indices)

    def load_partial_file_per_line(self, inpath, header=0, index_col=None, sep="\t", nrows=None,
                               skiprows=None, usecols=None, top_entries_cols=None,
                               specific_entries_cols=None, effects=None):
        print("Loading partial file per line...")
        if inpath.endswith(".gz"):
            fhi = gzip.open(inpath, 'rt')
        else:
            fhi = open(inpath, 'r')

        columns = None
        indices = None if index_col is None else []
        lines = []
        header_row = header if skiprows is None else skiprows + header
        max_rows = nrows
        n_values = None
        top_entries = {}
        if skiprows is not None and max_rows is not None:
            max_rows += skiprows
        i = 0
        for i, line in enumerate(fhi):
            if i % 1e6 == 0:
                print("  Parsed {:,} lines".format(i), end='\r')

            if skiprows is not None and i < skiprows:
                continue
            if max_rows is not None and i > max_rows:
                break

            values = line.strip("\n").split(sep)
            if n_values is None:
                n_values = len(values)
            if len(values) != n_values:
                print("  Error, unequal number of columns in the input file, skip line")
                continue

            if index_col is not None:
                indices.append(values[index_col])
                del values[index_col]

            if usecols is not None:
                values = [values[i] for i in usecols]

            if i == header_row:
                columns = values
                continue
            if specific_entries_cols is not None and effects is not None:
                key = self.extract_info(data=values, query=specific_entries_cols["key"])
                if key not in effects:
                    continue

                if "value" in specific_entries_cols.keys():
                    value = self.extract_info(data=values, query=specific_entries_cols["value"])
                    if effects[key] != value:
                        continue

            if top_entries_cols is None:
                lines.append(values)
            else:
                key = values[top_entries_cols["key"]]
                value = values[top_entries_cols["value"]]
                if key not in top_entries:
                    top_entries[key] = (value, values)
                elif key in top_entries and value < top_entries[key][0]:
                    top_entries[key] = (value, values)
                else:
                    pass
        fhi.close()
        print("  Parsed {:,} lines".format(i))

        if top_entries_cols is not None:
            for _, (_, values) in top_entries.items():
                lines.append(values)

        return pd.DataFrame(lines, columns=columns, index=indices)

    def add_missing_info(self, df):
        return df

    @staticmethod
    def get_other_allele(row, alleles_column, pattern, effect_allele):
        match = re.match(pattern, row[alleles_column])
        allele1 = match.group(1)
        allele2 = match.group(2)
        if row[effect_allele] != allele1 and row[effect_allele] != allele2:
            return np.nan
        elif row[effect_allele] == allele1:
            return allele2
        else:
            return allele1

    def standardize_format(self, df):
        if df is None:
            return None

        df_info = []
        columns = []
        for i, (_, row) in enumerate(df.iterrows()):
            row_info = []
            for argument, query in self.columns.items():
                row_info.append(self.extract_info(data=row, query=query))
                if i == 0:
                    columns.append(self.name + " " + argument)
            df_info.append(row_info)
        standard_df = pd.DataFrame(df_info, columns=columns)
        standard_df = standard_df.replace("NA", np.nan)
        standard_df.dropna(axis=1, how='all', inplace=True)

        dtypes = {
            self.name + " gene": str,
            self.name + " SNP": str,
            self.name + " EA": str,
            self.name + " OA": str,
            self.name + " effect_size": float,
            self.name + " se": float,
            self.name + " pvalue": float,
            self.name + " FDR": float,
            self.name + " N": int,
            self.name + " MAF": float,
            self.name + " label": str
        }
        standard_df = standard_df.astype({key:value for key, value in dtypes.items() if key in standard_df.columns})
        standard_df.index = standard_df[self.name + " gene"] + "_" + standard_df[self.name + " SNP"]
        return standard_df


##############################################################################################################

class Bryois(Dataset):
    def __init__(self, *args, **kwargs):
        super(Bryois, self).__init__(*args, **kwargs)
        self.class_name = "Bryois"

        # Set the class variables.
        self.hgnc_gene = [("symbol_ensembl", "([a-zA-Z0-9-.]+)_ENSG[0-9]+", None)]
        self.ensembl_gene = [("symbol_ensembl", "[a-zA-Z0-9-.]+_(ENSG[0-9]+)", None)]
        self.rsid_snp = [("SNP", None, None)]
        self.chr_pos_snp = [("SNP_id_hg38", "chr((([0-9]{1,2}|X|Y|MT):[0-9]+))", None)]
        self.beta_effect_size = [("beta", None, None)]
        # self.zscore_effect_size = [(None, None, None)]
        self.nominal_pvalue = [("pval", None, None)]
        # self.bonferroni_pvalue = [(None, None, None)]
        self.permuted_pvalue = [("bpval", None, None)]

        # Define the uniform column info.
        self.columns.update({
            "gene": self.get_gene_label(),
            "SNP": self.get_snp_label(),
            "EA": [("effect_allele", None, None)],
            "OA": [("other_allele", None, None)],
            "effect_size": self.get_effect_size_label(),
            # "se": [(None, None, None)],
            "pvalue": self.get_pvalue_label(),
            "FDR": [("adj_p", None, None)] if self.type == "discovery" else [(None, None, None)],
            "N": [("N", None, None)],
            "MAF": [("MAF", None, None)],
            "label": self.hgnc_gene
        })

        # File paths.
        self.snp_pos_path = os.path.join(self.path, "snp_pos.txt.gz")
        self.top_effects_path = os.path.join(self.path, "41593_2022_1128_MOESM3_ESM.xlsx")
        self.all_effects_path = os.path.join(self.path, self.cell_type + ".<CHR>.gz")

        # Other variables.
        self.excel_gene_trans = {
            "ENSG00000099785": "MARCHF2",
            "ENSG00000100167": "SEPTIN3",
            "ENSG00000108387": "SEPTIN4",
            "ENSG00000117791": "MTARC2",
            "ENSG00000122545": "SEPTIN7",
            "ENSG00000136536": "MARCHF7",
            "ENSG00000138758": "SEPTIN11",
            "ENSG00000139266": "MARCHF9",
            "ENSG00000140623": "SEPTIN12",
            "ENSG00000144583": "MARCHF4",
            "ENSG00000145416": "MARCHF1",
            "ENSG00000145495": "MARCHF6",
            "ENSG00000154997": "SEPTIN14",
            "ENSG00000164402": "SEPTIN8",
            "ENSG00000165406": "MARCHF8",
            "ENSG00000168385": "SEPTIN2",
            "ENSG00000173077": "DELEC1",
            "ENSG00000173838": "MARCHF10",
            "ENSG00000173926": "MARCHF3",
            "ENSG00000180096": "SEPTIN1",
            "ENSG00000183654": "MARCHF11",
            "ENSG00000184640": "SEPTIN9",
            "ENSG00000184702": "SEPTIN5",
            "ENSG00000186205": "MTARC1",
            "ENSG00000186522": "SEPTIN10",
            "ENSG00000198060": "MARCHF5"
        }
        self.all_effects_columns = ["symbol_ensembl", "SNP", "dist_TSS", "pval", "beta"]
        self.all_effects_dtypes = {"symbol_ensembl": str, "SNP": str, "dist_TSS": int, "pval": float, "beta": float}
        self.n = {
            "Astrocytes": 192,
            "Endothelial.cells": 154,
            "Excitatory.neurons": 191,
            "Inhibitory.neurons": 173,
            "Microglia": 190,
            "Oligodendrocytes": 192,
            "OPCs...COPs": 188,
            "Pericytes": 165
        }[self.cell_type]

    def load_snp_pos_data(self, effects, specific_entries_cols):
        specific_entries_cols = self.trans_label_to_index(
            inpath=self.snp_pos_path,
            label_dict=specific_entries_cols
        )
        snp_pos_df =  self.load_partial_file(
            inpath=self.snp_pos_path,
            effects=effects,
            specific_entries_cols=specific_entries_cols
        )
        return snp_pos_df

    def get_partial_file_column_index(self):
        column_indices = {}
        for index, value in enumerate(self.all_effects_columns):
            column_indices[value] = index

        return column_indices

    def get_top_effects(self):
        df = self.load_excel(inpath=self.top_effects_path, sheet_name="Table S2", skiprows=3, dtype=str)
        df = df.loc[df["cell_type"] == self.cell_type.replace("...", " / ").replace(".", " "), :]
        for key, value in self.excel_gene_trans.items():
            df.loc[df["ensembl"] == key, "symbol"] = value
        df = df.astype({
            "cell_type": str,
            "symbol": str,
            "ensembl": str,
            "SNP": str,
            "effect_allele": str,
            "other_allele": str,
            "dist_TSS": int,
            "beta": float,
            "bpval": float,
            "adj_p": float,
            "beta_metabrain": float,
            "p_metabrain": float,
            "Replication": str}
        )
        df["symbol_ensembl"] = df["symbol"] + "_" + df["ensembl"]
        if df.shape[0] == 0:
            print("Error, no effects found")
            return None

        snp_pos_df = self.load_snp_pos_data(effects=set(df["SNP"]), specific_entries_cols={"key": [("SNP", None, None)]})
        df = df.merge(snp_pos_df, on=["SNP", "effect_allele", "other_allele"], how="left")
        df = self.add_missing_info(df=df)
        return df

    def get_effects(self, effects=None, mode="top"):
        if effects is None:
            snp_pos_df = self.load_file(self.snp_pos_path)
        else:
            snp_pos_df = self.load_snp_pos_data(effects=set(effects.values()), specific_entries_cols={"key": self.columns["SNP"]})

        specific_entries_cols = None
        if mode == "specific":
            label_dict = {"key": self.columns["gene"], "value": self.columns["SNP"]}
            if self.columns["SNP"][0] not in self.all_effects_columns:
                extra_info = self.create_dict_from_df(
                    df=snp_pos_df,
                    key=[("SNP", None, None)],
                    value=self.columns["SNP"]
                )
                label_dict = {"key": self.columns["gene"], "value": [("SNP", extra_info, None)]}

            specific_entries_cols = self.update_label_dict(
                label_dict=label_dict,
                colname_to_index_dict=self.get_partial_file_column_index()
            )

        df_list = []
        for chromosome in CHROMOSOMES:
            inpath = os.path.join(self.path, self.all_effects_path.replace("<CHR>", chromosome))

            df = None
            if mode == "all":
                df = self.load_file(inpath, header=None, sep=" ")
            elif mode == "specific":
                df = self.load_partial_file(
                    inpath,
                    header=None,
                    sep=" ",
                    effects=effects,
                    specific_entries_cols=specific_entries_cols
                )
            else:
                print("Error, mode '{}' not implemented in get_effects() - {}".format(mode, self.class_name))
                exit()

            df_list.append(df)

        if len(df_list) == 0:
            return None

        df = pd.concat(df_list, axis=0)
        df.columns = self.all_effects_columns
        df = df.astype(self.all_effects_dtypes)
        df = df.merge(snp_pos_df, on="SNP", how="left")
        return df

    def add_missing_info(self, df):
        if df.empty:
            return df

        # Add missing sample size. This is not always accurate but okay.
        df["N"] = self.n
        return df


##############################################################################################################


class LIMIX(Dataset):
    def __init__(self, *args, **kwargs):
        super(LIMIX, self).__init__(*args, **kwargs)
        self.class_name = "LIMIX"

        # Set the class variables.
        self.hgnc_gene = [("feature_id", None, None)]
        self.ensembl_gene = [("ENSG", None, None)]
        self.rsid_snp = [("snp_id", None, None)]
        # self.chr_pos_snp = [(None, None, None)]
        self.beta_effect_size = [("beta", None, None)]
        # self.zscore_effect_size = [(None, None, None)]
        self.nominal_pvalue = [("p_value", None, None)]
        self.bonferroni_pvalue = [("p_value_bonf_corr", None, None)]
        self.permuted_pvalue = [("empirical_feature_p_value", None, None)]

        # Define the uniform column info.
        self.columns.update({
            "gene": self.get_gene_label(),
            "SNP": self.get_snp_label(),
            "EA": [("assessed_allele", None, None)],
            "OA": [("other_allele", None, None)],
            "effect_size": self.get_effect_size_label(),
            "se": [("beta_se", None, None)],
            "pvalue": self.get_pvalue_label(),
            # "FDR": [(None, None, None)],
            "N": [("n_samples", None, None)],
            "MAF": [("maf", None, None)],
            "label": self.hgnc_gene
        })

        # File paths.
        # self.genotype_inpath = os.path.join(self.path, "genotype_input", self.ancestry, self.ancestry + "_imputed_hg38_stats_filtered.vars.gz")
        self.top_effects_path = os.path.join(self.path.replace("<CT>", self.cell_type), "top_qtl_results_all.txt")
        self.all_effects_path = os.path.join(self.path.replace("<CT>", self.cell_type), "qtl_results_all.txt")
        self.feature_metadata_file = os.path.join(self.path.replace("<CT>", self.cell_type), "qtl", "feature_metadata_<CHUNK>.txt")
        self.snp_metadata_file = os.path.join(self.path.replace("<CT>", self.cell_type), "qtl", "snp_metadata_<CHUNK>.txt")
        self.h5_file = os.path.join(self.path.replace("<CT>", self.cell_type), "qtl", "qtl_results_<CHUNK>.h5")
        self.snp_qc_metrics = os.path.join(self.path.replace("<CT>", self.cell_type), "qtl", "snp_qc_metrics_naContaining_feature_<FEATURE>.txt")

        # Other class variables.
        self.qtl_results_pattern = "qtl_results_([0-9]{1,2}|X|Y|MT)_([0-9]+)_([0-9]+).h5"
        self.qtl_chunks = None

    def get_qtl_chunks(self):
        if self.qtl_chunks is None:
            self.qtl_chunks = self.load_qtl_chunks()

        return self.qtl_chunks

    def load_qtl_chunks(self):
        qtl_results = glob.glob(self.h5_file.replace("<CHUNK>", "*"))
        qtl_results_data = []
        for filepath in qtl_results:
            basename = os.path.basename(filepath)
            match = re.match(self.qtl_results_pattern, basename)
            qtl_results_data.append([basename, match.group(1), int(match.group(2)), int(match.group(3))])

        df = pd.DataFrame(qtl_results_data, columns=["Filename", "Chromosome", "Start", "End"])
        df.insert(1, "Chunk", df["Chromosome"] + "_" + df["Start"].astype(str) + "_" + df["End"].astype(str))
        df = df.astype({"Filename": str,
                        "Chunk": str,
                        "Chromosome": str,
                        "Start": int,
                        "End": int})
        df.sort_values(by=["Chromosome", "Start"], key=natsort_keygen(), ascending=True, inplace=True)
        return df

    def get_effects(self, effects=None, mode="top"):
        all_effects_path = self.all_effects_path
        if all_effects_path and not os.path.exists(all_effects_path) and os.path.exists(all_effects_path + ".gz"):
            all_effects_path = all_effects_path + ".gz"

        top_effects_path = self.top_effects_path
        if top_effects_path and not os.path.exists(top_effects_path) and os.path.exists(top_effects_path + ".gz"):
            top_effects_path = top_effects_path + ".gz"

        df = None
        if mode == "all":
            df = self.load_file(all_effects_path)
        elif mode == "top":
            if os.path.exists(top_effects_path):
                df = self.load_file(top_effects_path)
            else:
                df = self.load_h5_files(top=True)
        elif mode == "specific":
            if os.path.exists(all_effects_path):
                specific_entries_cols = self.trans_label_to_index(
                    inpath=all_effects_path,
                    label_dict={"key": self.columns["gene"], "value": self.columns["SNP"]}
                )
                df = self.load_partial_file(
                    all_effects_path,
                    effects=effects,
                    specific_entries_cols=specific_entries_cols
                )
            else:
                df = self.load_h5_files(effects=effects)
        else:
            print("Error, mode '{}' not implemented in get_effects() - {}".format(mode, self.class_name))
            exit()

        return df

    def load_h5_files(self, effects=None, top=False):
        qlt_chunks_df = self.get_qtl_chunks()
        df_list = []
        for index, row in qlt_chunks_df.iterrows():
            df = self.load_h5_file(chunk=row["Chunk"], effects=effects, top=top)
            if df is not None:
                df_list.append(df)

        if len(df_list) == 0:
            return None

        return pd.concat(df_list, axis=0)

    def load_h5_file(self, chunk, effects=None, top=False):
        print("Loading chunk...")
        feature_metadata_file = self.feature_metadata_file.replace("<CHUNK>", chunk)
        if not os.path.exists(feature_metadata_file) and os.path.exists(feature_metadata_file + ".gz"):
            feature_metadata_file = feature_metadata_file + ".gz"

        snp_metadata_file = self.snp_metadata_file.replace("<CHUNK>", chunk)
        if not os.path.exists(snp_metadata_file) and os.path.exists(snp_metadata_file + ".gz"):
            snp_metadata_file = snp_metadata_file + ".gz"

        h5_file = self.h5_file.replace("<CHUNK>", chunk)

        try:
            if os.path.exists(feature_metadata_file):
                ffea_df = pd.read_table(feature_metadata_file, sep="\t")
            else:
                print("  Warning, skipping: '{}' missing feature metadata file.".format(chunk))
                return None
        except pd.errors.EmptyDataError:
            print("  Error, issue in feature annotation. Skipping '{}'.".format(chunk))
            return None

        feature_id_to_snp_dict = self.create_dict_from_df(df=ffea_df, key=[("feature_id", None, None)], value=self.columns["gene"])
        if effects is not None:
            if len(set(effects.keys()).intersection(set(feature_id_to_snp_dict.values()))) == 0:
                print("  Warning, skipping: '{}' no overlapping features.".format(chunk))
                return None

        try:
            if os.path.exists(snp_metadata_file):
                fsnp_df = pd.read_table(snp_metadata_file, sep="\t")
            else:
                print("  Warning, skipping: '{}' missing SNP metadata file.".format(chunk))
                return None
        except pd.errors.EmptyDataError:
            print("  Error, issue in snp annotation. Skipping '{}'.".format(chunk))
            return None

        snp_to_snp_id_dict = self.create_dict_from_df(df=fsnp_df, key=self.columns["SNP"], value=[("snp_id", None, None)])
        if effects is not None:
            if len(set(effects.values()).intersection(set(snp_to_snp_id_dict.keys()))) == 0:
                print("  Warning, skipping: '{}' no overlapping SNPs.".format(chunk))
                return None

        ffea_df = ffea_df.rename(index=str,
                                 columns={"chromosome": "feature_chromosome",
                                          "start": "feature_start",
                                          "end": "feature_end"})
        fsnp_df = fsnp_df.rename(index=str,
                                 columns={"chromosome": "snp_chromosome",
                                          "position": "snp_position"})

        frez = h5py.File(h5_file, 'r')
        frez_keys = [k.replace('_i_', '') for k in list(frez.keys())]

        df_list = []
        for frez_key in frez_keys:
            if effects is not None and feature_id_to_snp_dict[frez_key] not in effects.keys():
                continue

            frez_df = pd.DataFrame(np.array(frez[frez_key]))
            frez_df['feature_id'] = frez_key
            frez_df = frez_df.astype({"beta": float,
                                      "beta_se": float,
                                      "empirical_feature_p_value": float,
                                      "p_value": float,
                                      "snp_id": str,
                                      "feature_id": str})


            if effects is not None:
                snp = effects[feature_id_to_snp_dict[frez_key]]
                if (snp not in snp_to_snp_id_dict) or (snp_to_snp_id_dict[snp] not in set(frez_df["snp_id"])):
                    continue
                frez_df = frez_df.loc[frez_df["snp_id"] == snp_to_snp_id_dict[snp] ,:]

            df_list.append(frez_df)
            del frez_df

        if len(df_list) == 0:
            print("  Warning, skipping: '{}' no overlapping feature-SNP combinations.".format(chunk))
            return None

        df = pd.concat(df_list, axis=0)
        del df_list

        df = pd.merge(df, ffea_df, on='feature_id', how='left')

        if(len(glob.glob(self.snp_qc_metrics.replace("<FEATURE>", "*"))) > 0):
            print("  Error, code not implemented yet")
            exit()
            # tmp_df = pd.DataFrame(columns=df.columns)
            # for key in frez_keys:
            #     qc_metrics_inpath = self.snp_qc_metrics.replace("<FEATURE>", key)
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

        if effects is None:
            df['nTotalTestsPerFeat'] = np.nan
            for index, value in df['feature_id'].value_counts().iteritems():
                df.loc[df['feature_id'] == index, 'nTotalTestsPerFeat'] = value

        if top:
            df = df.sort_values(by=['empirical_feature_p_value', 'p_value'], ascending=[True, True])
            df = df.groupby(df['feature_id']).first()
            df.reset_index(drop=False, inplace=True)

        del ffea_df, fsnp_df

        print("Loaded chunk: {} "
              "with shape: {}".format(chunk,
                                      df.shape))
        return df

    def add_missing_info(self, df):
        if df.empty:
            return df

        # Add Bonferroni corrected p-value.
        if "nTotalTestsPerFeat" in df:
            df["p_value_bonf_corr"] = df[self.nominal_pvalue].astype(float) * df["nTotalTestsPerFeat"].astype(float)
            df.loc[df["p_value_bonf_corr"] > 1., "p_value_bonf_corr"] = 1.

        # Add other allele for check later.
        df["other_allele"] = df.apply(lambda row: self.get_other_allele(row=row, alleles_column="snp_id", pattern="(?:[0-9]{1,2}|X|Y|MT):[0-9]+:([A-Z]+):([A-Z]+)", effect_allele="assessed_allele"), axis=1)

        return df

##############################################################################################################

class LIMIX_REDUCED(LIMIX):
    def __init__(self, *args, **kwargs):
        super(LIMIX, self).__init__(*args, **kwargs)
        self.class_name = "LIMIX_REDUCED"

        # Set the class variables.
        self.hgnc_gene = [("feature_id", None, None)]
        self.ensembl_gene = [("ENSG", None, None)]
        # self.rsid_snp = [(None, None, None)]
        self.chr_pos_snp = [("snp_id", "(([0-9]{1,2}|X|Y|MT):([0-9]+))", None)]
        self.beta_effect_size = [("beta", None, None)]
        self.zscore_effect_size = [("z_score", None, None)]
        self.nominal_pvalue = [("p_value", None, None)]
        self.bonferroni_pvalue = [("p_value_bonf_corr", None, None)]
        self.permuted_pvalue = [("empirical_feature_p_value", None, None)]

        # Define the uniform column info.
        self.columns.update({
            "gene": self.get_gene_label(),
            "SNP": self.get_snp_label(),
            "EA": [("assessed_allele", None, None)],
            "OA": [("other_allele", None, None)],
            "effect_size": self.get_effect_size_label(),
            "se": [("beta_se", None, None)],
            "pvalue": self.get_pvalue_label(),
            "FDR": [("global_pvalue", None, None)] if self.type == "discovery" else [(None, None, None)],
            "N": [("n_samples", None, None)],
            "MAF": [("maf", None, None)],
            "label": self.hgnc_gene
        })

        appendix = ".wg3_Ye_wg3_wijst2018_wg3_sawcer_wg3_oneK1k_wg3_okada_wg3_Li_wg3_Franke_split_v3_wg3_Franke_split_v2_"
        if self.cell_type == "CD4_T":
            appendix = ".wg3_Ye_wg3_wijst2018_wg3_Trynka_wg3_sawcer_wg3_oneK1k_wg3_okada_wg3_Nawijn_wg3_Li_wg3_Franke_split_v3_wg3_Franke_split_v2_"

        # File paths.
        # self.genotype_inpath = None
        self.top_effects_path = os.path.join(self.path, self.cell_type + appendix + "top.txt")
        self.all_effects_path = os.path.join(self.path, self.cell_type + appendix + "all.txt")
        self.feature_metadata_file = None
        self.snp_metadata_file = None
        self.h5_file = None
        self.snp_qc_metrics = None

        # Other class variables.
        self.qtl_results_pattern = None
        self.qtl_chunks = None

##############################################################################################################

class mbQTL(Dataset):
    def __init__(self, *args, **kwargs):
        super(mbQTL, self).__init__(*args, **kwargs)
        self.class_name = "mbQTL"

        # Set the class variables.
        self.hgnc_gene = [("Gene", None, None)]
        # self.ensembl_gene = [(None, None, None)]
        self.rsid_snp = [("SNP", None, None)]
        self.chr_pos_snp = [("SNPChr", None, ":"), ("SNPPos", None, None)]
        self.beta_effect_size = [("MetaBeta", None, None)]
        self.zscore_effect_size = [("MetaPZ", None, None)]
        self.nominal_pvalue = [("MetaP", None, None)]
        self.bonferroni_pvalue = [("MetaPBonfCorr", None, None)]
        self.permuted_pvalue = [("BetaAdjustedMetaP", None, None)]

        # Define the uniform column info.
        self.columns.update({
            "gene": self.get_gene_label(),
            "SNP": self.get_snp_label(),
            "EA": [("SNPEffectAllele", None, None)],
            "OA": [("SNPOtherAllele", None, None)],
            "effect_size": self.get_effect_size_label(),
            "se": [("MetaSE", None, None)],
            "pvalue": self.get_pvalue_label(),
            # "FDR": [(None, None, None)],
            "N": [("MetaPN", None, None)],
            "MAF": [("MAF", None, None)],
            "label": self.hgnc_gene
        })

        # File paths.
        self.top_effects_path = os.path.join(self.path.replace("<CT>", self.cell_type), self.cell_type + "-TopEffects.txt")
        self.all_effects_path = os.path.join(self.path.replace("<CT>", self.cell_type), self.cell_type + "-AllEffects.txt")

    def add_missing_info(self, df):
        if df.empty:
            return df

        # Add the other allele column.
        df["SNPOtherAllele"] = df.apply(lambda row: self.get_other_allele(row=row, alleles_column="SNPAlleles", pattern="([A-Z]+)\/([A-Z]+)", effect_allele="SNPEffectAllele"), axis=1)

        # Update type of N.
        df["MetaPN"] = df["MetaPN"].astype(float).astype(int)

        # Add Bonferroni corrected p-value. Only exists in top_inpath.
        if "NrTestedSNPs" in df:
            df["MetaPBonfCorr"] = df["MetaP"].astype(float) * df["NrTestedSNPs"].astype(float)
            df.loc[df["MetaPBonfCorr"] > 1., "MetaPBonfCorr"] = 1.

        # Add MAF.
        df["MAF"] = df["SNPEffectAlleleFreq"].astype(float)
        df.loc[df["MAF"] > 0.5, "MAF"] = 1 - df.loc[df["MAF"] > 0.5, "MAF"]

        return df


##############################################################################################################


class mbQTL_MetaBrain(mbQTL):
    def __init__(self, *args, **kwargs):
        super(mbQTL, self).__init__(*args, **kwargs)
        self.class_name = "mbQTL_WG3"

        # Set the class variables.
        self.hgnc_gene = [("GeneSymbol", None, None)]
        self.ensembl_gene = [("Gene", "(ENSG[0-9]+).[0-9]+", None)]
        self.rsid_snp = [("SNP", "(?:[0-9]{1,2}|X|Y|MT):[0-9]+:(rs[0-9]+):[A-Z]+_[A-Z]+", None)]
        self.chr_pos_snp = [("SNPChr", None, ":"), ("SNPPos", None, None)]
        self.beta_effect_size = [("MetaBeta", None, None)]
        self.zscore_effect_size = [("MetaPZ", None, None)]
        self.nominal_pvalue = [("MetaP", None, None)]
        self.bonferroni_pvalue = [("MetaPBonfCorr", None, None)]
        self.permuted_pvalue = [("BetaAdjustedMetaP", None, None)]

        # Define the uniform column info.
        self.columns.update({
            "gene": self.get_gene_label(),
            "SNP": self.get_snp_label(),
            "EA": [("SNPEffectAllele", None, None)],
            "OA": [("SNPOtherAllele", None, None)],
            "effect_size": self.get_effect_size_label(),
            "se": [("MetaSE", None, None)],
            "pvalue": self.get_pvalue_label(),
            "FDR": [("qval", None, None)] if self.type == "discovery" else [(None, None, None)],
            "N": [("MetaPN", None, None)],
            "MAF": [("MAF", None, None)],
            "label": self.hgnc_gene
        })

        # File paths.
        self.top_effects_path = os.path.join(self.path, "merged-withqval.txt.gz")
        self.all_effects_path = os.path.join(self.path, "chr<CHR>-AllEffects-sorted.bgz.txt.gz")

    def get_top_effects(self):
        df = self.load_file(self.top_effects_path)
        return self.add_missing_info(df=df)

    def get_effects(self, effects=None, mode="top"):
        df_list = []
        for chromosome in CHROMOSOMES:
            inpath = os.path.join(self.path, self.all_effects_path.replace("<CHR>", chromosome))

            df = None
            if mode == "all":
                df = self.load_file(inpath, header=None, sep="\t")
            elif mode == "specific":
                specific_entries_cols = self.trans_label_to_index(
                    inpath=inpath,
                    label_dict={"key": self.columns["gene"], "value": self.columns["SNP"]}
                )
                df = self.load_partial_file(
                    inpath,
                    effects=effects,
                    specific_entries_cols=specific_entries_cols,
                    ignore_tabix=True
                )
            else:
                print("Error, mode '{}' not implemented in get_effects() - {}".format(mode, self.class_name))
                exit()

            df_list.append(df)

        if len(df_list) == 0:
            return None

        return pd.concat(df_list, axis=0)


##############################################################################################################


class eQTLMappingPipeline(Dataset):
    def __init__(self, *args, **kwargs):
        super(eQTLMappingPipeline, self).__init__(*args, **kwargs)
        self.class_name = "eQTLMappingPipeline"

        # Set the class variables.
        self.hgnc_gene = [("GeneSymbol", None, None)]
        self.ensembl_gene = [("Gene", None, None)]
        self.rsid_snp = [("SNP", None, None)]
        self.chr_pos_snp = [("SNPChr", None, ":"), ("SNPPos", None, None)]
        # self.beta_effect_size = [(None, None, None)]
        self.zscore_effect_size = [("Zscore", None, None)]
        self.nominal_pvalue = [("Pvalue", None, None)]
        self.bonferroni_pvalue = [("BonferroniP", None, None)]
        self.permuted_pvalue = [("FDR", None, None)]

        # Define the uniform column info.
        self.columns.update({
            "gene": self.get_gene_label(),
            "SNP": self.get_snp_label(),
            "EA": [("AssessedAllele", None, None)],
            "OA": [("OtherAllele", None, None)],
            "effect_size": self.get_effect_size_label(),
            # "se": [(None, None, None)],
            "pvalue": self.get_pvalue_label(),
            # "FDR": [(None, None, None)],
            "N": [("NrSamples", None, None)],
            # "MAF": [(None, None, None)],
            "label": self.hgnc_gene
        })

        # File paths.
        self.top_effects_path = os.path.join(self.path, "2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")
        self.all_effects_path = os.path.join(self.path, "cis-eQTLs_full_20180905.txt.gz")


##############################################################################################################


class eQTLgenPhase2(Dataset):
    def __init__(self, *args, **kwargs):
        super(eQTLgenPhase2, self).__init__(*args, **kwargs)
        self.class_name = "eQTLgenPhase2"

        # Set the class variables.
        # self.hgnc_gene = [(None, None, None)]
        self.ensembl_gene = [("phenotype", None, None)]
        self.rsid_snp = [("variant", None, None)]
        self.chr_pos_snp = [("chromosome_variant", None, ":"), ("bp_variant", None, None)]
        self.beta_effect_size = [("beta", None, None)]
        self.zscore_effect_size = [("z_score", None, None)]
        self.nominal_pvalue = [("p_value", None, None)]
        # self.bonferroni_pvalue = [(None, None, None)]
        # self.permuted_pvalue = [(None, None, None)]

        # Define the uniform column info.
        self.columns.update({
            "gene": self.get_gene_label(),
            "SNP": self.get_snp_label(),
            "EA": [("allele_eff", None, None)],
            "OA": [("allele_ref", None, None)],
            "effect_size": self.get_effect_size_label(),
            "se": [("standard_error", None, None)],
            "pvalue": self.get_pvalue_label(),
            # "FDR": [(None, None, None)],
            "N": [("sample_size", None, None)],
            "MAF": [("MAF", None, None)],
            # "label": self.hgnc_gene
        })

        # File paths.
        self.all_effects_path = os.path.join(self.path, "output_empirical_4GenPC20ExpPC_2023-05-27_interestingGeneSnpCombi_fixedConcatenation.csv") # Not really all effects, just the sc-eQTLgen ones
        # self.top_effects_path = None

    def get_top_effects(self):
        all_effects_path = self.all_effects_path
        if not os.path.exists(all_effects_path) and os.path.exists(all_effects_path + ".gz"):
            all_effects_path = all_effects_path + ".gz"

        top_entries_cols = self.trans_label_to_index(
            inpath=all_effects_path,
            label_dict={"key": self.columns["gene"], "value": self.columns["Pvalue"]}
        )
        df = self.load_partial_file(
            all_effects_path,
            top_entries_cols=top_entries_cols
        )
        return self.add_missing_info(df=df)

    def add_missing_info(self, df):
        if df.empty:
            return df

        # Add MAF.
        df["MAF"] = df["allele_eff_freq"].astype(float)
        df.loc[df["MAF"] > 0.5, "MAF"] = 1 - df.loc[df["MAF"] > 0.5, "MAF"]

        return df


if __name__ == '__main__':
    m = main()
    m.start()
