#!/usr/bin/env python3

"""
File:         replication.py
Created:      2024/02/28
Last Changed: 2024/06/22
Author:       M.Vochteloo

Copyright (C) 2024 M.Vochteloo
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
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adjustText import adjust_text
import cProfile
import pstats

# R imports. This is needed to prevent loading of the RProfile file.
import rpy2.rinterface as rinterface
rinterface.embedded.set_initoptions(options=["rpy2", "--quiet", "--no-save", "--no-init-file"])
rinterface.initr()
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from rpy2.rinterface_lib.embedded import RRuntimeError

# Local application imports.

"""
Syntax:
./replication.py -h
"""

# Metadata
__program__ = "Replication Plot"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@umcg.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)

CHROMOSOMES = [str(chromosome) for chromosome in range(1, 23)] + ["X", "Y", "MT"]
BATCHES = [str(batch) for batch in range(0, 100000)]
METHODS = ["CUSTOM", "LIMIX", "mbQTL", "eQTLMappingPipeline", "eQTLgenPhase2", "Bryois", "Bryois_REDUCED", "Fujita", "DeconQTL", "PICALO", "GTEx"]
EFFECTS = ["zscore", "beta"]


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.discovery_method = getattr(arguments, 'discovery_method')
        self.discovery_path = getattr(arguments, 'discovery_path')
        self.discovery_all_filename = getattr(arguments, 'discovery_all_filename')
        self.discovery_top_filename = getattr(arguments, 'discovery_top_filename')
        self.discovery_name = getattr(arguments, 'discovery_name')
        self.discovery_cell_type = getattr(arguments, 'discovery_cell_type')
        self.discovery_class_settings = getattr(arguments, 'discovery_class_settings')
        self.discovery_rm_dupl = getattr(arguments, 'discovery_rm_dupl')

        self.replication_method = getattr(arguments, 'replication_method')
        self.replication_path = getattr(arguments, 'replication_path')
        self.replication_all_filename = getattr(arguments, 'replication_all_filename')
        self.replication_top_filename = getattr(arguments, 'replication_top_filename')
        self.replication_name = getattr(arguments, 'replication_name')
        self.replication_cell_type = getattr(arguments, 'replication_cell_type')
        self.replication_class_settings = getattr(arguments, 'replication_class_settings')
        self.replication_rm_dupl = getattr(arguments, 'replication_rm_dupl')

        self.gene = "gene_" + getattr(arguments, 'gene')
        self.snp = "SNP_" + getattr(arguments, 'snp')
        self.pvalue = getattr(arguments, 'pvalue') + "_pvalue"
        self.effect = getattr(arguments, 'effect')
        self.allow_infer = getattr(arguments, 'allow_infer')
        self.alpha = getattr(arguments, 'alpha')
        self.fdr_calc_method = getattr(arguments, 'fdr_calc_method')
        self.log_modulus = getattr(arguments, 'log_modulus')
        self.gene_translate_path = getattr(arguments, 'gene_translate')
        self.cell_type_names_path = getattr(arguments, 'cell_type_names')
        self.palette_path = getattr(arguments, 'palette')
        self.extensions = getattr(arguments, 'extensions')
        outdir = getattr(arguments, 'outdir')
        self.force = getattr(arguments, 'force')
        self.save = getattr(arguments, 'save')
        self.verbose = getattr(arguments, 'verbose')
        self.dev = getattr(arguments, 'dev')
        self.qvalue_truncp_script = getattr(arguments, 'qvalue_truncp')
        self.rb_script = getattr(arguments, 'rb')

        for name, method, class_setting in [
            ("discovery", self.discovery_method, self.discovery_class_settings),
            ("replication", self.replication_method, self.replication_class_settings)]:
            if method == "CUSTOM" and (class_setting is None or not os.path.exists(class_setting)):
                print("Error, --{name}_class_settings needs to exist if --{name}_method is CUSTOM.".format(name=name))
                exit()
            if method != "CUSTOM" and class_setting is not None and os.path.exists(class_setting):
                print("Warning, ignoring --{name}_class_settings if --{name}_method is not CUSTOM.".format(name=name))

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

        # Matching cell types. Used to match different names for the same cell type. The value
        # name will be looked up in the palette if they match.
        self.cell_type_names = None
        if self.cell_type_names_path is not None:
            with open(self.cell_type_names_path) as f:
                self.cell_type_names = json.load(f)
            f.close()

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
                            help="Method of the discovery summary statistics.")
        parser.add_argument("--discovery_path",
                            type=str,
                            required=True,
                            default=None,
                            help="Basedir of the discovery summary statistics.")
        parser.add_argument("--discovery_all_filename",
                            type=str,
                            required=False,
                            default=None,
                            help="Update the default all effects filename of the discovery class.")
        parser.add_argument("--discovery_top_filename",
                            type=str,
                            required=False,
                            default=None,
                            help="Update the default top effects filename of the discovery class.")
        parser.add_argument("--discovery_name",
                            type=str,
                            required=True,
                            default=None,
                            help="Name of the discovery summary statistics.")
        parser.add_argument("--discovery_cell_type",
                            type=str,
                            required=False,
                            default="",
                            help="Cell type of the discovery summary statistics. Default: ''.")
        parser.add_argument("--discovery_class_settings",
                            required=False,
                            default=None,
                            help="JSON file containing custom discovery class settings.")
        parser.add_argument("--discovery_rm_dupl",
                            type=str,
                            choices=["none", "all", "allbutfirst"],
                            default="none",
                            help="How to deal with duplicates in the discovery summary statistics. Options: none) throw error and exit, all) remove all duplicates, allbutfirst) removed all except the first effect. Default: 'none'.")
        parser.add_argument("--replication_method",
                            type=str,
                            required=True,
                            choices=METHODS,
                            help="Method of the replication summary statistics.")
        parser.add_argument("--replication_path",
                            required=True,
                            default=None,
                            help="Basedir of the replication summary statistics.")
        parser.add_argument("--replication_all_filename",
                            type=str,
                            required=False,
                            default=None,
                            help="Update the default all effects filename of the replication class.")
        parser.add_argument("--replication_top_filename",
                            type=str,
                            required=False,
                            default=None,
                            help="Update the default top effects filename of the replication class.")
        parser.add_argument("--replication_name",
                            type=str,
                            required=True,
                            default=None,
                            help="Name of the replication summary statistics.")
        parser.add_argument("--replication_cell_type",
                            type=str,
                            required=False,
                            default="",
                            help="Cell type of the replication summary statistics. Default: ''.")
        parser.add_argument("--replication_class_settings",
                            required=False,
                            default=None,
                            help="JSON file containing custom discovery class settings.")
        parser.add_argument("--replication_rm_dupl",
                            type=str,
                            choices=["none", "all", "mismatched"],
                            default="none",
                            help="How to deal with duplicates in the replication summary statistics. Options: none) throw error and exit, all) remove all duplicates, mismatches) removed duplicates for which the effect allele does not match. Default: 'none'.")
        parser.add_argument("--gene",
                            type=str,
                            choices=["ensembl", "hgnc"],
                            default="ensembl",
                            help="Which gene format to select on. Default: 'ensembl'.")
        parser.add_argument("--snp",
                            type=str,
                            choices=["chr:pos", "rsid"],
                            default="chr:pos",
                            help="Which variant format to select on. Default: 'chr:pos'.")
        parser.add_argument("--pvalue",
                            type=str,
                            choices=["permuted", "bonferroni", "nominal"],
                            default="permuted",
                            help="Which pvalue to use. Default: 'permuted'.")
        parser.add_argument("--effect",
                            type=str,
                            choices=EFFECTS,
                            default="zscore",
                            help="What to consider as the effect column. Default: 'zscore'.")
        parser.add_argument("--allow_infer",
                            action='store_true',
                            help="Allow for inferring summary stats information. Note that inferred info are approximations and not exact. Default: False.")
        parser.add_argument("--alpha",
                            type=float,
                            required=False,
                            default=0.05,
                            help="The significance threshold to use. Default: 0.05.")
        parser.add_argument("--fdr_calc_method",
                            type=str,
                            choices=["qvalues", "bh_fdr", "none"],
                            default="qvalues",
                            help="The multiple testing correction method to use. Default: 'qvalues'.")
        parser.add_argument("--log_modulus",
                            action='store_true',
                            help="Transfer the effect column into log space while maintaining effect direction. Default: False.")
        parser.add_argument("--gene_translate",
                            type=str,
                            required=False,
                            default=None,
                            help="A tab seperated file containing ensembl - HGNC gene mappings. Default: None.")
        parser.add_argument("--cell_type_names",
                            type=str,
                            required=False,
                            default=None,
                            help="A JSON file containing method specific to normalised cell type names. Default: None.")
        parser.add_argument("--palette",
                            type=str,
                            required=False,
                            default=None,
                            help="A color palette file. Default: None.")
        parser.add_argument("--extensions",
                            nargs="+",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default=["png"],
                            help="The figure file extension. Default: 'png'.")
        parser.add_argument("--outdir",
                            type=str,
                            default=None,
                            help="The output directory. Default: current work directory.")
        parser.add_argument("--force",
                            action='store_true',
                            help="Whether to ignore previously loaded summary statistics. Default: False.")
        parser.add_argument("--save",
                            action='store_true',
                            help="Whether to store loaded summary statistics. Default: False.")
        parser.add_argument("--verbose",
                            action='store_true',
                            help="Print additional log messages. Default: False.")
        parser.add_argument("--dev",
                            action='store_true',
                            help="Print additional log messages. Default: False.")

        # Required external scripts.
        parser.add_argument("--qvalue_truncp",
                            type=str,
                            default="qvalue_truncp.R",
                            help="The path to the qvalues script. Default: 'qvalue_truncp.R'.")
        parser.add_argument("--rb",
                            type=str,
                            default="Rb.R",
                            help="The path to the Rb script. Default: 'Rb.R'.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        # Start dev.
        profiler = None
        if self.dev:
            profiler = self.enable_profiler()

        print("### Loading data classes and validating input arguments... ###")
        disc = self.get_method_class(method=self.discovery_method)(
            type="discovery",
            name=self.discovery_name,
            path=self.discovery_path,
            all_filename=self.discovery_all_filename,
            top_filename=self.discovery_top_filename,
            gene_translate_path=self.gene_translate_path,
            cell_type=self.discovery_cell_type,
            class_settings=self.discovery_class_settings,
            allow_infer=self.allow_infer,
            verbose=self.verbose
        )
        print(disc)

        repl = self.get_method_class(method=self.replication_method)(
            type="replication",
            name=self.replication_name,
            path=self.replication_path,
            all_filename=self.replication_all_filename,
            top_filename=self.replication_top_filename,
            gene_translate_path=self.gene_translate_path,
            cell_type=self.replication_cell_type,
            class_settings=self.replication_class_settings,
            allow_infer=self.allow_infer,
            verbose=self.verbose
        )
        print(repl)

        # Making sure the snp, gene, and effect data are actually present before we start doing anything.
        self.check_data_availability(class_object=disc, class_type="discovery")
        self.check_data_availability(class_object=repl, class_type="replication")

        # Compile the regex for speed-up.
        disc.compile_columns()
        repl.compile_columns()

        # First get the top effect per gene in the discovery dataset.
        print("### Loading discovery data... ###")
        discovery_top_filepath = os.path.join(self.data_outdir, disc.get_id() + "_TopEffects.txt.gz")
        if os.path.exists(discovery_top_filepath) and not self.force:
            discovery_top_df = self.load_file(discovery_top_filepath, index_col=None)
        else:
            discovery_top_df = disc.get_top_effects(gene=self.gene, snp=self.snp)
            if self.save:
                self.save_file(df=discovery_top_df, outpath=discovery_top_filepath, index=False)
        discovery_top_df = disc.postprocess(df=discovery_top_df)
        print(discovery_top_df)
        print("\n")

        if discovery_top_df is None or discovery_top_df.shape[0] == 0:
            print("Error, discovery dataframe is empty.")
            exit()

        # Create a dict of top effects with keys being gene names and values being SNPs.
        discovery_eqtls, n_disc_effects, n_disc_duplicated = disc.get_eqtl_effects(df=discovery_top_df, gene=self.gene, snp=self.snp)
        print("\tDiscovery has {:,} effects, {:,} were excluded due to duplication\n".format(n_disc_effects, n_disc_duplicated))
        if n_disc_effects == 0:
            print("Error, no effects in discovery dataframe to query.")
            exit()

        # Select the specific top effects from the discovery datasets in the replication dataset.
        print("### Loading replication data... ###")
        replication_filepath = os.path.join(self.data_outdir, disc.get_id() + "_TopEffects_LookupIn_" + repl.get_id() + ".txt.gz")
        if os.path.exists(replication_filepath) and not self.force:
            replication_df = self.load_file(replication_filepath, index_col=None)
        else:
            replication_df = repl.get_specific_effects(effects=discovery_eqtls, gene=self.gene, snp=self.snp)
            if self.save:
                self.save_file(df=replication_df, outpath=replication_filepath, index=False)
        replication_df = repl.postprocess(df=replication_df)
        print(replication_df)
        print("\n")

        if replication_df is None or replication_df.shape[0] == 0:
            print("Error, replication dataframe is empty.")
            exit()

        # Overlap the data frames and harmonise direction of effect.
        print("### Overlapping data... ###")
        overlapped_df = self.overlap_summary_stats(
            disc_df=disc.standardize_format(df=discovery_top_df, gene=self.gene, snp=self.snp),
            repl_df=repl.standardize_format(df=replication_df, gene=self.gene, snp=self.snp, include_fdr=False),
            disc_name=disc.get_name(),
            repl_name=repl.get_name()
        )
        self.save_file(df=overlapped_df, outpath=os.path.join(self.data_outdir, disc.get_id() + "_Disc_" + repl.get_id() + "_Repl_MergedEffects.txt.gz"))
        print(overlapped_df)
        print("\n")

        if overlapped_df is None or overlapped_df.shape[0] == 0:
            print("Error, merged dataframe is empty.")
            exit()

        # Remove missingness.
        overlapped_df = overlapped_df.loc[~overlapped_df[repl.get_name() + " " + self.effect].isna(), :]

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
        print("### Plotting data... ###")
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

        # Start dev.
        if self.dev:
            self.disable_profiler(profiler=profiler)

    @staticmethod
    def enable_profiler():
        profiler = cProfile.Profile()
        profiler.enable()
        return profiler

    @staticmethod
    def disable_profiler(profiler, sort='cumtime', n_stats=40):
        profiler.disable()
        stats = pstats.Stats(profiler).sort_stats(sort)
        stats.print_stats(n_stats)

    def print_arguments(self):
        print("Arguments:")
        print("  > Discovery method: {}".format(self.discovery_method))
        print("  > Discovery path: {}".format(self.discovery_path))
        print("  > Discovery all filename: {}".format(self.discovery_all_filename))
        print("  > Discovery top filename: {}".format(self.discovery_top_filename))
        print("  > Discovery name: {}".format(self.discovery_name))
        print("  > Discovery cell type: {}".format(self.discovery_cell_type))
        print("  > Discovery class settings: {}".format(self.discovery_class_settings))
        print("  > Discovery remove duplicates: {}".format(self.discovery_rm_dupl))
        print("  > Replication method: {}".format(self.replication_method))
        print("  > Replication path: {}".format(self.replication_path))
        print("  > Replication all filename: {}".format(self.replication_all_filename))
        print("  > Replication top filename: {}".format(self.replication_top_filename))
        print("  > Replication name: {}".format(self.replication_name))
        print("  > Replication cell type: {}".format(self.replication_cell_type))
        print("  > Replication class settings: {}".format(self.replication_class_settings))
        print("  > Replication remove duplicates: {}".format(self.replication_rm_dupl))
        print("  > Cell type names path: {}".format(self.cell_type_names_path))
        print("  > Palette path: {}".format(self.palette_path))
        print("  > Gene: {}".format(self.gene))
        print("  > SNP: {}".format(self.snp))
        print("  > P-value: {}".format(self.pvalue))
        print("  > Effect: {}".format(self.effect))
        print("  > Allow infer: {}".format(self.allow_infer))
        print("  > Alpha: {}".format(self.alpha))
        print("  > FDR calc. method: {}".format(self.fdr_calc_method))
        print("  > Log modulus transform: {}".format(self.log_modulus))
        print("  > Data directory: {}".format(self.data_outdir))
        print("  > Plot directory: {}".format(self.plot_outdir))
        print("  > Extensions: {}".format(self.extensions))
        print("  > Force: {}".format(self.force))
        print("  > Save: {}".format(self.save))
        print("  > Verbose: {}".format(self.verbose))
        print("  > Dev: {}".format(self.dev))
        print("  > qvalues truncp script: {}".format(self.qvalue_truncp_script))
        print("  > Rb script: {}".format(self.rb_script))
        print("")

    @staticmethod
    def get_method_class(method):
        if method == "CUSTOM":
            return CUSTOM
        elif method == "LIMIX":
            return LIMIX
        elif method == "mbQTL":
            return mbQTL
        elif method == "eQTLMappingPipeline":
            return eQTLMappingPipeline
        elif method == "eQTLgenPhase2":
            return eQTLgenPhase2
        elif method == "Bryois":
            return Bryois
        elif method == "Bryois_REDUCED":
            return Bryois_REDUCED
        elif method == "Fujita":
            return Fujita
        elif method == "DeconQTL":
            return DeconQTL
        elif method == "PICALO":
            return PICALO
        elif method == "GTEx":
            return GTEx
        else:
            print("Error, unexpected method '{}'".format(method))
            exit()

    def check_data_availability(self, class_object, class_type):
        columns_of_interest = [self.gene, self.snp, self.effect, "EA"]
        if class_type == "discovery":
            columns_of_interest.append(self.pvalue)
        elif class_type == "replication":
            columns_of_interest.append("nominal_pvalue")
        else:
            print("Error in check_data_availability(), class_type '{}' was unexpected.".format(class_type))
            exit()

        for column in columns_of_interest:
            if not class_object.contains_data(column):
                print("Error, {} class '{}' is missing '{}' column.".format(class_type, class_object.get_class_name(), column))
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
        # No need to save stuff if it is empty anyway.
        if df is None or df.shape[0] == 0:
            return

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

        # Make sure we have the EA for both datasets.
        if not disc_name + " EA" in disc_df.columns:
            print("Error, discovery is missing the effect allele column")
            exit()
        if not repl_name + " EA" in repl_df.columns:
            print("Error, replication is missing the effect allele column")
            exit()

        # Set index, prefer to match on alphabetically sorted alleles as well to ensure we are looking at the same variant.
        if disc_name + " alleles" in disc_df.columns and repl_name + " alleles" in repl_df.columns:
            disc_df.index = disc_df.index + "_" + self.reorder_alleles(df=disc_df, column=disc_name + " alleles")
            repl_df.index = repl_df.index + "_" + self.reorder_alleles(df=repl_df, column=repl_name + " alleles")
        else:
            print("Warning, could not verify that the alleles match. Assuming it is fine; use with caution!")

        # Resolve duplicates in discovery data.
        if len(set(disc_df.index)) != disc_df.shape[0]:
            if self.discovery_rm_dupl == "none":
                print("Error, discovery contains duplicate indices")
                self.print_duplicate_indices(df=disc_df)
                exit()
            elif self.discovery_rm_dupl == "all":
                print("Warning, discovery contains duplicate indices. Resolving this by removing"
                      " all duplicates.")
                disc_df = self.remove_duplicates(df=disc_df)
            elif self.discovery_rm_dupl == "allbutfirst":
                print("Warning, discovery contains duplicate indices. Attempting to resolve duplicates"
                      " by removing all effects except the first.")
                disc_df = self.remove_allbutfirst_duplicates(df=disc_df)
            else:
                print("Error, unexpected argument for --discovery_rm_dupl.")
                exit()

            if len(set(disc_df.index)) != disc_df.shape[0]:
                print("Error, discovery still contains duplicate indices.")
                self.print_duplicate_indices(df=disc_df)
                exit()

        # Resolve duplicates in replication data.
        if len(set(repl_df.index)) != repl_df.shape[0]:
            if self.replication_rm_dupl == "none":
                print("Error, replication contains duplicate indices")
                self.print_duplicate_indices(df=repl_df)
                exit()
            elif self.replication_rm_dupl == "all":
                print("Warning, replication contains duplicate indices. Resolving this by removing"
                      " all duplicates.")
                repl_df = self.remove_duplicates(df=repl_df)
            elif self.replication_rm_dupl == "mismatched":
                print("Warning, replication contains duplicate indices. Attempting to resolve duplicates"
                      " by removing effects that require flipping.")
                repl_df = self.remove_mismatched_duplicates(disc_df=disc_df, disc_name=disc_name, repl_df=repl_df, repl_name=repl_name)
            else:
                print("Error, unexpected argument for --replication_rm_dupl.")
                exit()

            if len(set(repl_df.index)) != repl_df.shape[0]:
                print("Error, replication still contains duplicate indices.")
                self.print_duplicate_indices(df=repl_df)
                exit()

        print("Discovery dataframe:")
        print(disc_df)

        print("\nReplication dataframe:")
        print(repl_df)
        print("")

        n_overlap = len(set(disc_df.index).intersection(set(repl_df.index)))
        print("Number of effects overlapping: {:,}".format(n_overlap))
        if n_overlap == 0:
            exit()

        # Merge the data.
        df = disc_df.merge(repl_df, left_index=True, right_index=True, how='left')
        print("Merged dataframe has shape: {}\n".format(df.shape))

        # Check which effects are missing in the replication data frame.
        na_mask = (df[repl_name + " " + self.effect].isna()).to_numpy()

        # Flip the effects. Make sure we flip both the beta and z-score and whatever else for effect columns there might be.
        df[repl_name + " flip"] = False
        df.loc[~na_mask, repl_name + " flip"] = df.loc[~na_mask, repl_name + " EA"] != df.loc[~na_mask, disc_name + " EA"]
        for effect in EFFECTS:
            if repl_name + " " + effect in df.columns:
                df[repl_name + " " + effect] = df[repl_name + " " + effect] * df[repl_name + " flip"].map({True: -1, False: 1})
        if repl_name + " AF" in df.columns:
            df.loc[df[repl_name + " flip"], repl_name + " AF"] = 1 - df.loc[df[repl_name + " flip"], repl_name + " AF"]

        # Fix the effect allele and other allele columns.
        if disc_name + " OA" in df.columns:
            # If the discovery has the other allele we can just copy all the info from the discovery.
            df[repl_name + " OA"] = df[disc_name + " OA"]
            df[repl_name + " EA"] = df[disc_name + " EA"]
        elif repl_name + " OA" in df.columns:
            # If the discovery does not have the other allele but the replication does. We then first need to update
            # the other allele column and then we can just copy the effect allele from the discovery.
            df.loc[df[repl_name + " flip"], repl_name + " OA"] = df.loc[df[repl_name + " flip"], repl_name + " EA"]
            df[repl_name + " EA"] = df[disc_name + " EA"]
        else:
            # If neither one has the other allele we can technically move the flipped allees to the other allele column
            # but since we aren't sure about the alleles I rather just keep them empty.
            pass
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

        print("")
        if disc_name + " FDR" not in df.columns:
            print("Warning, discovery summary statistics did not include a dataset level multiple testing corrected p-value. "
                  "Calculating '{}' over '{}' p-values instead.".format(self.fdr_calc_method, self.pvalue))
            df[disc_name + " FDR"] = fdr_calc_function(p=df[disc_name + " " + self.pvalue])
        discovery_mask = (df[disc_name + " FDR"] <= self.alpha).to_numpy()
        print("Number of discovery significant effects: {:,}".format(np.sum(discovery_mask)))

        mask = np.logical_and(discovery_mask, ~na_mask)
        df[repl_name + " FDR"] = np.nan
        if np.sum(mask) > 1:
            df.loc[mask, repl_name + " FDR"] = fdr_calc_function(p=df.loc[mask, repl_name + " nominal_pvalue"])

        replication_mask = (df[repl_name + " FDR"] <= self.alpha).to_numpy()
        print("Number of replication significant effects: {:,}".format(np.sum(replication_mask)))

        return df

    @staticmethod
    def reorder_alleles(df, column):
        return df[column].apply(lambda value: "/".join(sorted(value.split("/"))))

    @staticmethod
    def remove_duplicates(df):
        return df.loc[~df.index.duplicated(keep=False), :]

    @staticmethod
    def remove_allbutfirst_duplicates(df):
        return df.groupby(df.index).first()

    @staticmethod
    def remove_mismatched_duplicates(disc_df, disc_name, repl_df, repl_name):
        # Create a set of effects just considering gene_snp_EA. We consider this the reference.
        disc_effects = set(disc_df.index + "_" + disc_df[disc_name + " EA"])

        # Construct a mask for the effects in the replication data frame that are duplicated.
        # Also, construct a mask that checks if the gene_snp_EA is matching the reference set (discovery effects).
        duplicated = repl_df.index.duplicated(keep=False)
        matched_ea = (repl_df.index + "_" + repl_df[repl_name + " EA"]).isin(disc_effects)

        # Only keep effects that are unique or are matching the discovery effects. Note that
        # this can still contain duplicates if the gene_snp_EA combinations occurs more than once
        # in the replication data frame. In that case I have no clue how to resolve this.
        return repl_df.loc[~duplicated | matched_ea, :]

    @staticmethod
    def print_duplicate_indices(df):
        dups = df.index[df.index.duplicated()]
        for index in dups:
            print(df.loc[[index], :].T)
        print("\tFound {:,} duplicates: {}".format(len(dups), ", ".join(dups)))

    def calculate_replication_stats(self, df, disc_name, repl_name):
        repl_stats_df = df.copy()
        df.sort_values(by=disc_name + " " + self.pvalue, inplace=True)

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
            rb_se = np.nan
            rb_pvalue = np.nan

            if n > 0:
                ac = self.ac(
                    df=repl_stats_df.loc[mask, :],
                    x=disc_name + " " + self.effect,
                    y=repl_name + " " + self.effect
                )
            if n > 1:
                coef, coefp = self.pearson_r(
                    df=repl_stats_df.loc[mask, :],
                    x=disc_name + " " + self.effect,
                    y=repl_name + " " + self.effect
                )

                if disc_signif and not repl_signif:
                    pi1 = self.pi1(
                        df=repl_stats_df.loc[mask, :],
                        x=repl_name + " nominal_pvalue"
                    )

                    rb, rb_se, rb_pvalue = self.rb(
                        df=repl_stats_df.loc[mask, :],
                        b1=disc_name + " beta",
                        se1=disc_name + " beta_se",
                        b2=repl_name + " beta",
                        se2=repl_name + " beta_se",
                    )

            replication_stats.append([disc_signif, repl_signif, n, coef, coefp, ac, pi1, rb, rb_se, rb_pvalue])
        replication_stats_df = pd.DataFrame(replication_stats,
                                            columns=["Disc significant",
                                                     "Repl significant",
                                                     "N",
                                                     "Coef",
                                                     "CoefP",
                                                     "AC",
                                                     "pi1",
                                                     "Rb",
                                                     "RbSE",
                                                     "RbP"])
        replication_stats_df.index = replication_stats_df["Disc significant"].astype(str) + "_" + replication_stats_df["Repl significant"].astype(str)
        return replication_stats_df

    @staticmethod
    def bh_fdr(p):
        return multitest.multipletests(p, method='fdr_bh')[1]

    @staticmethod
    def qvalues(p):
        qvalues = pd.Series(np.nan, index=p)
        try:
            qvalue = importr("qvalue")
            pvals = robjects.FloatVector(p)
            qobj = robjects.r['qvalue'](pvals)
            qvalues = np.array(qobj.rx2('qvalues'))
        except RRuntimeError as e:
            print("Error, could not calculate qvalues.")
            print(e)
            exit()

        return qvalues

    @staticmethod
    def pearson_r(df, x, y):
        return stats.pearsonr(df[x], df[y])

    @staticmethod
    def ac(df, x, y):
        lower_quadrant = ((df[x] < 0) & (df[y] < 0)).sum()
        upper_quadrant = ((df[x] > 0) & (df[y] > 0)).sum()
        n_concordant = lower_quadrant + upper_quadrant
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

        pi1 = 1 - np.array(pi0)
        return min(pi1, 1)

    def rb(self, df, b1, se1, b2, se2, theta=0):
        valid = True
        for column in [b1, se1, b2, se2]:
            if column not in df.columns:
                print("Warning, could not calculate Rb due to missing column '{}'.".format(column))
                valid = False
                continue

            if df[column].std() == 0:
                print("Warning, could not calculate Rb due to zero variance in column '{}'.".format(column))
                valid = False
        if not valid:
            return np.nan, np.nan, np.nan

        robjects.r("source('{}')".format(self.rb_script))
        b1 = robjects.FloatVector(df[b1])
        se1 = robjects.FloatVector(df[se1])
        b2 = robjects.FloatVector(df[b2])
        se2 = robjects.FloatVector(df[se2])
        calcu_cor_true = robjects.globalenv['calcu_cor_true']

        try:
            result = calcu_cor_true(b1, se1, b2, se2, theta)
        except RRuntimeError as e:
            print("Warning, could not calculate Rb.")
            print(e)
            return np.nan, np.nan, np.nan

        rb, rb_se, rb_pvalue = np.array(result)[0]
        if not np.isnan(rb):
            rb = max(-1., min(rb, 1.))

        return rb, rb_se, rb_pvalue

    def translate_cell_type(self, cell_type):
        if self.cell_type_names is not None and cell_type in self.cell_type_names:
            return self.cell_type_names[cell_type]
        return cell_type

    def plot_replication(self, df, replication_stats_df, disc_name, repl_name, disc_cell_type, repl_cell_type, log_modulus=False, title=""):
        # Plot the HGNC symbols of the top effects. By default take the discovery HGNC symbol unless
        # it does not exist, then take the replication or else take the discovery ensembl IDs.
        if disc_name + " gene_hgnc" in df.columns:
            label = disc_name + " gene_hgnc"
        elif repl_name + " gene_hgnc" in df.columns:
            label = repl_name + " gene_hgnc"
        else:
            label = disc_name + " gene_ensembl"

        disc_effect = disc_name + " " + self.effect
        disc_signif = disc_name + " FDR"
        repl_effect = repl_name + " " + self.effect
        repl_signif = repl_name + " FDR"

        plot_df = df[[label, disc_effect, disc_signif, repl_effect, repl_signif]].copy()

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

        if log_modulus:
            old_disc_effect = disc_effect
            old_repl_effect = repl_effect

            disc_effect = "{} log {}".format(disc_name, self.effect)
            repl_effect = "{} log {}".format(repl_name, self.effect)

            plot_df[disc_effect] = self.calc_log_modulus(df=plot_df, column=old_disc_effect)
            plot_df[repl_effect] = self.calc_log_modulus(df=plot_df, column=old_repl_effect)

            del old_disc_effect, old_repl_effect

        plot_df["facecolors"] = "#808080"
        plot_df.loc[(plot_df[disc_signif] <= self.alpha) & (plot_df[repl_signif] <= self.alpha), "facecolors"] = color

        print("Plotting column 1.")
        xlim, ylim = self.scatterplot(
            df=plot_df,
            fig=fig,
            ax=axes[0, 0],
            x=disc_effect,
            y=repl_effect,
            title="All overlapping",
            color=color,
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
            col=0
        )

        print("Plotting column 2.")
        xlim, ylim = self.scatterplot(
            df=plot_df.loc[plot_df[disc_signif] <= self.alpha, :],
            fig=fig,
            ax=axes[0, 1],
            x=disc_effect,
            y=repl_effect,
            facecolors="facecolors",
            ylabel="",
            title=disc_name + " signif.",
            color=color,
            annotations=["N = {:,}".format(replication_stats_df.loc["True_False", "N"]),
                         "r = {:.2f}".format(replication_stats_df.loc["True_False", "Coef"]),
                         "AC = {:.0f}%".format(replication_stats_df.loc["True_False", "AC"]),
                         "\u03C01 = {:.2f}".format(replication_stats_df.loc["True_False", "pi1"]),
                         "Rb = {:.2f} (±{:.2f})".format(replication_stats_df.loc["True_False", "Rb"], replication_stats_df.loc["True_False", "RbSE"])]
        )
        shared_xlim, shared_ylim = self.update_limits(
            shared_xlim=shared_xlim,
            shared_ylim=shared_ylim,
            xlim=xlim,
            ylim=ylim,
            row=0,
            col=1
        )

        print("Plotting column 3.")
        xlim, ylim = self.scatterplot(
            df=plot_df.loc[(plot_df[disc_signif] <= self.alpha) & (plot_df[repl_signif] <= self.alpha), :],
            fig=fig,
            ax=axes[0, 2],
            x=disc_effect,
            y=repl_effect,
            label=label,
            ylabel="",
            title="Both signif.",
            color=color,
            annotations=["N = {:,}".format(replication_stats_df.loc["True_True", "N"]),
                         "r = {:.2f}".format(replication_stats_df.loc["True_True", "Coef"]),
                         "AC = {:.0f}%".format(replication_stats_df.loc["True_True", "AC"]),
                         "%repl = {:.0f}%".format((100 / replication_stats_df.loc["True_False", "N"]) * replication_stats_df.loc["True_True", "N"])]
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
                    label=None, max_labels=15, xlabel=None, ylabel=None, title="",
                    color="#000000", ci=95, annotations=None):
        sns.despine(fig=fig, ax=ax)

        if xlabel is None:
            xlabel = x
        if ylabel is None:
            ylabel = y

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

        low_x, high_x = ax.get_xlim()
        low_y, high_y = ax.get_ylim()
        low = max(low_x, low_y)
        high = min(high_x, high_y)
        ax.plot([low, high], [low, high], ls='--', color="#808080", alpha=0.3, zorder=-1)

        y_pos = 0.95
        for annotation in annotations:
            # Doing a little check since the replication stats are calculated outside this function.
            # If the sample size don't match there might be an issue.
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
    def __init__(self, type, name, path, all_filename, top_filename, gene_translate_path, cell_type, class_settings, allow_infer, verbose):
        self.class_name = None
        self.type = type
        self.name = name
        self.path = path
        self.all_filename = all_filename
        self.top_filename = top_filename
        self.cell_type = cell_type
        self.gene_translate_path = gene_translate_path if (gene_translate_path is not None and os.path.exists(gene_translate_path)) else None
        self.class_settings = class_settings
        self.allow_infer = allow_infer
        self.verbose = verbose

        # Default input files.
        self.all_effects_path = None
        self.top_effects_path = None

        # Default effects info.
        self.effects_header = None
        self.effects_line = None

        # Gene translate dicts.
        self.ensembl_to_hgnc = None
        self.hgnc_to_ensembl = None

        # Default sample size, assume equal for all effects.
        self.n = None

        # Default settings.
        self.ignore_tabix = True

        # Option 0: data does not exist [(None, None, suffix)]
        # Option 1: data is a full singular column [("A", None, suffix)]
        # Option 2: data is part of a singular column [("A", "(a-zA-Z]+)_", suffix)], make sure that group 1 of the regex captures the info of interest
        # Option 3: data is two or more full columns [("A", None, suffix), ("B", None, suffix)]
        # Option 4: data is two or more full columns with a symbol inbetween `[("A", null, ":"), ("B", null, null)]`
        # Option 5: data is a combination of option 2 and option 3 [("A", "(a-zA-Z]+)_", suffix), ("B", None, suffix)]
        # Option 6: data needs info from other file [("A", {"A": "a"}, suffix)], prepare info as a translate dictionary
        # TODO; one thing I did not consider yet is if you want to use a regex pattern and then apply a translate dict. Since this
        #  did not occur yet I will just deal with it once I encounter it.
        self.na = [(None, None, None)]
        self.columns = {
            "gene_hgnc": self.na,
            "gene_ensembl": self.na,
            "SNP_rsid": self.na,
            "SNP_chr:pos": self.na,
            "alleles": self.na,
            "EA": self.na,
            "OA": self.na,
            "beta": self.na,
            "beta_se": self.na,
            "n_tests": self.na,
            "nominal_pvalue": self.na,
            "permuted_pvalue": self.na,
            "bonf_pvalue": self.na,
            "zscore": self.na,
            "FDR": self.na,
            "N": self.na,
            "AF": self.na,
            "MAN": self.na,
            "MAF": self.na
        }

    def get_class_name(self):
        return self.class_name

    def get_type(self):
        return self.type

    def get_name(self):
        return self.name

    def get_path(self):
        return self.path

    def update_all_filename(self, filename):
        if self.all_filename is None and filename != "":
            self.all_filename = filename

    def get_all_filename(self):
        return self.all_filename

    def update_top_filename(self, filename):
        if self.top_filename is None and filename != "":
            self.top_filename = filename

    def get_top_filename(self):
        return self.top_filename

    def set_all_effects_path(self):
        if self.all_filename is None:
            return
        self.all_effects_path = self.get_fpath(fpath=os.path.join(self.path, self.all_filename))

    def get_all_effects_path(self):
        return self.all_effects_path

    def set_top_effects_path(self):
        if self.top_filename is None:
            return
        self.top_effects_path = self.get_fpath(fpath=os.path.join(self.path, self.top_filename))

    def get_top_effects_path(self):
        return self.top_effects_path

    def get_example_file(self):
        # Function that grabs a file preferring the all effects file. I assume
        # that the format of both files is the same so it doesn't matter which one I take.
        if self.type == "discovery":
            first_file = self.top_effects_path
            second_file = self.all_effects_path
        elif self.type == "replication":
            first_file = self.all_effects_path
            second_file = self.top_effects_path
        else:
            return None

        if first_file is not None:
            fpath = first_file
        elif second_file is not None:
            fpath = second_file
        else:
            return None

        wildcards, tmp_fpath = self.replace_wildcards(fpath=fpath)
        if "<CHR>" in wildcards or "<BATCH>" in wildcards:
            example_fpath, _, _ = self.search_chr_batch_fpath(fpath=tmp_fpath)
        else:
            example_fpath, _, _ = self.get_existing_fpath(fpath=tmp_fpath)

        return example_fpath

    def set_effects_info(self):
        # Extract the header and first line of an example file to automatically determine the format
        # of certain columns. This is needed since methods like mbQTL have their 'Gene' and 'SNP' based
        # on the input files meaning they can be different depending what input data you used to run the
        # method. To make this more dynamic I try to determine if the format on the fly based on the first
        # line of a given file.
        example_file = self.get_example_file()
        if example_file is None:
            print("Error, could not find any data file(s).")
            exit()
        self.effects_header, self.effects_line = self.get_file_header_and_first_line(example_file)

    def get_cell_type(self):
        return self.cell_type

    def load_gene_translate_dicts(self):
        if not os.path.exists(self.gene_translate_path):
            print("Error, gene translate file does not exist.")
            return

        df = pd.read_csv(self.gene_translate_path, sep="\t", header=None, index_col=None)
        if df.iloc[:, 0].nunique() != df.shape[0]:
            print("Error, ensembl IDs should be unique.")
            exit()
        self.ensembl_to_hgnc = dict(zip(df.iloc[:, 0], df.iloc[:, 1]))

        if df.iloc[:, 1].nunique() != df.shape[0]:
            mask = df.iloc[:, 1].duplicated()
            print("Warning, gene translate had {} duplicate HGNC symbols. Removing all duplicates.".format(mask.sum()))
            df = df.loc[~mask, :]

        self.hgnc_to_ensembl = dict(zip(df.iloc[:, 1], df.iloc[:, 0]))

        del df

    def get_ensembl_to_hgnc(self):
        if self.ensembl_to_hgnc is None:
            self.load_gene_translate_dicts()

        return self.ensembl_to_hgnc

    def get_hgnc_to_ensembl(self):
        if self.hgnc_to_ensembl is None:
            self.load_gene_translate_dicts()

        return self.hgnc_to_ensembl

    def get_id(self):
        if self.cell_type is None or self.cell_type == "":
            return self.get_name()
        return self.name + "_" + self.cell_type

    def get_allow_infer(self):
        return self.allow_infer

    def get_n(self):
        return self.n

    def get_ignore_tabix(self):
        return self.ignore_tabix

    def get_na(self):
        return self.na

    def select_columns(self):
        # TODO: would be good to build this logic in for the get_top_effects and get_all_effects functions as well
        if self.type == "discovery":
            return self.get_discovery_columns()
        elif self.type == "replication":
            return self.get_replication_columns()
        else:
            print("Error, unexpected type '{}' in {}".format(self.type, self.name))
            exit()

    def get_discovery_columns(self):
        return self.columns

    def get_replication_columns(self):
        return self.columns

    def get_columns(self):
        return self.columns

    def compile_columns(self):
        columns = self.get_columns()
        compiled_columns = {}
        for column, query in columns.items():
            compiled_columns[column] = self.compile_regex(query=query)

        self.columns = compiled_columns

    @staticmethod
    def compile_regex(query):
        compiled_query = []
        for (column, pattern, suffix) in query:
            if isinstance(pattern, str):
                pattern = re.compile(pattern)
            compiled_query.append((column, pattern, suffix))
        return compiled_query

    def get_column(self, column):
        if column not in self.columns:
            return self.na
        return self.columns[column]

    def has_column(self, column):
        return self.get_column(column) != self.na

    def contains_data(self, label):
        # This method determines if a method has certain info column taking into
        # account that certain info columns can be deduced or inferred from others.
        if self.has_column(label):
            return True

        # We should only translate the discovery gene names to prevent issues later on.
        if label == "gene_ensembl" and self.contains_data("gene_hgnc") and self.gene_translate_path is not None:
            if self.type == "replication":
                print("Warning, translating HGNC to Ensembl gene names is not possible on the fly for the replication dataset, please use 'hgnc' for --gene.")
                return False
            return True
        if label == "gene_hgnc" and self.contains_data("gene_ensembl") and self.gene_translate_path is not None:
            if self.type == "replication":
                print("Warning, translating Ensembl to HGNC gene names is not possible on the fly for the replication dataset, please use 'ensembl' for --gene.")
                return False
            return True

        # In some cases the data might not have a certain column but it can be deduced from other columns.
        if label == "OA" and self.contains_data("alleles") and self.contains_data("EA"):
            return True

        if label == "bonferroni_pvalue" and self.contains_data("nominal_pvalue") and self.contains_data("n_tests"):
            return True

        if label == "zscore":
            if self.contains_data("beta") and self.contains_data("nominal_pvalue"):
                return True
            elif self.contains_data("beta") and self.contains_data("MAF") and self.contains_data("N"):
                if not self.allow_infer:
                    print("\tWarning, zscore could be inferred from beta, MAF, and N info. Turn on --allow_infer to allow this estimation.")
                    return False
                return True

        if label == "N" and (self.n is not None or (self.contains_data("MAN") and self.contains_data("MAF"))):
            return True

        if label == "MAF" and self.contains_data("AF"):
            return True

        if label == "beta_se" and self.contains_data("zscore") and self.contains_data("MAF") and self.contains_data("N"):
            if not self.allow_infer:
                print("\tWarning, beta_se could be inferred from zscore, MAF, and N info. Turn on --allow_infer to allow this estimation.")
                return False
            return True

        if label == "nominal_pvalue":
            if self.contains_data("zscore"):
                return True
            elif self.contains_data("beta") and self.contains_data("MAF") and self.contains_data("N"):
                if not self.allow_infer:
                    print("\tWarning, nominal_pvalue could be inferred from beta, MAF, and N info. Turn on --allow_infer to allow this estimation.")
                    return False
                return True

        return False

    def get_fpath(self, fpath):
        # TODO: see if this can be faster

        # This function searches for a file to see if it exists. If he gzipped / not gzipped version
        # exists of the file then it will return that one automatically.
        if fpath is None:
            return None

        wildcards, tmp_fpath = self.replace_wildcards(fpath=fpath)

        # If CHR or BATCH is in the input file that means we have multiple files
        # that need to be combined to create the full 'all' or 'top' file.
        if "<CHR>" in wildcards or "<BATCH>" in wildcards:
            existing_fpath, fix, n_fpaths = self.search_chr_batch_fpath(fpath=tmp_fpath)
        else:
            existing_fpath, fix, n_fpaths = self.get_existing_fpath(fpath=tmp_fpath)

            # Handling edge case here where if there was an asterix in the fpath
            # which matches exactly one file we will just assume the user
            # wanted to use that file alone and therefore return it.
            if n_fpaths == 1:
                return existing_fpath

        if existing_fpath is None:
            return None

        # Applying to fix before returning the original file path to deal with
        # file that were expected to be gzipped or not.
        return (fpath + fix if not fix.startswith("-") else fpath.rstrip(fix.lstrip("-"))).replace("<CT>", self.cell_type)

    def replace_wildcards(self, fpath):
        wildcards = re.findall("<[A-Z]+>", fpath)
        tmp_fpath = fpath
        for wildcard in wildcards:
            if wildcard == "<CT>":
                # Cell type wildcards can simply be replaced.
                tmp_fpath = fpath.replace("<CT>", self.cell_type)
            elif wildcard == "<CHR>" or wildcard == "<BATCH>":
                # These cannot be handled here but need to be iterated over.
                continue
            else:
                # All other wildcards I simply replace with * to
                # make the input more generic.
                tmp_fpath = tmp_fpath.replace(wildcard, "*")

        return wildcards, tmp_fpath


    def search_chr_batch_fpath(self, fpath):
        # Loop through the chromosomes and batches and return the first existing file.
        # I do not check if all chromosomes are there since I cannot assume which ones
        # were supposed to be there.
        for chromosome in CHROMOSOMES:
            chr_fpath = fpath
            if "<CHR>" in chr_fpath:
                chr_fpath = chr_fpath.replace("<CHR>", chromosome)

            for batch in BATCHES:
                chr_batch_fpath = chr_fpath
                if "<BATCH>" in chr_batch_fpath:
                    chr_batch_fpath = chr_batch_fpath.replace("<BATCH>", batch)

                existing_fpath_info = self.get_existing_fpath(fpath=chr_batch_fpath)
                if existing_fpath_info[0] is not None:
                     return existing_fpath_info

        return None, None, None


    def get_existing_fpath(self, fpath):
        if "*" in fpath:
            return self.get_existing_fpath_wildcard(fpath=fpath)
        return (*self.get_existing_fpath_default(fpath=fpath), None)

    @staticmethod
    def get_existing_fpath_wildcard(fpath):
        # Use glob to find all files matching the asterix pattern. Only return the first.
        # Also, handle adding or removing .gz suffix.
        fpaths = glob.glob(fpath)
        if len(fpaths) > 0:
            return fpaths[0], "", len(fpaths)

        fpaths = glob.glob(fpath + ".gz")
        if len(fpaths) > 0:
            return fpaths[0], ".gz", len(fpaths)

        fpaths = glob.glob(fpath.rstrip(".gz"))
        if len(fpaths) > 0:
            return fpaths[0], "-.gz", len(fpaths)

        return None, None, None

    @staticmethod
    def get_existing_fpath_default(fpath):
        # Handle adding or removing .gz suffix.
        if os.path.exists(fpath):
            return fpath, ""

        if os.path.exists(fpath + ".gz"):
            return fpath + ".gz", ".gz"

        if os.path.exists(fpath.rstrip(".gz")):
            return fpath.rstrip(".gz"), "-.gz"

        return None, None

    def get_eqtl_effects(self, df, gene, snp):
        genes = df.apply(lambda row: self.extract_info(data=row, query=self.get_column(gene)), axis=1)
        snps = df.apply(lambda row: self.extract_info(data=row, query=self.get_column(snp)), axis=1)

        # Save this as {gene: {snp1, snp2}} to allow for multiple snps per gene to be extracted.
        effects = {}
        n_effects = 0
        n_duplicates = 0
        for gene, snp in zip(genes, snps):
            if gene is None or snp is None:
                continue

            if gene not in effects:
                effects[gene] = set()

            if snp in effects[gene]:
                print("Warning, {} - {} is duplicated.".format(gene, snp))
                n_duplicates += 1
                continue

            effects[gene].add(snp)
            n_effects += 1

        return effects, n_effects, n_duplicates

    @staticmethod
    def get_effect_snps(effects):
        all_snps = set()
        for snps in effects.values():
            all_snps.update(snps)
        return all_snps

    def create_dict_from_df(self, df, key, value):
        keys = df.apply(lambda row: self.extract_info(data=row, query=key), axis=1)
        values = df.apply(lambda row: self.extract_info(data=row, query=value), axis=1)

        # To deal with failed regex patterns.
        dictionary = {}
        for key, value in zip(keys, values):
            if key is None or value is None:
                continue
            dictionary[key] = value

        return dictionary

    def extract_info(self, data, query):
        info = ""
        for (column, pattern, suffix) in query:
            # Option 0: Skip if missing.
            if column is None:
                continue

            # Edge case where the index is renamed to 'Unnamed: 0' by pandas.
            if column == '' and 'Unnamed: 0' in data:
                column = 'Unnamed: 0'

            # Skip if column does not exist.
            if isinstance(column, str) and column not in data:
                continue

            # Save full column value.
            if pattern is None:
                info += str(data[column])
            elif isinstance(pattern, str) or isinstance(pattern, re.Pattern):
                # Extract part of column value based on regex pattern.
                try:
                    info += re.match(pattern, str(data[column])).group(1)
                except AttributeError:
                    if self.verbose:
                        print("Warning, pattern did not match in extract_info():\tre.match({}, {}).group(1). Ignoring this row.".format(pattern, str(data[column])))
                    return None
            elif isinstance(pattern, dict):
                # Translate column value based on some dictionary. Skip if it is not in the dictionary.
                if not data[column] in pattern:
                    return None
                info += pattern[data[column]]
            else:
                print("Error, unexpected pattern '{}' ({}) in extract_info()".format(pattern, type(pattern)))
                exit()

            # Add the suffix.
            if suffix is not None:
                info += suffix

        # Replace empty string with None.
        if info == "":
            info = None

        return info

    def trans_label_to_index(self, inpath, label_dict):
        # Function to translate column labels to column indices.
        colname_to_index_dict = self.get_file_colname_to_index(inpath)
        return self.update_label_dict(label_dict=label_dict, colname_to_index_dict=colname_to_index_dict)

    def get_file_colname_to_index(self, inpath, header=0, sep="\t"):
        # Function to create a label to index dict from a file header.
        colnames = self.get_file_header(inpath=inpath, header=header, sep=sep)
        return dict(zip(colnames, range(len(colnames))))

    @staticmethod
    def get_file_head(inpath, header=0, nrows=0, sep="\t"):
        # Function that loads in a file and stored the header and/or the
        # first n number of lines.
        if inpath is None or not os.path.exists(inpath):
            raise FileNotFoundError("Warning, '{}' file not found".format(inpath))

        if inpath.endswith(".gz"):
            fhin = gzip.open(inpath, 'rt')
        else:
            fhin = open(inpath, 'r')

        colnames = None
        data = []
        for i, line in enumerate(fhin):
            if (i > header) and (i > nrows):
                break
            values = line.rstrip("\n").split(sep)

            if i == header:
                colnames = values
                continue

            data.append(values)
        fhin.close()

        if colnames is None and len(data) == 0:
            print("Error, empty file.")
            exit()

        return colnames, data

    def get_file_header(self, inpath, header=0, sep="\t"):
        header, _ = self.get_file_head(inpath=inpath, header=header, sep=sep)
        return header

    def get_first_line(self, inpath, header=0, sep="\t"):
        header, line = self.get_file_head(inpath=inpath, header=header, nrows=1, sep=sep)
        return dict(zip(header, line[0]))

    def get_file_header_and_first_line(self, inpath, header=0, sep="\t"):
        header, line = self.get_file_head(inpath=inpath, header=header, nrows=1, sep=sep)
        return header, dict(zip(header, line[0]))

    @staticmethod
    def update_label_dict(label_dict, colname_to_index_dict):
        # Update the keys of a {colname: [None, None, None} dict based on another dict.
        indices_dict = {}
        for argument, query in label_dict.items():
            index_query = []
            for (column, pattern, sep) in query:
                index_query.append((colname_to_index_dict[column], pattern, sep))
            indices_dict[argument] = index_query
            del index_query

        return indices_dict

    def get_all_effects(self):
        # Function that loads in all effects from the 'all_effects_path' file.
        try:
            return self.get_effects_wrapper(
                inpath=self.get_all_effects_path(),
                func=self.extract_all_effects
            )
        except FileNotFoundError as e:
            print(e)
            print("Failed to load 'all_effects_path'. Selecting top effects from 'all_effects_path' instead.")

        try:
            return self.get_effects_wrapper(
                inpath=self.get_top_effects_path(),
                func=self.extract_all_effects
            )
        except FileNotFoundError as e:
            print(e)
            print("Failed to load 'top_effects_path'.")

        return None

    def extract_all_effects(self, inpath, **kwargs):
        # inpath = all_effects_path
        return self.load_file(inpath, **kwargs)

    def get_top_effects(self, gene=None, snp=None):
        # Function that loads in the top effects from the 'top_effects_path'. If this does not exist,
        # select the top effects from the 'all_effects_path' based on nominal p-value.
        # TODO: this will fail if a dataset does not have nominal p-values in their 'all_effects_path'.. would anyone do that?
        #  Probably not right?
        try:
            return self.get_effects_wrapper(
                inpath=self.get_top_effects_path(),
                func=self.extract_all_effects
            )
        except FileNotFoundError as e:
            print(e)
            print("Failed to load 'top_effects_path'. Selecting top effects from 'all_effects_path' instead.")

        # Making sure the required columns are available.
        # TODO: potentially this could still work if a dataset has no nominal p-values but it
        #   does have z-scores. However, I need to rewrite the load_partial_file function to
        #   pick the effects with the highest zscore instead of the lowest p-value...
        for label in [gene, "nominal_pvalue"]:
            if not self.contains_data(label=label):
                print("Error, get_top_effects() from the all effects file is unavailable for {} since"
                      " there is no '{}' data available.".format(self.class_name, label))
                exit()

        try:
            return self.get_effects_wrapper(
                inpath=self.get_all_effects_path(),
                label_dict={"key": self.get_column(gene), "value": self.get_column("nominal_pvalue")},
                func=self.extract_top_effects
            )
        except FileNotFoundError as e:
            print(e)
            print("Failed to load 'all_effects_path'.")

        return None

    def extract_top_effects(self, inpath, label_dict=None, top_entries_cols=None, **kwargs):
        # inpath = all_effects_path
        if top_entries_cols is None and label_dict is not None:
            # Since this function works on column indices instead of column names
            # we first need to translate the labels to indices.
            top_entries_cols = self.trans_label_to_index(
                inpath=inpath,
                label_dict=label_dict
            )

        # Next we can start loading the file.
        df = self.load_partial_file(
            inpath,
            top_entries_cols=top_entries_cols,
            ignore_tabix=self.ignore_tabix,
            **kwargs
        )
        return df

    def get_specific_effects(self, effects=None, gene=None, snp=None):
        # Making sure the required columns are available.
        for label in [gene, snp]:
            if not self.contains_data(label=label):
                print("Error, get_specific_effects() from the all effects file is unavailable for {} since"
                      " there is no '{}' data available.".format(self.class_name, label))
                exit()

        try:
            # This function loads in specific gene-SNP combos from the 'all_effects_path' file.
            return self.get_effects_wrapper(
                inpath=self.get_all_effects_path(),
                label_dict={"key": self.get_column(gene), "value": self.get_column(snp)},
                func=self.extract_specific_effects,
                effects=effects
            )
        except FileNotFoundError as e:
            print(e)
            print("Failed to load 'all_effects_path'. Selecting specific effects from 'top_effects_path' instead.")

        try:
            # Using the top file instead.
            return self.get_effects_wrapper(
                inpath=self.get_top_effects_path(),
                label_dict={"key": self.get_column(gene), "value": self.get_column(snp)},
                func=self.extract_specific_effects,
                effects=effects
            )
        except FileNotFoundError as e:
            print(e)
            print("Failed to load 'top_effects_path'.")

        return None

    def extract_specific_effects(self, inpath, label_dict=None, specific_entries_cols=None, effects=None, **kwargs):
        # inpath = all_effects_inpath
        if specific_entries_cols is None and label_dict is not None:
            # Since this function works on column indices instead of column names
            # we first need to translate the labels to indices.
            specific_entries_cols = self.trans_label_to_index(
                inpath=inpath,
                label_dict=label_dict
            )

        # Next we can start loading the file.
        df = self.load_partial_file(
            inpath,
            effects=effects,
            specific_entries_cols=specific_entries_cols,
            ignore_tabix=self.ignore_tabix,
            **kwargs
        )
        return df

    def get_effects_wrapper(self, inpath, func, *args, **kwargs):
        # This function wraps around extract_(all/top/specific)effects() to allow for
        # processing of these files per chromosome and/or per batch.
        if inpath is None:
            raise FileNotFoundError("Warning, '{}' file not found".format(inpath))

        has_multiple_chr = "<CHR>" in inpath
        has_multiple_batches = "<BATCH>" in inpath

        # First check if the input files consists of multiple
        # chromosomes and / or batches. If so we can keep it simple.
        if not has_multiple_chr and not has_multiple_batches:
            return func(inpath=inpath, *args, **kwargs)

        # Combine the multiple files into one big data frame.
        df_list = []
        for chromosome in CHROMOSOMES:
            batch_start = False
            for batch in BATCHES:
                full_inpath, _ = self.get_existing_fpath_default(os.path.join(self.path, inpath.replace("<CHR>", chromosome).replace("<BATCH>", batch)))

                if full_inpath is None:
                    if not batch_start:
                        # Here we allow batches to start counting at 0 but also at 1.
                        # as a result we might iterate through all N max batches
                        # if a chromosome has 0 batches in total but that doesn't
                        # take any time anyway.
                        continue
                    else:
                        break

                df_list.append(func(inpath=full_inpath, *args, **kwargs))
                batch_start = True

                if not has_multiple_batches:
                    break

            if not has_multiple_chr:
                break

        if len(df_list) == 0:
            return None

        return pd.concat(df_list, axis=0)

    @staticmethod
    def load_file(inpath, header=0, index_col=None, sep="\t", low_memory=True,
                  nrows=None, skiprows=None, usecols=None):
        if inpath is None or not os.path.exists(inpath):
            raise FileNotFoundError("Warning, '{}' file not found".format(inpath))

        print("Loading file...")
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows,
                         usecols=usecols)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def load_h5(inpath):
        if inpath is None or not os.path.exists(inpath):
            raise FileNotFoundError("Warning, '{}' file not found".format(inpath))

        print("Loading file...")
        f = h5py.File(inpath, 'r')
        print("\tLoaded h5 file: {} "
              "with {} keys".format(os.path.basename(inpath),
                                      len(f.keys())))
        return f

    @staticmethod
    def load_excel(inpath, sheet_name, skiprows=None, dtype=None):
        if inpath is None or not os.path.exists(inpath):
            raise FileNotFoundError("Warning, '{}' file not found".format(inpath))

        print("Loading file...")
        df = pd.read_excel(inpath, sheet_name=sheet_name, skiprows=skiprows, dtype=dtype)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    def load_partial_file(self, inpath, header=0, index_col=None, sep="\t", nrows=None,
                          skiprows=None, usecols=None, top_entries_cols=None,
                          specific_entries_cols=None, effects=None, ignore_tabix=False):
        if inpath is None or not os.path.exists(inpath):
            raise FileNotFoundError("Warning, '{}' file not found".format(inpath))

        if os.path.exists(inpath + ".tbi") and not ignore_tabix:
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

        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    def load_partial_file_w_tabix(self, inpath, header=0, index_col=None, sep="\t", usecols=None,
                                  top_entries_cols=None, specific_entries_cols=None, effects=None):
        # TODO: this class now only works by tabix query on the SNP position instead of the gene position.
        #  Not sure how I am supposed to know the gene chr:pos from the HGNC symbol / ensembl ID without using API's
        #  or first parsing the file line by line anyway (that info is in the line). So unless you indexed on SNP pos
        #  this function is a lot slower than the line by line one and you would be better off using that one.
        #  Note to self, if I do somehow get the gene position info then I should could the number of SNPs within a cis
        #  window and add that as NEntries to make it in line with the line by line code again.
        raise NotImplementedError("This code sucks, just use the line by line function by setting ignore_tabix=True")

        print("Loading partial file with tabix...")
        tb = tabix.open(inpath)
        columns = self.get_file_header(inpath=inpath, header=header, sep=sep)
        lines = []
        indices = None if index_col is None else []
        top_entries = {}
        n_effects = len(effects)
        i = 0
        keys_found = set()
        pairs_found = set()
        for i, snp in enumerate(effects.values()):
            if i % 10 == 0:
                print("\tParsed {:,} / {:,} effects".format(i, n_effects), end='\r')
            chr, pos = snp.split(":")
            try:
                records = tb.query(chr, int(pos), int(pos))
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
                    if key is None or key not in effects:
                        continue
                    keys_found.add(key)

                    if "value" in specific_entries_cols.keys():
                        value = self.extract_info(data=values, query=specific_entries_cols["value"])
                        if value is None or value not in effects[key]:
                            continue
                        pairs_found.add(key + "_" + value)

                if top_entries_cols is None:
                    lines.append(values)
                else:
                    key = self.extract_info(data=values, query=top_entries_cols["key"])
                    value = float(self.extract_info(data=values, query=top_entries_cols["value"]))
                    if key is None or value is None:
                        continue

                    if key not in top_entries:
                        top_entries[key] = [value, values]
                    elif key in top_entries and value < top_entries[key][0]:
                        top_entries[key] = [value, values]
                    else:
                        pass
        print("\tParsed {:,} / {:,} effects".format(i, n_effects))

        if specific_entries_cols is not None and isinstance(effects, dict):
            n_pairs_found = len(pairs_found)
            n_query_pairs = sum([len(value) for value in effects.values()])
            pairs_missing = n_query_pairs - n_pairs_found
            n_keys_found = len(keys_found)
            keys_missing = len(effects) - n_keys_found
            values_missing = pairs_missing - keys_missing
            print("\tFound {:,} / {:,} key-value pairs. Of the {:,} missing pairs {:,} did not match on key and {:,} did not match on key-value".format(n_pairs_found, n_query_pairs, pairs_missing, keys_missing, values_missing))

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
        keys_found = set()
        pairs_found = set()
        for i, line in enumerate(fhi):
            if i % 1e6 == 0:
                print("\tParsed {:,} lines".format(i), end='\r')

            if skiprows is not None and i < skiprows:
                continue
            if max_rows is not None and i > max_rows:
                break

            values = line.strip("\n").split(sep)
            if n_values is None:
                n_values = len(values)
            if len(values) != n_values:
                print("\tError, unequal number of columns in the input file, skip line")
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
                if key is None or key not in effects:
                    continue
                keys_found.add(key)

                if "value" in specific_entries_cols.keys():
                    value = self.extract_info(data=values, query=specific_entries_cols["value"])
                    if value is None or value not in effects[key]:
                        continue
                    pairs_found.add(key + "_" + value)

            if top_entries_cols is None:
                lines.append(values)
            else:
                key = self.extract_info(data=values, query=top_entries_cols["key"])
                value = float(self.extract_info(data=values, query=top_entries_cols["value"]))
                if key is None or value is None:
                    continue

                # TODO: could be nice to add an option here to only select certain type of
                #   variants to be picked as top effect. For example, exclude INDELS.
                if key not in top_entries:
                    top_entries[key] = [value, values, 1]
                elif key in top_entries:
                    top_entries[key][2] += 1
                    if value < top_entries[key][0]:
                        top_entries[key][0:2] = [value, values]
                else:
                    pass
        fhi.close()
        print("\tParsed {:,} lines".format(i))

        if specific_entries_cols is not None and isinstance(effects, dict):
            n_pairs_found = len(pairs_found)
            n_query_pairs = sum([len(value) for value in effects.values()])
            pairs_missing = n_query_pairs - n_pairs_found
            n_keys_found = len(keys_found)
            keys_missing = len(effects) - n_keys_found
            values_missing = pairs_missing - keys_missing
            print("\tFound {:,} / {:,} key-value pairs. Of the {:,} missing pairs {:,} did not match on key and {:,} did not match on key-value".format(n_pairs_found, n_query_pairs, pairs_missing, keys_missing, values_missing))

        if top_entries_cols is not None:
            columns += ["NEntries"]
            for _, (_, values, nentries) in top_entries.items():
                values += [nentries]
                lines.append(values)

        return pd.DataFrame(lines, columns=columns, index=indices)

    def postprocess(self, df):
        if df is None or df.shape[0] == 0:
            df = None
        df = self.update_df(df=df)
        df = self.add_missing_info(df=df)
        return df

    def update_df(self, df):
        return df

    def add_missing_info(self, df):
        if df is None:
            return df
        print("Adding missing info to '{}':".format(self.class_name))

        # Fill in missing gene info.
        if not self.has_column("gene_ensembl") and self.contains_data("gene_hgnc") and self.gene_translate_path is not None:
            print("\tAdding ensembl column")
            df["gene_ensembl"] = self.translate_column(df=df, key_column=self.get_column("gene_hgnc"), translate_dict=self.get_hgnc_to_ensembl())
            self.columns["gene_ensembl"] = [("gene_ensembl", None, None)]
        if not self.has_column("gene_hgnc") and self.contains_data("gene_ensembl") and self.gene_translate_path is not None:
            print("\tAdding HGNC column")
            df["gene_hgnc"] = self.translate_column(df=df, key_column=self.get_column("gene_ensembl"), translate_dict=self.get_ensembl_to_hgnc())
            self.columns["gene_hgnc"] = [("gene_hgnc", None, None)]

        # Fill in the other allele.
        if not self.has_column("OA") and self.contains_data("OA"):
            print("\tAdding OA column")
            df["OA"] = self.impute_other_allele(df=df, alleles_column=self.get_column("alleles"), effect_allele=self.get_column("EA"))
            self.columns["OA"] = [("OA", None, None)]

        # Fill in the bonforroni corrected pvalue. Only usefull for the discovery.
        if self.type == "discovery" and not self.has_column("bonferroni_pvalue") and self.contains_data("bonferroni_pvalue"):
            print("\tAdding bonferroni_pvalue column")
            df["bonferroni_pvalue"] = self.calc_bonferroni_pvalue(df=df, pvalue_column=self.get_column("nominal_pvalue"), ntests_column=self.get_column("n_tests"))
            self.columns["bonferroni_pvalue"] = [("bonferroni_pvalue", None, None)]

        # Fill in N.
        if not self.has_column("N") and self.contains_data("N"):
            if self.n is not None:
                print("\tAdding N column")
                df["N"] = self.n
            else:
                print("\tAdding N column")
                df["N"] = self.calc_n_from_man(df=df, man_column=self.get_column("MAN"), maf_column=self.get_column("MAF"))
            self.columns["N"] = [("N", None, None)]

        # Fill in the MAF.
        if not self.has_column("MAF") and self.contains_data("MAF"):
            print("\tAdding MAF column")
            df["MAF"] = self.calc_maf(df=df, af_column=self.get_column("AF"))
            self.columns["MAF"] = [("MAF", None, None)]

        # Fill in the zscore.
        if not self.has_column("zscore") and self.contains_data("zscore"):
            if self.contains_data("nominal_pvalue"):
                print("\tAdding zscore column")
                df["zscore"] = self.calc_zscore_from_pvalue(df=df, beta_column=self.get_column("beta"), pvalue_column=self.get_column("nominal_pvalue"))
            else:
                print("\tAdding zscore column (infered)")
                # TODO: this code works but is not perfect. It would be better to look up the nominal p-values to calculate
                #  the z-score but if that info isn't available or takes too long... best I can do
                df["zscore"] = self.calc_zscore_from_beta(df=df, beta_column=self.get_column("beta"), maf_column=self.get_column("MAF"), n_column=self.get_column("N"))
            self.columns["zscore"] = [("zscore", None, None)]

        # Fill in the beta standard error.
        if not self.has_column("beta_se") and self.contains_data("beta_se"):
            print("\tAdding beta_se column (inferred)")
            df["beta_se"] = self.calc_beta_se(df=df, zscore_column=self.get_column("zscore"), maf_column=self.get_column("MAF"), n_column=self.get_column("N"))
            self.columns["beta_se"] = [("beta_se", None, None)]

        # Fill in the nominal p-value.
        if not self.has_column("nominal_pvalue") and self.contains_data("nominal_pvalue"):
            # We know zscore is there because it is checked above.
            print("\tAdding nominal_pvalue column")
            df["nominal_pvalue"] = self.calc_pvalue_from_zscore(df=df, zscore_column=self.get_column("zscore"))
            self.columns["nominal_pvalue"] = [("nominal_pvalue", None, None)]

        print("")
        return df

    def translate_column(self, df, key_column, translate_dict):
        return df.apply(lambda row: self.translate_value(row=row, key_column=key_column, translate_dict=translate_dict), axis=1)

    def translate_value(self, row, key_column, translate_dict):
        value = self.extract_info(data=row, query=key_column)
        if value is None or value not in translate_dict:
            return np.nan

        return translate_dict[value]

    def impute_other_allele(self, df, alleles_column, effect_allele):
        return df.apply(lambda row: self.extract_other_allele(row=row, alleles_column=alleles_column, effect_allele_column=effect_allele), axis=1)

    def extract_other_allele(self, row, alleles_column, effect_allele_column):
        alleles = self.extract_info(data=row, query=alleles_column)
        if alleles is None:
            return np.nan
        allele1, allele2 = alleles.split("/")

        effect_allele = self.extract_info(data=row, query=effect_allele_column)
        if effect_allele is None or (effect_allele != allele1 and effect_allele != allele2):
            return np.nan
        elif effect_allele == allele1:
            return allele2
        else:
            return allele1

    def calc_bonferroni_pvalue(self, df, pvalue_column, ntests_column):
        return df.apply(lambda row: self.bonferroni_pvalue(row=row, pvalue_column=pvalue_column, ntests_column=ntests_column), axis=1)

    def bonferroni_pvalue(self, row, pvalue_column, ntests_column):
        try:
            nominal_pvalue = float(self.extract_info(data=row, query=pvalue_column))
            ntests = float(self.extract_info(data=row, query=ntests_column))
        except ValueError:
            return np.nan

        return min(nominal_pvalue * ntests, 1)

    def calc_n_from_man(self, df, man_column, maf_column):
        return df.apply(lambda row: self.man_to_n(row=row, man_column=man_column, maf_column=maf_column), axis=1)


    def man_to_n(self, row, man_column, maf_column):
        try:
            man = int(self.extract_info(data=row, query=man_column))
            maf = float(self.extract_info(data=row, query=maf_column))
        except ValueError:
            return np.nan

        n = int(man / maf)

        return n

    def calc_zscore_from_pvalue(self, df, beta_column, pvalue_column):
        return df.apply(lambda row: self.beta_pvalue_to_zscore(row=row, beta_column=beta_column, pvalue_column=pvalue_column), axis=1)

    def beta_pvalue_to_zscore(self, row, beta_column, pvalue_column):
        try:
            nominal_pvalue = float(self.extract_info(data=row, query=pvalue_column))
            beta = float(self.extract_info(data=row, query=beta_column))
        except ValueError:
            return np.nan

        if nominal_pvalue < 0 or nominal_pvalue > 1:
            print("Error in calc_zscore(), p-value outside range.")
            exit()

        if nominal_pvalue >= 0.9999999999999999:
            zscore = -1.3914582123358836e-16
        elif nominal_pvalue <= 1e-323:
            zscore = -38.467405617144344
        else:
            zscore = stats.norm.ppf(nominal_pvalue / 2)

        if beta > 0:
            zscore *= -1.

        return zscore

    def calc_pvalue_from_zscore(self, df, zscore_column):
        return df.apply(lambda row: self.zscore_to_pvalue(row=row, zscore_column=zscore_column), axis=1)

    def zscore_to_pvalue(self, row, zscore_column):
        """
        non-zero output range is between -37.677121 and 37.677121.
        """
        try:
            zscore = float(self.extract_info(data=row, query=zscore_column))
        except ValueError:
            return np.nan

        pvalue = stats.norm.cdf(-abs(zscore)) * 2

        return pvalue

    def calc_zscore_from_beta(self, df, beta_column, maf_column, n_column):
        return df.apply(lambda row: self.beta_to_zscore(row=row, beta_column=beta_column, maf_column=maf_column, n_column=n_column), axis=1)

    def beta_to_zscore(self, row, beta_column, maf_column, n_column):
        # This is a rewritten form of the zscore_to_beta_se() function
        # to calculate the zscore from a given beta, maf and n.
        try:
            beta = float(self.extract_info(data=row, query=beta_column))
            maf = float(self.extract_info(data=row, query=maf_column))
            n = float(self.extract_info(data=row, query=n_column))
        except ValueError:
            return np.nan

        return ((-1j * 2 ** (1 / 2) * beta * (-1 + maf) ** (1 / 2) * maf ** (1 / 2) * n ** (1 / 2)) / (1 - 2 * beta ** 2 * maf + 2 * beta ** 2 * maf ** 2) ** (1 / 2)).real

    def calc_maf(self, df, af_column):
        return df.apply(lambda row: self.af_to_maf(row=row, af_column=af_column), axis=1)

    def af_to_maf(self, row, af_column):
        try:
            af = float(self.extract_info(data=row, query=af_column))
        except ValueError:
            return np.nan
        return min(af, 1 - af)

    def calc_beta_se(self, df, zscore_column, maf_column, n_column):
        return df.apply(lambda row: self.zscore_to_beta_se(row=row, zscore_column=zscore_column, maf_column=maf_column, n_column=n_column), axis=1)

    def zscore_to_beta_se(self, row, zscore_column, maf_column, n_column):
        try:
            zscore = float(self.extract_info(data=row, query=zscore_column))
            maf = float(self.extract_info(data=row, query=maf_column))
            n = float(self.extract_info(data=row, query=n_column))
        except ValueError:
            return np.nan

        a = 2 * maf * (1 - maf) * (n + (zscore * zscore))
        # derived_beta = zscore / a ** (1/2)
        derived_beta_se = 1 / a ** (1 / 2)
        return derived_beta_se

    def standardize_format(self, df, gene, snp, include_fdr=True):
        if df is None:
            return None

        # Transform the input data into a standardized format with the standard
        # column names for each type of information.
        df_info = []
        columns = []
        for i, (_, row) in enumerate(df.iterrows()):
            row_info = []
            for argument, query in self.columns.items():
                row_info.append(self.extract_info(data=row, query=query))
                if i == 0:
                    columns.append(self.name + " " + argument)
            df_info.append(row_info)
        standard_df = pd.DataFrame(df_info, columns=columns).replace("NA", np.nan).replace("nan", np.nan).replace("-", np.nan).dropna(axis=1, how='all').dropna(subset=[self.name + " " + gene, self.name + " " + snp])

        dtypes = {
            self.name + " gene_hgnc": str,
            self.name + " gene_ensembl": str,
            self.name + " SNP_rsid": str,
            self.name + " SNP_chr:pos": str,
            self.name + " alleles": str,
            self.name + " EA": str,
            self.name + " OA": str,
            self.name + " beta": float,
            self.name + " beta_se": float,
            self.name + " n_tests": int,
            self.name + " nominal_pvalue": float,
            self.name + " permuted_pvalue": float,
            self.name + " bonferroni_pvalue": float,
            self.name + " zscore": float,
            self.name + " FDR": float,
            self.name + " N": int,
            self.name + " AF": float,
            self.name + " MAN": int,
            self.name + " MAF": float
        }

        # Change the data types.
        for column, dtype in dtypes.items():
            if column not in standard_df.columns:
                continue
            # this needs to be done in two steps to deal with things like '180.0'.
            if dtype == int:
                standard_df[column] = standard_df[column].astype(float)

                # This only works if there are no NaN values.
                if standard_df[column].isna().sum() == 0:
                    standard_df[column] = standard_df[column].astype(int)
            else:
                standard_df[column] = standard_df[column].astype(dtype)

        # Check for p-values out of range.
        for pval_column in ["nominal_pvalue", "permuted_pvalue", "bonferroni_pvalue"]:
            if self.name + " " + pval_column not in standard_df.columns:
                continue
            invalid_pval_range_mask = (standard_df[self.name + " " + pval_column] < 0) | (standard_df[self.name + " " + pval_column] > 1)
            if invalid_pval_range_mask.sum() > 0:
                standard_df = standard_df.loc[~invalid_pval_range_mask, :]
                print(f"Warning: {invalid_pval_range_mask.sum():,} rows were removed since the value in column {pval_column} is outside of a normal p-value range.")

        # Remove FDR if we do not want it (in the case of replication).
        if self.name + " FDR" in standard_df and not include_fdr:
            standard_df.drop([self.name + " FDR"], axis=1, inplace=True)

        standard_df.index = standard_df[self.name + " " + gene] + "_" + standard_df[self.name + " " + snp]
        return standard_df

    def __str__(self):
        return "{} dataset:\n  class_name = {}\n  name = {}\n  path = {}\n  " \
               "all_filename = {}\n  top_filename = {}\n  cell_type = {}\n  " \
               "class_settings = {}\n  allow_infer = {}\n  verbose = {}\n  " \
               "all_effects_path = {}\n  top_effects_path = {}\n  n = {}\n  " \
               "ignore_tabix = {}\n  na = {}\n  columns = {}\n".format(
            self.type.title(),
            self.class_name,
            self.name,
            self.path,
            self.all_filename,
            self.top_filename,
            self.cell_type,
            self.class_settings,
            self.allow_infer,
            self.verbose,
            self.all_effects_path,
            self.top_effects_path,
            self.n,
            self.ignore_tabix,
            self.na,
            "\n    " + "\n    ".join(["{} = {}".format(key, value) for key, value in self.columns.items() if value != self.na]),
        )

##############################################################################################################

class CUSTOM(Dataset):
    def __init__(self, *args, **kwargs):
        super(CUSTOM, self).__init__(*args, **kwargs)
        class_settings = self.load_json()
        self.class_name = class_settings["class_name"]

        # Update filenames.
        self.update_all_filename(filename=class_settings["all_effects_filename"])
        self.update_top_filename(filename=class_settings["top_effects_filename"])

        # Set file paths.
        self.set_all_effects_path()
        self.set_top_effects_path()

        # Default sample size.
        self.n = class_settings["n"]

        # Columns that are in the original file.
        self.columns.update(class_settings["columns"])

    def load_json(self):
        if not os.path.exists(self.class_settings):
            raise FileNotFoundError("Error, '{}' class settings does not exist.".format(self.class_settings))

        with open(self.class_settings) as f:
            class_settings = json.load(f)
        f.close()

        if "class_name" not in class_settings:
            class_settings["class_name"] = "CUSTOM"

        if "n" not in class_settings:
            class_settings["n"] = None

        for filenames in ["all_effects_filename", "top_effects_filename"]:
            if filenames not in class_settings:
                class_settings[filenames] = ""

        if "columns" not in class_settings:
            class_settings["columns"] = {}

        valid = True
        for column_name in class_settings["columns"].keys():
            column_value = class_settings["columns"][column_name]
            if isinstance(column_value, str):
                class_settings["columns"][column_name] = [(column_value, None, None)]
            elif isinstance(column_value, list):
                for part_column_value in column_value:
                    # Make sure each value in the list has 3 elements.
                    if len(part_column_value) != 3:
                        print("Error in {}, part CUSTOM class setting {} is invalid. Expecting list with length 3 but '{}' has length '{}'.".format(self.name,str(column_value), str(part_column_value), len(part_column_value)))
                        valid = False

                    # Make sure the first element is a string.
                    if not isinstance(part_column_value[0], str):
                        print("Error in {}, part of CUSTOM class setting '{}' is invalid. Expecting string on the first index but '{}' has value '{}'.".format(self.name,str(column_value), str(part_column_value), str(part_column_value[0])))
                        valid = False

                    # Make sure the second element is a string or a dict.
                    if not part_column_value[1] is None and not isinstance(part_column_value[1], str) and not isinstance(part_column_value[1], dict):
                        print("Error in {}, part of CUSTOM class setting '{}' is invalid. Expecting null or a string on the second index but '{}' has value '{}'.".format(self.name,str(column_value), str(part_column_value), str(part_column_value[1])))
                        valid = False

                    # Make sure the third element is a string.
                    if not part_column_value[2] is None and not isinstance(part_column_value[2], str):
                        print("Error in {}, part of CUSTOM class setting '{}' is invalid. Expecting null or a string on the third index but '{}' has value '{}'.".format(self.name, str(column_value), str(part_column_value), str(part_column_value[2])))
                        valid = False
            else:
                print("Error in {}, CUSTOM class setting {} is invalid. Expecting string or list.".format(self.name, str(column_value)))
                valid = False

        if not valid:
            exit()

        return class_settings


##############################################################################################################


class LIMIX(Dataset):
    def __init__(self, *args, **kwargs):
        super(LIMIX, self).__init__(*args, **kwargs)
        self.class_name = "LIMIX"

        # Update filenames.
        self.update_all_filename(filename="qtl_results_all.txt")
        self.update_top_filename(filename="top_qtl_results_all.txt")

        # Set file paths.
        self.set_all_effects_path()
        self.set_top_effects_path()

        # Set effects info.
        self.set_effects_info()

        # Other file paths. Note that <CHUNK> and <FEATURE> are not filled in are remain 'wildcards'.
        self.feature_metadata_file = self.get_fpath(os.path.join(self.path, "qtl", "feature_metadata_<CHUNK>.txt"))
        self.snp_metadata_file = self.get_fpath(os.path.join(self.path, "qtl", "snp_metadata_<CHUNK>.txt"))
        self.qtl_results_file = self.get_fpath(os.path.join(self.path, "qtl", "qtl_results_<CHUNK>.h5"))
        self.snp_qc_metrics = self.get_fpath(os.path.join(self.path, "qtl", "snp_qc_metrics_naContaining_feature_<FEATURE>.txt"))

        # Other class variables.
        self.chunk_pattern = "(([0-9]{1,2}|X|Y|MT)_([0-9]+)_([0-9]+))"
        self.qtl_chunks = None

        # Columns that are in the original file.
        self.columns.update({
            "gene_hgnc": [("feature_id", None, None)],
            "gene_ensembl": [("ENSG", "(ENSG[0-9]+)", None)],
            "SNP_rsid": self.get_snp_rsid_column(),
            "SNP_chr:pos": [("snp_chromosome", "([0-9]{1,2}|X|Y|MT)", ":"), ("snp_position", "([0-9]+)", None)],
            "alleles": self.get_alleles_column(),
            "EA": [("assessed_allele", None, None)],
            # "OA": [(None, None, None)],
            "beta": [("beta", None, None)],
            "beta_se": [("beta_se", None, None)],
            "n_tests": [("nTotalTestsPerFeat", None, None)],
            "nominal_pvalue": [("p_value", None, None)],
            "permuted_pvalue": [("empirical_feature_p_value", None, None)],
            "bonferroni_pvalue": [("p_value_bonf_corr", None, None)],
            # "zscore": [(None, None, None)],
            "FDR": [("global_pvalue", None, None)],
            "N": [("n_samples", None, None)],
            # "AF": [(None, None, None)],
            # "MAN": [(None", None, None)],
            "MAF": [("maf", None, None)]
        })

    def get_snp_rsid_column(self):
        if re.match("rs[0-9]+$", self.effects_line["snp_id"]):
            return [("snp_id", "(rs[0-9]+)", None)]
        else:
            return [(None, None, None)]

    def get_alleles_column(self):
        if re.match("(?:[0-9]{1,2}|X|Y|MT):[0-9]+:[A-Z]+:[A-Z]+", self.effects_line["snp_id"]):
            return [("snp_id", "(?:[0-9]{1,2}|X|Y|MT):[0-9]+:([A-Z]+):[A-Z]+", "/"), ("snp_id", "(?:[0-9]{1,2}|X|Y|MT):[0-9]+:[A-Z]+:([A-Z]+)", None)]
        else:
            return [(None, None, None)]

    def get_all_effects(self, gene=None, snp=None):
        try:
            return self.get_effects_wrapper(
                inpath=self.get_all_effects_path(),
                func=self.extract_all_effects
            )
        except FileNotFoundError as e:
            print(e)
            print("Failed to load 'all_effects_path'. Selecting all effects from 'qtl chunks' instead.")

        try:
            return self.load_qtl_data(gene=gene, snp=snp)
        except FileNotFoundError as e:
            print(e)
            print("Failed to load 'qtl chunks'. Selecting all effects from 'top_effects_path' insetad.")

        try:
            return self.get_effects_wrapper(
                inpath=self.get_top_effects_path(),
                func=self.extract_all_effects
            )
        except FileNotFoundError as e:
            print(e)
            print("Failed to load 'top_effects_path'.")

        return None

    def get_top_effects(self, gene=None, snp=None):
        df = super(LIMIX, self).get_top_effects(gene=gene, snp=snp)
        if df is not None:
            return df

        print("Selecting all effects from 'qtl chunks' instead.")
        try:
            self.load_qtl_data(gene=gene, snp=snp, top=True)
        except FileNotFoundError as e:
            print(e)
            print("Failed to load 'qtl chunks'.")

        return None

    def get_specific_effects(self, effects=None, gene=None, snp=None):
        # Making sure the required columns are available.
        for label in [gene, snp]:
            if not self.contains_data(label=label):
                print("Error, get_specific_effects() from the all effects file is unavailable for {} since"
                      " there is no '{}' data available.".format(self.class_name, label))
                exit()

        try:
            # This function loads in specific gene-SNP combos from the 'all_effects_path' file.
            return self.get_effects_wrapper(
                inpath=self.get_all_effects_path(),
                label_dict={"key": self.get_column(gene), "value": self.get_column(snp)},
                func=self.extract_specific_effects,
                effects=effects
            )
        except FileNotFoundError as e:
            print(e)
            print("Failed to load 'all_effects_path'. Selecting all effects from 'qtl chunks' instead.")

        try:
            return self.load_qtl_data(effects=effects, gene=gene, snp=snp)
        except FileNotFoundError as e:
            print(e)
            print("Failed to load chunk files. Selecting specific effects from 'top_effects_path' instead.")

        try:
            # Using the top file instead.
            return self.get_effects_wrapper(
                inpath=self.get_top_effects_path(),
                label_dict={"key": self.get_column(gene), "value": self.get_column(snp)},
                func=self.extract_specific_effects,
                effects=effects
            )
        except FileNotFoundError as e:
            print(e)
            print("Failed to load 'top_effects_path'.")

        return None

    def load_qtl_data(self, effects=None, gene=None, snp=None, top=False):
        qlt_chunks_df = self.get_qtl_chunks()
        if qlt_chunks_df.shape[0] == 0:
            raise FileNotFoundError("Warning, '{}' chunk files not found".format(self.qtl_results_file.replace("<CHUNK>", "*")))

        df_list = []
        for index, row in qlt_chunks_df.iterrows():
            df = self.load_qtl_chunk(chunk=row["Chunk"], effects=effects, gene=gene, snp=snp, top=top)
            if df is not None:
                df_list.append(df)

        if len(df_list) == 0:
            return None

        return pd.concat(df_list, axis=0)

    def get_qtl_chunks(self):
        if self.qtl_chunks is None:
            self.qtl_chunks = self.load_qtl_chunks()

        return self.qtl_chunks

    def load_qtl_chunks(self):
        qtl_results = glob.glob(self.qtl_results_file.replace("<CHUNK>", "*"))

        # Extract the chunk info: chr_start_end
        qtl_results_data = []
        for filepath in qtl_results:
            match = re.match(self.qtl_results_file.replace("<CHUNK>", self.chunk_pattern), filepath)
            qtl_results_data.append([os.path.basename(filepath), match.group(1), match.group(2), int(match.group(3)), int(match.group(4))])

        df = pd.DataFrame(qtl_results_data, columns=["Filename", "Chunk", "Chromosome", "Start", "End"])
        df.sort_values(by=["Chromosome", "Start"], key=natsort_keygen(), ascending=True, inplace=True)
        return df

    def load_qtl_chunk(self, chunk, effects=None, gene=None, snp=None, top=False):
        print("Loading chunk...")
        ffea_df = self.load_file(self.feature_metadata_file.replace("<CHUNK>", chunk))
        ffea_df = ffea_df.rename(index=str,
                                 columns={"chromosome": "feature_chromosome",
                                          "start": "feature_start",
                                          "end": "feature_end"})

        feature_id_to_gene_dict = self.create_dict_from_df(df=ffea_df, key=[("feature_id", None, None)], value=self.get_column(gene))
        if effects is not None:
            if len(set(effects.keys()).intersection(set(feature_id_to_gene_dict.values()))) == 0:
                print("\tWarning, skipping: '{}' no overlapping features.".format(chunk))
                return None

        fsnp_df = self.load_file(self.snp_metadata_file.replace("<CHUNK>", chunk))
        fsnp_df = fsnp_df.rename(index=str,
                                 columns={"chromosome": "snp_chromosome",
                                          "position": "snp_position"})

        snp_to_snp_id_dict = self.create_dict_from_df(df=fsnp_df, key=self.get_column(snp), value=[("snp_id", None, None)])
        if effects is not None:
            if len(set(self.get_effect_snps(effects=effects)).intersection(set(snp_to_snp_id_dict.keys()))) == 0:
                print("\tWarning, skipping: '{}' no overlapping SNPs.".format(chunk))
                return None

        frez = self.load_h5(inpath=self.qtl_results_file.replace("<CHUNK>", chunk))
        frez_keys = [k.replace('_i_', '') for k in list(frez.keys())]

        df_list = []
        for frez_key in frez_keys:
            if effects is not None and feature_id_to_gene_dict[frez_key] not in effects.keys():
                continue

            frez_df = pd.DataFrame(np.array(frez[frez_key]))
            frez_df['snp_id'] = frez_df['snp_id'].str.decode("utf-8")
            frez_df['feature_id'] = frez_key
            frez_df = frez_df.astype({"beta": float,
                                      "beta_se": float,
                                      "empirical_feature_p_value": float,
                                      "p_value": float,
                                      "snp_id": str,
                                      "feature_id": str})

            if effects is not None:
                snp_ids = []

                # Translate the feature_id (should be HGNC symbol) to --gene parameter.
                # Note that it might just translate hgnc to hgnc if --gene = hgnc.
                # Then, look up that gene in the effects dict and see which --snp's
                # we are looking for.
                for snp in effects[feature_id_to_gene_dict[frez_key]]:
                    # For each of the snp's we query, translate the --snp value
                    # to snp_id (should be rsid) and check if that is present in
                    # the frex_df.
                    if (snp not in snp_to_snp_id_dict) or (snp_to_snp_id_dict[snp] not in set(frez_df["snp_id"])):
                        continue

                    # If it is present, we will keep the snp_id value.
                    snp_ids.append(snp_to_snp_id_dict[snp])

                # Subset on the snps from the effects dict.
                frez_df = frez_df.loc[frez_df["snp_id"].isin(snp_ids), :]

            if frez_df.shape[0] == 0:
                continue

            df_list.append(frez_df)
            del frez_df

        if len(df_list) == 0:
            print("\tWarning, skipping: '{}' no overlapping feature-SNP combinations.".format(chunk))
            return None

        df = pd.concat(df_list, axis=0)
        del df_list

        df = pd.merge(df, ffea_df, on='feature_id', how='left')

        if (self.snp_qc_metrics is not None and len(glob.glob(self.snp_qc_metrics.replace("<FEATURE>", "*"))) > 0):
            print("\tError, code not implemented yet")
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

        # If we are keeping all effects we can calculate how many
        # variants we had per feature.
        if effects is None:
            df['nTotalTestsPerFeat'] = np.nan
            for index, value in df['feature_id'].value_counts().iteritems():
                df.loc[df['feature_id'] == index, 'nTotalTestsPerFeat'] = value

        # Keep the variant with the lowest permuted p-value per feature.
        if top:
            df = df.sort_values(by=['empirical_feature_p_value', 'p_value'], ascending=[True, True])
            df = df.groupby(df['feature_id']).first()
            df.reset_index(drop=False, inplace=True)

        del ffea_df, fsnp_df

        print("\tLoaded chunk: {} "
              "with shape: {}".format(chunk,
                                      df.shape))
        return df

    def __str__(self):
        class_str = super(LIMIX, self).__str__()
        return (class_str.rstrip("\n") +
                "\n  feature_metadata_file = {}\n  snp_metadata_file = {}\n  " \
                "qtl_results_file = {}\n  snp_qc_metrics = {}\n  chunk_pattern = {}\n".format(
                    self.feature_metadata_file,
                    self.snp_metadata_file,
                    self.qtl_results_file,
                    self.snp_qc_metrics,
                    self.chunk_pattern
                ))

##############################################################################################################


class mbQTL(Dataset):
    def __init__(self, *args, **kwargs):
        super(mbQTL, self).__init__(*args, **kwargs)
        self.class_name = "mbQTL"

        # Update filenames.
        self.update_all_filename(filename="<CT>-AllEffects.txt")
        self.update_top_filename(filename="<CT>-TopEffects.txt")

        # Set file paths.
        self.set_all_effects_path()
        self.set_top_effects_path()

        # Set effects info.
        self.set_effects_info()

        # Columns that are in the original file.
        self.columns.update({
            "gene_hgnc": [("GeneSymbol", None, None)],
            "gene_ensembl": self.get_gene_ensembl_column(),
            "SNP_rsid": self.get_snp_rsid_column(),
            "SNP_chr:pos": [("SNPChr", "([0-9]{1,2}|X|Y|MT)", ":"), ("SNPPos", "([0-9]+)", None)],
            "alleles": [("SNPAlleles", None, None)],
            "EA": [("SNPEffectAllele", None, None)],
            # "OA": [(None, None, None)],
            "beta": self.get_beta_column(),
            "beta_se": [("MetaSE", None, None)],
            "n_tests": self.get_n_tests_column(),
            "nominal_pvalue": [("MetaP", None, None)],
            "permuted_pvalue": self.get_permuted_pvalue_column(),
            # "bonferroni_pvalue": [(None, None, None)],
            "zscore": [("MetaPZ", None, None)],
            "FDR": self.get_fdr_column(),
            "N": [("MetaPN", None, None)],
            "AF": [("SNPEffectAlleleFreq", None, None)],
            # "MAN": [(None", None, None)],
            # "MAF": [(None, None, None)]
        })

    def get_beta_column(self):
        if self.effects_header is None:
            return [(None, None, None)]
        
        if "MetaBeta" in self.effects_header:
            return [("MetaBeta", None, None)]
        elif "MetaR" in self.effects_header:
            return [("MetaR", None, None)]
        else:
            return [(None, None, None)]

    def get_n_tests_column(self):
        if self.effects_header is None:
            return [(None, None, None)]

        if "NrTestedSNPs" in self.effects_header:
            return [("NrTestedSNPs", None, None)]
        else:
            return [(None, None, None)]

    def get_permuted_pvalue_column(self):
        if self.effects_header is None:
            return [(None, None, None)]

        if "BetaAdjustedMetaP" in self.effects_header:
            return [("BetaAdjustedMetaP", None, None)]
        else:
            return [(None, None, None)]

    def get_gene_ensembl_column(self):
        if self.effects_line is None:
            return [(None, None, None)]
        
        if re.match("ENSG[0-9]+$", self.effects_line["Gene"]):
            return [("Gene", "(ENSG[0-9]+)", None)]
        elif re.match("(ENSG[0-9]+).[0-9]+", self.effects_line["Gene"]):
            return [("Gene", "(ENSG[0-9]+).[0-9]+", None)]
        else:
            return [(None, None, None)]

    def get_snp_rsid_column(self):
        if self.effects_line is None:
            return [(None, None, None)]
        
        if re.match("rs[0-9]+$", self.effects_line["SNP"]):
            return [("SNP", "(rs[0-9]+)", None)]
        elif re.match("(?:[0-9]{1,2}|X|Y|MT):[0-9]+:(rs[0-9]+):[A-Z]+_[A-Z]+", self.effects_line["SNP"]):
            return [("SNP", "(?:[0-9]{1,2}|X|Y|MT):[0-9]+:(rs[0-9]+):[A-Z]+_[A-Z]+", None)]
        else:
            return [(None, None, None)]

    def get_fdr_column(self):
        if self.effects_header is None:
            return [(None, None, None)]
        
        if "qval" in self.effects_header:
            return [("qval", None, None)]
        else:
            return [(None, None, None)]


##############################################################################################################


class eQTLMappingPipeline(Dataset):
    def __init__(self, *args, **kwargs):
        super(eQTLMappingPipeline, self).__init__(*args, **kwargs)
        self.class_name = "eQTLMappingPipeline"

        # Update filenames.
        self.update_all_filename(filename="cis-eQTLs_full_20180905.txt")
        self.update_top_filename(filename="2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt") # TODO: for some reason this still has multiple snps per gene in there.

        # Set file paths.
        self.set_all_effects_path()
        self.set_top_effects_path()

        # Columns that are in the original file.
        self.columns.update({
            "gene_hgnc": [("GeneSymbol", None, None)],
            "gene_ensembl": [("Gene", "(ENSG[0-9]+)", None)],
            "SNP_rsid": [("SNP", "(rs[0-9]+)", None)],
            "SNP_chr:pos": [("SNPChr", "([0-9]{1,2}|X|Y|MT)", ":"), ("SNPPos", "([0-9]+)", None)],
            "alleles": [("AssessedAllele", None, "/"), ("OtherAllele", None, None)],
            "EA": [("AssessedAllele", None, None)],
            "OA": [("OtherAllele", None, None)],
            # "beta": [(None, None, None)],
            # "beta_se": [(None, None, None)],
            "n_tests": [("NEntries", None, None)], # TODO: only because the top file has multiple snps per gene in there
            "nominal_pvalue": [("Pvalue", None, None)],
            # "permuted_pvalue": [(None, None, None)],
            "bonferroni_pvalue": [("BonferroniP", None, None)],
            "zscore": [("Zscore", None, None)],
            # "FDR": [("FDR", None, None)], # TODO: not sure what kind of FDR this is
            # "N": [("NrSamples", None, None)],
            # "AF": [(None, None, None)],
            # "MAN": [(None", None, None)],
            # "MAF": [(None, None, None)],
        })

    def get_top_effects(self, gene=None, snp=None):
        # Making sure the required columns are available.
        for label in [gene, "nominal_pvalue"]:
            if not self.contains_data(label=label):
                print("Error, get_top_effects() from the all effects file is unavailable for {} since"
                      " there is no '{}' data available.".format(self.class_name, label))
                exit()

        try:
            return self.get_effects_wrapper(
                inpath=self.get_top_effects_path(),
                label_dict={"key": self.get_column(gene), "value": self.get_column("nominal_pvalue")},
                func=self.extract_top_effects
            )
        except FileNotFoundError as e:
            print(e)
            print("Failed to load 'top_effects_path'. Selecting top effects from 'all_effects_path' instead.")

        try:
            return self.get_effects_wrapper(
                inpath=self.get_all_effects_path(),
                label_dict={"key": self.get_column(gene), "value": self.get_column("nominal_pvalue")},
                func=self.extract_top_effects
            )
        except FileNotFoundError as e:
            print(e)
            print("Failed to load 'all_effects_path'.")

        return None

##############################################################################################################


class eQTLgenPhase2(Dataset):
    def __init__(self, *args, **kwargs):
        super(eQTLgenPhase2, self).__init__(*args, **kwargs)
        self.class_name = "eQTLgenPhase2"

        # Update filenames.
        self.update_all_filename(filename="output_empirical_4GenPC20ExpPC_2023-05-27_interestingGeneSnpCombi_fixedConcatenation.csv")  # Not really all effects, just the sc-eQTLgen ones
        # self.update_top_filename(filename=None)

        # Set file paths.
        self.set_all_effects_path()
        self.set_top_effects_path()

        # Columns that are in the original file.
        self.columns.update({
            # "gene_hgnc": [(None, None, None)],
            "gene_ensembl": [("phenotype", "(ENSG[0-9]+)", None)],
            "SNP_rsid": [("variant", "(rs[0-9]+)", None)],
            "SNP_chr:pos": [("chromosome_variant", "([0-9]{1,2}|X|Y|MT)", ":"), ("bp_variant", "([0-9]+)", None)],
            "alleles": [("allele_ref", None, "/"), ("allele_eff", None, None)],
            "EA": [("allele_eff", None, None)],
            # "OA": [(None, None, None)],
            "beta": [("beta", None, None)],
            "beta_se": [("standard_error", None, None)],
            # "n_tests": [(None, None, None)],
            "nominal_pvalue": [("p_value", None, None)],
            # "permuted_pvalue": [(None, None, None)],
            # "bonferroni_pvalue": [(None, None, None)],
            "zscore": [("z_score", None, None)],
            # "FDR": [(None, None, None)],
            "N": [("sample_size", None, None)],
            "AF": [("allele_eff_freq", None, None)],
            # "MAN": [(None", None, None)],
            # "MAF": [(None, None, None)]
        })


##############################################################################################################

class Bryois(Dataset):
    def __init__(self, *args, **kwargs):
        super(Bryois, self).__init__(*args, **kwargs)
        self.class_name = "Bryois"

        # IMPORTANT: Bryois et al. published their significant results in supplementary table
        # '41593_2022_1128_MOESM3_ESM.xlsx'. However, this table does not contain the nominal
        # p-values. If --pvalue is set to nominal or bonferroni I require these values even
        # if Bryois is usd as discovery dataset. Since finding the nominal p-value requires
        # me to parse the excel file as well as all <CT>.<CHR>.gz files I figured it would
        # be better to add this column once using 'bryois2022/get_bryois_top_effects.py'
        # and use that file as top_effects_path instead.

        # Update filenames.
        self.update_all_filename(filename="<CT>.<CHR>.gz")
        self.update_top_filename(filename="<CT>.top.gz")

        # Set file paths.
        self.set_all_effects_path()
        self.set_top_effects_path()
        self.snp_pos_path = self.get_fpath(fpath=os.path.join(self.path, "snp_pos.txt.gz"))

        # Other variables.
        self.all_effects_columns = ["symbol_ensembl", "SNP", "dist_TSS", "pval", "beta"]

        # Number of samples depends on the cell type.
        self.n = {
            "Astrocytes": 189,
            "Endothelial.cells": 144,
            "Excitatory.neurons": 187,
            "Inhibitory.neurons": 178,
            "Microglia": 187,
            "Oligodendrocytes": 192,
            "OPCs...COPs": 187,
            "Pericytes": 150
        }[self.cell_type]

        # Columns that are in the original file.
        self.columns.update({
            "gene_hgnc": [("symbol_ensembl", "([a-zA-Z0-9-.]+)_ENSG[0-9]+", None)],
            "gene_ensembl": [("symbol_ensembl", "[a-zA-Z0-9-.]+_(ENSG[0-9]+)", None)],
            "SNP_rsid": [("SNP", "(rs[0-9]+)", None)],
            "SNP_chr:pos": [("SNP_id_hg38", "chr((([0-9]{1,2}|X|Y|MT):[0-9]+))", None)],
            "alleles": [("effect_allele", None, "/"), ("other_allele", None, None)],
            "EA": [("effect_allele", None, None)],
            "OA": [("other_allele", None, None)],
            "beta": [("beta", None, None)],
            # "beta_se": [(None, None, None)],
            # "n_tests": [(None, None, None)],
            "nominal_pvalue": [("pval", None, None)],
            "permuted_pvalue": [("bpval", None, None)],
            # "bonferroni_pvalue": [(None, None, None)],
            # "zscore": [(None, None, None)],
            "FDR": [("adj_p", None, None)],
            # "N": [(None, None, None)],
            # "AF": [(None, None, None)],
            # "MAN": [(None", None, None)],
            "MAF": [("MAF", None, None)]
        })

    def load_snp_pos_data(self, effects, specific_entries_cols):
        specific_entries_cols = self.trans_label_to_index(
            inpath=self.snp_pos_path,
            label_dict=specific_entries_cols
        )
        snp_pos_df = self.load_partial_file(
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

    def get_all_effects(self):
        df = super(Bryois, self).get_all_effects()
        if df is None:
            return None

        df.columns = self.all_effects_columns

        snp_pos_df = self.load_snp_pos_data(effects=set(df["SNP"]), specific_entries_cols={"key": [("SNP", None, None)]})
        df = df.merge(snp_pos_df, how="left")
        return df

    # def extract_all_effects(self, inpath, **kwargs):
    #     # TODO: should this actually be overwritten? Top effects is tab seperated with header I thought.
    #     # Need to overwrite the default here since Bryois files do not have a header and are space seperated.
    #     return super(Bryois, self).extract_all_effects(inpath=inpath, header=None, sep=" ")

    def get_top_effects(self, gene=None, snp=None):
        df = super(Bryois, self).get_top_effects(gene=gene, snp=snp)
        if df is None:
            return None

        snp_pos_df = self.load_snp_pos_data(effects=set(df["SNP"]), specific_entries_cols={"key": [("SNP", None, None)]})
        df = df.merge(snp_pos_df, how="left")
        return df

    def extract_top_effects(self, inpath, label_dict=None, top_entries_cols=None, **kwargs):
        return super(Bryois, self).extract_top_effects(inpath=inpath, label_dict=label_dict, top_entries_cols=top_entries_cols, header=None, sep=" ")

    def get_specific_effects(self, effects=None, gene=None, snp=None):
        # TODO: implement alternative option to load from
        #  'top_effects_path' if 'all_effects_path' is missing.

        snp_pos_df = self.load_snp_pos_data(effects=self.get_effect_snps(effects=effects), specific_entries_cols={"key": self.get_column(snp)})
        label_dict = {"key": self.get_column(gene), "value": self.get_column(snp)}

        if self.get_column(snp)[0][0] not in self.all_effects_columns:
            extra_info = self.create_dict_from_df(
                df=snp_pos_df,
                key=[("SNP", None, None)],
                value=label_dict["value"]
            )
            label_dict["value"] = [("SNP", extra_info, None)]

        specific_entries_cols = self.update_label_dict(
            label_dict=label_dict,
            colname_to_index_dict=self.get_partial_file_column_index()
        )

        all_effects_path = self.get_all_effects_path()
        df = self.get_effects_wrapper(
            inpath=all_effects_path,
            effects=effects,
            func=self.extract_specific_effects,
            specific_entries_cols=specific_entries_cols
        )
        df.columns = self.all_effects_columns
        df = df.merge(snp_pos_df, how="left")
        return df

    def extract_specific_effects(self, inpath, label_dict=None, specific_entries_cols=None, effects=None, **kwargs):
        # Need to overwrite the default here since Bryois files do not have a header and are space seperated.
        return super(Bryois, self).extract_specific_effects(inpath=inpath, label_dict=label_dict, specific_entries_cols=specific_entries_cols, effects=effects, header=None, sep=" ")

    def __str__(self):
        class_str = super(Bryois, self).__str__()
        return (class_str.rstrip("\n") +
                "\n  snp_pos_path = {}\n  all_effects_columns = {}\n".format(
                    self.snp_pos_path,
                    self.all_effects_columns
                ))



##############################################################################################################


class Bryois_REDUCED(Dataset):
    def __init__(self, *args, **kwargs):
        super(Bryois_REDUCED, self).__init__(*args, **kwargs)
        self.class_name = "Bryois_REDUCED"

        # Update filenames.
        # self.update_all_filename(filename=None)
        self.update_top_filename(filename="JulienBryois2021SummaryStats.txt.gz")

        # Set file paths.
        self.set_all_effects_path()
        self.set_top_effects_path()

        # Columns that are in the original file.
        self.columns.update({
            # "gene_hgnc": [(None, None, None)],
            "gene_ensembl": [("", "(ENSG[0-9]+)_rs[0-9]+", None)],
            "SNP_rsid": [("", "ENSG[0-9]+_(rs[0-9]+)", None)],
            # "SNP_chr:pos": [(None, None, None)],
            # "alleles": [(None, None, None)],
            "EA": [("effect_allele", None, None)],
            # "OA": [(None, None, None)],
            "beta": [(self.cell_type + " beta", None, None)],
            # "beta_se": [(None, None, None)],
            # "n_tests": [(None, None, None)],
            "nominal_pvalue": [(self.cell_type + " p-value", None, None)],
            # "permuted_pvalue": [(None", None, None)],
            # "bonferroni_pvalue": [(None, None, None)],
            # "zscore": [(None, None, None)],
            # "FDR": [(None, None, None)],
            # "N": [(None, None, None)],
            # "AF": [(None, None, None)],
            # "MAN": [(None", None, None)],
            # "MAF": [(None, None, None)],
        })

    def update_df(self, df):
        # Remove empty rows and unwanted columns.
        df = df.loc[df[self.cell_type + " beta"] != "", [col for col in df.columns if (not col.endswith("p-value") and not col.endswith("beta")) or col.startswith(self.cell_type)]]

        # TODO: Could still merge snp_pos file into this.
        return df

##############################################################################################################

class Fujita(Dataset):
    def __init__(self, *args, **kwargs):
        super(Fujita, self).__init__(*args, **kwargs)
        self.class_name = "Fujita"

        # Update filenames.
        self.update_all_filename(filename="celltype-eqtl-sumstats.<CT>.tsv")
        # self.update_top_filename(filename=None)

        # Set file paths.
        self.set_all_effects_path()
        self.set_top_effects_path()

        # Default sample size.
        self.n = 424 # TODO: might not be the same for all cell types?

        # Columns that are in the original file.
        self.columns.update({
            "gene_hgnc": [("gene_symbol", None, None)],
            "gene_ensembl": [("gene_id", "(ENSG[0-9]+)", None)],
            "SNP_rsid": [("snps", "(rs[0-9]+)", None)],
            "SNP_chr:pos": [("chr38", "chr([0-9]{1,2}|X|Y|MT)", ":"), ("pos38", "([0-9]+)", None)],
            "alleles": [("ALT", None, "/"), ("REF", None, None)],
            "EA": [("ALT", None, None)],  # Don't ask me why but for some reason the alternative allele is the effect allele
            "OA": [("REF", None, None)],
            "beta": [("beta", None, None)],
            "beta_se": [("se", None, None)],
            "n_tests": [("NEntries", None, None)],
            "nominal_pvalue": [("pvalue", None, None)],
            # "permuted_pvalue": [(None, None, None)],
            # "bonferroni_pvalue": [(None, None, None)],
            # "zscore": [(None, None, None)],
            # "FDR": [(None, None, None)],
            # "N": [(None, None, None)],
            "AF": [("ALT_AF", None, None)],
            # "MAN": [(None", None, None)],
            # "MAF": [(None, None, None)]
        })

##############################################################################################################


class DeconQTL(Dataset):
    def __init__(self, *args, **kwargs):
        super(DeconQTL, self).__init__(*args, **kwargs)
        self.class_name = "DeconQTL"

        # TODO: work in progress, untested code

        # Update filenames.
        # self.update_all_filename(filename=None)
        self.update_top_filename(filename="deconvolutionResults.txt.gz")

        # Set file paths.
        self.set_all_effects_path()
        self.set_top_effects_path()
        self.alleles_path = self.get_fpath(fpath="/groups/umcg-biogen/prm03/projects/2022-DeKleinEtAl/output/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-10-18-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_alleles.txt.gz") # TODO: this is hard-coded
        self.genotypes_stats_path = self.get_fpath(fpath=os.path.join(self.path, "geno_stats.txt.gz"))

        # Other variables.
        self.beta_pattern = "Beta[0-9]+_<CT>:GT"
        self.gene_pattern = [("Unnamed: 0", "(ENSG[0-9]+.[0-9]+)_(?:[0-9]{1,2}|X|Y|MT):[0-9]+:rs[0-9]+:[A-Z]+_[A-Z]+", None)]
        self.snp_pattern = [("Unnamed: 0", "ENSG[0-9]+.[0-9]+_((?:[0-9]{1,2}|X|Y|MT):[0-9]+:rs[0-9]+:[A-Z]+_[A-Z]+)", None)]
        self.columns["beta"] = self.get_beta_column()

        # Columns that are in the original file.
        self.columns.update({
            # "gene_hgnc": [(None, None, None)],
            "gene_ensembl": [("", "(ENSG[0-9]+).[0-9]+_(?:[0-9]{1,2}|X|Y|MT):[0-9]+:rs[0-9]+:[A-Z]+_[A-Z]+", None)],
            "SNP_rsid": [("", "ENSG[0-9]+.[0-9]+_(?:[0-9]{1,2}|X|Y|MT):[0-9]+:(rs[0-9]+):[A-Z]+_[A-Z]+", None)],
            "SNP_chr:pos": [("", "ENSG[0-9]+.[0-9]+_(([0-9]{1,2}|X|Y|MT):[0-9]+):rs[0-9]+:[A-Z]+_[A-Z]+", None)],
            "alleles": [("", "ENSG[0-9]+.[0-9]+_(?:[0-9]{1,2}|X|Y|MT):[0-9]+:rs[0-9]+:([A-Z]+)_[A-Z]+", "/"), ("", "ENSG[0-9]+.[0-9]+_(?:[0-9]{1,2}|X|Y|MT):[0-9]+:rs[0-9]+:[A-Z]+_([A-Z]+)", None)],
            "EA": [("EA", None, None)],
            "OA": [("OA", None, None)],
            # "beta": [(None, None, None)],
            # "beta_se": [(None, None, None)],
            # "n_tests": [(None, None, None)],
            "nominal_pvalue": [(self.cell_type + "_pvalue", None, None)],
            # "permuted_pvalue": [(None", None, None)],
            # "bonferroni_pvalue": [(None, None, None)],
            # "zscore": [(None, None, None)],
            # "FDR": [(None, None, None)],
            # "N": [("N", None, None)],
            # "AF": [(None, None, None)],
            # "MAN": [(None", None, None)],
            # "MAF": [(None, None, None)],
        })

    def get_beta_column(self):
        header = self.get_file_header(inpath=self.top_effects_path)

        pattern = self.beta_pattern.replace("<CT>", self.cell_type)
        column_of_interest = []
        for column in header:
            if re.match(pattern, column):
                column_of_interest.append(column)

        if len(column_of_interest) == 1:
            return [(column_of_interest[0], None, None)]
        elif len(column_of_interest) == 0:
            print("Error, could not find any columns matching '{}'.".format(pattern))
            exit()
        else:
            print("Error, found multiple columns matching '{}': {}.".format(pattern, ", ".join(column_of_interest)))
            exit()

    def update_df(self, df):
        df["SNP"] = df.apply(lambda row: self.extract_info(data=row, query=self.snp_pattern), axis=1)

        # Load the alleles path.
        # We assume here that the alleles are encoded as 0/2 thus 2 is the effect allele (should be the case).
        alleles_df = self.load_file(inpath=self.alleles_path)
        alleles_df.rename(columns={"Unnamed: 0": "SNP"}, inplace=True)
        alleles_df[["OA", "EA"]] = alleles_df["Alleles"].str.split("/", n=1, expand=True)
        alleles_df.drop(["Alleles"], axis=1, inplace=True)

        # Load the genotype stats.
        genotype_stats_df = self.load_file(inpath=self.genotypes_stats_path)
        genotype_stats_df.rename(columns={"Unnamed: 0": "SNP"}, inplace=True)

        if not alleles_df["SNP"].equals(genotype_stats_df["SNP"]):
            print("Error, the SNP column in the alleles dataframe and genotype stats dataframe do not match.")
            exit()

        # Merge the info stats.
        genotype_stats_df.drop(["SNP"], axis=1, inplace=True)
        info_df = alleles_df.merge(genotype_stats_df, left_index=True, right_index=True, how="left")
        info_df = info_df.loc[info_df["mask"] == 1, [col for col in info_df.columns if col != "mask"]]
        info_df.drop_duplicates(inplace=True)

        if len(set(info_df["SNP"])) != info_df.shape[0]:
            print("Error, the info data frame has duplicate SNP values. Unable to merge into the summary statistics dataframe.")
            exit()

        return info_df.merge(df, on="SNP", how="right")



##############################################################################################################


class PICALO(Dataset):
    def __init__(self, *args, **kwargs):
        super(PICALO, self).__init__(*args, **kwargs)
        self.class_name = "PICALO"

        # Update filenames.
        # self.update_all_filename(filename=None)
        self.update_top_filename(filename="<CT>.txt.gz")

        # Set file paths.
        self.set_all_effects_path()
        self.set_top_effects_path()
        self.alleles_path = self.get_fpath(fpath="/groups/umcg-biogen/prm03/projects/2022-VochtelooEtAl-PICALO/preprocess_scripts/prepare_picalo_files/2022-03-24-MetaBrain_CortexEUR_NoENA_NoRNAseqAlignmentMetrics_GT1AvgExprFilter_PrimaryeQTLs_UncenteredPCA/genotype_alleles_table.txt.gz") # TODO: this is hard-coded
        self.genotypes_stats_path = self.get_fpath(fpath=os.path.join(self.path, "genotype_stats.txt.gz"))

        # Columns that are in the original file.
        self.columns.update({
            # "gene_hgnc": [(None, None, None)],
            "gene_ensembl": [("gene", "(ENSG[0-9]+).[0-9]+", None)],
            "SNP_rsid": [("SNP", "(?:[0-9]{1,2}|X|Y|MT):[0-9]+:(rs[0-9]+):[A-Z]+_[A-Z]+", None)],
            "SNP_chr:pos": [("SNP", "(([0-9]{1,2}|X|Y|MT):([0-9]+)):rs[0-9]+:[A-Z]+_[A-Z]+", None)],
            "alleles": [("SNP", "(?:[0-9]{1,2}|X|Y|MT):[0-9]+:rs[0-9]+:([A-Z]+)_[A-Z]+", "/"), ("SNP", "(?:[0-9]{1,2}|X|Y|MT):[0-9]+:rs[0-9]+:[A-Z]+_([A-Z]+)", "")],
            "EA": [("EA", None, None)],
            "OA": [("OA", None, None)],
            "beta": [("beta-interaction", None, None)],
            "beta_se": [("std-interaction", None, None)],
            # "n_tests": [(None, None, None)],
            "nominal_pvalue": [("p-value", None, None)],
            # "permuted_pvalue": [(None, None, None)],
            # "bonferroni_pvalue": [(None, None, None)],
            # "zscore": [(None, None, None)],
            "FDR": [("FDR", None, None)],
            "N": [("N", None, None)],
            # "AF": [(None, None, None)],
            # "MAN": [(None", None, None)],
            "MAF": [("MAF", None, None)]
        })

    def update_df(self, df):
        # Load the alleles path.
        # We assume here that the alleles are encoded as 0/2 thus 2 is the effect allele (should be the case).
        alleles_df = self.load_file(inpath=self.alleles_path)
        alleles_df[["OA", "EA"]] = alleles_df["Alleles"].str.split("/", n=1, expand=True)
        alleles_df.drop(["Alleles"], axis=1, inplace=True)

        # Load the genotype stats.
        genotype_stats_df = self.load_file(inpath=self.genotypes_stats_path)
        genotype_stats_df.rename(columns={"N": "N_genotype"}, inplace=True)

        if not alleles_df["SNP"].equals(genotype_stats_df["SNP"]):
            print("Error, the SNP column in the alleles dataframe and genotype stats dataframe do not match.")
            exit()

        # Merge the info stats.
        genotype_stats_df.drop(["SNP"], axis=1, inplace=True)
        info_df = alleles_df.merge(genotype_stats_df, left_index=True, right_index=True, how="left")
        info_df = info_df.loc[info_df["mask"] == 1, [col for col in info_df.columns if col != "mask"]]
        info_df.reset_index(inplace=True, drop=True)

        if not info_df["SNP"].equals(df["SNP"]):
            print("Error, the SNP column in the info dataframe and summary statistics dataframe do not match.")
            exit()

        df.drop(["SNP"], axis=1, inplace=True)
        return info_df.merge(df, left_index=True, right_index=True, how="right")


##############################################################################################################


class GTEx(Dataset):
    def __init__(self, *args, **kwargs):
        super(GTEx, self).__init__(*args, **kwargs)
        self.class_name = "GTEx"

        # Update filenames.
        self.update_all_filename(filename="<CT>.v8.signif_variant_gene_pairs.txt.gz")
        self.update_top_filename(filename="<CT>.v8.egenes.txt.gz")

        # Set file paths.
        self.set_all_effects_path()
        self.set_top_effects_path()

        # Columns that are in the original file.
        self.columns.update(self.select_columns())

    def get_discovery_columns(self):
        return {
            "gene_hgnc": [("gene_name", None, None)],
            "gene_ensembl": [("gene_id", "(ENSG[0-9]+).[0-9]+", None)],
            "SNP_rsid": [("rs_id_dbSNP151_GRCh38p7", None, None)],
            "SNP_chr:pos": [("chr", "chr([0-9]{1,2}|X|Y|MT)", ":"), ("variant_pos", None, None)],
            "alleles": [("ref", None, "/"), ("alt", None, "")],
            "EA": [("variant_id", "chr(?:[0-9]{1,2}|X|Y|MT)_[0-9]+_[A-Z]_([A-Z]+)_b38", None)],
            # "OA": [(None, None, None)],
            "beta": [("slope", None, None)],
            "beta_se": [("slope_se", None, None)],
            # "n_tests": [(None, None, None)],
            "nominal_pvalue": [("pval_nominal", None, None)],
            "permuted_pvalue": [("pval_perm", None, None)],
            # "bonferroni_pvalue": [(None, None, None)],
            # "zscore": [(None, None, None)],
            "FDR": [("qval", None, None)],
            # "N": [(None, None, None)],
            # "AF": [(None, None, None)],
            "MAN": [("minor_allele_samples", None, None)],
            "MAF": [("maf", None, None)]
        }

    def get_replication_columns(self):
        return {
            # "gene_hgnc": [(None, None, None)],
            "gene_ensembl": [("gene_id", "(ENSG[0-9]+).[0-9]+", None)],
            # "SNP_rsid": [(None, None, None)],
            "SNP_chr:pos": [("variant_id", "chr([0-9]{1,2}|X|Y|MT)_[0-9]+_[A-Z]_[A-Z]+_b38", ":"), ("variant_id", "chr(?:[0-9]{1,2}|X|Y|MT)_([0-9]+)_[A-Z]_[A-Z]+_b38", None)],
            "alleles": [("variant_id", "chr(?:[0-9]{1,2}|X|Y|MT)_[0-9]+_([A-Z])_[A-Z]+_b38", "/"), ("variant_id", "chr(?:[0-9]{1,2}|X|Y|MT)_[0-9]+_[A-Z]_([A-Z]+)_b38", "")],
            "EA": [("variant_id", "chr(?:[0-9]{1,2}|X|Y|MT)_[0-9]+_[A-Z]_([A-Z]+)_b38", None)],
            # "OA": [(None, None, None)],
            "beta": [("slope", None, None)],
            "beta_se": [("slope_se", None, None)],
            # "n_tests": [(None, None, None)],
            "nominal_pvalue": [("pval_nominal", None, None)],
            # "permuted_pvalue": [(None, None, None)],
            # "bonferroni_pvalue": [(None, None, None)],
            # "zscore": [(None, None, None)],
            # "FDR": [(None, None, None)],
            # "N": [(None, None, None)],
            # "AF": [(None, None, None)],
            "MAN": [("ma_samples", None, None)],
            "MAF": [("maf", None, None)]
        }

##############################################################################################################



if __name__ == '__main__':
    m = main()
    m.start()
