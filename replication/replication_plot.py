#!/usr/bin/env python3

"""
File:         replication_plot.py
Created:      2024/02/28
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
import argparse
import gzip
import glob
import os
import re
import pickle
import json
from abc import abstractmethod

# Third party imports.
import h5py
import numpy as np
import pandas as pd
from natsort import natsort_keygen
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
./replication_plot.py -h
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

class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.discovery_path = getattr(arguments, 'discovery_path')
        self.discovery_name = getattr(arguments, 'discovery_name')
        self.discovery_cell_type = getattr(arguments, 'discovery_cell_type')
        self.replication_path = getattr(arguments, 'replication_path')
        self.replication_name = getattr(arguments, 'replication_name')
        self.replication_cell_type = getattr(arguments, 'replication_cell_type')
        self.ancestry = getattr(arguments, 'ancestry')
        self.cell_level = getattr(arguments, 'cell_level')
        self.alpha = getattr(arguments, 'alpha')
        self.signif_column = getattr(arguments, 'signif_column')
        self.palette_path = getattr(arguments, 'palette')
        self.extensions = getattr(arguments, 'extension')
        self.force = getattr(arguments, 'force')
        self.qvalues_script = getattr(arguments, 'qvalues')
        self.rb_script = getattr(arguments, 'rb')

        self.disc = self.get_class(name=self.discovery_name, path=self.discovery_path)
        self.repl = self.get_class(name=self.replication_name, path=self.replication_path)

        # Set variables.
        self.outdir = os.path.join(str(os.path.dirname(os.path.abspath(__file__))), "replication_plot", "{}discovery_{}replication".format(self.discovery_name, self.replication_name))
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Loading palette.
        self.palette = None
        if self.palette_path is not None:
            with open(self.palette_path) as f:
                self.palette = json.load(f)
            f.close()

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

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
                            required=True,
                            default=None,
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
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--ancestry",
                            type=str,
                            required=False,
                            default="EUR",
                            help="")
        parser.add_argument("--cell_level",
                            type=str,
                            required=False,
                            default="L1",
                            help="")
        parser.add_argument("--alpha",
                            type=str,
                            required=False,
                            default=0.05,
                            help="")
        parser.add_argument("--signif_column",
                            type=str,
                            required=False,
                            default="pvalue",
                            help="")
        parser.add_argument("--palette",
                            type=str,
                            required=False,
                            default=None,
                            help="A color palette file.")
        parser.add_argument("-e",
                            "--extension",
                            nargs="+",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")
        parser.add_argument("--force",
                            action='store_true',
                            help="")

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

    def get_class(self, name, path):
        if "Bryois" in name:
            return Bryois(name=name, path=path)
        elif "scMetaBrain" in name:
            return scMetaBrain(name=name, path=path, ancestry=self.ancestry, cell_level=self.cell_level)
        elif "mbQTL" in name:
            return mbQTL(name=name, path=path)
        else:
            print("Error, unexpected name '{}'".format(name))
            exit()

    def start(self):
        self.print_arguments()

        # # Determine overlapping variants.
        # overlapping_variants = self.get_overlapping_variants()
        #
        # # Set the overlapping variants.
        # self.disc.set_overlapping_variants(overlapping_variants)
        #

        # First get the top effect per gene in the discovery dataset.
        discovery_top_filepath = os.path.join(self.outdir, self.discovery_name + self.discovery_cell_type + "TopEffects.txt.gz")
        if os.path.exists(discovery_top_filepath) and not self.force:
            discovery_top_df = self.load_file(discovery_top_filepath)
        else:
            discovery_top_df = self.disc.get_top_effects(cell_type=self.discovery_cell_type)
            self.save_file(df=discovery_top_df, outpath=discovery_top_filepath)
        print(discovery_top_df)

        # Create a dict of top effects with keys being HGNC names and values being RS IDs.
        discovery_eqtls = self.disc.get_eqtls_from_summary_stats(df=discovery_top_df)

        # Select the specific top effects from the discovery datasets in the replication dataset.
        replication_filepath = os.path.join(self.outdir, self.replication_name + self.replication_cell_type + "SpecificEffects.txt.gz")
        if os.path.exists(replication_filepath) and not self.force:
            replication_df = self.load_file(replication_filepath)
        else:
            replication_df = self.repl.get_specific_effects(cell_type=self.replication_cell_type, effects=discovery_eqtls)
            self.save_file(df=replication_df, outpath=replication_filepath)
        print(replication_df)

        # Overlap the data frames and harmonise direction of effect.
        overlapped_df = self.overlap_summary_stats(
            disc_df=discovery_top_df,
            repl_df=replication_df,
            disc_name=self.disc.get_name(),
            repl_name=self.repl.get_name()
        )
        print(overlapped_df)

        if overlapped_df is None or overlapped_df.shape[0] == 0:
            print("No rows to plot.")
            exit()

        # Compare.
        self.plot(
            df=overlapped_df,
            disc_name=self.disc.get_name(),
            repl_name=self.repl.get_name(),
            disc_cell_type=self.discovery_cell_type,
            repl_cell_type=self.replication_cell_type,
        )

    def print_arguments(self):
        print("Arguments:")
        print("  > Discovery path: {}".format(self.discovery_path))
        print("  > Discovery name: {}".format(self.discovery_name))
        print("  > Discovery cell type: {}".format(", ".join(self.discovery_cell_type)))
        print("  > Replication path: {}".format(self.replication_path))
        print("  > Replication name: {}".format(self.replication_name))
        print("  > Replication cell type: {}".format(", ".join(self.replication_cell_type)))
        print("  > Palette path: {}".format(self.palette_path))
        print("  > Ancestry: {}".format(self.ancestry))
        print("  > Cell level: {}".format(self.cell_level))
        print("  > Alpha: {}".format(self.alpha))
        print("  > Signif. column: {}".format(self.signif_column))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Extensions: {}".format(self.extensions))
        print("  > Force: {}".format(self.force))
        print("  > qvalues script: {}".format(self.qvalues_script))
        print("  > Rb script: {}".format(self.rb_script))
        print("")


    def get_overlapping_variants(self):
        print("Finding overlapping variants...")
        outpath = os.path.join(self.outdir, "overlapping_variants.pkl")

        overlapping_variants = None
        if os.path.exists(outpath) and not self.force:
            with gzip.open(outpath, 'rb') as fp:
                overlapping_variants = pickle.load(fp)
            fp.close()
        else:
            discovery_variants = self.disc.get_variants()
            replication_variants = self.repl.get_variants()
            overlapping_variants = discovery_variants.intersection(replication_variants)

            with gzip.open(outpath, 'wb') as fp:
                pickle.dump(overlapping_variants, fp)
            fp.close()

            del discovery_variants, replication_variants

        print("\tloaded {:,} variants".format(len(overlapping_variants)))

        return overlapping_variants

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
        print("Overlapping summary statistics...")
        print("\tDiscovery dataframe '{}' has shape: {}".format(disc_name, disc_df.shape))
        print("\tReplication dataframe '{}' has shape: {}".format(repl_name, repl_df.shape))

        df = disc_df.merge(repl_df, left_index=True, right_index=True, how='left')
        print(df)

        na_mask = (df[repl_name + " pvalue"].isna()).to_numpy()
        print("\tReplication N-effects overlapping: {:,}".format(np.sum(na_mask)))
        if np.sum(na_mask) == df.shape[0]:
            return None

        if repl_name + " OA" in df.columns:
            df = df.loc[~na_mask & ((df[disc_name + " EA"] == df[repl_name + " EA"]) | (df[disc_name + " EA"] == df[repl_name + " OA"])), :]
            print("\tMismatched variants removed dataframe has shape: {}".format(df.shape))
            if df.shape[0] == 0:
                return None
        else:
            print("\tWarning, could not verify that the alleles match. Assuming it is fine.")

        df[repl_name + " flip"] = df[repl_name + " EA"] != df[disc_name + " EA"]
        df[repl_name + " beta"] = df[repl_name + " beta"] * df[repl_name + " flip"].map({True: -1, False: 1})
        df[repl_name + " EA"] = df[disc_name + " EA"]
        if repl_name + " OA" in df.columns:
            df[repl_name + " OA"] = df[disc_name + " OA"]
        print("\t{:,} betas where flipped".format(df[repl_name + " flip"].sum()))

        df[disc_name + " FDR"] = multitest.multipletests(df[disc_name + " pvalue"], method='fdr_bh')[1]
        discovery_mask = (df[disc_name + " FDR"] <= self.alpha).to_numpy()
        print("\tDiscovery N-signif: {:,}".format(np.sum(discovery_mask)))

        na_mask = (df[repl_name + " pvalue"].isna()).to_numpy()
        mask = np.logical_and(discovery_mask, na_mask)
        df[repl_name + " FDR"] = np.nan
        if np.sum(mask) > 1:
            df.loc[mask, repl_name + " FDR"] = multitest.multipletests(df.loc[mask, repl_name + " pvalue"], method='fdr_bh')[1]

        replication_mask = (df[repl_name + " FDR"] <= self.alpha).to_numpy()
        print("\tReplication N-signif: {:,}".format(np.sum(replication_mask)))

        return df


    def plot(self, df, disc_name, repl_name, disc_cell_type, repl_cell_type, title=""):
        print("Plotting {}...".format(title))

        color = "#808080"
        if (disc_cell_type == repl_cell_type) and (disc_cell_type in self.palette.keys()):
            color = self.palette[disc_cell_type]

        plot_df = df.loc[~df[repl_name + " pvalue"].isna(), :].copy()
        plot_df.sort_values(by=disc_name + " pvalue", inplace=True)

        nrows = 1
        ncols = 3

        self.shared_ylim = {i: (0, 1) for i in range(nrows)}
        self.shared_xlim = {i: (0, 1) for i in range(ncols)}

        sns.set(rc={'figure.figsize': (ncols * 8, nrows * 6)})
        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharex='none',
                                 sharey='none')
        if nrows == 1:
            axes = axes[np.newaxis, ...]

        disc_rb_beta, disc_rb_se = self.prepare_rb_info(df=plot_df, name=disc_name)
        repl_rb_beta, repl_rb_se = self.prepare_rb_info(df=plot_df, name=repl_name)

        plot_df[disc_name + " log beta"] = self.log_modulus_beta(plot_df[disc_name + " beta"])
        plot_df[repl_name + " log beta"] = self.log_modulus_beta(plot_df[repl_name + " beta"])

        print("\tPlotting column 1.")
        xlim, ylim, stats1 = self.scatterplot(
            df=plot_df,
            fig=fig,
            ax=axes[0, 0],
            x=disc_name + " log beta",
            y=repl_name + " log beta",
            xlabel=disc_name + " log eQTL beta",
            ylabel=repl_name + " log eQTL beta",
            title="All",
            color=color,
            include_ylabel=True
        )
        self.update_limits(xlim, ylim, 0, 0)

        print("\tPlotting column 2.")
        xlim, ylim, stats2 = self.scatterplot(
            df=plot_df.loc[plot_df[disc_name + " " + self.signif_column] <= self.alpha, :],
            fig=fig,
            ax=axes[0, 1],
            x=disc_name + " log beta",
            y=repl_name + " log beta",
            xlabel=disc_name + " log eQTL beta",
            ylabel="",
            title=disc_name + " signif.",
            color=color,
            include_ylabel=False,
            pi1_column=repl_name + " pvalue",
            rb_columns=[(disc_rb_beta, disc_rb_se), (repl_rb_beta, repl_rb_se)]
        )
        self.update_limits(xlim, ylim, 0, 1)

        print("\tPlotting column 3.")
        xlim, ylim, stats3 = self.scatterplot(
            df=plot_df.loc[(plot_df[disc_name + " " + self.signif_column] <= self.alpha) & (plot_df[repl_name + " " + self.signif_column] <= self.alpha), :],
            fig=fig,
            ax=axes[0, 2],
            x=disc_name + " log beta",
            y=repl_name + " log beta",
            label=disc_name + " label",
            xlabel=disc_name + " log eQTL beta",
            ylabel="",
            title="Both signif.",
            color=color,
            include_ylabel=False
        )
        self.update_limits(xlim, ylim, 0, 2)

        for (m, n), ax in np.ndenumerate(axes):
            (xmin, xmax) = self.shared_xlim[n]
            (ymin, ymax) = self.shared_ylim[m]

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
            filename = "{}_{}_replication_in_{}_{}.{}".format(disc_name, disc_cell_type, repl_name, repl_cell_type, extension)
            fig.savefig(os.path.join(self.outdir, filename))
            print("\tSaved {}".format(filename))
        plt.close()

    def prepare_rb_info(self, df, name):
        se_column = name + " se"
        if se_column not in df.columns:
            self.pvalue_to_zscore(df=df,
                                  beta_col=name + " beta",
                                  p_col=name + " pvalue",
                                  zscore_col=name + " derived z-score")
            self.zscore_to_beta(df=df,
                                z_col=name + " derived z-score",
                                maf_col=name + " MAF",
                                n_col=name + " N",
                                beta_col=name + " derived beta",
                                se_col=name + " derived se")

            return name + " derived beta", name + " derived se"
        return name + " beta", se_column

    @staticmethod
    def pvalue_to_zscore(df, beta_col, p_col, zscore_col="derived z-score"):
        p_values = df[p_col].to_numpy()
        zscores = stats.norm.ppf(p_values / 2)
        mask = np.ones_like(p_values)
        mask[df[beta_col] > 0] = -1
        df[zscore_col] = zscores * mask
        df.loc[df[p_col] == 1, zscore_col] = 0
        df.loc[df[p_col] == 0, zscore_col] = -40.

    @staticmethod
    def zscore_to_beta(df, z_col, maf_col, n_col, beta_col="derived beta", se_col="derived se"):
        chi = df[z_col] * df[z_col]
        a = 2 * df[maf_col] * (1 - df[maf_col]) * (df[n_col] + chi)
        df[beta_col] = df[z_col] / a ** (1/2)
        df[se_col] = 1 / a ** (1/2)

    @staticmethod
    def log_modulus_beta(series):
        s = series.copy()
        data = []
        for index, beta in s.T.items():
            data.append(np.log(abs(beta) + 1) * np.sign(beta))
        new_df = pd.Series(data, index=s.index)

        return new_df

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

    def calculate_p1(self, p):
        robjects.r("source('{}')".format(self.qvalues_script))
        p = robjects.FloatVector(p)
        qvalue_truncp = robjects.globalenv['qvalue_truncp']
        pi0 = qvalue_truncp(p)[0]
        return 1 - np.array(pi0)

    def calculate_rb(self, b1, se1, b2, se2, theta=0):
        robjects.r("source('{}')".format(self.rb_script))
        b1 = robjects.FloatVector(b1)
        se1 = robjects.FloatVector(se1)
        b2 = robjects.FloatVector(b2)
        se2 = robjects.FloatVector(se2)
        calcu_cor_true = robjects.globalenv['calcu_cor_true']
        rb = calcu_cor_true(b1, se1, b2, se2, theta)
        return np.array(rb)[0]

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

##############################################################################################################

class Dataset:
    def __init__(self, name, path):
        self.class_name = None
        self.name = name
        self.path = path
        self.columns = {
            "gene": None,
            "SNP": None,
            "EA": None,
            "OA": None,
            "beta": None,
            "se": None,
            "pvalue": None,
            "N": None,
            "MAF": None,
            "label": None
        }

        # Set the class variables.
        self.snp_pos_df = None

        # Class variables.
        self.overlapping_variants = None

    def get_class_name(self):
        return self.class_name

    def get_name(self):
        return self.name

    def get_path(self):
        return self.path

    @staticmethod
    def translate_cell_type(cell_type):
        return cell_type

    @abstractmethod
    def get_snp_pos_df(self):
        print("Error, function 'get_snp_pos_df()' not implemented in {}".format(self.class_name))
        exit()

    def get_variants(self):
        df = self.get_snp_pos_df()
        return set(df.index)

    @abstractmethod
    def set_overlapping_variants(self, overlapping_variants):
        print("Error, function 'set_overlapping_variants()' not implemented in {}".format(self.class_name))
        exit()

    def get_all_effects(self, cell_type):
        df = self.get_effects(cell_type, mode="all")
        return self.standardize_format(df, cell_type)

    def get_all_overlapping_effects(self, cell_type):
        if self.overlapping_variants is None:
            print("Error, overlapping variants not set.")
            exit()
        df = self.get_effects(cell_type, mode="all-overlapping")
        return self.standardize_format(df, cell_type)

    def get_top_effects(self, cell_type):
        df = self.get_effects(cell_type, mode="top")
        return self.standardize_format(df, cell_type)

    def get_top_overlapping_effects(self, cell_type):
        if self.overlapping_variants is None:
            print("Error, overlapping variants not set.")
            exit()
        df = self.get_effects(cell_type, mode="top-overlapping")
        return self.standardize_format(df, cell_type)

    def get_specific_effects(self, cell_type, effects):
        df = self.get_effects(cell_type, effects=effects, mode="specific")
        return self.standardize_format(df, cell_type)

    @abstractmethod
    def get_effects(self, df, effects=None, mode="top"):
        print("Error, function 'get_effects()' not implemented in {}".format(self.class_name))
        exit()

    def standardize_format(self, df, cell_type):
        if df is None:
            print("Error, empty data frame.")
            exit()

        columns_of_interest = []
        new_column_names = []
        for key, value in self.columns.items():
            if value is not None:
                if "CELL_TYPE" in value:
                    value = value.replace("CELL_TYPE", self.translate_cell_type(cell_type))
                columns_of_interest.append(value)
                new_column_names.append("{} {}".format(self.name, key))

        df = df[columns_of_interest].copy()
        df.columns = new_column_names

        dtypes = {
            self.name + " gene": str,
            self.name + " SNP": str,
            self.name + " EA": str,
            self.name + " OA": str,
            self.name + " beta": float,
            self.name + " se": float,
            self.name + " pvalue": float,
            self.name + " N": int,
            self.name + " MAF": float,
            self.name + " label": str
        }
        df = df.astype({key:value for key, value in dtypes.items() if key in df.columns})
        return df

    def get_eqtls_from_summary_stats(self, df):
        return dict(zip(df[self.name + " gene"], df[self.name + " SNP"]))

    @staticmethod
    def load_file(inpath, header=0, index_col=None, sep="\t", low_memory=True,
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
    def load_partial_file(inpath, header=0, index_col=None, sep="\t", nrows=None,
                          skiprows=None, usecols=None, split_cols=None, filter=None, top_effect_cols=None,
                          specific_effect_cols=None, effects=None):
        print("Loading partial file...")
        if inpath.endswith(".gz"):
            fhi = gzip.open(inpath, 'rt')
        else:
            fhi = open(inpath, 'r')

        top_effects = {}

        columns = None
        indices = None if index_col is None else []
        lines = []
        header_row = header if skiprows is None else skiprows + header
        max_rows = nrows
        if skiprows is not None and max_rows is not None:
            max_rows += skiprows
        for i, line in enumerate(fhi):
            if skiprows is not None and i < skiprows:
                continue
            if max_rows is not None and i > max_rows:
                break

            values = line.strip("\n").split(sep)
            if split_cols is not None:
                splitted_values = []
                for col_index, value in enumerate(values):
                    if col_index in split_cols.keys():
                        splitted_values.extend(values[col_index].split(split_cols[col_index]))
                    else:
                        splitted_values.append(values[col_index])
                values = splitted_values

            if index_col is not None:
                indices.append(values[index_col])
                del values[index_col]

            if usecols is not None:
                values = [values[i] for i in usecols]

            if i == header_row:
                columns = values
                continue
            if filter is not None and not filter(values):
                continue
            if specific_effect_cols is not None and effects is not None:
                gene = values[specific_effect_cols["gene"]]
                snp = values[specific_effect_cols["snp"]]
                if gene not in effects or effects[gene] != snp:
                    continue

            if top_effect_cols is None:
                lines.append(values)
            else:
                key = values[top_effect_cols["key"]]
                value = values[top_effect_cols["value"]]
                if key not in top_effects:
                    top_effects[key] = (value, values)
                elif key in top_effects and value < top_effects[key][0]:
                    top_effects[key] = (value, values)
                else:
                    pass

        if top_effect_cols is not None:
            for _, (_, values) in top_effects.items():
                lines.append(values)

        df = pd.DataFrame(lines, columns=columns, index=indices)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df


##############################################################################################################

class Bryois(Dataset):
    def __init__(self, name, path):
        super().__init__(name=name, path=path)
        self.class_name = "Bryois"

        # Define the important columns.
        self.columns = {
            "gene": "HGNC",
            "SNP": "SNP",
            "EA": "effect_allele",
            "OA": "other_allele",
            "beta": "Beta",
            "se": None,
            "pvalue": "Nominal p-value",
            "N": "N",
            "MAF": "MAF",
            "label": "HGNC"
        }

    def get_snp_pos_df(self):
        if self.snp_pos_df is None:
            snp_pos_df = self.load_snp_pos_df()
            snp_pos_df.index = snp_pos_df["SNP_id_hg38"]
            snp_pos_df.index.name = None
            self.snp_pos_df = snp_pos_df

        return self.snp_pos_df

    def load_snp_pos_df(self):
        snp_pos_inpath = os.path.join(self.path, "snp_pos.txt.gz")
        if not os.path.exists(snp_pos_inpath):
            print("Error, '{}' file not found".format(os.path.basename(snp_pos_inpath)))
            exit()

        return self.load_file(inpath=snp_pos_inpath)

    def set_overlapping_variants(self, overlapping_variants):
        snp_pos_df = self.get_snp_pos_df()
        trans_dict = dict(zip(snp_pos_df["SNP_id_hg38"], snp_pos_df["SNP"]))
        self.overlapping_variants = set([trans_dict[snp_id] for snp_id in overlapping_variants])

    @staticmethod
    def translate_cell_type(cell_type):
        ct_link = {
            "AST": "Astrocytes",
            "END": "Endothelial",
            "EX": "Excitatory.neurons",
            "IN": "Inhibitory.neurons",
            "MIC": "Microglia",
            "OPC": "OPCs...COPs",
            "OLI": "Oligodendrocytes",
            "PER": "Pericytes"
        }
        return ct_link[cell_type]

    def get_effects(self, cell_type, effects=None, mode="top"):
        snp_pos_df = self.get_snp_pos_df()

        df_list = []
        for chromosome in CHROMOSOMES:
            inpath = os.path.join(self.path, "{}.{}.gz".format(self.translate_cell_type(cell_type), chromosome))
            if not os.path.exists(inpath):
                print("Error, '{}' file not found".format(os.path.basename(inpath)))
                exit()

            df = None
            if mode == "all":
                df = self.load_file(inpath, header=None, sep=" ")
            elif mode == "all-overlapping":
                df = self.load_partial_file(inpath, header=None, sep=" ", split_cols={0: "_"}, filter=lambda x: x[1] in self.overlapping_variants)
            elif mode == "top":
                df = self.load_partial_file(inpath, header=None, sep=" ", split_cols={0: "_"}, top_effect_cols={"key": 0, "value": 3})
            elif mode == "top-overlapping":
                df = self.load_partial_file(inpath, header=None, sep=" ", split_cols={0: "_"}, top_effect_cols={"key": 0, "value": 3}, filter=lambda x: x[1] in self.overlapping_variants)
            elif mode == "specific":
                df = self.load_partial_file(inpath, header=None, sep=" ", split_cols={0: "_"}, effects=effects, specific_effect_cols={"gene": 0, "snp": 2})
            else:
                print("Error, mode '{}' not implemented in {}".format(mode, self.class_name))
                exit()

            df_list.append(df)

        if len(df_list) == 0:
            return None

        df = pd.concat(df_list, axis=0)
        df.columns = ["HGNC", "Gene", "SNP", "Distance to TSS", "Nominal p-value", "Beta"]
        df = df.astype({"HGNC": str, "Gene": str, "SNP": str, "Distance to TSS": int, "Nominal p-value": float, "Beta": float})
        df = df.merge(snp_pos_df, on="SNP", how="left")
        df["N"] = 196
        df.index = df[self.columns["gene"]] + "_" + df[self.columns["SNP"]]
        return df


##############################################################################################################


class scMetaBrain(Dataset):
    def __init__(self, name, path, ancestry, cell_level):
        super().__init__(name=name, path=path)
        self.class_name = "scMetaBrain"
        self.ancestry = ancestry
        self.cell_level = cell_level

        # Set the class variables.
        self.qtl_chunks = None

        # Define the important columns.
        self.columns = {
            "gene": "feature_id",
            "SNP": "snp_id",
            "EA": "assessed_allele",
            "OA": None,
            "beta": "beta",
            "se": "beta_se",
            "pvalue": "p_value",
            "N": "n_samples",
            "MAF": "maf",
            "label": "feature_id"
        }

    def get_ancestry(self):
        return self.ancestry

    def get_snp_pos_df(self):
        if self.snp_pos_df is None:
            snp_pos_df = self.load_snp_pos_df()
            snp_pos_df.index = "chr" + snp_pos_df["CHR"].astype(str) + ":" + snp_pos_df["POS"].astype(str)
            snp_pos_df.index.name = None
            self.snp_pos_df = snp_pos_df

        return self.snp_pos_df

    def load_snp_pos_df(self):
        snp_pos_inpath = os.path.join(self.path, "genotype_input/{ancestry}/{ancestry}_imputed_hg38_stats_filtered.vars.gz".format(ancestry=self.ancestry))
        if not os.path.exists(snp_pos_inpath):
            print("Error, '{}' file not found".format(os.path.basename(snp_pos_inpath)))
            exit()

        return self.load_partial_file(inpath=snp_pos_inpath, header=0, usecols=[0, 1, 2, 3], filter=lambda x: len(x[3]) == 3)

    def set_overlapping_variants(self, overlapping_variants):
        snp_pos_df = self.get_snp_pos_df()
        trans_dict = dict(zip(snp_pos_df.index, snp_pos_df["ID"]))
        self.overlapping_variants = set([trans_dict[snp_id] for snp_id in overlapping_variants])

    def get_qtl_chunks(self, cell_type):
        if self.qtl_chunks is None:
            self.qtl_chunks = self.load_qtl_chunks(cell_type)

        return self.qtl_chunks

    def load_qtl_chunks(self, cell_type):
        qtl_results = glob.glob(os.path.join(self.path, "output", self.ancestry, self.cell_level, cell_type, "qtl", "qtl_results_*.h5"))
        qtl_results_data = []
        for filepath in qtl_results:
            basename = os.path.basename(filepath)
            match = re.match("qtl_results_([0-9]{1,2}|X|Y|MT)_([0-9]+)_([0-9]+).h5", basename)
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

    @staticmethod
    def natural_keys(text):
        return [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', text)]

    def get_top_effects(self, cell_type):
        inpath = os.path.join(self.path, "output", self.ancestry, self.cell_level, cell_type, "top_qtl_results_all.txt.gz")
        if os.path.exists(inpath):
            df = self.load_file(inpath, index_col=0)
            df.index = df[self.columns["gene"]] + "_" + df[self.columns["SNP"]]
        else:
            df = self.get_effects(cell_type, mode="top")
        return self.standardize_format(df, cell_type)

    def get_effects(self, cell_type, effects=None, mode="top"):
        qlt_chunks_df = self.get_qtl_chunks(cell_type=cell_type)
        print(qlt_chunks_df)

        df_list = []
        for index, row in qlt_chunks_df.iterrows():
            df = None
            if mode == "top":
                df = self.load_h5_file(cell_type=cell_type, chunk=row["Chunk"], top=True)
            elif mode == "specific":
                df = self.load_h5_file(cell_type=cell_type, chunk=row["Chunk"], effects=effects)
            else:
                print("Error, mode '{}' not implemented in {}".format(mode, self.class_name))
                exit()

            if df is not None:
                df_list.append(df)

        if len(df_list) == 0:
            return None

        df =  pd.concat(df_list, axis=0)
        df.index = df[self.columns["gene"]] + "_" + df[self.columns["SNP"]]
        return df

    def load_h5_file(self, cell_type, chunk, effects=None, top=False):
        print("Loading chunk...")
        feature_metadata_file = os.path.join(self.path, "output", self.ancestry, self.cell_level, cell_type, "qtl", "feature_metadata_{}.txt.gz".format(chunk))
        snp_metadata_file = os.path.join(self.path, "output", self.ancestry, self.cell_level, cell_type, "qtl", "snp_metadata_{}.txt.gz".format(chunk))
        h5_file = os.path.join(self.path, "output", self.ancestry, self.cell_level, cell_type, "qtl", "qtl_results_{}.h5".format(chunk))

        try:
            if os.path.exists(feature_metadata_file):
                ffea_df = pd.read_table(feature_metadata_file, sep="\t")
            elif os.path.exists(feature_metadata_file + ".gz"):
                ffea_df = pd.read_table(feature_metadata_file + ".gz", sep="\t")
            else:
                print("  Warning, skipping: '{}' missing feature metadata file.".format(chunk))
                return None
        except pd.errors.EmptyDataError:
            print("  Error, issue in feature annotation. Skipping '{}'.".format(chunk))
            return None

        if effects is not None:
            chunk_overlapping_features = set(effects.keys()).intersection(set(ffea_df["feature_id"].values.tolist()))
            if len(chunk_overlapping_features) == 0:
                print("  Warning, skipping: '{}' no overlapping features.".format(chunk))
                return None

        try:
            if os.path.exists(snp_metadata_file):
                fsnp_df = pd.read_table(snp_metadata_file, sep="\t")
            elif os.path.exists(snp_metadata_file + ".gz"):
                fsnp_df = pd.read_table(snp_metadata_file + ".gz", sep="\t")
            else:
                print("  Warning, skipping: '{}' missing SNP metadata file.".format(chunk))
                return None
        except pd.errors.EmptyDataError:
            print("  Error, issue in snp annotation. Skipping '{}'.".format(chunk))
            return None

        if effects is not None:
            chunk_overlapping_snps = set(effects.values()).intersection(set(fsnp_df["snp_id"].values.tolist()))
            if len(chunk_overlapping_snps) == 0:
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
            if effects is not None and frez_key not in effects.keys():
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
                snp_id = effects[frez_key]
                if snp_id not in frez_df["snp_id"].values:
                    continue
                frez_df = frez_df.loc[frez_df["snp_id"] == snp_id ,:]

            df_list.append(frez_df)
            del frez_df

        if len(df_list) == 0:
            print("  Warning, skipping: '{}' no overlapping feature-SNP combinations.".format(chunk))
            return None

        df = pd.concat(df_list, axis=0)
        del df_list

        df = pd.merge(df, ffea_df, on='feature_id', how='left')

        if(len(glob.glob(os.path.join(self.path, "output", self.ancestry, self.cell_level, cell_type, "qtl", "snp_qc_metrics_naContaining_feature_*.txt"))) > 0):
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

        if top:
            df = df.sort_values(by=['empirical_feature_p_value', "p_value"], ascending=[True, True])
            df = df.groupby(df['feature_id']).first()
            df.reset_index(drop=False, inplace=True)

        del ffea_df, fsnp_df

        print("\tLoaded chunk: {} "
              "with shape: {}".format(chunk,
                                      df.shape))
        return df

##############################################################################################################


class mbQTL(Dataset):
    def __init__(self, name, path):
        super().__init__(name=name, path=path)
        self.class_name = "mbQTL"

        # Define the important columns.
        self.columns = {
            "gene": "Gene",
            "SNP": "SNP",
            "EA": "SNPEffectAllele",
            "OA": "SNPOtherAllele",
            "beta": "MetaBeta",
            "se": "MetaSE",
            "pvalue": "MetaP",
            "N": None,
            "MAF": None,
            "label": "Gene"
        }

    @staticmethod
    def translate_cell_type(cell_type):
        return cell_type.capitalize()

    def get_effects(self, cell_type, effects=None, mode="top"):
        inpath = os.path.join(self.path, "{}-TopEffects.txt".format(cell_type.lower()))
        if not os.path.exists(inpath):
            print("Error, '{}' file not found".format(os.path.basename(inpath)))
            exit()

        df = None
        if mode == "all":
            df = self.load_file(inpath)
        elif mode == "top":
            df = self.load_partial_file(inpath, top_effect_cols={"key": 0, "value": 11})
        elif mode == "specific":
            df = self.load_partial_file(inpath, effects=effects, specific_effect_cols={"gene": 0, "snp": 5})
        else:
            print("Error, mode '{}' not implemented in {}".format(mode, self.class_name))
            exit()

        # Add the other allele column.
        df[["AlleleA", "AlleleB"]] = df["SNPAlleles"].str.split("/", n=1, expand=True)
        df.insert(10, "SNPOtherAllele", df["AlleleA"])
        mask = df["SNPEffectAllele"] == df["SNPOtherAllele"]
        df.loc[mask, "SNPOtherAllele"] = df.loc[mask, "AlleleB"]
        df.drop(["AlleleA", "AlleleB"], axis=1, inplace=True)

        df.index = df[self.columns["gene"]] + "_" + df[self.columns["SNP"]]
        return df


if __name__ == '__main__':
    m = main()
    m.start()
