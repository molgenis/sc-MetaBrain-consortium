#!/usr/bin/env python3

"""
File:         compare_wg3_to_mbqtl.py
Created:      2023/04/27
Last Changed:
Author:       M.Vochteloo

Copyright (C) 2022 University Medical Center Groningen.
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
# from rpy2.robjects.packages import importr
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adjustText import adjust_text

# Local application imports.

"""
Syntax:
./compare_wg3_to_mbqtl.py \
    --work_dir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-25-mbQTLBryoisRepl \
    --dataset_outdir 2023-04-27-Mathys2019WG3_vs_mbQTL \
    --mbqtl_folder /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-25-mbQTLBryoisRepl/2023-04-25-BryoisRepl \
    --wg3_folder /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-11-WorkGroup3eQTLAndDEA/2023-04-12-Mathys2019 \
    --qvalues /groups/umcg-biogen/tmp02/users/umcg-mvochteloo/utils/qvalue_truncp.R \
    --rb /groups/umcg-biogen/tmp02/users/umcg-mvochteloo/utils/Rb.R
    
"""

# Metadata
__program__ = "Compare Workgroup 3 to mbQTL"
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
        self.mbqtl_folder = getattr(arguments, 'mbqtl_folder')
        self.extensions = getattr(arguments, 'extensions')
        self.verbose = getattr(arguments, 'verbose')
        self.qvalues_script = getattr(arguments, 'qvalues')
        self.rb_script = getattr(arguments, 'rb')

        # Pre-process the dataset output directory.
        date_str = datetime.now().strftime("%Y-%m-%d")
        if not re.match("^[0-9]{4}-(0[1-9]|1[0-2])-(0[1-9]|[1-2][0-9]|3[0-1])", dataset_outdir):
            dataset_outdir = "{}-{}".format(date_str, dataset_outdir)

        self.dataset_outdir = dataset_outdir
        self.dataset_plot_outdir = os.path.join(self.dataset_outdir, "plot")

        for dir in [dataset_outdir, self.dataset_outdir, self.dataset_plot_outdir]:
            if not os.path.exists(dir):
                os.makedirs(dir)

        ########################################################################

        self.cell_types = ["AST", "EX", "IN", "MIC", "OLI", "OPC", "PER"]

        self.palette = {
            "AST": "#D55E00",
            "EX": "#0072B2",
            "IN": "#56B4E9",
            "MIC": "#E69F00",
            "OLI": "#009E73",
            "OPC": "#F0E442",
            "PER": "#808080"
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
                            help="The working directory.")
        parser.add_argument("--dataset_outdir",
                            type=str,
                            required=True,
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
        parser.add_argument("--mbqtl_folder",
                            type=str,
                            required=True,
                            help="The path to the outputs from Meta-Beta-QTL "
                                 "output")
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

        plot_data = {}
        for cell_type in self.cell_types:
            print("  Working on '{}'".format(cell_type))

            # Load workgroup 3 data.
            wg3_file = os.path.join(self.wg3_folder, "output", self.annotation_level, cell_type, "qtl_results_BryoiseQTLs.txt.gz")
            if not os.path.exists(wg3_file):
                continue
            wg3_df = self.load_file(wg3_file, header=0, index_col=None)

            # Load the mbQTL data.
            mbqtl_file = os.path.join(self.mbqtl_folder, "{}-TopEffects.txt".format(cell_type))
            if not os.path.exists(mbqtl_file):
                continue
            mbqtl_df = self.load_file(mbqtl_file, header=0, index_col=None)

            # merge.
            df = wg3_df.merge(mbqtl_df, left_on="feature_id", right_on="Gene", how="inner")

            # flip.
            df["flip"] = df["SNPEffectAllele"] != df["assessed_allele"]
            df["MetaBeta"] = df["MetaBeta"] * df["flip"].map({True: -1, False: 1})

            # save.
            # self.save_file(df=df, outpath=os.path.join(self.dataset_outdir, "{}.txt.gz".format(cell_type)), header=True, index=False)

            plot_data[cell_type] = df

        print("\nVisualizing comparison")
        _ = self.visualise_data(plot_data=plot_data, interest=" log beta")
        _ = self.visualise_data(plot_data=plot_data, interest=" log se")
        _ = self.visualise_data(plot_data=plot_data, interest=" log10 pvalue")


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

    def save_file(self, df, outpath, header=True, index=False, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)

        if self.verbose:
            print("\tSaved dataframe: {} "
                  "with shape: {}".format(os.path.basename(outpath),
                                          df.shape))

    def visualise_data(self, plot_data, interest):
        cell_types = list(plot_data.keys())
        cell_types.sort()

        nrows = 4
        ncols = len(cell_types)

        self.shared_ylim = {i: (0, 0.5) for i in range(nrows)}
        self.shared_xlim = {i: (0, 0.5) for i in range(ncols)}

        replication_stats = []

        sns.set(rc={'figure.figsize': (ncols * 8, nrows * 5.3553)})
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
            plot_df = df.loc[:, ["feature_id",
                                 "beta",
                                 "beta_se",
                                 "p_value",
                                 "empirical_feature_p_value",
                                 "MetaBeta",
                                 "MetaSE",
                                 "MetaP",
                                 ]].copy()
            plot_df.columns = ["Gene symbol",
                               "WG3 beta",
                               "WG3 se",
                               "WG3 pvalue",
                               "WG3 perm pvalue",
                               "mbQTL beta",
                               "mbQTL se",
                               "mbQTL pvalue",
                               ]
            plot_df = plot_df.loc[(~plot_df["WG3 pvalue"].isna()) & (~plot_df["mbQTL pvalue"].isna()), :]
            plot_df.sort_values(by="WG3 pvalue", inplace=True)

            # Convert the interaction beta to log scale.
            plot_df["WG3 log beta"] = [np.log(abs(beta) + 1) * np.sign(beta) for beta in plot_df["WG3 beta"]]
            plot_df["mbQTL log beta"] = [np.log(abs(beta) + 1) * np.sign(beta) for beta in plot_df["mbQTL beta"]]

            plot_df["WG3 log se"] = [np.log(abs(beta) + 1) * np.sign(beta) for beta in plot_df["WG3 se"]]
            plot_df["mbQTL log se"] = [np.log(abs(beta) + 1) * np.sign(beta) for beta in plot_df["mbQTL se"]]

            plot_df["WG3 log10 pvalue"] = -np.log10(plot_df["WG3 pvalue"])
            plot_df["mbQTL log10 pvalue"] = -np.log10(plot_df["mbQTL pvalue"])

            include_ylabel = False
            if col_index == 0:
                include_ylabel = True

            if col_index == 0:
                for row_index, panel in enumerate(["A", "B", "C", "D"]):
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
                x="WG3{}".format(interest),
                y="mbQTL{}".format(interest),
                xlabel="",
                ylabel="mbQTL{}".format(interest),
                title=cell_type,
                color=self.palette[cell_type],
                include_ylabel=include_ylabel
            )
            self.update_limits(xlim, ylim, 0, col_index)

            print("\tPlotting row 2.")
            xlim, ylim, stats2 = self.scatterplot(
                df=plot_df.loc[plot_df["WG3 pvalue"] <= 0.05, :],
                fig=fig,
                ax=axes[1, col_index],
                x="WG3{}".format(interest),
                y="mbQTL{}".format(interest),
                xlabel="",
                ylabel="mbQTL{}\nWG3 pvalue <= 0.05".format(interest),
                title="",
                color=self.palette[cell_type],
                include_ylabel=include_ylabel,
                pi1_column="mbQTL pvalue",
                rb_columns=[("mbQTL beta", "mbQTL se"), ("WG3 beta", "WG3 se")]
            )
            self.update_limits(xlim, ylim, 1, col_index)

            print("\tPlotting row 3.")
            xlim, ylim, stats2 = self.scatterplot(
                df=plot_df.loc[plot_df["mbQTL pvalue"] <= 0.05, :],
                fig=fig,
                ax=axes[2, col_index],
                x="WG3{}".format(interest),
                y="mbQTL{}".format(interest),
                xlabel="",
                ylabel="mbQTL{}\nmbQTL pvalue <= 0.05".format(interest),
                title="",
                color=self.palette[cell_type],
                include_ylabel=include_ylabel,
                pi1_column="WG3 pvalue",
                rb_columns=[("WG3 beta", "WG3 se"), ("mbQTL beta", "mbQTL se")]
            )
            self.update_limits(xlim, ylim, 2, col_index)

            print("\tPlotting row 4.")
            xlim, ylim, stats3 = self.scatterplot(
                df=plot_df.loc[(plot_df["WG3 pvalue"] <= 0.05) & (plot_df["mbQTL pvalue"] <= 0.05), :],
                fig=fig,
                ax=axes[3, col_index],
                x="WG3{}".format(interest),
                y="mbQTL{}".format(interest),
                label="Gene symbol",
                xlabel="WG3{}".format(interest),
                ylabel="mbQTL{}".format(interest),
                title="",
                color=self.palette[cell_type],
                include_ylabel=include_ylabel
            )
            self.update_limits(xlim, ylim, 3, col_index)
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

            ax.set_xlim(xmin - xmargin - 0.2, xmax + xmargin)
            ax.set_ylim(ymin - ymargin, ymax + ymargin)

        # Add the main title.
        fig.suptitle("Workgroup 3 replication with mbQTL",
                     fontsize=40,
                     color="#000000",
                     weight='bold')

        for extension in self.extensions:
            fig.savefig(os.path.join(self.dataset_plot_outdir, "replication_plot_{}.{}".format(interest.replace(" ", ""), extension)))
        plt.close()

        # Construct the replication stats data frame.
        replication_stats_df = pd.concat(replication_stats, axis=0)
        replication_stats_df.dropna(inplace=True)

        return replication_stats_df


    def plot_heatmap(self, df, xlabel="", ylabel=""):
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
                            square=True, annot=df.round(2), fmt='',
                            cbar=False, annot_kws={"size": 16, "color": "#000000"},
                            ax=ax)

                plt.setp(ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=20, rotation=0))
                plt.setp(ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=20, rotation=90))

                ax.set_xlabel(xlabel, fontsize=14)
                ax.xaxis.set_label_position('top')

                ax.set_ylabel(ylabel, fontsize=14)
                ax.yaxis.set_label_position('right')
            else:
                ax.set_axis_off()

            col_index += 1
            if col_index > 1:
                col_index = 0
                row_index += 1

        for extension in self.extensions:
            fig.savefig(os.path.join(self.dataset_plot_outdir, "replication_stats.{}".format(extension)))
        plt.close()

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
                    rb = min(rb_est[0], 1)


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

    # @staticmethod
    # def qvalues(p):
    #     qvalue = importr("qvalue")
    #     pvals = robjects.FloatVector(p)
    #     qobj = robjects.r['qvalue'](pvals)
    #     return np.array(qobj.rx2('qvalues'))

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

    def print_arguments(self):
        print("Arguments:")
        print("  > Working directory:     {}".format(self.work_dir))
        print("  > Workgroup 3 folder:    {}".format(self.wg3_folder))
        print("  > Annotation level:      {}".format(self.annotation_level))
        print("  > mbQTL data folder:     {}".format(self.dataset_outdir))
        print("  > Plot extensions:       {}".format(", ".join(self.extensions)))
        print("  > Verbose:               {}".format(self.verbose))
        print("  > Qvalues script:        {}".format(self.qvalues_script))
        print("  > Rb script:             {}".format(self.rb_script))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
