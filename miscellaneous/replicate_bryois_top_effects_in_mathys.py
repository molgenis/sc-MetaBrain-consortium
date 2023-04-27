#!/usr/bin/env python3

"""
File:         replicate_bryois_top_effects_in_mathys.py
Created:      2023/04/25
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
./replicate_bryois_top_effects_in_mathys.py \
    --work_dir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-25-mbQTLBryoisRepl \
    --dataset_outdir 2023-04-25-Mathys2019 \
    --mbqtl_folder /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-25-mbQTLBryoisRepl/2023-04-25-BryoisRepl \
    --bryois_folder /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-20-ReplicateInBryois/Bryois2022 \
    --qvalues /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-20-ReplicateInBryois/qvalue_truncp.R \
    --rb /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-20-ReplicateInBryois/Rb.R
    
"""

# Metadata
__program__ = "Replicate Bryois et al. 2022 in Mathys et al. 2019"
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
        self.mbqtl_folder = getattr(arguments, 'mbqtl_folder')
        self.bryois_folder = getattr(arguments, 'bryois_folder')
        self.bryois_n = 196
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

        self.bryois_ct_dict = {
            "Astrocytes": "AST",
            "Excitatory neurons": "EX",
            "Inhibitory neurons": "IN",
            "Microglia": "MIC",
            "Oligodendrocytes": "OLI",
            "OPCs / COPs": "OPC",
            "Pericytes": "PER"
        }

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
        parser.add_argument("--mbqtl_folder",
                            type=str,
                            required=True,
                            help="The path to the outputs from Meta-Beta-QTL "
                                 "output")
        parser.add_argument("--bryois_folder",
                            type=str,
                            required=True,
                            help="The path to the Bryois et al. 2022 data")
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

        print("Loading data")
        bryois_top_effects = self.load_file(os.path.join(self.bryois_folder, "BryoisTopEffects.txt.gz"), header=0, index_col=None)

        plot_data = {}
        for bryois_ct, cell_type in self.bryois_ct_dict.items():
            print("  Working on '{}'".format(cell_type))

            # Load Bryois et al. 2022 data.
            bryois_df = self.load_file(os.path.join(self.bryois_folder, "merged", "{}.txt.gz".format(cell_type)), header=0, index_col=0)
            bryois_df.index = bryois_df["ENSG"] + "_" + bryois_df["SNP_id"]

            top_effects_df = bryois_top_effects.loc[bryois_top_effects["cell_type"] == bryois_ct, :].copy()
            top_effects_df.index = top_effects_df["ensembl"] + "_" + top_effects_df["SNP"]

            bryois_df = bryois_df.merge(top_effects_df[["adj_p"]], left_index=True, right_index=True, how="inner")

            # Load the Mathys et al. 2019 data.
            mathys_df = self.load_file(os.path.join(self.mbqtl_folder, "{}-TopEffects.txt".format(cell_type)), header=0, index_col=None)

            # merge.
            df = bryois_df.merge(mathys_df, left_on="HGNC", right_on="Gene", how="inner")

            # flip.
            df["flip"] = df["SNPEffectAllele"] != df["effect_allele"]
            df["MetaBeta"] = df["MetaBeta"] * df["flip"].map({True: -1, False: 1})

            df["BH-FDR"] = np.nan
            ieqtl_repl_mask = np.logical_and((df["adj_p"] <= 0.05).to_numpy(), (~df["MetaP"].isna()).to_numpy())
            n_overlap = np.sum(ieqtl_repl_mask)
            if n_overlap > 1:
                df.loc[ieqtl_repl_mask, "BH-FDR"] = multitest.multipletests(df.loc[ieqtl_repl_mask, "MetaP"], method='fdr_bh')[1]

            plot_data[cell_type] = df

        print("\nVisualizing comparison")
        replication_stats_df = self.visualise_data(plot_data=plot_data)

        plot_repl_stats_df = replication_stats_df.loc[replication_stats_df["label"] == "discovery significant", :].pivot(index="variable", values="value", columns="cell type")
        plot_repl_stats_df = plot_repl_stats_df.iloc[2:, :]
        plot_repl_stats_df.index.name = None
        plot_repl_stats_df.columns.name = None
        plot_repl_stats_df.loc["concordance", :] = plot_repl_stats_df.loc["concordance", :] / 100
        self.plot_heatmap(df=plot_repl_stats_df)
        del plot_repl_stats_df

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
            df["bryois_n_samples"] = self.bryois_n
            plot_df = df.loc[:, ["HGNC",
                                 "Nominal p-value",
                                 "Beta",
                                 "adj_p",
                                 "bryois_n_samples",
                                 "MetaP",
                                 "MetaBeta",
                                 "MetaSE",
                                 "BH-FDR",
                                 "SNPEffectAlleleFreq",
                                 ]].copy()
            plot_df.columns = ["Gene symbol",
                               "Bryois pvalue",
                               "Bryois beta",
                               "Bryois FDR",
                               "Bryois N",
                               "Mathys pvalue",
                               "Mathys beta",
                               "Mathys se",
                               "Mathys FDR",
                               "Mathys MAF"
                               ]
            plot_df = plot_df.loc[(~plot_df["Bryois pvalue"].isna()) & (~plot_df["Mathys pvalue"].isna()), :]
            plot_df.sort_values(by="Bryois pvalue", inplace=True)

            # Calculate the replication standard error.
            self.pvalue_to_zscore(df=plot_df,
                                  beta_col="Bryois beta",
                                  p_col="Bryois pvalue",
                                  prefix="Bryois ")
            self.zscore_to_beta(df=plot_df,
                                z_col="Bryois z-score",
                                maf_col="Mathys MAF",
                                n_col="Bryois N",
                                prefix="Bryois zscore-to-")

            # Convert the interaction beta to log scale.
            plot_df["Mathys log beta"] = [np.log(abs(beta) + 1) * np.sign(beta) for beta in plot_df["Mathys beta"]]
            plot_df["Bryois log beta"] = [np.log(abs(beta) + 1) * np.sign(beta) for beta in plot_df["Bryois beta"]]

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
                y="Mathys log beta",
                xlabel="",
                ylabel="Mathys log beta",
                title=cell_type,
                color=self.palette[cell_type],
                include_ylabel=include_ylabel
            )
            self.update_limits(xlim, ylim, 0, col_index)

            print("\tPlotting row 2.")
            xlim, ylim, stats2 = self.scatterplot(
                df=plot_df.loc[plot_df["Bryois FDR"] <= 0.05, :],
                fig=fig,
                ax=axes[1, col_index],
                x="Bryois log beta",
                y="Mathys log beta",
                xlabel="",
                ylabel="Mathys log beta",
                title="",
                color=self.palette[cell_type],
                include_ylabel=include_ylabel,
                pi1_column="Mathys pvalue",
                rb_columns=[("Bryois zscore-to-beta", "Bryois zscore-to-se"), ("Mathys beta", "Mathys se")]
            )
            self.update_limits(xlim, ylim, 1, col_index)

            print("\tPlotting row 3.")
            xlim, ylim, stats3 = self.scatterplot(
                df=plot_df.loc[(plot_df["Bryois FDR"] <= 0.05) & (plot_df["Mathys FDR"] <= 0.05), :],
                fig=fig,
                ax=axes[2, col_index],
                x="Bryois log beta",
                y="Mathys log beta",
                label="Gene symbol",
                xlabel="Bryois log beta",
                ylabel="Mathys log beta",
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
        fig.suptitle("Bryois et al. 2022 replication in Mathys et al. 2019",
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
        print("  > Working directory:        {}".format(self.work_dir))
        print("  > MBeta-QTL data folder:    {}".format(self.dataset_outdir))
        print("  > Bryois data folder:       {}".format(self.bryois_folder))
        print("  > Plot extensions:          {}".format(", ".join(self.extensions)))
        print("  > Verbose:                  {}".format(self.verbose))
        print("  > Qvalues script:           {}".format(self.qvalues_script))
        print("  > Rb script:                {}".format(self.rb_script))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
