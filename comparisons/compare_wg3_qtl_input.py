#!/usr/bin/env python3

"""
File:         compare_wg3_qtl_input.py
Created:      2024/03/04
Last Changed:
Author:       M.Vochteloo

Copyright (C) 2020 University Medical Center Groningen.

A copy of the BSD 3-Clause "New" or "Revised" License can be found in the
LICENSE file in the root directory of this source tree.
"""

# Standard imports.
from __future__ import print_function
import argparse
import os

# Third party imports.
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Compare WG3 QTL input"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@rug.nl"
__license__ = "BSD (3-Clause)"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)

"""
Syntax: 
./compare_wg3_qtl_input.py -h
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.data_path1 = getattr(arguments, 'data1')
        self.name1 = " ".join(getattr(arguments, 'name1'))
        self.data_path2 = getattr(arguments, 'data2')
        self.name2 = " ".join(getattr(arguments, 'name2'))
        self.cell_type = getattr(arguments, 'cell_type')
        self.extensions = getattr(arguments, 'extensions')

        # Set variables.
        self.outdir = os.path.join(str(os.path.dirname(os.path.abspath(__file__))), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

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
                            help="show program's version number and exit.")
        parser.add_argument("-d1",
                            "--data1",
                            type=str,
                            required=True,
                            help="The path to the -d1 / --data1 directory.")
        parser.add_argument("-n1",
                            "--name1",
                            nargs="*",
                            type=str,
                            required=False,
                            default="",
                            help="The name of -d1 / --data2.")
        parser.add_argument("-d2",
                            "--data2",
                            type=str,
                            required=True,
                            help="The path to the -d2 / --data2 directory.")
        parser.add_argument("-n2",
                            "--name2",
                            nargs="*",
                            type=str,
                            required=False,
                            default="",
                            help="The name of -d2 / --data2.")
        parser.add_argument("-c",
                            "--cell_type",
                            type=str,
                            required=True,
                            help="")
        parser.add_argument("-e",
                            "--extensions",
                            type=str,
                            nargs="+",
                            default=["png"],
                            choices=["eps", "pdf", "pgf", "png", "ps", "raw", "rgba", "svg", "svgz"],
                            help="The output file format(s), default: ['png']")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        qtl_workdir1, reformat1 = self.get_qtl_workdir(inpath=self.data_path1)
        qtl_workdir2, reformat2 = self.get_qtl_workdir(inpath=self.data_path2)

        # CELL_TYPE.Exp.txt
        # CELL_TYPE.Qced.Normalized.SCs.Rds
        # CELL_TYPE.covariates.txt
        # CELL_TYPE.qtlInput.Pcs.txt
        # CELL_TYPE.qtlInput.txt

        for fillename_suffix in [".Qced.Normalized.SCs.Rds", ".covariates.txt", ".Exp.txt", ".qtlInput.txt", ".qtlInput.Pcs.txt", ".covariates.txt"]:
            filename_brief = fillename_suffix.rstrip(".txt").rstrip(".Rds").replace(".", "_")

            if fillename_suffix == ".Qced.Normalized.SCs.Rds":
                continue

            print("Loading '*{}'.".format(fillename_suffix))
            df1 = self.load_data(workdir=qtl_workdir1, fillename_suffix=fillename_suffix, index_col=0)
            print(df1)
            print(df1.columns.tolist())
            df2 = self.load_data(workdir=qtl_workdir2, fillename_suffix=fillename_suffix, index_col=0)
            print(df2)
            print(df2.columns.tolist())

            print(df1.equals(df2))

            if fillename_suffix in [".Exp.txt", ".qtlInput.Pcs.txt", ".qtlInput.txt"]:
                if fillename_suffix == ".qtlInput.Pcs.txt":
                    if reformat1:
                        df1.index = [index.split(";;")[0] + "_" + index.split(";;")[1].split("_")[0] for index in df1.index]
                    if reformat2:
                        df2.index = [index.split(";;")[0] + "_" + index.split(";;")[1].split("_")[0] for index in df2.index]

                corr_df, annot_df = self.correlate_dataframes(df1=df1, df2=df2)
                print(corr_df)
                print(annot_df)

                self.plot_heatmap(df=corr_df,
                                  annot_df=annot_df,
                                  xlabel=self.name1,
                                  ylabel=self.name2,
                                  filename="{}{}_Spearman_corr_heatmap".format(self.cell_type, filename_brief))

            break

    def get_qtl_workdir(self, inpath):
        old_covariates_indir = os.path.join(inpath, "input", "L1")
        if os.path.exists(old_covariates_indir):
            return old_covariates_indir, True

        new_covariates_indir = os.path.join(inpath, "expression_input", "EUR", "L1", self.cell_type, "PreQC")
        if os.path.exists(new_covariates_indir):
            return new_covariates_indir, False

    def load_data(self, workdir, fillename_suffix, index_col=None):
        inpath = os.path.join(workdir, self.cell_type + fillename_suffix)
        if not os.path.exists(inpath):
            print("Error, covariates file could not been found.")
            exit()

        if fillename_suffix.endswith(".txt"):
            return pd.read_csv(inpath, sep="\t", index_col=index_col)
        else:
            print("Error, file not supported.")
            exit()

    @staticmethod
    def correlate_dataframes(df1, df2, triangle=False, drop_paired_zero=False, method="Spearman"):
        corr_df = pd.DataFrame(np.nan, index=df1.columns, columns=df2.columns)
        annot_df = pd.DataFrame("", index=df1.columns, columns=df2.columns)

        for i, df1_colname in enumerate(df1.columns):
            for j, df2_colname in enumerate(df2.columns):
                if triangle and i < j:
                    continue
                corr_data = df1[[df1_colname]].merge(df2[[df2_colname]], left_index=True, right_index=True)
                if drop_paired_zero:
                    corr_data.loc[corr_data.sum(axis=1) == 0, :] = np.nan
                corr_data.dropna(inplace=True)

                coef = np.nan
                pvalue = np.nan
                n = corr_data.shape[0]
                if np.min(corr_data.std(axis=0)) > 0:
                    if method == "Pearson":
                        coef, pvalue = stats.pearsonr(corr_data.iloc[:, 1], corr_data.iloc[:, 0])
                    elif method == "Spearman":
                        coef, pvalue = stats.spearmanr(corr_data.iloc[:, 1], corr_data.iloc[:, 0])

                corr_df.loc[df1_colname, df2_colname] = coef
                annot_df.loc[df1_colname, df2_colname] = "{:.2f}\np={:.2e}\nn={:,}".format(coef, pvalue, n)

        return corr_df, annot_df

    def plot_heatmap(self, df, annot_df, xlabel="", ylabel="",
                     vmin=-1, vmax=1, filename="plot"):
        cmap = sns.diverging_palette(246, 24, as_cmap=True)

        fig, axes = plt.subplots(nrows=2,
                                 ncols=2,
                                 figsize=(1 * df.shape[1] + 10, 1 * df.shape[0] + 10),
                                 gridspec_kw={"width_ratios": [0.2, 0.8],
                                              "height_ratios": [0.8, 0.2]})
        sns.set(color_codes=True)

        annot_df.fillna("", inplace=True)

        row_index = 0
        col_index = 0
        for _ in range(4):
            ax = axes[row_index, col_index]
            if row_index == 0 and col_index == 1:
                sns.heatmap(df, cmap=cmap, vmin=vmin, vmax=vmax, center=0,
                            square=True, annot=annot_df, fmt='',
                            cbar=False, annot_kws={"size": 10, "color": "#000000"},
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

        plt.tight_layout()
        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "{}.{}".format(filename, extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Data path 1: {}".format(self.data_path1))
        print("  > Name 1: {}".format(self.name1))
        print("  > Data path 2: {}".format(self.data_path2))
        print("  > Name 2: {}".format(self.name2))
        print("  > Cell type: {}".format(self.cell_type))
        print("  > Outpath {}".format(self.outdir))
        print("  > Extensions: {}".format(self.extensions))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
