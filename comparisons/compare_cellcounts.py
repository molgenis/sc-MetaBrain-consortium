#!/usr/bin/env python3

"""
File:         compare_cellcounts.py
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
import random
import glob
import os

# Third party imports.
import numpy as np
import pandas as pd
import pyreadr
from scipy import stats
from scipy import sparse
import scanpy
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Compare Cellcounts"
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
        self.input_dir = getattr(arguments, 'input_dir')
        self.extensions = getattr(arguments, 'extensions')

        # Set variables.
        self.outdir = os.path.join(str(os.path.dirname(os.path.abspath(__file__))), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.n_cells = 100
        self.n_genes = 100

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
        parser.add_argument("-i",
                            "--input_dir",
                            type=str,
                            required=True,
                            help="The input directory.")
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

        df1 = self.get_count_matrices(type="CellRanger")
        df2 = self.get_count_matrices(type="CellBender")
        df = df1.merge(df2, on="Pool")
        del df1, df2
        print(df)

        for _, row in df.iterrows():
            print("Processing Pool '{}'".format(row["Pool"]))
            cellranger_df = self.load_h5_data(inpath=row["CellRanger"])
            print(cellranger_df)
            cellbender_df = self.load_h5_data(inpath=row["CellBender"])
            print(cellbender_df)

            row_verlap = set(cellranger_df.index).intersection(set(cellbender_df.index))
            nrows = len(row_verlap)
            print("\tN overlapping rows: {:,}".format(nrows))

            col_verlap = set(cellranger_df.columns).intersection(set(cellbender_df.columns))
            ncols = len(col_verlap)
            print("\tN overlapping columns: {:,}".format(ncols))

            cellranger_df = cellranger_df.loc[list(row_verlap), list(col_verlap)]
            cellbender_df = cellbender_df.loc[list(row_verlap), list(col_verlap)]

            print(cellranger_df)
            print(cellbender_df)

            delta_df = self.calculate_delta_dataframe(df1=cellranger_df, df2=cellbender_df)
            print(delta_df)

            col_filtered_cellranger_df = cellranger_df.copy()
            col_filtered_cellbender_df = cellbender_df.copy()
            filtered_delta_df = delta_df.copy()
            if self.n_genes is not None:
                random_columns = random.sample(range(0, ncols), self.n_genes)
                col_filtered_cellranger_df = cellranger_df.iloc[:, random_columns].copy()
                col_filtered_cellbender_df = cellbender_df.iloc[:, random_columns].copy()
                filtered_delta_df = delta_df.iloc[:, :self.n_genes].copy()

            corr_df, annot_df = self.correlate_dataframes(df1=col_filtered_cellranger_df, df2=col_filtered_cellbender_df)
            self.plot_heatmap(df=corr_df,
                              annot_df=annot_df,
                              xlabel="CellRanger",
                              ylabel="CellBender",
                              filename="{}Spearman_corr_heatmap_perGene".format(row["Pool"]))

            row_filtered_cellranger_df = cellranger_df.copy()
            row_filtered_cellbender_df = cellbender_df.copy()
            if self.n_genes is not None:
                random_indices = random.sample(range(0, nrows), self.n_cells)
                row_filtered_cellranger_df = cellranger_df.iloc[random_indices, :].copy()
                row_filtered_cellbender_df = cellbender_df.iloc[random_indices, :].copy()
                filtered_delta_df = filtered_delta_df.iloc[:self.n_cells, :].copy()

            corr_df, annot_df = self.correlate_dataframes(df1=row_filtered_cellranger_df.T, df2=row_filtered_cellbender_df.T)
            self.plot_heatmap(df=corr_df,
                              annot_df=annot_df,
                              xlabel="CellRanger",
                              ylabel="CellBender",
                              filename="{}Spearman_corr_heatmap_perCell".format(row["Pool"]))

            self.plot_heatmap(df=filtered_delta_df,
                              annot_df=filtered_delta_df,
                              xlabel="Genes",
                              ylabel="Cells",
                              filename="{}delta_heatmap_perCell".format(row["Pool"]),
                              vmax = filtered_delta_df.max(axis=0).max(),
                              cbar=True)
            exit()
        exit()


    def get_count_matrices(self, type):
        pool_paths = glob.glob(os.path.join(self.input_dir, type, "*"))

        data = []
        for pool_path in pool_paths:
            pool = os.path.basename(pool_path)

            if type == "CellRanger":
                inpath = os.path.join(self.input_dir, type, pool, "outs", "filtered_feature_bc_matrix.h5")
                if os.path.exists(inpath):
                    data.append([pool, inpath])
            elif type == "CellBender":
                inpath = os.path.join(self.input_dir, type, pool, "cellbender_remove_background_output_filtered.h5")
                if os.path.exists(inpath):
                    data.append([pool, inpath])
            else:
                print("Error, unknown type")
                exit()

        return pd.DataFrame(data, columns=["Pool", type])

    def load_h5_data(self, inpath):
        print(inpath)
        data = scanpy.read_10x_h5(inpath)
        df = pd.DataFrame(sparse.csr_matrix.todense(data.X), index=data.obs_names, columns=data.var_names)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def calculate_delta_dataframe(df1, df2):
        delta_df = df1.subtract(df2)

        col_order = delta_df.sum(axis=0)
        col_order.sort_values(inplace=True, ascending=False)

        row_order = delta_df.sum(axis=1)
        row_order.sort_values(inplace=True, ascending=False)

        return delta_df.loc[row_order.index, col_order.index]


    @staticmethod
    def correlate_dataframes(df1, df2, triangle=False, drop_paired_zero=True, method="Pearson"):
        print(df1)
        print(df2)

        corr_df = pd.DataFrame(np.nan, index=df1.columns, columns=df2.columns)
        annot_df = pd.DataFrame("{:.2f}\np={:.2e}\nn={:,}".format(np.nan, np.nan, 0), index=df1.columns, columns=df2.columns)

        for i, df1_colname in enumerate(df1.columns):
            for j, df2_colname in enumerate(df2.columns):
                if triangle and i < j:
                    continue

                corr_data = df1[[df1_colname]].merge(df2[[df2_colname]], left_index=True, right_index=True)
                corr_data.dropna(inplace=True)
                if corr_df.shape[0] == 0:
                    continue

                if drop_paired_zero:
                    corr_data.loc[corr_data.sum(axis=1) == 0, :] = np.nan
                    corr_data.dropna(inplace=True)

                coef = np.nan
                pvalue = np.nan
                n = corr_data.shape[0]
                if n == 0:
                    coef = 1
                    pvalue = 0
                else:
                    if n >= 2 and np.min(corr_data.std(axis=0)) > 0:
                        if method == "Pearson":
                            coef, pvalue = stats.pearsonr(corr_data.iloc[:, 1], corr_data.iloc[:, 0])
                        elif method == "Spearman":
                            coef, pvalue = stats.spearmanr(corr_data.iloc[:, 1], corr_data.iloc[:, 0])

                corr_df.loc[df1_colname, df2_colname] = coef
                annot_df.loc[df1_colname, df2_colname] = "{:.2f}\np={:.2e}\nn={:,}".format(coef, pvalue, n)

        return corr_df, annot_df

    def plot_heatmap(self, df, annot_df, xlabel="", ylabel="",
                     vmin=-1, vmax=1, cbar=False, filename="plot"):
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
                            cbar=cbar, annot_kws={"size": 10, "color": "#000000"},
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
        print("  > Inpout directory: {}".format(self.input_dir))
        print("  > Outpath {}".format(self.outdir))
        print("  > Extensions: {}".format(self.extensions))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
