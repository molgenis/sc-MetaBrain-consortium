#!/usr/bin/env python3

"""
File:         plot_confusion_matrix.py
Created:      2023/04/06
Last Changed: 2023/04/07
Author:       M.Vochteloo

Copyright (C) 2022 M.Vochteloo
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
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Plot Convusion Matrix"
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

"""
Syntax: 

./plot_confusion_matrix.py \
    --workdir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-02-02-WorkGroup2CellType/2023-04-05-Mathys2019/all
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.workdir = getattr(arguments, 'workdir')
        self.extensions = getattr(arguments, 'extension')

        self.query_cell_types = ["broad.cell.type", "Subcluster"]
        self.reference_cell_types = ["predicted.cluster", "predicted.major_subclass", "predicted.minor_subclass", "predicted.cross_species_cluster"]

        self.matches = (
            ("broad.cell.type", "predicted.major_subclass"),
            ("Subcluster", "predicted.minor_subclass"),
        )

        self.major_subclass_map = {
            "astrocyte": "Ast",
            "endothelial cell": "End",
            "excitatory neuron": "Ex",
            "inhibitory neuron": "In",
            "perivascular macrophage": "Mic",
            "oligodendrocyte": "Oli",
            "oligodendrocyte precursor cell": "Opc",
            "pericyte": "Per"
        }

        # Creating output directory
        self.outdir = os.path.join(self.workdir, "plots", "confusion_matrices")
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add other arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit")
        parser.add_argument("--workdir",
                            type=str,
                            required=True,
                            help="The working directory.")
        parser.add_argument("--extension",
                            nargs="+",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()
        print("Loading data")
        meta_data_df = self.load_file(os.path.join(self.workdir, "data", "seurat_metadata.csv"))
        meta_data_df["predicted.major_subclass"] = meta_data_df["predicted.major_subclass"].map(self.major_subclass_map)

        print("\tPlotting.")
        for query_cell_type, reference_cell_type in self.matches:
            subset = meta_data_df[[query_cell_type, reference_cell_type]].copy()
            subset.columns = ["real", "assigned"]
            subset["real"] = subset["real"].astype(str)
            subset["assigned"] = subset["assigned"].astype(str)
            confusion_df, annotation_df = self.create_confusion_matrix(df=subset)

            self.plot_heatmap(df=confusion_df,
                              annot_df=annotation_df,
                              xlabel="assigned",
                              ylabel="real",
                              outfile="{}_vs_{}_confusion_matrix".format(query_cell_type.replace(".", "_"), reference_cell_type.replace(".", "_")))

    @staticmethod
    def load_file(inpath, header=0, index_col=None, sep=",", skiprows=None,
                  nrows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         skiprows=skiprows, nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def create_confusion_matrix(df):
        row_counts = list(zip(*np.unique(df["real"], return_counts=True)))
        row_counts.sort(key=lambda x: x[0])
        row_labels = ["{} [n={:,.0f}]".format(label, size) for label, size in row_counts]

        col_counts = list(zip(*np.unique(df["assigned"], return_counts=True)))
        col_counts.sort(key=lambda x: x[0])
        col_labels = ["{} [n={:,.0f}]".format(label, size) for label, size in col_counts]

        confusion_df = pd.DataFrame(np.nan, index=row_labels, columns=col_labels)
        annotation_df = pd.DataFrame("", index=row_labels, columns=col_labels)
        for row_label, (row_value, row_count) in zip(row_labels, row_counts):
            for col_label, (col_value, _) in zip(col_labels, col_counts):
                n_overlap = df.loc[(df["real"] == row_value) & (df["assigned"] == col_value), :].shape[0]
                frac_overlap = np.nan
                if n_overlap > 0:
                    frac_overlap = n_overlap / row_count
                    annotation_df.loc[row_label, col_label] = "{:.2f}\nn={:,.0f}".format(frac_overlap, n_overlap)
                confusion_df.loc[row_label, col_label] = frac_overlap

        return confusion_df, annotation_df

    def plot_heatmap(self, df, annot_df, xlabel="", ylabel="",
                     outfile="confusion_matrix"):
        cmap = sns.diverging_palette(246, 24, as_cmap=True)

        fig, axes = plt.subplots(nrows=2,
                                 ncols=2,
                                 figsize=(1 * df.shape[1] + 10, 1 * df.shape[0] + 10),
                                 gridspec_kw={"width_ratios": [0.2, 0.8],
                                              "height_ratios": [0.8, 0.2]})
        sns.set(color_codes=True)

        # annot_df.fillna("", inplace=True)

        row_index = 0
        col_index = 0
        for _ in range(4):
            ax = axes[row_index, col_index]
            if row_index == 0 and col_index == 1:
                sns.heatmap(df, cmap=cmap, vmin=-1, vmax=1, center=0,
                            square=True, annot=annot_df, fmt='',
                            cbar=False, annot_kws={"size": 14, "color": "#000000"},
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
            fig.savefig(os.path.join(self.outdir, "{}.{}".format(outfile, extension)))
        plt.close()


    def print_arguments(self):
        print("Arguments:")
        print("  > Working directory: {}".format(self.workdir))
        print("  > Extensions: {}".format(", ".join(self.extensions)))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()