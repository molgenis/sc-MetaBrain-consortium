#!/usr/bin/env python3

"""
File:         compare_confusion_matrices.py
Created:      2023/04/06
Last Changed: 2023/04/10
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
import itertools
import argparse
import glob
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
__program__ = "Plot Confusion Matrices"
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

./compare_confusion_matrices.py \
    --workdir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-02-02-WorkGroup2CellType/2023-04-05-Mathys2019 \
    --truths broad.cell.type Subcluster \
    --predictions predicted.major_subclass predicted.minor_subclass \
    --extension pdf
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.workdir = getattr(arguments, 'workdir')
        self.truths = getattr(arguments, 'truths')
        self.predictions = getattr(arguments, 'predictions')
        self.extensions = getattr(arguments, 'extension')

        self.ct_abbr_map = {
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

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

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
        parser.add_argument("--truths",
                            nargs="+",
                            type=str,
                            required=True,
                            help="The meta data column considered the truth.")
        parser.add_argument("--predictions",
                            nargs="+",
                            type=str,
                            required=True,
                            help="The meta data column containing the "
                                 "azimuth predictions.")
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
        data = []
        folders = []
        for folder in glob.glob(os.path.join(self.workdir, "*")):
            folder = os.path.basename(folder)
            metadata_path = os.path.join(self.workdir, folder, "data", "seurat_metadata.csv")
            if not os.path.exists(metadata_path):
                continue

            folder_numbers = np.inf
            if folder.startswith("downsample"):
                try:
                    folder_numbers = int(folder.replace("downsample", ""))
                except ValueError:
                    pass

            print("\tloading '{}' results".format(folder))
            meta_data_df = self.load_file(metadata_path)
            meta_data_df["folder"] = folder
            for column in self.truths + self.predictions:
                meta_data_df[column] = [self.ct_abbr_map[value] if value in self.ct_abbr_map else str(value) for value in meta_data_df[column]]
            meta_data_df = meta_data_df[["folder"] + self.truths + self.predictions]
            data.append(meta_data_df)
            folders.append((folder, folder_numbers))

        df = pd.concat(data, axis=0)
        folders.sort(key=lambda x: -x[1])
        folders = [folder[0] for folder in folders]

        print("Plotting.")
        truth_combinations = list(itertools.combinations(self.truths, 2))
        for x, y in truth_combinations:
            print("\tcomparing truth '{}' vs '{}'".format(x, y))
            self.create_multi_confusion_matrix(df=df,
                                               x=x,
                                               y=y,
                                               col="folder",
                                               col_order=folders)

        prediction_combinations = list(itertools.combinations(self.predictions, 2))
        for x, y in prediction_combinations:
            print("\tcomparing predictions '{}' vs '{}'".format(x, y))
            self.create_multi_confusion_matrix(df=df,
                                               x=x,
                                               y=y,
                                               col="folder",
                                               col_order=folders)

        for x in self.predictions:
            for y in self.truths:
                print("\tcomparing predictions '{}' vs '{}'".format(x, y))
                self.create_multi_confusion_matrix(df=df,
                                                   x=x,
                                                   y=y,
                                                   col="folder",
                                                   col_order=folders)

    @staticmethod
    def load_file(inpath, header=0, index_col=None, sep=",", skiprows=None,
                  nrows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         skiprows=skiprows, nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    def create_multi_confusion_matrix(self, df, x, y, col, col_order):
        confusion_data = {}
        annotation_data = {}
        tpr_data = {}
        for column in df[col].unique():
            confusion_df, annotation_df, tpr = self.create_confusion_matrix(
                df=df.loc[df[col] == column, [x, y]].copy(),
                x=x,
                y=y)
            confusion_data[column] = confusion_df
            annotation_data[column] = annotation_df
            tpr_data[column] = tpr

        self.plot_multi_heatmap(
            confusion_data=confusion_data,
            annotation_data=annotation_data,
            tpr_data=tpr_data,
            col_order=col_order,
            outfile="{}_vs_{}_confusion_matrix".format(x.replace(".", "_"), y.replace(".", "_")))

    @staticmethod
    def create_confusion_matrix(df, x, y):
        row_counts = list(zip(*np.unique(df[y], return_counts=True)))
        row_counts.sort(key=lambda x: x[0])
        row_labels = ["{} [n={:,.0f}]".format(label, size) for label, size in row_counts]

        col_counts = list(zip(*np.unique(df[x], return_counts=True)))
        col_counts.sort(key=lambda x: x[0])
        col_labels = ["{} [n={:,.0f}]".format(label, size) for label, size in col_counts]

        confusion_df = pd.DataFrame(np.nan, index=row_labels, columns=col_labels)
        annotation_df = pd.DataFrame("", index=row_labels, columns=col_labels)
        tp_count = 0
        tp_total = 0
        for row_label, (row_value, row_count) in zip(row_labels, row_counts):
            for col_label, (col_value, _) in zip(col_labels, col_counts):
                n_overlap = df.loc[(df[y] == row_value) & (df[x] == col_value), :].shape[0]
                frac_overlap = np.nan
                if n_overlap > 0:
                    frac_overlap = n_overlap / row_count
                    annotation_df.loc[row_label, col_label] = "{:.2f}\nn={:,.0f}".format(frac_overlap, n_overlap)

                if row_value == col_value:
                    tp_count += n_overlap
                    tp_total += row_count
                confusion_df.loc[row_label, col_label] = frac_overlap

        tpr = np.nan
        if tp_total > 0 :
            tpr = tp_count / tp_total

        return confusion_df, annotation_df, tpr

    def plot_multi_heatmap(self, confusion_data, annotation_data, tpr_data,
                           col_order, outfile="confusion_matrix"):
        cmap = sns.diverging_palette(246, 24, as_cmap=True)

        nplots = len(col_order)
        figsize_width = sum([confusion_df.shape[1] for confusion_df in confusion_data.values()])
        figsize_length = max([confusion_df.shape[0] for confusion_df in confusion_data.values()])
        width_ratio_text = 0.2 / nplots
        width_ratio_plot = 0.8 / nplots
        fig, axes = plt.subplots(nrows=2,
                                 ncols=nplots * 2,
                                 sharex='none',
                                 sharey='none',
                                 figsize=(figsize_width + 5 * nplots, figsize_length + 5 * nplots),
                                 gridspec_kw={"width_ratios": [width_ratio_text, width_ratio_plot] * nplots,
                                              "height_ratios": [0.8, 0.2]})
        sns.set(color_codes=True)

        for col_index, col_value in enumerate(col_order):
            axes[1, col_index * 2].set_axis_off()
            axes[1, col_index * 2 + 1].set_axis_off()
            axes[0, col_index * 2].set_axis_off()

            ax = axes[0, col_index * 2 + 1]
            data = confusion_data[col_value]
            sns.heatmap(data,
                        cmap=cmap,
                        vmin=-1,
                        vmax=1,
                        center=0,
                        square=True,
                        yticklabels=True,
                        annot=annotation_data[col_value],
                        fmt='',
                        cbar=False,
                        annot_kws={"size": 10},
                        ax=ax)

            # tmp_ylabels = [""] * data.shape[0]
            # if col_index == 1:
            #     tmp_ylabels = ax.get_ymajorticklabels()

            # if ylabels is None:
            #     ylabels = data.index.tolist()
            # elif data.index.tolist() != ylabels:
            #     print("Error, not identical row labels.")
            #     plt.close()
            #     exit()

            plt.setp(ax.set_yticklabels(ax.get_ymajorticklabels(),
                                        fontsize=25,
                                        rotation=0))
            plt.setp(ax.set_xticklabels(ax.get_xmajorticklabels(),
                                        fontsize=25,
                                        rotation=90))

            tmp_title = col_value
            if not np.isnan(tpr_data[col_value]):
                tmp_title += " [TPR: {:.2f}]".format(tpr_data[col_value])
            ax.set_title(tmp_title, fontsize=40)

        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "{}.{}".format(outfile, extension)))
            print("Saved '{}.{}'.".format(outfile, extension))
        plt.close()


    def print_arguments(self):
        print("Arguments:")
        print("  > Working directory: {}".format(self.workdir))
        print("  > Truth column: {}".format(self.truths))
        print("  > Predictions columns: {}".format(self.predictions))
        print("  > Extensions: {}".format(", ".join(self.extensions)))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()