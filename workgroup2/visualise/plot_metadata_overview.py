#!/usr/bin/env python3

"""
File:         plot_metadata_overview.py
Created:      2024/01/17
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
import argparse
import glob
import json
import os
import re

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

# Local application imports.

# Metadata
__program__ = "Plot Metadata Overview"
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
./plot_metadata_overview.py -h
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.workdir = getattr(arguments, 'workdir')
        self.mapping = getattr(arguments, 'mapping')
        self.mapping_key = getattr(arguments, 'key')
        self.ct_columns = getattr(arguments, 'columns')
        self.palette = getattr(arguments, 'palette')
        self.extensions = getattr(arguments, 'extension')

        # Creating output directory
        self.outdir = os.path.join(self.workdir, "plots", "metadata_overview")
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
        parser.add_argument("--mapping",
                            type=str,
                            required=False,
                            default=None,
                            help="A cell mapping file.")
        parser.add_argument("--key",
                            type=str,
                            required=False,
                            default=None,
                            help="A cell mapping key.")
        parser.add_argument("--columns",
                            nargs="*",
                            type=str,
                            required=False,
                            default=None,
                            help="")
        parser.add_argument("--palette",
                            type=str,
                            required=False,
                            default=None,
                            help="A color palette file.")
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
        meta_data_list = []
        for path in glob.glob(os.path.join(self.workdir, "*")):
            dataset = os.path.basename(path)
            metadata_path1 = os.path.join(self.workdir, dataset, "map", "azimuth.metadata.tsv.gz")
            metadata_path2 = os.path.join(self.workdir, dataset, "map", "azimuth_all.metadata.tsv.gz")
            if not os.path.exists(metadata_path1) and not os.path.exists(metadata_path2):
                continue

            dataset_label = re.match("[0-9]{4}-[0-9]{2}-[0-9]{2}-(.+)", dataset).group(1)
            if "-sequencing-round-1" in dataset_label:
                dataset_label = dataset_label.replace("-sequencing-round-1", "-SR1")
            if "-sequencing-round-2" in dataset_label:
                dataset_label = dataset_label.replace("-sequencing-round-2", "-SR2")

            metadata_path = None
            if os.path.exists(metadata_path1):
                metadata_path = metadata_path1
            else:
                metadata_path = metadata_path2

            print("  loading '{}' results".format(dataset))
            meta_data_df = self.load_file(metadata_path)
            meta_data_df["dataset"] = dataset_label
            meta_data_list.append(meta_data_df)

        meta_data_df = pd.concat(meta_data_list, axis=0)
        ct_columns = [col for col in meta_data_df.columns if col.startswith("predicted.") and not col.endswith(".score")]

        if self.mapping is not None:
            print("Adding cell mapping")
            mapping_df = self.load_file(self.mapping, sep=";", header=0)
            if self.mapping_key is None:
                self.mapping_key = mapping_df.columns[0]
            ct_columns.append(*[col for col in mapping_df if col != self.mapping_key])
            meta_data_df = meta_data_df.merge(mapping_df, on=self.mapping_key, how="left")

        palette = None
        if self.palette is not None:
            print("Loading palette")
            with open(self.palette) as f:
                palette = json.load(f)
            f.close()

        if self.ct_columns is not None:
            ct_columns = set(ct_columns).intersection(set(self.ct_columns))

        print("Plotting.")
        for ct_col in ct_columns:
            print("\t" + ct_col)
            counts = meta_data_df[["Barcode", ct_col]].groupby(ct_col).count()
            counts.reset_index(inplace=True, drop=False)
            print(counts)

            self.plot_barplot(
                df=counts,
                x=ct_col,
                y="Barcode",
                xlabel="Cell Type",
                ylabel="Counts",
                title=ct_col,
                palette=palette,
                filename=ct_col.replace("predicted.", "") + "_combined_cell_counts_barplot"
            )

            counts = meta_data_df[["Barcode", ct_col, "dataset"]].groupby([ct_col, "dataset"]).count()
            counts.reset_index(inplace=True, drop=False)
            print(counts)

            self.plot_violinplot(
                df=counts,
                x=ct_col,
                y="Barcode",
                hue=ct_col,
                xlabel="Cell Type",
                ylabel="Counts",
                title=ct_col,
                palette=palette,
                filename=ct_col.replace("predicted.", "") + "_combined_cell_counts_violin"
            )

            counts = meta_data_df[["Barcode", ct_col, "dataset", "Pool"]].groupby([ct_col, "dataset", "Pool"]).count()
            counts.reset_index(inplace=True, drop=False)
            print(counts)

            self.plot_violinplot_per_group(
                df=counts,
                x="dataset",
                y="Barcode",
                group=ct_col,
                xlabel="Cell Type",
                ylabel="Counts",
                title=ct_col,
                palette=palette,
                filename=ct_col.replace("predicted.", "") + "_combined_cell_counts_violin_per_ct"
            )

            total_counts = meta_data_df[["Barcode", "dataset", "Pool"]].groupby(["dataset", "Pool"]).count()
            total_counts.columns = ["Total"]
            total_counts.reset_index(inplace=True, drop=False)
            fraction_counts = counts.merge(total_counts, on=["dataset", "Pool"], how="left")
            fraction_counts["fraction"] = (fraction_counts["Barcode"] / fraction_counts["Total"]) * 100
            print(fraction_counts)

            self.plot_violinplot_per_group(
                df=fraction_counts,
                x="dataset",
                y="fraction",
                group=ct_col,
                xlabel="Cell Type",
                ylabel="%cells",
                title=ct_col,
                palette=palette,
                filename=ct_col.replace("predicted.", "") + "_combined_cell_fractions_violin_per_ct"
            )

    @staticmethod
    def load_file(inpath, header=0, index_col=None, sep="\t", skiprows=None,
                  nrows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         skiprows=skiprows, nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    def plot_barplot(self, df, x="x", y="y", xlabel="", ylabel="", title="",
                     palette=None, filename=""):
        sns.set_style("ticks")
        fig, ax = plt.subplots(figsize=(12, 12))

        sns.despine(fig=fig, ax=ax)

        color = None
        if palette is None:
            color = "#808080"
        else:
            for key in df[x].unique():
                if key not in palette:
                    color = "#808080"
                    palette = None
                    break

        df.sort_values(by=y, inplace=True, ascending=False)

        g = sns.barplot(x=x,
                        y=y,
                        data=df,
                        color=color,
                        palette=palette,
                        dodge=False,
                        order=df[x],
                        ax=ax)

        rotation = 0
        if df.shape[0] > 10:
            rotation = 30
        if df.shape[0] > 20:
            rotation = 90

        g.set_xticklabels(ax.get_xticklabels(), rotation=rotation, fontsize=14)
        for index, row in df.iterrows():
            g.text(row[x],
                   row[y],
                   "{:,}\n{:.2f}%\n ".format(row[y], (100 / df[y].sum()) * row[y]),
                   fontsize=14,
                   color='black',
                   ha="center")

        ax.set_title(title,
                     fontsize=22,
                     fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        plt.tight_layout()
        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "{}.{}".format(filename, extension)))
        plt.close()

    def plot_violinplot(self, df, x="x", y="y", hue=None, xlabel="",
             ylabel="", title="", palette=None, filename=""):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        self.violinplot(
            ax=ax,
            df=df,
            x=x,
            y=y,
            hue=hue,
            xlabel=xlabel,
            ylabel=ylabel,
            title=title,
            palette=palette)

        plt.tight_layout()
        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "{}.{}".format(filename, extension)))
        plt.close()

    def plot_violinplot_per_group(self, df, group, x="x", y="y", hue=None, xlabel="",
                        ylabel="", title="", palette=None, filename=""):
        groups = df[group].unique()

        for group_value in groups:
            if group_value not in palette:
                palette = None
                break

        ngroups = len(groups)
        ncols = int(np.ceil(np.sqrt(ngroups)))
        nrows = int(np.ceil(ngroups / ncols))

        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharex="none",
                                 sharey="none",
                                 figsize=(4 * ncols, 4 * nrows))
        sns.set(color_codes=True)

        row_index = 0
        col_index = 0
        for i in range(ncols * nrows):
            if nrows == 1 and ncols == 1:
                ax = axes
            elif nrows == 1 and ncols > 1:
                ax = axes[col_index]
            elif nrows > 1 and ncols == 1:
                ax = axes[row_index]
            else:
                ax = axes[row_index, col_index]

            if i < ngroups:
                color = "#808080"
                if palette is not None:
                    color = palette[groups[i]]

                self.violinplot(
                    ax=ax,
                    df=df.loc[df[group] == groups[i]],
                    x=x,
                    y=y,
                    hue=hue,
                    xlabel=xlabel if row_index == nrows - 1 else "",
                    xticks=row_index == nrows - 1,
                    ylabel=ylabel if col_index == 0 else "",
                    title=groups[i],
                    color=color,
                    palette=None)
            else:
                ax.set_axis_off()

            col_index += 1
            if col_index > (ncols - 1):
                col_index = 0
                row_index += 1

        fig.suptitle(title,
                     fontsize=28,
                     fontweight='bold')

        plt.tight_layout()
        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "{}.{}".format(filename, extension)))
        plt.close()

    @staticmethod
    def violinplot(ax, df, x="x", y="y", hue=None, xlabel="", xticks=True,
             ylabel="", title="", color=None, palette=None):

        if color is None:
            if palette is None:
                color = "#808080"
            else:
                for key in df[x].unique():
                    if key not in palette:
                        color = "#808080"
                        palette = None
                        break

        sns.violinplot(x=x,
                       y=y,
                       hue=hue,
                       data=df,
                       color=color,
                       palette=palette,
                       cut=0,
                       dodge=False,
                       ax=ax)

        plt.setp(ax.collections, alpha=.75)

        sns.boxplot(x=x,
                    y=y,
                    hue=hue,
                    data=df,
                    whis=np.inf,
                    color=color,
                    palette=palette,
                    dodge=False,
                    ax=ax)

        if ax.get_legend() is not None:
            ax.get_legend().remove()

        rotation = 0
        n_categories = len(df[x].unique())
        total_length = sum([len(x) for x in df[x].unique()])
        if n_categories > 10 or total_length > 50:
            rotation = 30
        if n_categories > 20 or total_length > 100:
            rotation = 90

        xticklabels = ax.get_xticklabels()
        if not xticks:
            xticklabels = ["" for _ in xticklabels]

        ax.set_xticklabels(xticklabels, rotation=rotation, horizontalalignment="right")

        ax.set_title(title,
                     fontsize=22,
                     fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')

    def print_arguments(self):
        print("Arguments:")
        print("  > Working directory: {}".format(self.workdir))
        print("  > Mapping file: {}".format(self.mapping))
        print("  > Mapping key: {}".format(self.mapping_key))
        print("  > Palette: {}".format(self.palette))
        print("  > Extensions: {}".format(", ".join(self.extensions)))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()