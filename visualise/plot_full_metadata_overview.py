#!/usr/bin/env python3

"""
File:         plot_full_metadata_overview.py
Created:      2024/01/24
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

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Plot Full Metadata Overview"
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
./plot_full_metadata_overview.py -h
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.wg1_workdir = getattr(arguments, 'wg1_workdir')
        self.wg1_datasets = getattr(arguments, 'wg1_datasets')
        self.wg2_workdir = getattr(arguments, 'wg2_workdir')
        self.wg2_datasets = getattr(arguments, 'wg2_datasets')
        self.cell_mapping = getattr(arguments, 'cell_mapping')
        self.cell_level = getattr(arguments, 'cell_level')
        self.palette = getattr(arguments, 'palette')
        outdir = getattr(arguments, 'outdir')
        self.extensions = getattr(arguments, 'extension')

        # Creating output directory
        self.base_outdir = os.path.join(outdir, "plots")
        self.outdir = os.path.join(outdir, "plots", "full_metadata_overview")
        for outdir in [self.base_outdir, self.outdir]:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

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
        parser.add_argument("--wg1_workdir",
                            type=str,
                            required=True,
                            help="The WG1 working directory.")
        parser.add_argument("--wg1_datasets",
                            nargs="*",
                            type=str,
                            required=True,
                            default=None,
                            help="The WG1 datasets to load.")
        parser.add_argument("--wg2_workdir",
                            type=str,
                            required=True,
                            help="The WG2 working directory.")
        parser.add_argument("--wg2_datasets",
                            nargs="*",
                            type=str,
                            required=True,
                            default=None,
                            help="The WG2 datasets to load.")
        parser.add_argument("--cell_mapping",
                            type=str,
                            required=False,
                            default=None,
                            help="A cell mapping file.")
        parser.add_argument("--cell_level",
                            type=str,
                            required=False,
                            default="L1",
                            help="The cell level to plot.")
        parser.add_argument("--palette",
                            type=str,
                            required=False,
                            default=None,
                            help="A color palette file.")
        parser.add_argument("--outdir",
                            type=str,
                            required=True,
                            help="The output directory.")
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

        print("Loading metadata")
        metadata = self.load_metadata()
        print(metadata)

        metadata["CellType"] = metadata["DropletType"]
        metadata.loc[metadata["DropletType"] == "singlet", "CellType"] = metadata[self.cell_level]

        palette = None
        if self.palette is not None:
            print("Loading palette")
            with open(self.palette) as f:
                palette = json.load(f)
            f.close()

        print("Plotting.")
        counts = metadata[["Barcode", "Dataset"]].groupby("Dataset").count()
        counts.reset_index(inplace=True, drop=False)
        counts.sort_values(by="Barcode", inplace=True, ascending=False)
        print(counts)
        dataset_order = counts["Dataset"].tolist()

        self.plot_barplot(
            df=counts,
            x="Dataset",
            y="Barcode",
            xlabel="Dataset",
            ylabel="Counts",
            title="#cells per dataset\nN={:,}".format(counts["Barcode"].sum()),
            palette=palette,
            filename="Dataset_barplot"
        )

        counts = metadata[["Barcode", "CellType"]].groupby("CellType").count()
        counts.reset_index(inplace=True, drop=False)
        counts.sort_values(by="Barcode", inplace=True, ascending=False)
        print(counts)
        ct_order = counts["CellType"].tolist()

        self.plot_barplot(
            df=counts,
            x="CellType",
            y="Barcode",
            xlabel="Cell Type",
            ylabel="Counts",
            title="#cells per cell type\nN={:,}".format(counts["Barcode"].sum()),
            palette=palette,
            filename=self.cell_level + "_CellType_barplot"
        )

        counts = metadata[["Barcode", "CellType", "Dataset"]].groupby(["CellType", "Dataset"]).count()
        counts.reset_index(inplace=True, drop=False)
        print(counts)

        self.plot_violinplot(
            df=counts,
            x="CellType",
            y="Barcode",
            order=ct_order,
            hue="CellType",
            xlabel="Cell Type",
            ylabel="Counts",
            title="#cells per cell type\nN={:,}".format(counts["Barcode"].sum()),
            palette=palette,
            filename=self.cell_level + "_CellType_violinplot"
        )

        counts = metadata[["Barcode", "CellType", "Dataset", "Pool"]].groupby(["CellType", "Dataset", "Pool"]).count()
        counts.reset_index(inplace=True, drop=False)
        print(counts)

        self.plot_violinplot_per_group(
            df=counts,
            x="Dataset",
            y="Barcode",
            plot_order=ct_order,
            hue_order=dataset_order,
            group="CellType",
            xlabel="Dataset",
            ylabel="Counts",
            title="CellType",
            palette=palette,
            filename=self.cell_level + "_CellType_violinplot_per_dataset"
        )


    def load_metadata(self):
        metadata_path = os.path.join(self.base_outdir, "metadata.txt.gz")
        if os.path.exists(metadata_path):
            return self.load_file(metadata_path)

        wg1_metadata = self.load_wg1_metadata(
            workdir=self.wg1_workdir,
            datasets=self.wg1_datasets
        )
        # wg1_metadata = self.load_wg_metadata(
        #     workdir=self.wg1_workdir,
        #     datasets=self.wg1_datasets,
        #     relative_filepaths=[("Step2-DemultiplexingAndDoubletRemoval", "CombinedResults", "Final_Assignments_demultiplexing_doublets.tsv.gz")]
        # )
        wg2_metadata = self.load_wg_metadata(
            workdir=self.wg2_workdir,
            datasets=self.wg2_datasets,
            relative_filepaths=[("map", "azimuth.metadata.tsv.gz"),
                                 ("map", "azimuth_all.metadata.tsv.gz")]
        )
        wg2_metadata = wg2_metadata.loc[:, [col for col in wg2_metadata.columns if not col.endswith(".score")]]

        wg1_metadata["Pool"] = wg1_metadata["Pool"].astype(str)
        wg2_metadata["Pool"] = wg2_metadata["Pool"].astype(str)

        if self.cell_mapping is not None:
            cell_mapping_df = self.load_file(self.cell_mapping, sep=";")
            wg2_metadata = wg2_metadata.merge(cell_mapping_df, how="left")

        metadata = wg1_metadata.merge(wg2_metadata, on=["Pool", "Barcode", "Dataset"], how="inner")
        self.save_file(df=metadata, outpath=metadata_path)
        return metadata

    def load_wg1_metadata(self, workdir, datasets):
        df_list = []
        for dataset in datasets:
            print("  loading '{}' results".format(dataset))
            for path in glob.glob(os.path.join(workdir, dataset, "Step2-DemultiplexingAndDoubletRemoval", "*")):
                pool = os.path.basename(path)
                if pool in ["CombinedResults", "manual_selection", "QC_figures", "log", "slurm_log", "genotypes"] or pool.endswith(".sh") or pool.endswith(".yaml") or pool.endswith(".yaml~"):
                    continue

                doubletfinder_path = os.path.join(workdir, dataset, "Step2-DemultiplexingAndDoubletRemoval", pool, "DoubletFinderRun1", "DoubletFinder_doublets_singlets.tsv.gz")
                if not os.path.exists(doubletfinder_path):
                    print("Error, could not find DoubletFinder input file: {}.".format(doubletfinder_path))
                    continue
                doubletfinder_df = self.load_file(doubletfinder_path)

                scdblfinder_path = os.path.join(workdir, dataset, "Step2-DemultiplexingAndDoubletRemoval", pool, "scDblFinderRun1", "scDblFinder_doublets_singlets.tsv.gz")
                if not os.path.exists(scdblfinder_path):
                    print("Error, could not find scDblFinder input file: {}.".format(scdblfinder_path))
                    continue
                scdblfinder_df = self.load_file(scdblfinder_path)


                df = doubletfinder_df.merge(scdblfinder_df, on="Barcode")
                df["Pool"] = pool
                df["Barcode"] = [barcode.split("-")[0] + "_" + pool for barcode in df["Barcode"]]
                df["Assignment"] = np.nan
                df["DropletType"] = "doublet"
                df.loc[(df["DoubletFinder_DropletType"] == "singlet") & (df["scDblFinder_DropletType"] == "singlet"), "DropletType"] = "singlet"
                df["Dataset"] = dataset.split("-")[-1]
                df_list.append(df[["Pool", "Barcode", "Assignment", "DropletType", "Dataset"]])

        return pd.concat(df_list, axis=0)

    def load_wg_metadata(self, workdir, datasets, relative_filepaths):
        df_list = []
        for dataset in datasets:
            print("  loading '{}' results".format(dataset))
            metadata_path = None
            for relative_filepath in relative_filepaths:
                tmp_metadata_path = os.path.join(workdir, dataset, *relative_filepath)
                if not os.path.exists(tmp_metadata_path):
                    continue
                if metadata_path is not None:
                    print("Error, two valid input files found.")
                    exit()
                metadata_path = tmp_metadata_path
            if metadata_path is None:
                print("Error, could not find input file.")
                continue

            df = self.load_file(metadata_path)
            df["Dataset"] = dataset.split("-")[-1]
            df_list.append(df)

        return pd.concat(df_list, axis=0)

    @staticmethod
    def load_file(inpath, header=0, index_col=None, sep="\t", skiprows=None,
                  nrows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         skiprows=skiprows, nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def save_file(df, outpath, header=True, index=False, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    def plot_barplot(self, df, x="x", y="y", order=None, xlabel="", ylabel="",
                     title="", palette=None, filename=""):
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

        g = sns.barplot(x=x,
                        y=y,
                        data=df,
                        color=color,
                        palette=palette,
                        dodge=False,
                        order=order,
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


    def plot_violinplot(self, df, x="x", y="y", order=None, hue=None,
                        xlabel="", ylabel="", title="", palette=None, filename=""):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        self.violinplot(
            ax=ax,
            df=df,
            x=x,
            y=y,
            order=order,
            hue=hue,
            xlabel=xlabel,
            ylabel=ylabel,
            title=title,
            palette=palette)

        plt.tight_layout()
        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "{}.{}".format(filename, extension)))
        plt.close()

    def plot_violinplot_per_group(self, df, group, x="x", y="y", plot_order=None, hue_order=None, hue=None,
                                  xlabel="", ylabel="", title="", palette=None, filename=""):
        groups = df[group].unique()
        if plot_order is not None and len(set(plot_order).symmetric_difference(set(groups))) == 0:
            groups = plot_order

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
                    order=hue_order,
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
    def violinplot(ax, df, x="x", y="y", order=None, hue=None, xlabel="", xticks=True,
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
                       order=order,
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
                    order=order,
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
        print("  > WG1 working directory: {}".format(self.wg1_workdir))
        print("  > WG1 datasets: {}".format(", ".join(self.wg1_datasets)))
        print("  > WG2 working directory: {}".format(self.wg2_workdir))
        print("  > WG2 datasets: {}".format(", ".join(self.wg2_datasets)))
        print("  > Cell mapping: {}".format(self.cell_mapping))
        print("  > Cell level: {}".format(self.cell_level))
        print("  > Palette: {}".format(self.palette))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Extensions: {}".format(", ".join(self.extensions)))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()