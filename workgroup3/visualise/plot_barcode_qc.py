#!/usr/bin/env python3

"""
File:         plot_barcode_qc.py
Created:      2024/04/10
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
__program__ = "Plot Barcode QC"
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
./plot_barcode_qc.py -h
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.metadata = getattr(arguments, 'metadata')
        self.qc_mad = getattr(arguments, 'qc_mad')
        self.palette = getattr(arguments, 'palette')
        outdir = getattr(arguments, 'outdir')
        self.extensions = getattr(arguments, 'extension')

        # Creating output directory
        self.base_outdir = os.path.join(outdir, "plots")
        self.outdir = os.path.join(outdir, "plots", "barcode_qc")
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
        parser.add_argument("--metadata",
                            type=str,
                            required=True,
                            help="")
        parser.add_argument("--qc_mad",
                            type=str,
                            required=True,
                            help="")
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

        print("Loading data")
        metadata_df = self.load_file(self.metadata)
        print(metadata_df)
        qc_mad_df = self.load_qc_mad(self.qc_mad)
        print(qc_mad_df)

        # metadata_df = metadata_df.loc[metadata_df["nCount_RNA"] > 500, :]

        palette = None
        if self.palette is not None:
            print("Loading palette")
            with open(self.palette) as f:
                palette = json.load(f)
            f.close()

        for qc_metric in ["nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"]:
            if qc_metric + ".tag" not in metadata_df:
                metadata_df[qc_metric + ".tag"] = "NotOutlier"

            df = metadata_df[[qc_metric, qc_metric + ".tag"]].copy()
            df.columns = ["metric", "tag"]

            lower_lim = None
            higher_lim = None
            if sum(df["tag"] == "NotOutlier") != df.shape[0]:
                lower_lim = df.loc[df["tag"] == "NotOutlier", "metric"].min().round(2)
                higher_lim = df.loc[df["tag"] == "NotOutlier", "metric"].max().round(2)
            print(qc_metric, df["metric"].min(), df["metric"].max(), df["metric"].median(), lower_lim, higher_lim)

            self.histplot(
                df=df,
                x="metric",
                min_vline=lower_lim,
                max_vline=higher_lim,
                xlabel=qc_metric,
                ylabel="Count",
                title=qc_metric,
                filename=qc_metric + "_histplot"
            )

            del df, lower_lim, higher_lim


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
    def load_qc_mad(inpath):
        qc_thresholds = {}
        with open(inpath, 'r') as f:
            for i, line in enumerate(f):
                if i == 0:
                    continue
                qc_metric, _, bound, _, mad_min, _, mad_max = line.rstrip("\n").split("\t")
                qc_thresholds[qc_metric] = {bound: (mad_min, mad_max)}
        return qc_thresholds


    def histplot(self, df, x="x", min_vline=None, max_vline=None, xlabel="", ylabel="", title="", filename="plot"):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        range = abs(df[x].max() - df[x].min())

        g = sns.histplot(data=df,
                         x=x,
                         kde=True,
                         binwidth=range / 100,
                         color="#000000",
                         ax=ax)

        if min_vline:
            ax.axvline(0, ls='--', color="#000000", alpha=0.5, zorder=-1,
                       linewidth=2)
        if max_vline:
            ax.axvline(0, ls='--', color="#000000", alpha=0.5, zorder=-1,
                       linewidth=2)

        ax.set_title(title,
                     fontsize=14,
                     fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=10,
                      fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=10,
                      fontweight='bold')

        ax.annotate(
            'N - NotOutlier = {:,}'.format(df.loc[df["tag"] == "NotOutlier"].shape[0]),
            xy=(0.7, 0.9),
            xycoords=ax.transAxes,
            color="#000000",
            fontsize=14,
            fontweight='bold'
        )
        ax.annotate(
            'N - Outlier = {:,}'.format(df.loc[df["tag"] == "Outlier"].shape[0]),
            xy=(0.7, 0.85),
            xycoords=ax.transAxes,
            color="#000000",
            fontsize=14,
            fontweight='bold'
        )
        ax.annotate(
            'Min: {:.2f}'.format(df["metric"].min()),
            xy=(0.7, 0.80),
            xycoords=ax.transAxes,
            color="#000000",
            fontsize=14,
            fontweight='bold'
        )
        ax.annotate(
            'Max: {:.2f}'.format(df["metric"].max()),
            xy=(0.7, 0.75),
            xycoords=ax.transAxes,
            color="#000000",
            fontsize=14,
            fontweight='bold'
        )
        ax.annotate(
            'Median: {:.2f}'.format(df["metric"].median()),
            xy=(0.7, 0.70),
            xycoords=ax.transAxes,
            color="#000000",
            fontsize=14,
            fontweight='bold'
        )
        ax.annotate(
            'Incl. lower limit = {}'.format(min_vline),
            xy=(0.7, 0.65),
            xycoords=ax.transAxes,
            color="#000000",
            fontsize=14,
            fontweight='bold'
        )
        ax.annotate(
            'Incl. upper limit = {}'.format(max_vline),
            xy=(0.7, 0.60),
            xycoords=ax.transAxes,
            color="#000000",
            fontsize=14,
            fontweight='bold'
        )

        plt.tight_layout()
        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "{}.{}".format(filename, extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Metadata file: {}".format(self.metadata))
        print("  > QC MAD file: {}".format(self.qc_mad))
        print("  > Palette: {}".format(self.palette))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Extensions: {}".format(", ".join(self.extensions)))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()