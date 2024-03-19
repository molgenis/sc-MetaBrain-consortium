#!/usr/bin/env python3

"""
File:         plot_barcode_overlap.py
Created:      2023/08/18
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
import itertools
import argparse
import glob
import os

# Third party imports.
import pandas as pd
import matplotlib.pyplot as plt
import upsetplot as up

# Local application imports.

# Metadata
__program__ = "Plot Barcode Overlap"
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
./plot_barcode_overlap.py -h
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.cellranger_workdir = getattr(arguments, 'cell_ranger')
        self.cellbender_workdirs = getattr(arguments, 'cell_bender')
        self.cellbender_workdir_labels = getattr(arguments, 'cell_bender_labels')
        self.out_filename = getattr(arguments, 'outfile')
        self.extensions = getattr(arguments, 'extension')

        if len(self.cellbender_workdirs) != len(self.cellbender_workdir_labels):
            print("Error, CellBender work directories and labels need to have the same length.")
            exit()

        # Set variables.
        self.outdir = os.path.join(str(os.path.dirname(os.path.abspath(__file__))), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

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
        parser.add_argument("-cr",
                            "--cell_ranger",
                            type=str,
                            required=True,
                            help="The path to the CellRanger input directory.")
        parser.add_argument("-cb",
                            "--cell_bender",
                            nargs="+",
                            type=str,
                            required=True,
                            help="The path to the CellBender input directories.")
        parser.add_argument("-cbl",
                            "--cell_bender_labels",
                            nargs="+",
                            type=str,
                            required=True,
                            help="The labels for the CellBender input directories.")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            required=True,
                            help="The name of the outfile.")
        parser.add_argument("-e",
                            "--extension",
                            nargs="+",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        data = {}
        print("Load CellRanger data.")
        for subir in glob.glob(os.path.join(self.cellranger_workdir, "*")):
            if os.path.isdir(subir):
                folder = os.path.basename(subir)
                barcodes_path = os.path.join(subir, "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
                if os.path.exists(barcodes_path):
                    barcodes = set(pd.read_csv(barcodes_path, sep="\t", header=None, index_col=None).iloc[:, 0].values)
                    if folder not in data:
                        data[folder] = {}
                    data[folder]["CellRanger"] = barcodes

        print("Load CellBender data.")
        for cellbender_workdir, cellbender_label in zip(self.cellbender_workdirs, self.cellbender_workdir_labels):
            for subir in glob.glob(os.path.join(cellbender_workdir, "*")):
                if os.path.isdir(subir):
                    folder = os.path.basename(subir)
                    barcodes_path = os.path.join(subir, "cellbender_remove_background_output_cell_barcodes.csv")
                    if os.path.exists(barcodes_path):
                        barcodes = set(pd.read_csv(barcodes_path, sep="\t", header=None, index_col=None).iloc[:, 0].values)
                        if folder not in data:
                            data[folder] = {}
                        data[folder][cellbender_label] = barcodes

        print("Finding overlap.")
        counts_data = []
        folders = []
        for folder, folder_data in data.items():
            counts = self.count(folder_data)

            counts_data.append(counts)
            folders.append(folder)
            counts = counts[counts > 0]

            print("Creating plot.")
            up.plot(counts, sort_by='cardinality', show_counts=True)
            plt.savefig(os.path.join(self.outdir, "{}_{}_barcodes_upsetplot.png".format(self.out_filename, folder)))
            plt.close()

        counts_df = pd.concat(counts_data, axis=1)
        counts_df.columns = folders
        print(counts_df)

        print("Creating combined plot.")
        summed_counts = counts_df.sum(axis=1)
        up.plot(summed_counts, sort_by='cardinality', show_counts=True)
        plt.savefig(os.path.join(self.outdir, "{}_barcodes_upsetplot.png".format(self.out_filename)))
        plt.close()
        exit()


    @staticmethod
    def count(input_data):
        combinations = []
        cols = list(input_data.keys())
        for i in range(1, len(cols) + 1):
            combinations.extend(list(itertools.combinations(cols, i)))

        indices = []
        data = []
        for combination in combinations:
            index = []
            for col in cols:
                if col in combination:
                    index.append(True)
                else:
                    index.append(False)

            background = set()
            for key in cols:
                if key not in combination:
                    work_set = input_data[key].copy()
                    background.update(work_set)

            overlap = None
            for key in combination:
                work_set = input_data[key].copy()
                if overlap is None:
                    overlap = work_set
                else:
                    overlap = overlap.intersection(work_set)

            duplicate_set = overlap.intersection(background)
            length = len(overlap) - len(duplicate_set)

            indices.append(index)
            data.append(length)

        s = pd.Series(data,
                      index=pd.MultiIndex.from_tuples(indices, names=cols))
        s.name = "value"
        return s


    def print_arguments(self):
        print("Arguments:")
        print("  > CellRanger work directory: {}".format(self.cellranger_workdir))
        for cellbender_workdir, cellbender_label in zip(self.cellbender_workdirs, self.cellbender_workdir_labels):
            print("  > CellBender work directory '{}': {}".format(cellbender_label, cellbender_workdir))
        print("  > Output filename: {}".format(self.out_filename))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Extensions: {}".format(self.extensions))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
