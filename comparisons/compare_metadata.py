#!/usr/bin/env python3

"""
File:         compare_metadata.py
Created:      2024/03/02
Last Changed:
Author:       M.Vochteloo

Copyright (C) 2020 University Medical Center Groningen.

A copy of the BSD 3-Clause "New" or "Revised" License can be found in the
LICENSE file in the root directory of this source tree.
"""

# Standard imports.
from __future__ import print_function
import argparse
import glob
import os

# Third party imports.
import pandas as pd
from scipy import stats


# Local application imports.

# Metadata
__program__ = "Compare Metadata"
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
./compare_metadata.py -h
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.row_data_path = getattr(arguments, 'row_data')
        self.row_name = " ".join(getattr(arguments, 'row_name'))
        self.col_data_path = getattr(arguments, 'col_data')
        self.col_name = " ".join(getattr(arguments, 'col_name'))
        self.type = getattr(arguments, 'type')

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
        parser.add_argument("-rd",
                            "--row_data",
                            type=str,
                            required=True,
                            help="The path to the data matrix.")
        parser.add_argument("-rn",
                            "--row_name",
                            nargs="*",
                            type=str,
                            required=False,
                            default="",
                            help="The name of -r / --row_data.")
        parser.add_argument("-cd",
                            "--col_data",
                            type=str,
                            required=True,
                            help="The path to the data matrix.")
        parser.add_argument("-cn",
                            "--col_name",
                            nargs="*",
                            type=str,
                            required=False,
                            default="",
                            help="The name of -c / --col_data.")
        parser.add_argument("-t",
                            "--type",
                            type=str,
                            required=True,
                            choices=["WG1", "WG2", "WG3"],
                            help="")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Searching row data.")
        row_df = self.get_pools(self.row_data_path)
        if row_df.shape[0] == 0:
            print("\tNo data found.")
            exit()

        print("Searching col data.")
        col_df = self.get_pools(self.col_data_path)
        if col_df.shape[0] == 0:
            print("\tNo data found.")
            exit()

        print("Loading and compare per folder.")
        row_df.columns = ["Folder", self.row_name + " Metadata"]
        col_df.columns = ["Folder", self.col_name + " Metadata"]
        df = self.load_and_compare(row_df=row_df, col_df=col_df)
        print(df)

        print("Summary per metadata feature:")
        summary_df = df[["Feature", "Match"]].groupby("Feature").mean()
        for index, row in summary_df.iterrows():
            print("\tColumn: {}\t{:.2f} match".format(index, row["Match"]))


    def get_pools(self, inpath):
        if self.type == "WG1":
            return self.get_wg1_metadata(inpath=inpath)
        elif self.type == "WG2":
            return self.get_wg2_metadata(inpath=inpath)
        elif self.type == "WG3":
            return self.get_wg3_metadata(inpath=inpath)
        else:
            return pd.DataFrame()

    @staticmethod
    def get_wg1_metadata(inpath):
        data = []
        fpaths = glob.glob(os.path.join(inpath, "*"))
        for fpath in fpaths:
            folder = os.path.basename(fpath)
            metadata_path = os.path.join(inpath, folder, "CombinedResults", "Final_Assignments_demultiplexing_doublets.tsv.gz")
            if os.path.exists(metadata_path):
                data.append([folder, metadata_path])

        return pd.DataFrame(data, columns=["Folder", "Metadata"])

    @staticmethod
    def get_wg2_metadata(inpath):
        data = []
        metadata_path = os.path.join(inpath, "map", "azimuth_all.metadata.tsv.gz")
        if os.path.exists(metadata_path):
            data.append(["ALL", metadata_path])

        return pd.DataFrame(data, columns=["Folder", "Metadata"])

    @staticmethod
    def get_wg3_metadata(inpath):
        data = []
        # fpaths = glob.glob(os.path.join(inpath, "expression_input", "pools", "*"))
        # for fpath in fpaths:
        #     folder = os.path.basename(fpath).split(".")[0]
        #     metadata_path = os.path.join(inpath, "expression_input/pools/" + folder + ".metadata.tsv.gz")
        #     if os.path.exists(metadata_path):
        #         data.append([folder, metadata_path])

        metadata_path = os.path.join(inpath, "expression_input", "metadata.tagged.tsv.gz")
        if os.path.exists(metadata_path):
            data.append(["ALL", metadata_path])

        return pd.DataFrame(data, columns=["Folder", "Metadata"])

    def load_and_compare(self, row_df, col_df):
        data = []

        df = row_df.merge(col_df, on="Folder")
        for _, row in df.iterrows():
            row_df = pd.read_csv(row[self.row_name + " Metadata"], sep="\t", dtype=str)
            row_df.fillna("NA", inplace=True)
            row_n = row_df.shape[0]

            col_df = pd.read_csv(row[self.col_name + " Metadata"], sep="\t", dtype=str)
            col_df.fillna("NA", inplace=True)
            col_n = col_df.shape[0]

            overlapping_columns = [column for column in row_df.columns if column in col_df.columns]

            row_df.set_index("Barcode", inplace=True)
            col_df.set_index("Barcode", inplace=True)
            overlapping_columns.remove("Barcode")
            if self.type == "WG2" and row["Folder"] == "ALL":
                overlapping_columns.remove("Pool")

            row_df.columns = ["{} {}".format(self.row_name, column) for column in row_df.columns]
            col_df.columns = ["{} {}".format(self.col_name, column) for column in col_df.columns]

            df = row_df.merge(col_df, left_index=True, right_index=True)
            overlap_n = df.shape[0]

            if row["Folder"] == "ALL":
                for pool in df[self.row_name + " Pool"].unique():
                    row_n = (row_df[self.row_name + " Pool"] == pool).sum()
                    col_n = (col_df[self.col_name + " Pool"] == pool).sum()

                    subset = df.loc[df[self.row_name + " Pool"] == pool, :]
                    overlap_n = subset.shape[0]

                    for column in overlapping_columns:
                        match = self.match(df=subset, column=column)
                        data.append([pool, row_n, col_n, overlap_n, column, match])
            else:
                for column in overlapping_columns:
                    match = self.match(df=df, column=column)
                    data.append([row["Folder"], row_n, col_n, overlap_n, column, match])

        return pd.DataFrame(data, columns=["Folder", self.row_name + " N", self.col_name + " N", "Overlap N", "Feature", "Match"])

    def match(self, df, column):
        tmp_df = df[["{} {}".format(self.row_name, column),
                     "{} {}".format(self.col_name, column)]].copy()
        tmp_df.columns = ["X", "Y"]

        try:
            tmp_df["X-float"] = tmp_df["X"].astype(float)
            tmp_df["Y-float"] = tmp_df["Y"].astype(float)

            if tmp_df["X-float"].std() == 0 or tmp_df["Y-float"].std() == 0:
                raise ValueError

            pearson_coef, _ = stats.pearsonr(tmp_df["Y-float"], tmp_df["X-float"])
            return pearson_coef
        except ValueError:
            return (tmp_df["X"] == tmp_df["Y"]).sum() / df.shape[0]


    def print_arguments(self):
        print("Arguments:")
        print("  > Row data path: {}".format(self.row_data_path))
        print("  > Row name: {}".format(self.row_name))
        print("  > Col data path: {}".format(self.col_data_path))
        print("  > Col name: {}".format(self.col_name))
        print("  > Type: {}".format(self.type))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
