#!/usr/bin/env python3

"""
File:         parse_ega_metadata.py
Created:      2023/07/13
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
import os

# Third party imports.
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Parse EGA Metadata"
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

./parse_ega_metadata.py

### Columbia2022 ###
./parse_ega_metadata.py \
    -i /groups/umcg-biogen/tmp02/input/rawdata/single-cell/Roche/Bryois_Columbia_EGAD00001009168/metadata \
    -o /groups/umcg-biogen/tmp02/input/processeddata/single-cell/Columbia2022/

### RocheAD2022 ###
./parse_ega_metadata.py \
    -i /groups/umcg-biogen/tmp02/input/rawdata/single-cell/Roche/Bryois_RocheAD_EGAD00001009166/metadata \
    -o /groups/umcg-biogen/tmp02/input/processeddata/single-cell/RocheAD2022/

### RocheMS2022 ###
./parse_ega_metadata.py \
    -i /groups/umcg-biogen/tmp02/input/rawdata/single-cell/Roche/Bryois_RocheMS_EGAD00001009169/metadata \
    -o /groups/umcg-biogen/tmp02/input/processeddata/single-cell/RocheMS2022/
"""

class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.input_dir = getattr(arguments, 'input_dir')
        self.output_dir = getattr(arguments, 'output_dir')

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
                            help="")
        parser.add_argument("-o",
                            "--output_dir",
                            type=str,
                            required=True,
                            help="")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        sample_file_path = os.path.join(self.input_dir, "delimited_maps", "Sample_File.map")
        if os.path.exists(sample_file_path):
            sample_file_df = self.load_file(sample_file_path)
            sample_file_df.columns = ["Ind_Sample", "individualID", "fastqPath", "sampleID"]
            print(sample_file_df)
            sample_file_df = sample_file_df.loc[["fastq" in path for path in sample_file_df["fastqPath"]], :]
            sample_file_df["fastqID"] = [x.split("/")[1].split("_S1")[0] for x in sample_file_df["fastqPath"]]
            sample_file_df = sample_file_df[["fastqID", "individualID"]].drop_duplicates()
            if len(sample_file_df["individualID"].unique()) != sample_file_df.shape[0]:
                print("Error, not all individualID's are unique")
                exit()
            sample_file_df.sort_values(by="individualID", inplace=True)
            print(sample_file_df)
            # for index, row in sample_file_df.iterrows():
            #     print("{},{}".format(row[0], row[1]))
            self.save_file(df=sample_file_df, outpath=os.path.join(self.output_dir, "link_table.csv"))

    @staticmethod
    def load_file(inpath, header=None, index_col=None, sep="\t", skiprows=None,
                  nrows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         skiprows=skiprows, nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def save_file(df, outpath, header=True, index=False, sep=","):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    def print_arguments(self):
        print("Arguments:")
        print("  > Input directory: {}".format(self.input_dir))
        print("  > Output directory: {}".format(self.output_dir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
