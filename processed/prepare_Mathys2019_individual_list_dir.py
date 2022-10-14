#!/usr/bin/env python3

"""
File:         prepare_Mathys2019_individual_list_dir.py
Created:      2022/10/14
Last Changed:
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
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Prepare Mathys2019 individual_list_dir"
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

./prepare_Mathys2019_individual_list_dir.py \
    --workdir /groups/umcg-biogen/tmp01/input/processeddata/single-cell/Mathys2019 \
    --samplesheet_filepath /groups/umcg-biogen/tmp01/input/processeddata/single-cell/Mathys2019/samplesheet.txt \
    --idkey /groups/umcg-biogen/tmp01/input/ROSMAP-scRNAseq/meta/ROSMAP_IDkey.csv \
    --gte /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR/GTE-EUR-AMPAD-ROSMAP-V2.txt.gz
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.workdir = getattr(arguments, 'workdir')
        self.samplesheet = getattr(arguments, 'samplesheet_filepath')
        self.idkey = getattr(arguments, 'idkey')
        self.gte = getattr(arguments, 'gte')

        # Creating output directory
        self.out_dir = os.path.join(self.workdir, "individual_list_dir")
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

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
        parser.add_argument("--samplesheet_filepath",
                            type=str,
                            required=True,
                            help="")
        parser.add_argument("--idkey",
                            type=str,
                            required=True,
                            help="")
        parser.add_argument("--gte",
                            type=str,
                            required=True,
                            help="")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data")
        samplesheet_df = self.load_file(self.samplesheet)
        idkey_df = self.load_file(self.idkey, sep=",")
        gte_df = self.load_file(self.gte, header=None)
        print(samplesheet_df)
        print(idkey_df)
        print(gte_df)

        print("Merging data")
        df = samplesheet_df.copy()
        wgs_dicht = dict(zip(idkey_df["projid"], idkey_df["wgs_id"]))
        df["wgs_id"] = df["Pool"].map(wgs_dicht)

        gte_dicht = dict(zip(idkey_df.iloc[:, 0], idkey_df.iloc[:, 1]))
        df["GTE"] = df["Pool"].map(gte_dicht)

        df = df[["Pool", "GTE"]].dropna()

        print("Saving data")
        for _, (individual_id, genotype_id) in df.iterrows():
            self.save_file(
                outpath=os.path.join(self.out_dir, "{}.txt".format(individual_id)),
                lines=["{}_{}".format(individual_id, genotype_id)]
            )

    @staticmethod
    def load_file(inpath, header=0, index_col=None, sep="\t"):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def save_file(lines, outpath):
        with open(outpath, "w") as f:
            for line in lines:
                f.write(line + "\n")
        f.close()
        print("\tSaved jobfile: {}".format(os.path.basename(outpath)))

    def print_arguments(self):
        print("Arguments:")
        print("  > Working directory: {}".format(self.workdir))
        print("  > Samplesheet file: {}".format(self.samplesheet))
        print("  > IDkey file: {}".format(self.idkey))
        print("  > GTE file: {}".format(self.gte))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()