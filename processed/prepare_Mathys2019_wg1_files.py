#!/usr/bin/env python3

"""
File:         prepare_Mathys2019_wg1_files.py
Created:      2022/10/14
Last Changed: 2022/10/19
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
__program__ = "Prepare Mathys2019 Workgroup 1 Files"
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

./prepare_Mathys2019_wg1_files.py \
    --workdir /groups/umcg-biogen/tmp01/input/processeddata/single-cell/Mathys2019 \
    --sample_id_table /groups/umcg-biogen/tmp01/input/rawdata/single-cell/Mathys2019/metadata/snRNAseqPFC_BA10_biospecimen_metadata.csv \
    --idkey /groups/umcg-biogen/tmp01/input/ROSMAP-scRNAseq/meta/ROSMAP_IDkey.csv \
    --gte /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR/GTE-EUR-AMPAD-ROSMAP-V2.txt.gz
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.workdir = getattr(arguments, 'workdir')
        self.sample_id_table = getattr(arguments, 'sample_id_table')
        self.n_individuals_per_pool = getattr(arguments, 'n_individuals_per_pool')
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
        parser.add_argument("--sample_id_table",
                            type=str,
                            required=True,
                            help="The sample-ID link matrix.")
        parser.add_argument("--n_individuals_per_pool",
                            type=int,
                            default=1,
                            help="The number of individuals per pool.")
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

        print("Loading sample-ID table.")
        sit_df = self.load_file(self.sample_id_table)
        print("\t  > Found {} samples".format(sit_df.shape[0]))

        print("Loading data")
        idkey_df = self.load_file(self.idkey)
        gte_df = self.load_file(self.gte, header=None, sep="\t")

        print("Merging data")
        wgs_dicht = dict(zip(idkey_df["projid"], idkey_df["wgs_id"]))
        sit_df["genotype_id"] = sit_df["projid"].map(wgs_dicht)

        gte_dicht = dict(zip(idkey_df.iloc[:, 0], idkey_df.iloc[:, 1]))
        sit_df["expression_id"] = sit_df["projid"].map(gte_dicht)

        sit_df = sit_df[["projid", "genotype_id", "expression_id"]].dropna()
        print("\t  > {} samples have matching GTE IDs".format(sit_df.shape[0]))

        self.save_file(df=sit_df,
                       outpath=os.path.join(self.workdir, "sample_id_table.txt.gz"),
                       sep="\t")

        print("Create sample sheet")
        sample_sheet_df = sit_df[["projid"]].copy()
        sample_sheet_df.columns = ["Pool"]
        sample_sheet_df["N_Individuals"] = self.n_individuals_per_pool
        self.save_file(df=sample_sheet_df,
                       outpath=os.path.join(self.workdir, "samplesheet.txt"),
                       sep="\t")

        print("Saving data")
        for _, (projid_id, genotype_id, _) in sit_df.iterrows():
            self.save_lines_to_file(
                outpath=os.path.join(self.out_dir, "{}.txt".format(projid_id)),
                lines=["{}_{}".format(projid_id, genotype_id)]
            )

    @staticmethod
    def load_file(inpath, header=0, index_col=None, sep=","):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def save_file(df, outpath, header=True, index=False, sep=","):
        df.to_csv(outpath, sep=sep, index=index, header=header)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    @staticmethod
    def save_lines_to_file(lines, outpath):
        with open(outpath, "w") as f:
            for line in lines:
                f.write(line + "\n")
        f.close()
        print("\tSaved jobfile: {}".format(os.path.basename(outpath)))

    def print_arguments(self):
        print("Arguments:")
        print("  > Working directory: {}".format(self.workdir))
        print("  > N-individuals per pool: {}".format(self.n_individuals_per_pool))
        print("  > Sample-ID file: {}".format(self.sample_id_table))
        print("  > IDkey file: {}".format(self.idkey))
        print("  > GTE file: {}".format(self.gte))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()