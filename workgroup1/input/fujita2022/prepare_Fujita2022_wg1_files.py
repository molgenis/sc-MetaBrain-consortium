#!/usr/bin/env python3

"""
File:         prepare_Fujita2022_wg1_files.py
Created:      2023/02/27
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
__program__ = "Prepare Fujita2022 Workgroup 1 Files"
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

./prepare_Fujita2022_wg1_files.py \
    --workdir /groups/umcg-biogen/tmp01/input/processeddata/single-cell/Fujita2022 \
    --sample_id_table /groups/umcg-biogen/tmp01/input/rawdata/single-cell/Fujita2022/unique_individualID.txt.gz \
    --idkey /groups/umcg-biogen/tmp01/input/ROSMAP-scRNAseq/meta/ROSMAP_IDkey.csv \
    --clinical /groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/data/ROSMAP_clinical.csv \
    --vcf /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2023-02-02-AMP_AD/Step1-Imputation/vcf_all_merged/imputed_hg38.vcf.gz
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.workdir = getattr(arguments, 'workdir')
        self.sample_id_table = getattr(arguments, 'sample_id_table')
        self.n_individuals_per_pool = getattr(arguments, 'n_individuals_per_pool')
        self.idkey = getattr(arguments, 'idkey')
        self.clinical = getattr(arguments, 'clinical')
        self.vcf = getattr(arguments, 'vcf')

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
        parser.add_argument("--clinical",
                            type=str,
                            required=True,
                            help="")
        parser.add_argument("--vcf",
                            type=str,
                            required=True,
                            help="")

        return parser.parse_args()

    def start(self):
        self.print_arguments()
        print("Loading data")
        sit_df = self.load_file(self.sample_id_table)
        print("\t  > Found {} samples".format(sit_df.shape[0]))
        idkey_df = self.load_file(self.idkey)
        clinical_df = self.load_file(self.clinical)

        print("Adding genotype column")
        projid_dict = dict(zip(clinical_df["individualID"], clinical_df["projid"]))
        wgs_dict = dict(zip(idkey_df["projid"], idkey_df["wgs_id"]))
        sit_df["projid"] = sit_df["individualID"].map(projid_dict)
        sit_df["genotype_id"] = sit_df["projid"].map(wgs_dict)

        print("Subsetting GTE columns")
        gte_df = sit_df[["projid", "genotype_id"]].copy()
        gte_df["projid"] = gte_df["projid"].astype(str)
        print(gte_df)
        del sit_df

        print("Checking if all data is available")
        geno_samples = set(self.load_file(inpath=self.vcf, sep="\t", skiprows=39, nrows=1).columns.tolist()[9:])

        mask = []
        for _, (expr_id, geno_id) in gte_df.iterrows():
            geno_found = True
            if geno_id not in geno_samples:
                geno_found = False

            if not geno_found:
                print("\tExpression ID: {} [NA]\tGenotype ID: {} [{}]".format(expr_id, geno_id, geno_found))
                mask.append(False)
            else:
                mask.append(True)
        print("")
        print("\t  > {} samples available".format(sum(mask)))

        print("Subsetting data")
        gte_df = gte_df.loc[mask, :].copy()

        self.save_file(df=gte_df,
                       outpath=os.path.join(self.workdir, "Fujita_GTE.txt.gz"),
                       sep="\t")
        self.save_file(df=gte_df[["genotype_id"]],
                       header=False,
                       outpath=os.path.join(self.workdir, "Fujita_genotype_samples.txt"),
                       sep="\t")

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
        print("  > Working directory: {}".format(self.workdir))
        print("  > N-individuals per pool: {}".format(self.n_individuals_per_pool))
        print("  > Sample-ID file: {}".format(self.sample_id_table))
        print("  > IDkey file: {}".format(self.idkey))
        print("  > Clinical file: {}".format(self.clinical))
        print("  > VCF file: {}".format(self.vcf))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()