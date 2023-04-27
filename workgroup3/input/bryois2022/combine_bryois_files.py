#!/usr/bin/env python3

"""
File:         combine_bryois_files.py
Created:      2023/04/24
Last Changed:
Author:       M.Vochteloo

Copyright (C) 2020 M.Vochteloo
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
import numpy as np
import pandas as pd

# Local application imports.

"""
Syntax:
./combine_bryois_files.py \
    --work_dir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-20-ReplicateInBryois/Bryois2022
    
"""

# Metadata
__program__ = "Combine Bryois files"
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


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.work_dir = getattr(arguments, 'work_dir')
        self.verbose = getattr(arguments, 'verbose')

        self.cell_type_dict = {
            "Astrocytes": "AST",
            "Endothelial.cells": "END",
            "Excitatory.neurons": "EX",
            "Inhibitory.neurons": "IN",
            "Microglia": "MIC",
            "Oligodendrocytes": "OLI",
            "OPCs...COPs": "OPC",
            "Pericytes": "PER"
        }

        self.outdir = os.path.join(self.work_dir, "merged")
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
                            help="show program's version number and exit")
        parser.add_argument("--work_dir",
                            type=str,
                            required=True,
                            default=None,
                            help="The working directory.")
        parser.add_argument("--verbose",
                            action='store_true',
                            help="Print all info.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("\tLoading SNP position file")
        snp_pos_df = self.load_file(inpath=os.path.join(self.work_dir, "snp_pos.txt"), header=0, index_col=None)
        snp_pos_df.dropna(axis=1, inplace=True)
        snp_pos_df["chromosome"] = [chr.replace("chr", "") for chr in snp_pos_df["chr"]]
        snp_pos_df['allele1'] = np.minimum(snp_pos_df['effect_allele'], snp_pos_df['other_allele'])
        snp_pos_df['allele2'] = np.maximum(snp_pos_df['effect_allele'], snp_pos_df['other_allele'])
        snp_pos_df['alleles'] = snp_pos_df['allele1'] + snp_pos_df['allele2']
        snp_pos_df.drop(['allele1', 'allele2'], axis=1, inplace=True)
        snp_pos_df.columns = ["SNP_id" if col == "SNP" else col for col in snp_pos_df.columns]

        print("Loading Bryois et al. 2022 data")
        for filename, cell_type in self.cell_type_dict.items():
            print("  Processing {} [{}]".format(filename, cell_type))

            dfs = []
            for chr in range(1, 23):
                print("\tLoading chromosome {}".format(chr))
                ct_chr_path = os.path.join(self.work_dir, "{}.{}.gz".format(filename, chr))
                if not os.path.exists(ct_chr_path):
                    continue

                ct_chr_df = self.load_file(ct_chr_path, sep=" ", header=None, index_col=None)
                ct_chr_df.columns = ["Gene_id", "SNP_id", "Distance to TSS", "Nominal p-value", "Beta"]
                ct_chr_df = ct_chr_df.loc[ct_chr_df.groupby('Gene_id')["Nominal p-value"].idxmin(), :]
                dfs.append(ct_chr_df)

                del ct_chr_df
                break

            if len(dfs) == 0:
                print("  Warning, no data loaded for this cell type")
                return None

            df = pd.concat(dfs, axis=0)

            if df.shape[0] == 0:
                print("  Warning, no data loaded for this cell type")
                return None

            gene_id_df = df["Gene_id"].str.split("_", n=None, expand=True)
            gene_id_df.columns = ["HGNC", "ENSG"]
            df = pd.concat([gene_id_df, df], axis=1)
            del gene_id_df

            print("  Merging with SNP position file")
            df = df.merge(snp_pos_df, on="SNP_id", how="left")

            print("  Saving file")
            outpath = os.path.join(self.outdir, cell_type)
            self.save_file(df=df, outpath=outpath + ".txt.gz", index=False)
            # self.save_file(df=df, outpath=outpath + ".pkl")
            # self.save_file(df=df, outpath=outpath + ".xlsx")

    def load_file(self, inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)

        if self.verbose:
            print("\tLoaded dataframe: {} "
                  "with shape: {}".format(os.path.basename(inpath),
                                          df.shape))
        return df

    def save_file(self, df, outpath, header=True, index=True, sep="\t", na_rep="NA",
                  sheet_name="Sheet1"):
        if outpath.endswith('pkl'):
            df.to_pickle(outpath)
        elif outpath.endswith('xlsx'):
            df.to_excel(outpath,
                        sheet_name=sheet_name,
                        na_rep=na_rep,
                        header=header,
                        index=index)
        else:
            compression = 'infer'
            if outpath.endswith('.gz'):
                compression = 'gzip'

            df.to_csv(outpath, sep=sep, index=index, header=header,
                      compression=compression)

        if self.verbose:
            print("\tSaved dataframe: {} "
                  "with shape: {}".format(os.path.basename(outpath),
                                          df.shape))

    def print_arguments(self):
        print("Arguments:")
        print("  > Working directory:    {}".format(self.work_dir))
        print("  > Verbose:              {}".format(self.verbose))
        print("  > Output directory:     {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
