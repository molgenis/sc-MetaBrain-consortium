#!/usr/bin/env python3

"""
File:         multiple_testing_correction.py
Created:      2023/04/24
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
import os

# Third party imports.
import numpy as np
import pandas as pd
from statsmodels.stats import multitest
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

# Local application imports.

"""
Syntax:
./multiple_testing_correction.py -h
"""

# Metadata
__program__ = "Multiple testing correction"
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
        self.wg3_folder = getattr(arguments, 'wg3_folder')
        self.n = getattr(arguments, 'n')

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
        parser.add_argument("--wg3_folder",
                            type=str,
                            default=None,
                            help="The path to the outputs from WG3")
        parser.add_argument("--n",
                            type=int,
                            default=5000,
                            help="Select top n genes")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Caculating FDR")
        for ancestry_path in glob.glob(os.path.join(self.wg3_folder, "output", "*")):
            ancestry = os.path.basename(ancestry_path)
            for cell_level_path in glob.glob(os.path.join(self.wg3_folder, "output", ancestry, "*")):
                cell_level = os.path.basename(cell_level_path)
                for cell_type_path in glob.glob(os.path.join(self.wg3_folder, "output", ancestry, cell_level, "*")):
                    cell_type = os.path.basename(cell_type_path)
                    eqtl_inpath = os.path.join(self.wg3_folder, "output", ancestry, cell_level, cell_type, "top_qtl_results_all.txt")
                    if not os.path.exists(eqtl_inpath):
                        print("Warning, 'top_qtl_results.txt.gz' file not found for {}/{}/{}".format(ancestry, cell_level, cell_type))
                        continue
                    expr_inpath = os.path.join(self.wg3_folder, "expression_input", ancestry, cell_level, cell_type, "PostQC", cell_type + ".Exp.txt")
                    if not os.path.exists(expr_inpath):
                        print("Warning, '{}.Exp.txt' file not found for {}/{}/{}".format(cell_type, ancestry, cell_level, cell_type))
                        continue

                    # Load the eQTL data.
                    eqtl_df = self.load_file(eqtl_inpath, header=0, index_col=None)

                    # Load the expression data.
                    expr_df = self.load_file(expr_inpath, index_col=0)

                    # Select the genes.
                    avg_expr = expr_df.mean(axis=1)
                    avg_expr.sort_values(inplace=True, ascending=False)
                    top_genes = set(avg_expr.iloc[:self.n].index)
                    eqtl_df = eqtl_df.loc[eqtl_df["feature_id"].isin(top_genes), :]

                    # Calculate the FDR.
                    eqtl_df["Global_FDR"] = self.qvalues(eqtl_df["empirical_feature_p_value"])
                    eqtl_df["BH-FDR"] = multitest.multipletests(eqtl_df["p_value"], method='fdr_bh')[1]
                    print(eqtl_df)

                    # Save the output.
                    self.save_file(df=eqtl_df, outpath=eqtl_inpath.replace(".txt.gz", "_top{}_FDR_added.txt.gz".format(self.n)))
                    # self.save_file(df=eqtl_df, outpath=eqtl_inpath.replace(".txt.gz", "_top{}_FDR_added.xlsx".format(self.n)))

                    print("\tAncestry: {}\tCell level: {}\tCell type: {}\tN-genes: {:,}\tN-eQTLs {:,}\tN-eQTLs (BH-FDR) {:,}".format(ancestry, cell_level, cell_type, eqtl_df.shape[0], np.sum(eqtl_df["Global_FDR"] < 0.05), np.sum(eqtl_df["BH-FDR"] < 0.05)))

    @staticmethod
    def load_file(inpath, header=0, index_col=None, sep="\t"):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col)

        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))

        return df

    @staticmethod
    def qvalues(p):
        qvalue = importr("qvalue")
        pvals = robjects.FloatVector(p)
        qobj = robjects.r['qvalue'](pvals)
        return np.array(qobj.rx2('qvalues'))

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t", na_rep="NA",
                  sheet_name="Sheet1"):
        if outpath.endswith('xlsx'):
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

        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    def print_arguments(self):
        print("Arguments:")
        print("  > Workgroup 3 folder:    {}".format(self.wg3_folder))
        print("  > Top N genes:           {}".format(self.n))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
