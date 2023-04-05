#!/usr/bin/env python3

"""
File:         check_ancestry_and_sex_update_table.py
Created:      2023/01/17
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
import os
import pandas as pd
import numpy as np

# Third party imports.

# Local application imports.

# Metadata
__program__ = "Check Ancesty and Sex Update Table"
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
  
./check_ancestry_and_sex_update_table.py
"""


class main():
    def __init__(self):
        self.ancestry_table = "/groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2023-01-20-AMP_AD/Step1-Imputation/pca_sex_checks/ancestry_update_remove.tsv"
        self.sex_table = "/groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2023-01-20-AMP_AD/Step1-Imputation/pca_sex_checks/check_sex_update_remove.tsv"
        self.bulk_sex_table = "/groups/umcg-biogen/prm03/projects/2022-DeKleinEtAl/output/2020-02-03-phenotype-table/2020-09-04.brain.phenotypes.withReannotatedDiagnosis.txt.gz"
        self.bulk_ancestry_table = "/groups/umcg-biogen/tmp01/input/AMP-AD/2022-01-18-SampleAssignment/allGTSamplesPopulationAssignment.txt.gz"

        self.sex_dict = {
            "M": 1,
            "F": 2,
            np.nan: 0,
            "no expression available": 0
        }

    def start(self):
        sc_ancestry_df = self.load_file(self.ancestry_table)
        print(sc_ancestry_df)

        sc_ancestry_df["UPDATE/REMOVE/KEEP"] = "UPDATE"
        self.save_file(df=sc_ancestry_df, outpath=self.ancestry_table)
        exit()

        # bulk_ancestry_df = self.load_file(self.bulk_ancestry_table)
        # bulk_ancestry_df.columns = ["IID", "Bulk ancestry"]
        # print(bulk_ancestry_df)
        #
        # ancestry_df = sc_ancestry_df.merge(bulk_ancestry_df, on="IID", how="left")
        # ancestry_df["UPDATE/REMOVE/KEEP"] = "UPDATE"
        # # ancestry_df.loc[ancestry_df["PCA_Assignment"] == ancestry_df["Bulk ancestry"], "UPDATE/REMOVE/KEEP"] = "UPDATE"
        # # ancestry_df.loc[ancestry_df["Provided_Ancestry"] == ancestry_df["Bulk ancestry"], "UPDATE/REMOVE/KEEP"] = "KEEP"
        # print(ancestry_df)
        #
        #
        # exit()

        sc_sex_df = self.load_file(self.sex_table)
        print(sc_sex_df)

        sc_sex_df["Round"] = 0
        sc_sex_df.loc[sc_sex_df["F"] <= 0.5, "Round"] = 2
        sc_sex_df.loc[sc_sex_df["F"] > 0.5, "Round"] = 1

        bulk_sex_df = self.load_file(self.bulk_sex_table)
        #  print(bulk_sex_df.loc[bulk_sex_df["individualID"] == "1950", :])
        bulk_sex_df.loc[bulk_sex_df["individualID"] == "1950", "sex.by.expression"] = "F"
        # print(bulk_sex_df.loc[bulk_sex_df["individualID"] == "1950", :])
        # exit()
        bulk_sex_df = bulk_sex_df.loc[:, ["individualID", "Gender", "sex.by.expression"]].copy()
        bulk_sex_df.drop_duplicates(inplace=True)

        bulk_sex_df["Bulk gender"] = 0
        bulk_sex_df.loc[bulk_sex_df["Gender"] == "M", "Bulk gender"] = 1
        bulk_sex_df.loc[bulk_sex_df["Gender"] == "F", "Bulk gender"] = 2

        bulk_sex_df["Bulk sex by expression"] = 0
        bulk_sex_df.loc[bulk_sex_df["sex.by.expression"] == "M", "Bulk sex by expression"] = 1
        bulk_sex_df.loc[bulk_sex_df["sex.by.expression"] == "F", "Bulk sex by expression"] = 2

        bulk_sex_df = bulk_sex_df.loc[:, ["individualID", "Bulk gender", "Bulk sex by expression"]]
        bulk_sex_df.columns = ["IID", "Bulk gender", "Bulk sex by expression"]
        print(bulk_sex_df)

        sex_df = sc_sex_df.merge(bulk_sex_df, on="IID", how="left")
        sex_df["Bulk gender"].fillna(0, inplace=True)
        sex_df["Bulk sex by expression"].fillna(0, inplace=True)
        sex_df["UPDATE/REMOVE/KEEP"] = "UPDATE"
        # sex_df.loc[sex_df["Bulk sex by expression"] == 0, "UPDATE/REMOVE/KEEP"] = "UPDATE"
        # sex_df.loc[sex_df["SNPSEX"] == sex_df["Bulk sex by expression"], "UPDATE/REMOVE/KEEP"] = "UPDATE"
        sex_df.loc[(sex_df["SNPSEX"] == 0) & (sex_df["PEDSEX"] == sex_df["Bulk sex by expression"]), "UPDATE/REMOVE/KEEP"] = "KEEP"
        print(sex_df)

        sex_df = sex_df.iloc[:, :7]
        print(sex_df)
        self.save_file(df=sex_df, outpath=self.sex_table)



    @staticmethod
    def load_file(inpath, header=0, index_col=None, sep="\t"):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col)
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

if __name__ == '__main__':
    m = main()
    m.start()