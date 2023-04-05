#!/usr/bin/env python3

"""
File:         overview_scmetabrain_samples.py
Created:      2023/01/19
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
import pandas as pd
import numpy as np
import os

# Third party imports.

# Local application imports.

# Metadata
__program__ = "Overview scMetaBrain Samples"
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

./overview_scmetabrain_samples.py
"""


class main():
    def __init__(self):
        self.imputation_workdir = "/groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2023-02-02-AMP_AD/"
        self.amp_ad_psam = "Step1-Imputation/pca_sex_checks/updated_psam.psam"
        self.amp_ad_vcf = "Step1-Imputation/vcf_all_merged/imputed_hg38.vcf.gz"
        self.rosmap_key = "/groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/data/ROSMAP_IDkey.csv"
        self.rosmap_bio = "/groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/data/ROSMAP_biospecimen_metadata.csv"
        self.rosmap_cli = "/groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/data/ROSMAP_clinical.csv"
        self.mathys = "/groups/umcg-biogen/tmp01/input/rawdata/single-cell/Mathys2019/metadata/snRNAseqPFC_BA10_biospecimen_metadata.csv"
        self.zhou = "/groups/umcg-biogen/tmp01/input/rawdata/single-cell/Zhou2020/snRNAseqAD_TREM2_biospecimen_metadata.csv"
        self.cain = "/groups/umcg-biogen/tmp01/input/rawdata/single-cell/Cain2020/unique_specimenID.txt.gz"
        self.fujita = "/groups/umcg-biogen/tmp01/input/rawdata/single-cell/Fujita2022/unique_individualID.txt.gz"

        # Mathys 2018 = 35 individuals
        # Zhou 2020 = 32 individuals
        # Cain 2020 = 24 individuals
        # Fujita 2022 = 465 individuals
        # total = 556 individuals

    def start(self):
        amp_ad_psam_df = self.load_file(inpath=os.path.join(self.imputation_workdir, self.amp_ad_psam))
        # print(amp_ad_psam_df)
        # print(amp_ad_psam_df["Provided_Ancestry"].value_counts())

        rosmap_key_df = self.load_file(inpath=self.rosmap_key, sep=",")
        # print(rosmap_key_df)

        rosmap_bio_df = self.load_file(inpath=self.rosmap_bio, sep=",")
        # print(rosmap_bio_df)
        rosmap_specimen_to_individual_dict = dict(zip(rosmap_bio_df["specimenID"], rosmap_bio_df["individualID"]))

        rosmap_cli_df = self.load_file(inpath=self.rosmap_cli, sep=",")
        # print(rosmap_cli_df)
        rosmap_individual_to_proj_dict = dict(zip(rosmap_cli_df["individualID"], rosmap_cli_df["projid"]))

        mathys_df = self.load_file(inpath=self.mathys, sep=",")
        mathys_df["Dataset"] = "Mathys"
        n_mathys = mathys_df.shape[0]
        print(mathys_df)

        zhou_df = self.load_file(inpath=self.zhou, sep=",")
        zhou_df["Dataset"] = "Zhou"
        n_zhou = zhou_df.shape[0]
        print(zhou_df)
        zhou_df["projid"] = zhou_df["individualID"].map(rosmap_individual_to_proj_dict)

        cain_df = self.load_file(inpath=self.cain, sep=",")
        cain_df["Dataset"] = "Cain"
        n_cain = cain_df.shape[0]
        print(cain_df)
        cain_df["individualID"] = cain_df["specimenID"].map(rosmap_specimen_to_individual_dict)
        cain_df["projid"] = cain_df["individualID"].map(rosmap_individual_to_proj_dict)

        fujita_df = self.load_file(inpath=self.fujita, sep=",")
        fujita_df["Dataset"] = "Fujita"
        n_fujita = fujita_df.shape[0]
        print(fujita_df)
        fujita_df["projid"] = fujita_df["individualID"].map(rosmap_individual_to_proj_dict)

        datasets_df = pd.concat([mathys_df[["projid", "specimenID", "individualID", "Dataset"]],
                                 zhou_df[["projid", "specimenID", "individualID", "Dataset"]],
                                 cain_df[["projid", "specimenID", "individualID", "Dataset"]],
                                 fujita_df[["projid", "individualID", "Dataset"]],
                                 ])
        print("######")
        print(datasets_df)
        print("")

        duplicates = datasets_df["projid"].value_counts()
        duplicates = duplicates[duplicates > 2]
        print("Duplicates:")
        for projid, occurances in duplicates.iteritems():
            print(projid, occurances)
            print(datasets_df.loc[datasets_df["projid"] == projid, :])
            print("")
        exit()

        datasets_df = datasets_df.merge(rosmap_key_df[["projid", "wgs_id"]], on="projid", how="left").merge(amp_ad_psam_df[["IID", "Provided_Ancestry"]], left_on="wgs_id", right_on="IID", how="left")
        # datasets_df = datasets_df.merge(rosmap_key_df[["projid", "wgs_id"]], on="projid", how="left")
        print(datasets_df)
        print(datasets_df["Provided_Ancestry"].value_counts())
        print(datasets_df["Dataset"].value_counts())

        amp_ad_samples = self.load_file(inpath=os.path.join(self.imputation_workdir, self.amp_ad_vcf), skiprows=39, nrows=1).columns.tolist()[9:]
        # print(amp_ad_samples[:10])

        # print("Checking samples:")
        # for _, (proj_id, specimen_id, individual_id, dataset, wgs_id, iid, ancestry) in datasets_df.iterrows():
        #     # if np.isnan(proj_id) or np.isnan(proj_id) or np.isnan(proj_id):
        #     #     exit()
        #     if wgs_id not in amp_ad_samples:
        #         print("\tprojid: {}\tspecimen_id: {}\tindividual_id: {}\tDataset: {}\twgs_id: {}\tIID: {}\tProvided_Ancestry: {}\t".format(proj_id, specimen_id, individual_id, dataset, wgs_id, iid, ancestry))

        for dataset, n in zip(["Mathys", "Zhou", "Cain", "Fujita"], [n_mathys, n_zhou, n_cain, n_fujita]):
            subset = datasets_df.loc[datasets_df["Dataset"] == dataset, :].copy()
            n_rna_seqsamples = subset.shape[0]
            subset = subset.loc[subset["wgs_id"].isin(amp_ad_samples), :]
            n_genotype_samples = subset.shape[0]
            print("Dataset: {}\tN: {}\tN-RNA-seq: {}\tN-Genotype: {}".format(dataset, n, n_rna_seqsamples, n_genotype_samples))
            self.save_file(
                df=subset[["wgs_id"]],
                index=False,
                header=False,
                outpath=os.path.join(self.imputation_workdir, "Step1-Imputation/split_per_dataset/{}_samples.txt".format(dataset))
            )

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