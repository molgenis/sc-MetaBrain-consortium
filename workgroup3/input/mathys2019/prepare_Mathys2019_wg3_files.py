#!/usr/bin/env python3

"""
File:         prepare_Mathys2019_wg3_files.py
Created:      2023/04/12
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
import gzip
import os

# Third party imports.
import numpy as np
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Prepare Mathys2019 Workgroup 3 Files"
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

./prepare_Mathys2019_wg3_files.py
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        self.workdir = "/groups/umcg-biogen/tmp01/input/processeddata/single-cell/datasets/Mathys2019"
        self.meta_data1 = "/groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-02-02-WorkGroup2CellType/2023-04-05-Mathys2019/all/data/seurat_metadata.csv"
        self.meta_data2 = "/groups/umcg-biogen/tmp01/input/rawdata/single-cell/Mathys2019/metadata/snRNAseqPFC_BA10_biospecimen_metadata.csv"
        self.meta_data3 = "/groups/umcg-biogen/tmp01/input/rawdata/single-cell/Mathys2019/metadata/snRNAseqPFC_BA10_assay_scRNAseq_metadata.csv"
        self.fastq_folder = "/groups/umcg-biogen/tmp01/input/rawdata/single-cell/Mathys2019/fastq"

    def start(self):
        print("Loading data")
        meta_data_df1 = self.load_file(self.meta_data1)
        meta_data_df1 = meta_data_df1[["projid"]].copy()
        meta_data_df1.index = ["{}_{}".format(index.split("-")[0], row["projid"]) for index, row in meta_data_df1.iterrows()]
        meta_data_df1.index.name = "Barcode"
        meta_data_df1["feature barcode"] = [barcode.split("-")[0] for barcode in meta_data_df1.index]
        meta_data_df1.reset_index(drop=False, inplace=True)
        meta_data_df1.drop_duplicates(inplace=True)
        print(meta_data_df1)
        # Pool in  this file???

        meta_data_df2 = self.load_file(self.meta_data2, index_col=None)
        meta_data_df2 = meta_data_df2[["projid", "specimenID", "tissue"]].copy()
        meta_data_df2.columns = ["projid", "specimenID", "biomaterial"]
        meta_data_df2.drop_duplicates(inplace=True)
        print(meta_data_df2)

        # df = meta_data_df1.merge(meta_data_df2, on="projid", how="left")
        #
        # df["lane"] = np.nan
        # for specimen_id in df["specimenID"].unique():
        #     for lane in np.arange(1, 5):
        #         fastq_file = os.path.join(self.fastq_folder, "{}_L00{}_R1_001.fastq.gz".format(specimen_id, lane))
        #         print("\tReading '{}'".format(fastq_file))
        #         feature_barcodes = self.parse_fastq(fastq_file)
        #         mask = (df["specimenID"] == specimen_id) & (df["feature barcode"].isin(feature_barcodes))
        #         print("\t  loaded {} unique feature barcodes, {} overlap".format(len(feature_barcodes), np.sum(mask)))
        #         df.loc[mask, "lane"] = lane
        #         print(df.loc[mask, :])
        # print(df)

        meta_data_df3 = self.load_file(self.meta_data3, index_col=None)
        meta_data_df3 = meta_data_df3[["specimenID", "platform", "sequencingBatch", "libraryPreparationMethod"]].copy()
        meta_data_df3.columns = ["specimenID", "sequencing_platform", "sequencingBatch", "scrna_platform"]
        meta_data_df3.drop_duplicates(inplace=True)
        print(meta_data_df3)

        df = meta_data_df1.merge(meta_data_df2, on="projid", how="left").merge(meta_data_df3, on="specimenID", how="left")
        # df = df.merge(meta_data_df3, on="specimenID", how="left")

        """
        Barcode                         sequencing_platform sequencing_run  sequencing_lane scrna_platform  plate_based umi_based   biomaterial sorting cell_treatment  sample_condition
        AAACCTGAGAAACCAT_180920_lane1   NovoSeq_100PE       180920_run1     180920_lane1    10x_3end_v3     N           Y           PBMC        none    UT              healthy
        TTTGTCAGTTGAGGTG_181231_lane2   NovoSeq_100PE       180920_run1     181231_lane2    10x_3end_v3     N           Y           PBMC        none    3hCA            healthy
        TTTGTCATCGCCTGTT_181231_lane2   BGISeq500_100PE     181231_run2     181231_lane2    10x_3end_v2     N           Y           CD4T        facs    UT              UC
        """

        print("Building dataframe")
        df["sequencing_run"] = df["projid"].astype(str) + "_run" + df["sequencingBatch"].astype(str)
        # df["sequencing_lane"] = df["projid"].astype(str) + "_lane" + df["lane"].astype(str)
        df["sequencing_lane"] = df["projid"].astype(str) + "_lane1"
        df["plate_based"] = np.nan
        df["umi_based"] = np.nan
        df["sorting"] = np.nan
        df["cell_treatment"] = "UT"
        df["sample_condition"] = "healthy"
        df = df[["Barcode", "sequencing_platform", "sequencing_run", "sequencing_lane", "scrna_platform", "plate_based", "umi_based", "biomaterial", "sorting", "cell_treatment", "sample_condition"]].copy()
        print(df)
        print(df.columns.tolist())

        print("Saving File")
        self.save_file(
            df=df,
            outpath=os.path.join(self.workdir, "Mathys_cell_annotations.tsv")
        )

    @staticmethod
    def load_file(inpath, header=0, index_col=0, sep=",", skiprows=None,
                  nrows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         skiprows=skiprows, nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def parse_fastq(inpath):
        feature_barcodes = set()
        with gzip.open(inpath, 'rt') as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    feature_barcodes.add(line[:16])
                # if i > 100:
                #     break
        f.close()

        return feature_barcodes

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