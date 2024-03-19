#!/usr/bin/env python3

"""
File:         prepare_Mathys2019_wg3_files.py
Created:      2023/04/12
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
import gzip
import os
import glob

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
./prepare_Mathys2019_wg3_files.py -h
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        self.workdir = "/groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2024-01-10-scMetaBrain-Workgroup3DownstreamAnalyses/2024-02-01-CellAnnotations"
        # self.barcodes_dir = "/groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-08-28-CellBender-v0.3.0/2023-09-10-Mathys2019-Default"
        self.barcodes_dir = "/groups/umcg-biogen/tmp02/input/processeddata/single-cell/Mathys2019/"
        self.meta_data1 = "/groups/umcg-biogen/tmp02/input/processeddata/single-cell/Mathys2019/GTE.tsv"
        self.meta_data2 = "/groups/umcg-biogen/tmp02/input/rawdata/single-cell/Mathys2019/metadata/snRNAseqPFC_BA10_biospecimen_metadata.csv"
        self.meta_data3 = "/groups/umcg-biogen/tmp02/input/rawdata/single-cell/Mathys2019/metadata/snRNAseqPFC_BA10_assay_scRNAseq_metadata.csv"
        self.meta_data4 = "/groups/umcg-biogen/prm03/projects/2022-DeKleinEtAl/output/2020-02-03-phenotype-table/2020-09-04.brain.phenotypes.withReannotatedDiagnosis.txt.gz"

    def start(self):
        print("Loading data")
        barcodes_list = []
        for fpath in glob.glob(os.path.join(self.barcodes_dir, "*")):
            pool = os.path.basename(fpath)
            # barcodes_fpath = os.path.join(fpath, "cellbender_remove_background_output_cell_barcodes.csv")
            barcodes_fpath = os.path.join(fpath, "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
            if not os.path.exists(barcodes_fpath):
                continue
            barcodes_df = pd.read_csv(barcodes_fpath, sep="\t", header=None, index_col=None)
            barcodes_df.columns = ["Barcode"]
            barcodes_df["projid"] = int(pool)
            barcodes_df["Barcode"] = [x.split("-")[0] + "_" + pool for x in barcodes_df["Barcode"]]
            barcodes_list.append(barcodes_df)
        barcodes_df = pd.concat(barcodes_list, axis=0)
        print(barcodes_df)

        meta_data_df1 = self.load_file(self.meta_data1, sep="\t", index_col=None)
        meta_data_df1.columns = ["genotype_id", "projid"]
        print(meta_data_df1)

        meta_data_df2 = self.load_file(self.meta_data2, index_col=None)
        meta_data_df2 = meta_data_df2[["projid", "specimenID", "tissue"]].copy()
        meta_data_df2.columns = ["projid", "specimenID", "biomaterial"]
        meta_data_df2.drop_duplicates(inplace=True)
        print(meta_data_df2)

        meta_data_df3 = self.load_file(self.meta_data3, index_col=None)
        meta_data_df3 = meta_data_df3[["specimenID", "platform", "sequencingBatch", "libraryPreparationMethod"]].copy()
        meta_data_df3.columns = ["specimenID", "sequencing_platform", "sequencingBatch", "scrna_platform"]
        meta_data_df3.drop_duplicates(inplace=True)
        print(meta_data_df3)

        meta_data_df4 = self.load_file(self.meta_data4, sep="\t", index_col=None)
        meta_data_df4 = meta_data_df4[["genotype_id", "reannotated_diangosis"]].copy()
        meta_data_df4.columns = ["genotype_id", "sample_condition"]
        meta_data_df4.dropna(inplace=True)
        print(meta_data_df4)

        print(meta_data_df1.merge(meta_data_df4, on="genotype_id", how="left"))

        df = barcodes_df.merge(meta_data_df1, on="projid", how="left").merge(meta_data_df2, on="projid", how="left").merge(meta_data_df3, on="specimenID", how="left").merge(meta_data_df4, on="genotype_id", how="left")
        # df = df.merge(meta_data_df3, on="specimenID", how="left")

        print("Building dataframe")
        df["sequencing_run"] = df["projid"].astype(str) + "_run" + df["sequencingBatch"].astype(str)
        # df["sequencing_lane"] = df["projid"].astype(str) + "_lane" + df["lane"].astype(str)
        df["sequencing_lane"] = df["projid"].astype(str) + "_lane1"
        df["plate_based"] = "N"
        df["umi_based"] = "Y"
        df["sorting"] = np.nan
        df["cell_treatment"] = "UT"
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