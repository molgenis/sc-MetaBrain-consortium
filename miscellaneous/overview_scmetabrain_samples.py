#!/usr/bin/env python3

"""
File:         overview_scmetabrain_samples.py
Created:      2023/01/19
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
from pathlib import Path
import itertools
import os

# Third party imports.
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import upsetplot as up
import matplotlib.pyplot as plt

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
./overview_scmetabrain_samples.py -h
"""


class main():
    def __init__(self):
        self.rosmap_bio = "/groups/umcg-biogen/tmp02/input/rawdata/single-cell/AMP-AD/metadata/ROSMAP_biospecimen_metadata.csv"
        self.rosmap_cli = "/groups/umcg-biogen/tmp02/input/rawdata/single-cell/AMP-AD/metadata/ROSMAP_clinical.csv"
        self.roche_ad = "/groups/umcg-biogen/tmp02/input/rawdata/single-cell/Roche/Bryois_RocheAD_EGAD00001009166/metadata/delimited_maps/Sample_File.map"
        self.roche_ms = "/groups/umcg-biogen/tmp02/input/rawdata/single-cell/Roche/Bryois_RocheMS_EGAD00001009169/metadata/delimited_maps/Sample_File.map"
        self.roche_columbia1 = "/groups/umcg-biogen/tmp02/input/rawdata/single-cell/Roche/Bryois_Columbia_EGAD00001009168/metadata/delimited_maps/Sample_File.map"
        self.roche_columbia2 = "/groups/umcg-biogen/tmp02/input/rawdata/single-cell/Roche/key_columbia_rosmap.txt"
        self.mathys = "/groups/umcg-biogen/tmp02/input/rawdata/single-cell/Mathys2019/metadata/snRNAseqPFC_BA10_biospecimen_metadata.csv"
        self.zhou = "/groups/umcg-biogen/tmp02/input/rawdata/single-cell/Zhou2020/metadata/snRNAseqAD_TREM2_biospecimen_metadata.csv"
        self.cain = "/groups/umcg-biogen/tmp02/input/rawdata/single-cell/Cain2023/processed/ROSMAP_Brain.snRNAseq_metadata_cells_20230420.csv"
        self.fujita = "/groups/umcg-biogen/tmp02/input/rawdata/single-cell/Fujita2022/Other/cell-annotation.csv"
        self.sea_ad = "/groups/umcg-biogen/tmp02/input/rawdata/single-cell/SEA-AD/Metadata/SEA-AD_biospecimen_metadata.csv"
        self.mc_brad = "/groups/umcg-biogen/tmp02/input/rawdata/single-cell/MC-BrAD/Metadata/MC-BrAD_biospecimen_metadata.csv"
        self.sun = "/groups/umcg-biogen/tmp02/input/rawdata/single-cell/Sun2023/Metadata/MIT_ROSMAP_Multiomics_biospecimen_metadata.csv"
        self.mam_body = "/groups/umcg-biogen/tmp02/input/rawdata/single-cell/MammillaryBody/Metadata/ROSMAP_MammillaryBody_biospecimen_metadata.csv"

        # Roche AD = 80
        # Roche MS = 246
        # Roche Columbia = 24
        # Mathys 2018 = 48
        # Zhou 2020 = 32
        # Cain 2020 = 24
        # Fujita 2022 = 465
        # total = 556

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent), 'plot')
        print(self.outdir)

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

    def start(self):
        rosmap_bio_df = self.load_file(inpath=self.rosmap_bio, sep=",")
        rosmap_specimen_to_individual_dict = dict(zip(rosmap_bio_df["specimenID"], rosmap_bio_df["individualID"]))

        rosmap_cli_df = self.load_file(inpath=self.rosmap_cli, sep=",")
        rosmap_individual_to_proj_dict = dict(zip(rosmap_cli_df["individualID"], rosmap_cli_df["projid"]))

        ################################################################################################

        print("RocheAD")
        roche_ad_df = self.load_file(inpath=self.roche_ad, header=None).iloc[:, [0]]
        roche_ad_df["projid"] = ["Ind" + x.split("_")[0] for x in roche_ad_df[0]]
        roche_ad_df["sampleID"] = ["Ind{}-Sample{}".format(*x.split("_"))  for x in roche_ad_df[0]]
        roche_ad_df["individualID"] = roche_ad_df["projid"]
        roche_ad_df.drop([0], axis=1, inplace=True)
        roche_ad_df.drop_duplicates(inplace=True)
        print(roche_ad_df)

        ################################################################################################

        print("RocheMS")
        roche_ms_df = self.load_file(inpath=self.roche_ms, header=None).iloc[:, [0, 2]]
        roche_ms_df = roche_ms_df.loc[[not x.startswith("output_genotype") for x in roche_ms_df[2]], :]
        roche_ms_df["projid"] = ["Ind" + x.split("_")[0] for x in roche_ms_df[0]]
        roche_ms_df["sampleID"] = ["Ind{}-Sample{}".format(*x.split("_"))  for x in roche_ms_df[0]]
        roche_ms_df["individualID"] = roche_ms_df["projid"]
        roche_ms_df.drop([0], axis=1, inplace=True)
        roche_ms_df.drop_duplicates(inplace=True)
        print(roche_ms_df)

        ################################################################################################

        print("RocheColumbia")
        # roche_columbia1_df = self.load_file(inpath=self.roche_columbia1, header=None).iloc[:, [0]]
        # roche_columbia1_df["projid"] = ["Ind" + x.split("_")[0] for x in roche_columbia1_df[0]]
        # roche_columbia1_df["sampleID"] = ["Ind{}-Sample{}".format(*x.split("_"))  for x in roche_columbia1_df[0]]
        # roche_columbia1_df["individualID"] = roche_columbia1_df["projid"]
        # roche_columbia1_df.drop([0], axis=1, inplace=True)
        # roche_columbia1_df.drop_duplicates(inplace=True)
        # print(roche_columbia1_df)

        roche_columbia_df = self.load_file(inpath=self.roche_columbia2)[["individual_id_anon", "sample_id_anon", "individual_id"]]
        roche_columbia_df.columns = ["projid", "sampleID", "specimenID"]
        roche_columbia_df["individualID"] = roche_columbia_df["specimenID"].map(rosmap_specimen_to_individual_dict)
        roche_columbia_df.loc[roche_columbia_df["specimenID"] == "ROS15738428", "individualID"] = "R4575675"
        roche_columbia_df.drop_duplicates(inplace=True)
        print(roche_columbia_df)

        ################################################################################################

        print("Mathys")
        mathys_df = self.load_file(inpath=self.mathys, sep=",")[["projid", "specimenID", "individualID"]]
        mathys_df.drop_duplicates(inplace=True)
        print(mathys_df)

        ################################################################################################

        print("Zhou")
        zhou_df = self.load_file(inpath=self.zhou, sep=",")[["individualID", "specimenID"]]
        zhou_df.drop_duplicates(inplace=True)
        zhou_df["projid"] = zhou_df["individualID"].map(rosmap_individual_to_proj_dict)
        print(zhou_df)

        ################################################################################################

        print("Cain")
        cain_df = self.load_file(inpath=self.cain, sep=",")[["specimenID"]]
        cain_df.drop_duplicates(inplace=True)
        cain_df["individualID"] = cain_df["specimenID"].map(rosmap_specimen_to_individual_dict)
        cain_df["projid"] = cain_df["individualID"].map(rosmap_individual_to_proj_dict)
        print(cain_df)

        ################################################################################################

        print("Fujita")
        fujita_df = self.load_file(inpath=self.fujita, sep=",")[["individualID"]]
        fujita_df.drop_duplicates(inplace=True)
        fujita_df["projid"] = fujita_df["individualID"].map(rosmap_individual_to_proj_dict)
        print(fujita_df)

        ################################################################################################

        print("SEA-AD")
        sea_ad_df = self.load_file(inpath=self.sea_ad, sep=",")[["individualID", "specimenID"]]
        sea_ad_df.drop_duplicates(inplace=True)
        print(sea_ad_df)

        ################################################################################################

        print("Sun")
        mc_brad_df = self.load_file(inpath=self.mc_brad, sep=",")[["individualID", "specimenID"]]
        mc_brad_df.drop_duplicates(inplace=True)
        print(mc_brad_df)

        ################################################################################################

        print("Sun")
        sun_ad_df = self.load_file(inpath=self.sun, sep=",")[["individualID", "specimenID"]]
        sun_ad_df.drop_duplicates(inplace=True)
        print(sun_ad_df)

        ################################################################################################

        print("SEA-MammillaryBody")
        mam_body_df = self.load_file(inpath=self.mam_body, sep=",")[["individualID", "specimenID"]]
        mam_body_df.drop_duplicates(inplace=True)
        print(mam_body_df)

        ################################################################################################

        print("Plotting samples.")
        for col in ["individualID"]:
            for filter in [False, True]:
                data = {
                    "RocheAD": set(roche_ad_df[col].values.tolist()),
                    "RocheMS": set(roche_ms_df[col].values.tolist()),
                    "RocheColumbia": set(roche_columbia_df[col].values.tolist()),
                    "MathysAD": set(mathys_df[col].values.tolist()),
                    "ZhouAD": set(zhou_df[col].values.tolist()),
                    "CainAD": set(cain_df[col].values.tolist()),
                    "Fujita2022": set(fujita_df[col].values.tolist()),
                    "SEA-AD": set(sea_ad_df[col].values.tolist()),
                    "MC-BrAD": set(mc_brad_df[col].values.tolist()),
                    "Sun": set(sun_ad_df[col].values.tolist()),
                    "MamBody": set(mam_body_df[col].values.tolist()),
                }
                if filter:
                    data = {
                        "RocheColumbia": set(roche_columbia_df[col].values.tolist()),
                        "MathysAD": set(mathys_df[col].values.tolist()),
                        "ZhouAD": set(zhou_df[col].values.tolist()),
                        "CainAD": set(cain_df[col].values.tolist()),
                        "Fujita2022": set(fujita_df[col].values.tolist()),
                        "Sun": set(sun_ad_df[col].values.tolist())
                    }
                counts = self.count(data)
                counts = counts[counts > 0]

                up.plot(counts, sort_by='cardinality', show_counts=True)
                plt.savefig(os.path.join(self.outdir, col + "_ROSMAP" if filter else "" + ".png"))
                plt.close()

                samples = set()
                for _, value in data.items():
                    samples.update(value)
                print(len(samples))

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
    def count(input_data):
        combinations = []
        cols = list(input_data.keys())
        for i in range(1, len(cols) + 1):
            combinations.extend(list(itertools.combinations(cols, i)))

        indices = []
        data = []
        for combination in combinations:
            index = []
            for col in cols:
                if col in combination:
                    index.append(True)
                else:
                    index.append(False)

            background = set()
            for key in cols:
                if key not in combination:
                    work_set = input_data[key].copy()
                    background.update(work_set)

            overlap = None
            for key in combination:
                work_set = input_data[key].copy()
                if overlap is None:
                    overlap = work_set
                else:
                    overlap = overlap.intersection(work_set)

            duplicate_set = overlap.intersection(background)
            length = len(overlap) - len(duplicate_set)

            indices.append(index)
            data.append(length)

        s = pd.Series(data,
                      index=pd.MultiIndex.from_tuples(indices, names=cols))
        s.name = "value"
        return s


if __name__ == '__main__':
    m = main()
    m.start()