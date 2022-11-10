#!/usr/bin/env python3

"""
File:         prepare_amp_ad_sex_file.py
Created:      2022/11/10
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

# Third party imports.
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Prepare AMP-AD Sex File"
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


def start():
    df1 = pd.read_csv(
        "/groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/2022-11-03-FilteredGenotypes/2-custom_vcffilter/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_22.recalibrated_variants_norm_vcffilter-filtered.vcf.gz",
        sep="\t", header=0, index_col=None, skiprows=171, nrows=1)
    samples = [x for x in df1.columns[9:]]
    del df1

    ### ROSMAP ###

    df2 = pd.read_csv("/groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/data/ROSMAP_IDkey.csv", sep=",", header=0, index_col=None)

    df3 = pd.read_csv("/groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/data/ROSMAP_clinical.csv", sep=",", header=0, index_col=None)

    df4 = df2.merge(df3, on="projid")[["wgs_id", "msex"]].drop_duplicates().dropna()
    df4["msex"] = df4["msex"].map({1.0: "M", 0.0: "F"})
    id_to_sex = dict(zip(df4["wgs_id"], df4["msex"]))
    del df2, df3, df4

    ### MAYO ###

    df5 = pd.read_csv("/groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/data/MayoRNAseq_individual_metadata_031422.csv", sep=",", header=0, index_col=None)
    df5 = df5[["individualID", "sex"]].drop_duplicates().dropna()
    df5["individualID"] = df5["individualID"].astype(str)
    df5["sex"] = df5["sex"].map({"male": "M", "female": "F"})
    id_to_sex.update(dict(zip(df5["individualID"], df5["sex"])))
    del df5

    ### MSBB ###

    df6 = pd.read_csv("/groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/data/MSBB_biospecimen_metadata.csv", sep=",", header=0, index_col=None)

    df7 = pd.read_csv("/groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/data/MSBB_individual_metadata.csv", sep=",", header=0, index_col=None)

    df8 = df6.merge(df7, on="individualID")[["specimenID", "sex"]].drop_duplicates().dropna()
    df8["sex"] = df8["sex"].map({"male": "M", "female": "F"})
    id_to_sex.update(dict(zip(df8["specimenID"], df8["sex"])))
    del df6, df7, df8

    ###

    data = []
    missing_samples = []
    for sample in samples:
        if sample in id_to_sex.keys():
            data.append([sample, id_to_sex[sample]])
        else:
            data.append([sample, "NA"])
            missing_samples.append(sample)

    df = pd.DataFrame(data, columns=["sample", "sex"])
    df.to_csv("AMP_AD_sexdata.csv", sep=",", header=False, index=False)

    print(df)
    print(df["sex"].value_counts())
    print("Missing: {}".format(", ".join(missing_samples)))

    # sex is NA in MayoRNAseq_individual_metadata_031422.csv: 1925, 1950, 1957
    # individualID is NA in MSBB_biospecimen_metadata.csv: 71729, 71823, 76354, 76655
    # NA in ROSMAP_clinical.csv: SM-CJEJI (projid = 90214403)


if __name__ == '__main__':
    start()