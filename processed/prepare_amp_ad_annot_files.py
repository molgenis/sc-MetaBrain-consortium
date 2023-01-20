#!/usr/bin/env python3

"""
File:         prepare_amp_ad_annot_files.py
Created:      2022/11/10
Last Changed: 2022/11/15
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
import numpy as np
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Prepare AMP-AD Annot Files"
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

./prepare_amp_ad_annot_files.py
"""


def start():
    psam_template = {
        '#FID': np.nan,
        'IID': np.nan,
        'PAT': 0,
        'MAT': 0,
        'SEX': 0,
        'Provided_Ancestry': np.nan,
        'genotyping_platform': np.nan,
        'array_available': "N",
        'wgs_available': "N",
        'wes_available': "N",
        'age': np.nan,
        'age_range': np.nan,
        'Study': np.nan,
        'smoking_status': np.nan,
        'hormonal_contraception_use_currently': np.nan,
        'menopause': np.nan,
        'pregnancy_status': np.nan
    }

    samples_df = pd.read_csv(
        "/groups/umcg-biogen/tmp01/input/AMP-AD/2017-12-08-joint-WGS/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_22.recalibrated_variants.vcf.gz",
        sep="\t", header=0, index_col=None, skiprows=198, nrows=1)
    samples = [str(x) for x in samples_df.columns[9:]]
    del samples_df

    print("Loading data")
    rosmap_df = construct_rosmap_df()
    mayo_df = construct_mayo_df()
    msbb_df = construct_msbb_df()

    # Merging data.
    coi = ["sex", "race", "platform", "wgs_available", "ageDeath", "Study"]
    df = pd.concat([rosmap_df[coi], mayo_df[coi], msbb_df[coi]], axis=0)

    # Construct PSAM data frame.
    psam_df = df.loc[[sample for sample in samples if sample in df.index], :].copy()
    psam_df.columns = ["SEX", "race", "genotyping_platform", "wgs_available", "ageDeath", "Study"]
    psam_df["SEX"] = round_float(psam_df["SEX"], na="0")
    psam_df["Provided_Ancestry"] = psam_df["race"].map(
        {"White": "EUR",
         "Black or African American": "AFR",
         "Hispanic or Latino": "AMR",
         "American Native or Alaskan Native": "AMR",
         "Black, Negro, African-American": "AFR",
         "NA": "EUR"})
    psam_df["Provided_Ancestry"].fillna("EUR", inplace=True)
    psam_df["age"] = round_float(psam_df["ageDeath"], na="NA")
    for col in psam_df.columns:
        if col not in psam_template:
            continue
        psam_df.loc[psam_df[col].isnull(), col] = psam_template[col]

    # Print overview.
    print("\nPSAM data:")
    print(psam_df)
    print("")
    for col in psam_df.columns:
        if col not in psam_template:
            continue
        print(col)
        print(psam_df[col].value_counts())
        print("")

    for column, default_value in psam_template.items():
        if column not in psam_df.columns:
            psam_df[column] = default_value

    psam_df = pd.concat([psam_df, pd.DataFrame(psam_template, index=[sample for sample in samples if sample not in psam_df.index])], axis=0)
    psam_df["#FID"] = psam_df.index
    psam_df["IID"] = psam_df.index
    psam_df["age_range"] = np.floor(psam_df["ageDeath"] / 10) * 10

    psam_df = psam_df.loc[samples, ['#FID', 'IID', 'PAT', 'MAT', 'SEX', 'Provided_Ancestry', 'genotyping_platform', 'array_available', 'wgs_available', 'wes_available', 'age', 'age_range', 'Study', 'smoking_status', 'hormonal_contraception_use_currently', 'menopause', 'pregnancy_status']]
    psam_df.fillna("NA", inplace=True)
    psam_df.to_csv("AMP_AD_annotdata.psam", sep="\t", index=False, header=True)

    sex_df = psam_df[["#FID", "SEX"]].copy()
    sex_df["SEX"] = sex_df["SEX"].map({np.nan: "M", 0: "M", 1: "M", 2: "F"})
    sex_df.to_csv("AMP_AD_sexdata.csv", sep=",", header=False, index=False)

    # sex is NA in MayoRNAseq_individual_metadata_031422.csv: 1925, 1950, 1957
    # individualID is NA in MSBB_biospecimen_metadata.csv: 71729, 71823, 76354, 76655
    # NA in ROSMAP_clinical.csv: SM-CJEJI (projid = 90214403)


def construct_rosmap_df():
    # Load data.
    idkey_df = pd.read_csv("/groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/data/ROSMAP_IDkey.csv", sep=",", header=0, index_col=None)

    # Load data.
    wgs_metadata_df = pd.read_csv("/groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/data/ROSMAP_assay_wholeGenomeSeq_metadata.csv", sep=",", header=0, index_col=None)
    wgs_metadata_df["specimenID"] = wgs_metadata_df["specimenID"].astype(str)
    wgs_metadata_df["wgs_available"] = "N"
    wgs_metadata_df.loc[wgs_metadata_df["assay"] == "wholeGenomeSeq", "wgs_available"] = "Y"

    biospecimen_metadata_df = pd.read_csv("/groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/data/ROSMAP_biospecimen_metadata.csv", sep=",", header=0, index_col=None)
    biospecimen_metadata_df["specimenID"] = biospecimen_metadata_df["specimenID"].astype(str)
    biospecimen_metadata_df["individualID"] = biospecimen_metadata_df["individualID"].astype(str)

    clinical_df = pd.read_csv("/groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/data/ROSMAP_clinical.csv", sep=",", header=0, index_col=None)
    clinical_df["individualID"] = clinical_df["individualID"].astype(str)

    # Merge together.
    df = idkey_df.merge(wgs_metadata_df, left_on="wgs_id", right_on="specimenID", how="left").merge(biospecimen_metadata_df, on="specimenID", how="left").merge(clinical_df, on="projid", how="left")
    df.index = df["wgs_id"]
    df.index.name = None

    # Reformat columns.
    df["sex"] = df["msex"].map({1.0: 1, 0.0: 2}, na_action='ignore')
    df["race"] = df["race"].map({1: "White",
                                 2: "Black, Negro, African-American",
                                 3: "Native American, Indian",
                                 4: "Eskimo",
                                 5: "Aleut",
                                 6: "Asian or Pacific Island",
                                 8: "NA",
                                 9: "NA"},
                                na_action='ignore')
    df["ageDeath"] = [90 if age == "90+" else float(age) for age in df["age_death"]]
    df["Study"] = df["Study"].fillna("ROSMAP")

    df = df.groupby(df.index).first()

    del wgs_metadata_df, biospecimen_metadata_df, clinical_df

    return df


def construct_mayo_df():
    # Load data.
    wgs_df = pd.read_csv("/groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/data/MayoRNAseq_assay_wholeGenomeSeq_metadata.csv", sep=",", header=0, index_col=None)
    wgs_df["specimenID"] = wgs_df["specimenID"].astype(str)
    wgs_df["wgs_available"] = "N"
    wgs_df.loc[wgs_df["assay"] == "wholeGenomeSeq", "wgs_available"] = "Y"

    biospecimen_df = pd.read_csv("/groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/data/MayoRNAseq_biospecimen_metadata.csv", sep=",", header=0, index_col=None)
    biospecimen_df["specimenID"] = biospecimen_df["specimenID"].astype(str)
    biospecimen_df["individualID"] = ["{:.0f}".format(x) for x in biospecimen_df["individualID"]]

    individual_df = pd.read_csv("/groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/data/MayoRNAseq_individual_metadata_031422.csv", sep=",", header=0, index_col=None)
    individual_df["individualID"] = individual_df["individualID"].astype(str)

    # Merge together.
    df = wgs_df.merge(biospecimen_df, on="specimenID", how="left").merge(individual_df, on="individualID", how="left")
    df.index = df["specimenID"]
    df.index.name = None

    # Reformat columns.
    df["sex"] = df["sex"].map({"male": 1, "female": 2}, na_action='ignore')
    df["ageDeath"] = [90 if age == "90_or_over" else float(age) for age in df["ageDeath"]]
    df["Study"] = "Mayo"

    df = df.groupby(df.index).first()

    del wgs_df, biospecimen_df, individual_df

    return df


def construct_msbb_df():
    # Load data.
    wgs_df = pd.read_csv("/groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/data/MSBB_assay_wholeGenomeSeq_metadata.csv", sep=",", header=0, index_col=None)
    wgs_df["specimenID"] = wgs_df["specimenID"].astype(str)
    wgs_df["wgs_available"] = "N"
    wgs_df.loc[wgs_df["assay"] == "wholeGenomeSeq", "wgs_available"] = "Y"

    biospecimen_df = pd.read_csv("/groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/data/MSBB_biospecimen_metadata.csv", sep=",", header=0, index_col=None)
    biospecimen_df["specimenID"] = biospecimen_df["specimenID"].astype(str)
    biospecimen_df["individualID"] = biospecimen_df["individualID"].astype(str)

    individual_df = pd.read_csv("/groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/data/MSBB_individual_metadata.csv", sep=",", header=0, index_col=None)
    individual_df["individualID"] = individual_df["individualID"].astype(str)

    # Merge together.
    df = wgs_df.merge(biospecimen_df, on="specimenID", how="left").merge(individual_df, on="individualID", how="left")
    df.index = df["specimenID"]
    df.index.name = None

    # Reformat columns.
    df["sex"] = df["sex"].map({"male": 1, "female": 2}, na_action='ignore')
    df["race"] = df["race"].map({"W": "White",
                                 "B": "Black or African American",
                                 "H": "Hispanic or Latino",
                                 "A": "American Native or Alaskan Native",
                                 "U": "NA"},
                                na_action='ignore')
    df["ageDeath"] = [90 if age == "90+" else float(age) for age in df["ageDeath"]]
    df["Study"] = "MSBB"

    df = df.groupby(df.index).first()

    del wgs_df, biospecimen_df, individual_df

    return df


def round_float(s, na):
    rounded_values = []
    for age in s:
        if np.isnan(age):
            rounded_values.append(na)
            continue

        if isinstance(age, float):
            rounded_values.append("{:.0f}".format(np.round(age)))
            continue

        rounded_values.append(age)

    return rounded_values


if __name__ == '__main__':
    start()