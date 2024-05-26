#!/usr/bin/env python3

"""
File:         parse_ega_metadata.py
Created:      2023/07/13
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
__program__ = "Parse EGA Metadata"
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
./parse_ega_metadata.py -h
"""

class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.input_dir = getattr(arguments, 'input_dir')
        self.vcf_path = getattr(arguments, 'vcf')
        self.output_dir = getattr(arguments, 'output_dir')

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
                            help="show program's version number and exit.")
        parser.add_argument("--input_dir",
                            type=str,
                            required=True,
                            help="")
        parser.add_argument("--vcf",
                            type=str,
                            required=False,
                            default="",
                            help="")
        parser.add_argument("--output_dir",
                            type=str,
                            required=True,
                            help="")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        sample_file_path = os.path.join(self.input_dir, "metadata", "delimited_maps", "Sample_File.map")
        ms_samples = []
        if os.path.exists(sample_file_path):
            sample_file_df = self.load_file(sample_file_path)
            sample_file_df.columns = ["Ind_Sample", "individualID", "fastqPath", "sampleID"]
            print(sample_file_df)
            sample_file_df = sample_file_df.loc[["fastq" in path for path in sample_file_df["fastqPath"]], :]
            sample_file_df["fastqID"] = [x.split("/")[1].split("_S1")[0] for x in sample_file_df["fastqPath"]]
            sample_file_df = sample_file_df[["fastqID", "individualID"]].drop_duplicates()
            if len(sample_file_df["individualID"].unique()) != sample_file_df.shape[0]:
                print("Error, not all individualID's are unique")
                exit()
            sample_file_df.sort_values(by="individualID", inplace=True)
            print(sample_file_df)
            ms_samples = [x.split("-")[0] for x in sample_file_df["fastqID"]]
            # for index, row in sample_file_df.iterrows():
            #     print("{},{}".format(row[0], row[1]))
            # self.save_file(df=sample_file_df, outpath=os.path.join(self.output_dir, "link_table.csv"))

        sample_meta_info_path = os.path.join(self.input_dir, "metadata", "delimited_maps", "Analysis_Sample_meta_info.map")
        if os.path.exists(sample_meta_info_path) and os.path.exists(self.vcf_path):
            sample_meta_info_df = self.load_file(sample_meta_info_path)
            info_df = self.split_info(sample_meta_info_df[2])
            del sample_meta_info_df

            info_df = info_df.groupby(info_df["subject_id"]).first()
            info_df.reset_index(drop=False, inplace=True)
            print(info_df)

            vcf_df = self.load_file(self.vcf_path, header=0, skiprows=31, nrows=2)
            vcf_sampels = [col for col in vcf_df.columns if col not in ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]]

            info_df = info_df[["subject_id", "gender", "region"]].copy()
            info_df.columns = ["IID", "SEX", "Provided_Ancestry"]
            info_df.insert(0, "#FID", info_df["IID"])
            info_df.insert(2, "PAT", 0)
            info_df.insert(3, "MAT", 0)
            info_df["SEX"] = info_df["SEX"].map({"male": "1", "female": "2"})
            info_df["Provided_Ancestry"] = info_df["Provided_Ancestry"].map({"Europe": "EUR"})
            info_df["genotyping_platform"] = "GSAv3IlluminaChIP"
            info_df["array_available"] = "Y"
            info_df["wgs_available"] = "N"
            info_df["wes_available"] = "N"
            info_df["age"] = "NA"
            info_df["age_range"] = "NA"
            info_df["Study"] = "RocheAD2022"
            info_df.loc[info_df["IID"].isin(ms_samples), "Study"] = "RocheMS2022"
            info_df["smoking_status"] = "NA"
            info_df["hormonal_contraception_use_currently"] = "NA"
            info_df["menopause"] = "NA"
            info_df["pregnancy_status"] = "NA"

            info_df.index = info_df["IID"]
            info_df = info_df.loc[vcf_sampels, :]

            print(info_df)
            self.save_file(df=info_df, outpath=os.path.join(self.output_dir, "Roche.psam"), sep="\t")

    @staticmethod
    def split_info(s):
        data = []
        for _, value in s.items():
            fields = value.split(";")

            row_data = {}
            for field in fields:
                splitted_fields = field.split("=")
                if len(splitted_fields) == 2:
                    try:
                        row_data[splitted_fields[0]] = float(splitted_fields[1])
                    except ValueError:
                        row_data[splitted_fields[0]] = splitted_fields[1]
                else:
                    row_data[splitted_fields[0]] = ""

            data.append(pd.Series(row_data))

        return pd.DataFrame(data)

    @staticmethod
    def load_file(inpath, header=None, index_col=None, sep="\t", skiprows=None,
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
        print("  > Input directory: {}".format(self.input_dir))
        print("  > VCF path: {}".format(self.vcf_path))
        print("  > Output directory: {}".format(self.output_dir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
