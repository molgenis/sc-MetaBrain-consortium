#!/usr/bin/env python3

"""
File:         prepare_roche_wg1_files.py
Created:      2024/03/06
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
import gzip
import re
import os

# Third party imports.
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Prepare Roche Workgroup 1 Files"
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
./prepare_roche_wg1_files.py -h
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.indir = getattr(arguments, 'indir')
        self.dataset = getattr(arguments, 'dataset')
        self.name = getattr(arguments, 'name')
        self.vcf = getattr(arguments, 'vcf')
        self.outdir = getattr(arguments, 'outdir')

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add other arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit")
        parser.add_argument("--indir",
                            type=str,
                            required=True,
                            help="")
        parser.add_argument("--dataset",
                            type=str,
                            default=None,
                            required=False,
                            help="")
        parser.add_argument("--name",
                            type=str,
                            required=True,
                            help="")
        parser.add_argument("--vcf",
                            type=str,
                            default=None,
                            required=False,
                            help="")
        parser.add_argument("--outdir",
                            type=str,
                            required=True,
                            help="The output directory.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading 'Sample_File.map'")
        sample_file_path = os.path.join(self.indir, "metadata", "delimited_maps", "Sample_File.map")
        print(sample_file_path)
        if not os.path.exists(sample_file_path):
            print("\tError, file does not exist.")
            exit()

        sample_file_df = self.load_file(inpath=sample_file_path, header=None)
        sample_file_df.columns = ["Ind_Sample", "SampleID", "fastqPath", "fileID"]
        sample_file_df.insert(0, "individual_id_anon", ["Ind" + id_value.split("_")[0] for id_value in sample_file_df["Ind_Sample"]])
        sample_file_df["fastqPath"] = [os.path.basename(fpath).rstrip(".cip") for fpath in sample_file_df["fastqPath"]]
        print(sample_file_df)

        print("Loading fastQ samples")
        fastq_df = self.load_fastq_files(inpath=self.indir)
        print(fastq_df)

        df = sample_file_df.merge(fastq_df, on="fastqPath", how="right")
        print(df)

        vcf_id_col = "individual_id_anon"
        if self.dataset is not None:
            print("Loading dataset biospecimen metadata")
            dataset_df = self.load_file(inpath=self.dataset)
            df = dataset_df.merge(df, how="right")
            df = df[["individual_id_anon", "individual_id", "SampleID", "sample_id_anon"]].drop_duplicates()
            vcf_id_col = "individual_id"
        else:
            df = df[["individual_id_anon", "SampleID", "sample_id_anon"]].drop_duplicates()
        print(df)

        if self.vcf is not None:
            print("Checking overlap with VCF")
            vcf_samples = self.load_vcf_samples(inpath=self.vcf)
            print(vcf_samples)
            df["FoundInVCF"] = [sample in vcf_samples for sample in df[vcf_id_col]]

            mask = []
            for index, row in df.iterrows():
                found = True
                if str(row[vcf_id_col]) != "nan" and not row["FoundInVCF"]:
                    print("  Error, {} not found in the VCF file.".format(row[vcf_id_col]))
                    found = False
                mask.append(found)

            # print("  Removed {:,} samples due to missing in VCF file.".format(len(mask) - sum(mask)))
            # df = df.loc[mask, : ]

        print("Merged data:")
        print(df)
        self.save_file(df=df,
                       outpath=os.path.join(self.outdir, self.name, "{}_full_link_table.tsv".format(self.name)),
                       sep="\t")

        print("Link table:")
        link_table = df[["sample_id_anon", "SampleID"]].copy()
        link_table.columns = ["fastqID", "SampleID"]
        print(link_table)
        self.save_file(df=link_table,
                       outpath=os.path.join(self.outdir, self.name, "{}_link_table.csv".format(self.name)),
                       sep=",")

        print("Individual coupling:")
        ind_coupling = df[["SampleID", vcf_id_col]].copy()
        ind_coupling.columns = ["Pool", "Assignment"]
        ind_coupling.dropna(inplace=True)
        print(ind_coupling)
        self.save_file(df=ind_coupling,
                       outpath=os.path.join(self.outdir, self.name, "{}_individual_coupling.tsv".format(self.name)))
        self.save_file(df=ind_coupling[["Assignment"]],
                       header=False,
                       outpath=os.path.join(self.outdir, self.name, "{}_genotype_samples.txt".format(self.name)))

        print("GTE:")
        gte = df[[vcf_id_col, "SampleID"]].copy()
        print(gte)
        self.save_file(df=gte,
                       outpath=os.path.join(self.outdir, self.name, "{}_GTE.tsv".format(self.name)),
                       header=False)


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
    def load_fastq_files(inpath):
        samples = []
        for root, dirs, files in os.walk(inpath):
            for file in files:
                match = re.match("(.+)_(S[0-9]+)_(L00[0-9])_(R[12])_001\.fastq\.gz$", file)
                if match is None:
                    continue
                samples.append([file, match.group(1), match.group(2), match.group(3), match.group(4)])

        return pd.DataFrame(samples, columns=["fastqPath", "sample_id_anon", "sample", "lane", "read_type"])

    @staticmethod
    def load_vcf_samples(inpath):
        header = None
        with gzip.open(inpath, 'rt') as f:
            for line in f:
                if line.startswith("##"):
                    continue
                header = line.strip("\n").split("\t")
                break
        f.close()
        return set([col for col in header if col not in ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]])

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

    def print_arguments(self):
        print("Arguments:")
        print("  > Input directory: {}".format(self.indir))
        print("  > Dataset metadata file: {}".format(self.dataset))
        print("  > Name: {}".format(self.name))
        print("  > VCF: {}".format(self.vcf))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()