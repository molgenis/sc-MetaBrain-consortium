#!/usr/bin/env python3

"""
File:         prepare_AMPAD_wg1_files.py
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
import glob
import re
import os

# Third party imports.
import numpy as np
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Prepare AMP-AD Workgroup 1 Files"
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
./prepare_AMPAD_wg1_files.py -h
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.metadata = getattr(arguments, 'metadata')
        self.fastqdir = getattr(arguments, 'fastqdir')
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
        parser.add_argument("--metadata",
                            type=str,
                            required=True,
                            help="")
        parser.add_argument("--fastqdir",
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

        print("Loading AMP-AD biospecimen metadata")
        metadata_df = self.load_file(inpath=os.path.join(self.metadata, "ROSMAP_biospecimen_metadata.csv"))
        # print(metadata_df)

        # scRNA-seq metadata.
        scrnaseq_metadata = metadata_df.loc[metadata_df["assay"] == "scrnaSeq", ["individualID", "specimenID"]].copy()
        scrnaseq_metadata.drop_duplicates(inplace=True)
        scrnaseq_metadata.columns = ["individualID", "scrnaSeqID"]
        print(scrnaseq_metadata)

        # WGS metadata.
        wholegenomeseq_metadata = metadata_df.loc[metadata_df["assay"] == "wholeGenomeSeq", ["individualID", "specimenID"]].copy()
        wholegenomeseq_metadata.drop_duplicates(inplace=True)
        wholegenomeseq_metadata.columns = ["individualID", "wholeGenomeSeqID"]
        print(wholegenomeseq_metadata)

        print("Loading fastQ samples")
        fastq_df = self.load_fastq_files(inpath=self.fastqdir)
        # print(fastq_df)

        scrnaseq_id = "scrnaSeqID"
        individual_id = "individualID"
        project_id = "individualID"
        genotype_id = "wholeGenomeSeqID"
        if self.name == "Mathys2019":
            scrnaseq_id = "specimenID"
            project_id = "projid"
            fastq_df[scrnaseq_id] = fastq_df["id"] + "_" + fastq_df["sample"]
        elif self.name == "Cain2023":
            fastq_df[scrnaseq_id] = fastq_df["id"]
        elif self.name == "Zhou2020":
            scrnaseq_id = "specimenID"
            fastq_df[scrnaseq_id] = [id.split("-")[1].split("_")[0] for id in fastq_df["id"]]

        fastq_df = fastq_df[[scrnaseq_id]].drop_duplicates()
        fastq_df.reset_index(drop=True, inplace=True)
        print(fastq_df)

        if self.dataset is None:
            fastq_df = fastq_df.merge(scrnaseq_metadata, on=scrnaseq_id, how="left")
        else:
            print("Loading dataset biospecimen metadata")
            dataset_df = self.load_file(inpath=self.dataset)
            dataset_df = dataset_df[list({scrnaseq_id, individual_id, project_id})]
            fastq_df = fastq_df.merge(dataset_df, on=scrnaseq_id, how="left")
        print(fastq_df)

        vcf_samples = None
        if self.vcf is not None:
            print("Checking overlap with VCF")
            vcf_samples = self.load_vcf_samples(inpath=self.vcf)

        # First whole genome seq ids.
        wholegenomeseq_dict = {}
        for _, row in fastq_df.iterrows():
            if row[individual_id] in wholegenomeseq_dict:
                continue
            wholegenomeseq_ids = list(wholegenomeseq_metadata.loc[wholegenomeseq_metadata[individual_id] == row[individual_id], genotype_id].values)
            if len(wholegenomeseq_ids) == 0:
                wholegenomeseq_dict[row[individual_id]] = "NA"
            elif len(wholegenomeseq_ids) == 1:
                wholegenomeseq_dict[row[individual_id]] = wholegenomeseq_ids[0]
            else:
                if vcf_samples is None:
                    print("Error, unable to match {} uniquely to a wholeGenomeSeqID. Options are {}.".format(row[individual_id], ", ".join(wholegenomeseq_ids)))
                    exit()
                else:
                    found_wholegenomeseq_ids = [wholegenomeseq_id for wholegenomeseq_id in wholegenomeseq_ids if wholegenomeseq_id in vcf_samples]
                    if len(found_wholegenomeseq_ids) == 0:
                        wholegenomeseq_dict[row[individual_id]] = "NA"
                    elif len(found_wholegenomeseq_ids) == 1:
                        wholegenomeseq_dict[row[individual_id]] = found_wholegenomeseq_ids[0]
                    else:
                        print("Error, unable to match {} uniquely to a wholeGenomeSeqID. Options are {}.".format(row[individual_id], ", ".join(wholegenomeseq_ids)))
                        exit()

        print(wholegenomeseq_dict)
        df = fastq_df.copy()
        fastq_df["wholeGenomeSeqID"] = fastq_df[individual_id].map(wholegenomeseq_dict)
        if "projid" in fastq_df.columns and vcf_samples is not None:
            for index, row in fastq_df.iterrows():
                if row[genotype_id] not in vcf_samples and "ROS" + str(row[project_id]) in vcf_samples:
                    fastq_df.loc[index, genotype_id] = "ROS" + str(row[project_id])
        del fastq_df
        print(df)

        if self.vcf is not None:
            print("Checking overlap with VCF")
            vcf_samples = self.load_vcf_samples(inpath=self.vcf)
            df["FoundInVCF"] = [sample in vcf_samples for sample in df[genotype_id]]

            mask = []
            for index, row in df.iterrows():
                found = True
                if str(row[genotype_id]) != "nan" and not row["FoundInVCF"]:
                    print("  Error, {} not found in the VCF file.".format(row[genotype_id]))
                    found = False

                    if "projid" in df.columns:
                        alt_individual_id = "ROS" + str(row[project_id])
                        if alt_individual_id in vcf_samples:
                            print("    Warning, alternative ID '{}' was found, replacing value. Use with caution.".format(alt_individual_id))
                            df.loc[index, genotype_id] = alt_individual_id
                            df.loc[index, "FoundInVCF"] = True
                            found = True
                mask.append(found)

            print("  Removed {:,} samples due to missing in VCF file.".format(len(mask) - sum(mask)))
            df.loc[mask, genotype_id] = np.nan

        print("Merged data:")
        print(df)
        self.save_file(df=df,
                       outpath=os.path.join(self.outdir, self.name, "{}_full_link_table.tsv".format(self.name)),
                       sep="\t")

        print("Link table:")
        link_table = df[[project_id, scrnaseq_id]].copy()
        link_table.columns = ["fastqID", "specimenID"]
        print(link_table)
        self.save_file(df=link_table,
                       outpath=os.path.join(self.outdir, self.name, "{}_link_table.csv".format(self.name)),
                       sep=",")

        print("Individual coupling:")
        ind_coupling = df[[project_id, genotype_id]].copy()
        ind_coupling.columns = ["Pool", "Assignment"]
        ind_coupling.dropna(inplace=True)
        print(ind_coupling)
        self.save_file(df=ind_coupling,
                       outpath=os.path.join(self.outdir, self.name, "{}_individual_coupling.tsv".format(self.name)))
        self.save_file(df=ind_coupling[["Assignment"]],
                       header=False,
                       outpath=os.path.join(self.outdir, self.name, "{}_genotype_samples.txt".format(self.name)))

        print("GTE:")
        gte = df[[genotype_id, project_id]].copy()
        print(gte)
        self.save_file(df=gte,
                       outpath=os.path.join(self.outdir, self.name, "{}_GTE.tsv".format(self.name)),
                       header=False)

    @staticmethod
    def load_file(inpath, header=0, index_col=None, sep=",", skiprows=None,
                  nrows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         skiprows=skiprows, nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def load_fastq_files(inpath):
        fpaths = glob.glob(inpath + "*.fastq.gz")

        samples = []
        for fpath in fpaths:
            match = re.match("(.+)_(S[0-9]+)_(L00[0-9])_(R[12])_001\.fastq\.gz$", os.path.basename(fpath))
            if match is None:
                continue
            samples.append([match.group(1), match.group(2), match.group(3), match.group(4)])

        return pd.DataFrame(samples, columns=["id", "sample", "lane", "read_type"])

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
        print("  > Metadata directory: {}".format(self.metadata))
        print("  > FastQ directory: {}".format(self.fastqdir))
        print("  > Dataset metadata file: {}".format(self.dataset))
        print("  > Name: {}".format(self.name))
        print("  > VCF: {}".format(self.vcf))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()