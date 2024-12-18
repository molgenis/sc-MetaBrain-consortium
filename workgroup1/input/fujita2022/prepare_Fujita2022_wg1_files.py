#!/usr/bin/env python3

"""
File:         prepare_Fujita2022_wg1_files.py
Created:      2023/02/27
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
__program__ = "Prepare Fujita2022 Workgroup 1 Files"
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
./prepare_Fujita2022_wg1_files.py -h
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.metadata = getattr(arguments, 'metadata')
        self.fastqdir = getattr(arguments, 'fastqdir')
        self.cell_annotation = getattr(arguments, 'cell_annotation')
        self.name = "Fujita2022"
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
        parser.add_argument("--cell_annotation",
                            type=str,
                            default=None,
                            required=False,
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
        print(metadata_df)

        # scRNA-seq metadata.
        scrnaseq_metadata = metadata_df.loc[metadata_df["assay"] == "scrnaSeq", ["individualID", "specimenID"]].copy()
        scrnaseq_metadata.drop_duplicates(inplace=True)
        scrnaseq_metadata.columns = ["individualID", "scrnaSeqID"]
        scrnaseq_metadata["PartialscrnaSeqID"] = [scrnaseq_id.split("_")[0] for scrnaseq_id in scrnaseq_metadata["scrnaSeqID"]]
        print(scrnaseq_metadata)
        print(len(scrnaseq_metadata["individualID"].unique()))
        print(scrnaseq_metadata.loc[scrnaseq_metadata["individualID"].duplicated(), :])

        # WGS metadata.
        wholegenomeseq_metadata = metadata_df.loc[metadata_df["assay"] == "wholeGenomeSeq", ["individualID", "specimenID", "organ", "tissue"]].copy()
        wholegenomeseq_metadata.columns = ["individualID", "wholeGenomeSeqID", "organ", "tissue"]
        print(wholegenomeseq_metadata)

        vcf_samples = None
        if self.vcf is not None:
            print("Checking overlap with VCF")
            vcf_samples = self.load_vcf_samples(inpath=self.vcf)
            wholegenomeseq_metadata = wholegenomeseq_metadata.loc[wholegenomeseq_metadata["wholeGenomeSeqID"].isin(vcf_samples), :]

        dups = set(wholegenomeseq_metadata.loc[wholegenomeseq_metadata["individualID"].duplicated(), "individualID"].values)
        clean_df = wholegenomeseq_metadata.loc[~wholegenomeseq_metadata["individualID"].isin(dups), :].copy()
        for dup_ind in dups:
            subset = wholegenomeseq_metadata.loc[wholegenomeseq_metadata["individualID"] == dup_ind, :]
            if sum(subset["organ"] == "brain") == 1:
                clean_df = pd.concat([clean_df, subset.loc[subset["organ"] == "brain", :]])
            elif sum(subset["tissue"] == "dorsolateral prefrontal cortex") == 1:
                clean_df = pd.concat([clean_df, subset.loc[subset["tissue"] == "dorsolateral prefrontal cortex", :]])
            else:
                print("Warning, could not resolve duplicates:")
                print(subset)
                clean_df = pd.concat([clean_df, subset])
        wholegenomeseq_metadata = clean_df
        del clean_df

        print(wholegenomeseq_metadata)
        # print(len(wholegenomeseq_metadata["individualID"].unique()))
        print(wholegenomeseq_metadata.loc[wholegenomeseq_metadata["individualID"].duplicated(), :])

        metadata = scrnaseq_metadata.merge(wholegenomeseq_metadata, how="left")
        print(metadata)
        del scrnaseq_metadata, wholegenomeseq_metadata

        print("Loading fastQ samples")
        fastq_df = self.load_fastq_files(inpath=self.fastqdir)
        fastq_df["PartialscrnaSeqID"] = [id_value.split("_")[0] for id_value in fastq_df["id"]]
        print(fastq_df)

        individuals_list_path = os.path.join(self.outdir, "Fujita2022", "individuals_list_dir")
        if not os.path.exists(individuals_list_path):
            os.makedirs(individuals_list_path)

        print("Loading cell annotations")
        individual_list_dict, pool_n_df = self.parse_fujita_cell_annotation(inpath=self.cell_annotation)

        print("Creating individual list files")
        pools_df = fastq_df[["id", "PartialscrnaSeqID"]].drop_duplicates().copy()
        print(pools_df)
        subset_list = []
        samplesheet = []
        for _, row in pools_df.iterrows():
            print("  Pool: {}".format(row["id"]))
            subset = metadata.loc[metadata["PartialscrnaSeqID"] == row["PartialscrnaSeqID"], :].copy()
            subset["id"] = row["id"]
            n = subset.shape[0]
            subset["N"] = n

            subset["cell_annot_match"] = False
            if not row["PartialscrnaSeqID"] in individual_list_dict:
                print("\tWarning, could not verify if individualIDs match")
            else:
                cell_annot_ind_ids = individual_list_dict[row["PartialscrnaSeqID"]]
                subset["cell_annot_match"] = [ind_id in cell_annot_ind_ids for ind_id in subset["individualID"]]

                n_matches = subset["cell_annot_match"].sum()
                if n_matches == n:
                    print("\tall individualIDs match")
                elif n_matches == len(cell_annot_ind_ids):
                    print("\tall overlapping individualIDs match")
                else:
                    print("\t{} / {} individualIDs match".format(n_matches, n))
                    print("\t\tMissing: {}".format(", ".join([cell_annot_ind_id for cell_annot_ind_id in cell_annot_ind_ids if cell_annot_ind_id not in subset["individualID"]])))
                    print(subset)

            interest = ["SM-CTED9", "SM-CJEHE", "SM-CJGGL", "SM-CJK4Y", "SM-CTEMN", "SM-CJGLY", "SM-CTDVR"]
            for bla in interest:
                if bla in subset["wholeGenomeSeqID"]:
                    print(subset)
                    exit()

            subset_list.append(subset)
            samplesheet.append([row["id"], n])

            self.save_file(pd.DataFrame(subset.loc[~subset["wholeGenomeSeqID"].isna(), "wholeGenomeSeqID"]), outpath=os.path.join(individuals_list_path, row["id"] + ".txt"), header=False, index=False)

        df = pd.concat(subset_list, axis=0)

        print("Merged data:")
        print(df)
        self.save_file(df=df,
                       outpath=os.path.join(self.outdir, self.name, "{}_full_link_table.tsv".format(self.name)),
                       sep="\t")

        samplesheet_df = pd.DataFrame(samplesheet, columns=["Pool", "N_Individuals"])
        print(samplesheet_df)
        self.save_file(df=samplesheet_df,
                       outpath=os.path.join(self.outdir, self.name, "{}_samplesheet.tsv".format(self.name)),
                       sep="\t")

        print("Link table:")
        link_table = df[["id", "id"]].copy()
        link_table.columns = ["fastqID", "specimenID"]
        print(link_table)
        self.save_file(df=link_table,
                       outpath=os.path.join(self.outdir, self.name, "{}_link_table.csv".format(self.name)),
                       sep=",")

        print("Genotype samples:")
        genotype_samples = df[["wholeGenomeSeqID"]].copy()
        genotype_samples.drop_duplicates(inplace=True)
        genotype_samples[genotype_samples == ""] = np.nan
        genotype_samples = genotype_samples.dropna()
        print(genotype_samples)
        self.save_file(df=genotype_samples,
                       header=False,
                       outpath=os.path.join(self.outdir, self.name, "{}_genotype_samples.txt".format(self.name)))

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
        if len(fpaths) == 0 and os.path.exists(inpath + "SYNAPSE_METADATA_MANIFEST.tsv"):
            fpaths = pd.read_csv(inpath + "SYNAPSE_METADATA_MANIFEST.tsv", sep="\t", header=0, index_col=None)["name"].values

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
    def parse_fujita_cell_annotation(inpath):
        individual_list_dict = {}
        with open(inpath, 'r') as f:
            for i, line in enumerate(f):
                if i == 0:
                    continue
                values = line.rstrip("\n").split(",")
                partial_scrna_seq_id = values[0].split("_")[0]
                individual_id = values[1]
                if partial_scrna_seq_id not in individual_list_dict:
                    individual_list_dict[partial_scrna_seq_id] = {individual_id}
                else:
                    individual_list_dict[partial_scrna_seq_id].add(individual_id)
        f.close()

        data = []
        for key, value in individual_list_dict.items():
            data.append([key, len(value)])

        return individual_list_dict, pd.DataFrame(data, columns=["PartialscrnaSeqID", "N"])

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
        print("  > Cell annotation file: {}".format(self.cell_annotation))
        print("  > Name: {}".format(self.name))
        print("  > VCF: {}".format(self.vcf))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()