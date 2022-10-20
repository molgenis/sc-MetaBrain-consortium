#!/usr/bin/env python3

"""
File:         print_WGS_VCF_filter_overview.py
Created:      2022/10/22
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
import argparse
import gzip
import os
import re

# Third party imports.
import pandas as pd
import numpy as np

# Local application imports.

# Metadata
__program__ = "Print WGS VCF Filter Overview"
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

### AMP-AD ###
./print_WGS_VCF_filter_overview.py \
    --workdir /groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/2022-10-20-FilteredGenotypes \
    --vcf_file_format NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_CHR.recalibrated_variants-filtered.log.gz \
    --exclude X Y others
"""


CHROMOSOMES = [str(chr) for chr in range(1, 23)] + ["X", "Y", "others"]

class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.workdir = getattr(arguments, 'workdir')
        self.vcf_file_format = getattr(arguments, 'vcf_file_format')
        self.exclude = getattr(arguments, 'exclude')

        self.chromosomes = [CHR for CHR in CHROMOSOMES if (self.exclude is None) or (CHR not in self.exclude)]

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
        parser.add_argument("--workdir",
                            type=str,
                            required=True,
                            help="The path to the working directory")
        parser.add_argument("--vcf_file_format",
                            type=str,
                            required=True,
                            help="The file format of the VCF files. CHR"
                                 "will be replaced with the chromosome number.")
        parser.add_argument("--exclude",
                            nargs="*",
                            type=str,
                            choices=CHROMOSOMES,
                            default=None,
                            help="Exclude certain chromosomes.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Printing overview:")
        data = []
        for chr in self.chromosomes:
            row = [chr, False, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, np.nan, False, "00:00:00"]

            job_logfile_path = os.path.join(self.workdir, "jobs", "output", "VCFFilter_CHR{}.out".format(chr))
            if os.path.exists(job_logfile_path):
                row[1] = True
                row[11:] = self.read_job_logfile(filepath=job_logfile_path)

            job_logfile_path = os.path.join(self.workdir, self.vcf_file_format.replace("CHR", chr))
            if os.path.exists(job_logfile_path):
                row[2:11] = self.read_filter_logfile(filepath=job_logfile_path)

            data.append(row)

        data.append([])
        df = pd.DataFrame(data, columns=["Chromosome", "Started", "PASSQC", "FailedVariantStats", "MultiAllelic", "IndelBelowVQSR", "IndelNonPass", "SNVBelowVQSR", "SNVNonPass", "IncorrectInbreedingCoeff", "BelowInbreedingCoeff", "Parsed", "Written", "PctKept", "Finished", "Elapsed"])
        df.iloc[df.shape[0] - 1, :2] = ["total", df["started"].all()]
        df.iloc[df.shape[0] - 1, 2:13] = df.iloc[:, 2:13].sum(axis=0)
        df.iloc[df.shape[0] - 1, 13:] = [np.round((df.iloc[df.shape[0] - 1, 12] / df.iloc[df.shape[0] - 1, 11]) * 100, 1), df["finished"].all(), ""]
        df.iloc[:, 2:13] = df.iloc[:, 2:13].astype(int)
        print(df)

        self.save_file(
            df=df,
            outpath=os.path.join(self.workdir, "VCFFilterSummaryStats.txt.gz")
        )

    @staticmethod
    def read_job_logfile(filepath):
        with open(filepath, 'r') as f:
            parsed = 0
            written = 0
            pct_kept = np.nan
            finished = False
            elapsed = ""
            for line in f:
                if re.match("([0-9]+) lines parsed, ([0-9]+) written", line):
                    match = re.search("([0-9]+) lines parsed, ([0-9]+) written", line)
                    parsed = int(match.group(1))
                    written = int(match.group(2))
                    pct_kept = np.round((written / parsed) * 100, 1)
                elif re.match("Done. How about that!", line):
                    finished = True
                elif re.match("[0-9]{7}.batch. ([0-9]{2}:[0-9]{2}:[0-9]{2})", line):
                    match = re.search("[0-9]{7}.batch. ([0-9]{2}:[0-9]{2}:[0-9]{2})", line)
                    elapsed = match.group(1)
                else:
                    pass
        f.close()

        return [parsed, written, pct_kept, finished, elapsed]

    @staticmethod
    def read_filter_logfile(filepath):
        pass_qc = 0
        failed_variant_stats = 0
        multi_allelic = 0
        indel_below_vqsr = 0
        indel_non_pass = 0
        snv_below_vqsr = 0
        snv_non_pass = 0
        incorrect_inbreeding_coeff = 0
        below_inbreeding_coeff = 0

        try:
            with gzip.open(filepath, 'rt') as f:
                for line in f:
                    reason = line.split("\t")[1]
                    if reason == "PASSQC":
                        pass_qc += 1
                    elif reason == "FailedVariantStats":
                        failed_variant_stats += 1
                    elif reason == "MultiAllelic":
                        multi_allelic += 1
                    elif reason == "IndelBelowVQSR":
                        indel_below_vqsr += 1
                    elif reason == "IndelNonPass":
                        indel_non_pass += 1
                    elif reason == "SNVBelowVQSR":
                        snv_below_vqsr += 1
                    elif reason == "SNVNonPass":
                        snv_non_pass += 1
                    elif reason.startswith("IncorrectInbreedingCoeff"):
                        incorrect_inbreeding_coeff += 1
                    elif reason.startswith("BelowInbreedingCoeff"):
                        below_inbreeding_coeff += 1
                    else:
                        pass
            f.close()
        except EOFError:
            pass

        return [pass_qc, failed_variant_stats, multi_allelic, indel_below_vqsr, indel_non_pass, snv_below_vqsr, snv_non_pass, incorrect_inbreeding_coeff, below_inbreeding_coeff]

    @staticmethod
    def save_file(df, outpath, header=True, index=False, sep=","):
        df.to_csv(outpath, sep=sep, index=index, header=header)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    def print_arguments(self):
        print("Arguments:")
        print("  > Working directory: {}".format(self.workdir))
        print("  > vcf_file_format: {}".format(self.vcf_file_format))
        print("  > Exclude: {}".format(self.exclude))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()