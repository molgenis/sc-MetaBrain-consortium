#!/usr/bin/env python3

"""
File:         filter_WGS_VCF_files.py
Created:      2022/10/19
Last Changed: 2022/10/20
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
from datetime import datetime
import subprocess
import argparse
import os
import re

# Third party imports.

# Local application imports.

# Metadata
__program__ = "Filter WGS VCF files"
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
./filter_WGS_VCF_files.py \
    --filter_script /groups/umcg-biogen/tmp01/input/rawdata/2017-GTExV8Genotypes/customVCFFilter.py \
    --vcf_indir /groups/umcg-biogen/tmp01/input/AMP-AD/2017-12-08-joint-WGS \
    --vcf_file_format NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_CHR.recalibrated_variants.vcf.gz \
    --outdir /groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD \
    --outfolder FilteredGenotypes
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.filter_script = getattr(arguments, 'filter_script')
        self.vcf_indir = getattr(arguments, 'vcf_indir')
        self.vcf_file_format = getattr(arguments, 'vcf_file_format')
        time = getattr(arguments, 'time')
        self.cpus_per_task = getattr(arguments, 'cpus_per_task')
        self.mem = getattr(arguments, 'mem')
        outdir = getattr(arguments, 'outdir')
        outfolder = getattr(arguments, 'outfolder')

        if not self.vcf_file_format.endswith(".vcf.gz"):
            print("Error: --vcf_file_format must end with '.vcf.gz'.")
            exit()
        if self.vcf_file_format.count('CHR') > 1:
            print("Error: --vcf_file_format must only contain 'CHR' once.")
            exit()

        time_dict = {
            "short": "05:59:59",
            "medium": "23:59:00",
            "long": "6-23:59:00"
        }
        self.time = time_dict[time]

        # Pre-process the dataset output directory.
        date_str = datetime.now().strftime("%Y-%m-%d")
        if not re.match("^[0-9]{4}-(0[1-9]|1[0-2])-(0[1-9]|[1-2][0-9]|3[0-1])", outfolder):
            outfolder = "{}-{}".format(date_str, outfolder)
        self.outdir = os.path.join(outdir, outfolder)
        self.jobdir = os.path.join(self.outdir, "jobs")
        self.jobs_outdir = os.path.join(self.jobdir, "output")
        for dir in [self.outdir, self.jobdir, self.jobs_outdir]:
            if not os.path.exists(dir):
                os.makedirs(dir)

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
        parser.add_argument("--filter_script",
                            type=str,
                            required=True,
                            help="The path to the customVCFFilter.py of "
                                 "Harm-Jan Westra.")
        parser.add_argument("--vcf_indir",
                            type=str,
                            required=True,
                            help="The directory containing the input VCF "
                                 "files.")
        parser.add_argument("--vcf_file_format",
                            type=str,
                            required=True,
                            help="The file format of the VCF files. CHR"
                                 "will be replaced with the chromosome number.")
        parser.add_argument("--time",
                            type=str,
                            default="medium",
                            choices=["short", "medium", "long"],
                            help="Sets the work allocation time a.k.a. "
                                 "walltime.")
        parser.add_argument("--cpus_per_task",
                            type=int,
                            default=1,
                            help="Requests X CPUs (cores) for your job.")
        parser.add_argument("--mem",
                            type=int,
                            default=1,
                            help="Requests X GB RAM total per job.")
        parser.add_argument("--outdir",
                            type=str,
                            required=True,
                            help="The path for the output directory.")
        parser.add_argument("--outfolder",
                            type=str,
                            required=True,
                            help="The name of the output directory.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Look for available files")
        files = []
        if "CHR" in self.vcf_file_format:
            for chr in [str(chr) for chr in range(1, 23)] + ["X", "Y", "others"]:
                files.append((chr, os.path.join(self.vcf_indir, self.vcf_file_format.replace("CHR", chr))))
        else:
            files = [("", os.path.join(self.vcf_indir, self.vcf_file_format))]
        files = [(chr, file) for chr, file in files if os.path.exists(file)]
        if len(files) > 0:
            for (chr, file) in files:
                print("  > CHR{} = {}".format(chr, os.path.basename(file)))
        else:
            print("  No files found.")
            exit()

        print("Generating job scripts")
        job_files = []
        for (chr, file) in files:
            joblog_path, jobfile_path = self.create_job_file(
                chr=chr,
                infile=file,
                outfile=os.path.join(self.outdir, os.path.basename(file).replace(".vcf.gz", ""))
            )
            job_files.append((joblog_path, jobfile_path))

        print("Starting job files.")
        for (joblog_path, jobfile_path) in job_files:
            completed = False
            if os.path.exists(joblog_path):
                for line in open(joblog_path, 'r'):
                    if "Done. How about that!" in line:
                        completed = True
                print("\tSkipping '{}'".format(os.path.basename(jobfile_path)))

            if not completed:
                self.run_command(command=['sbatch', jobfile_path])

    def create_job_file(self, chr, infile, outfile):
        job_name = "VCFFilter_CHR{}".format(chr)
        joblog_path = os.path.join(self.jobs_outdir, job_name + ".out")

        lines = ["#!/bin/bash",
                 "#SBATCH --job-name={}".format(job_name),
                 "#SBATCH --output={}".format(joblog_path),
                 "#SBATCH --error={}".format(joblog_path),
                 "#SBATCH --time={}".format(self.time),
                 "#SBATCH --cpus-per-task={}".format(self.cpus_per_task),
                 "#SBATCH --mem={}gb".format(self.mem),
                 "#SBATCH --nodes=1",
                 "#SBATCH --open-mode=append",
                 "#SBATCH --export=NONE",
                 "#SBATCH --get-user-env=L",
                 "",
                 "module load Python/3.7.4-GCCcore-7.3.0-bare",
                 "",
                 "python3 {} {} {}".format(self.filter_script, infile, outfile)]

        jobfile_path = os.path.join(self.jobdir, job_name + ".sh")
        with open(jobfile_path, "w") as f:
            for line in lines:
                f.write(line + "\n")
        f.close()
        print("\tSaved jobfile: {}".format(os.path.basename(jobfile_path)))
        return joblog_path, jobfile_path

    @staticmethod
    def run_command(command):
        print("\t" + " ".join(command))
        subprocess.call(command)

    @staticmethod
    def save_file(df, outpath, header=True, index=False, sep=","):
        df.to_csv(outpath, sep=sep, index=index, header=header)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    def print_arguments(self):
        print("Arguments:")
        print("  > filter_script: {}".format(self.filter_script))
        print("  > vcf_indir: {}".format(self.vcf_indir))
        print("  > vcf_file_format: {}".format(self.vcf_file_format))
        print("  > Time: {}".format(self.time))
        print("  > CPU-per-task: {}Gb".format(self.cpus_per_task))
        print("  > Memory: {}Gb".format(self.mem))
        print("  > outdir: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()