#!/usr/bin/env python3

"""
File:         filter_WGS_VCF_files.py
Created:      2022/10/19
Last Changed: 2022/10/26
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
    --filter_script /groups/umcg-biogen/tmp01/input/rawdata/2017-GTExV8Genotypes/customVCFFilterV3.py \
    --vcf_indir /groups/umcg-biogen/tmp01/input/AMP-AD/2017-12-08-joint-WGS \
    --vcf_file_format NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_CHR.recalibrated_variants.vcf.gz \
    --exclude X Y others \
    --annotations /groups/umcg-biogen/tmp01/annotation/GenomeReference/dbsnp/b37/All_b37_b151_20180418.vcf.gz \
    --outdir /groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD \
    --outfolder FilteredGenotypes \
    --dryrun
"""

CHROMOSOMES = [str(chr) for chr in range(1, 23)] + ["X", "Y", "others"]
TIME_DICT = {
    "short": "05:59:59",
    "medium": "23:59:00",
    "long": "6-23:59:00"
}

class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.filter_script = getattr(arguments, 'filter_script')
        self.vcf_indir = getattr(arguments, 'vcf_indir')
        self.vcf_file_format = getattr(arguments, 'vcf_file_format')
        self.exclude = getattr(arguments, 'exclude')
        self.annotations = getattr(arguments, 'annotations')
        self.time_per_chr = TIME_DICT[getattr(arguments, 'time_per_chr')]
        self.cpus_per_chr = getattr(arguments, 'cpus_per_chr')
        self.mem_per_chr = getattr(arguments, 'mem_per_chr')
        self.time_merge = TIME_DICT[getattr(arguments, 'time_merge')]
        self.cpus_merge = getattr(arguments, 'cpus_merge')
        self.mem_merge = getattr(arguments, 'mem_merge')
        self.naive_merge = getattr(arguments, 'naive_merge')
        outdir = getattr(arguments, 'outdir')
        outfolder = getattr(arguments, 'outfolder')
        self.dryrun = getattr(arguments, 'dryrun')

        if not self.vcf_file_format.endswith(".vcf.gz"):
            print("Error: --vcf_file_format must end with '.vcf.gz'.")
            exit()
        if self.vcf_file_format.count('CHR') > 1:
            print("Error: --vcf_file_format must only contain 'CHR' once.")
            exit()

        self.chromosomes = [CHR for CHR in CHROMOSOMES if (self.exclude is None) or (CHR not in self.exclude)]

        # Pre-process the dataset output directory.
        date_str = datetime.now().strftime("%Y-%m-%d")
        if not re.match("^[0-9]{4}-(0[1-9]|1[0-2])-(0[1-9]|[1-2][0-9]|3[0-1])", outfolder):
            outfolder = "{}-{}".format(date_str, outfolder)
        self.outdir = os.path.join(outdir, outfolder)
        self.norm_outdir = os.path.join(self.outdir, "1-bcftools_norm")
        self.filtered_outdir = os.path.join(self.outdir, "2-custom_vcffilter")
        self.bgzip_outdir = os.path.join(self.outdir, "3-bgzip")
        self.annotated_outdir = os.path.join(self.outdir, "4-bcftools_annotate")
        self.merged_outdir = os.path.join(self.outdir, "5-bcftools_concat")
        self.plink_outdir = os.path.join(self.outdir, "6-plink2_makepgen")
        self.jobdir = os.path.join(self.outdir, "jobs")
        self.jobs_outdir = os.path.join(self.jobdir, "output")
        for dir in [self.outdir, self.norm_outdir, self.filtered_outdir,
                    self.bgzip_outdir, self.annotated_outdir,
                    self.merged_outdir, self.plink_outdir, self.jobdir,
                    self.jobs_outdir]:
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
        parser.add_argument("--exclude",
                            nargs="*",
                            type=str,
                            choices=CHROMOSOMES,
                            default=None,
                            help="Exclude certain chromosomes.")
        parser.add_argument("--annotations",
                            type=str,
                            required=True,
                            help="Bgzip-compressed and tabix-indexed file "
                                 "with annotations..")
        parser.add_argument("--time_per_chr",
                            type=str,
                            default="medium",
                            choices=["short", "medium", "long"],
                            help="Sets the work allocation time a.k.a. "
                                 "walltime per chromosome job.")
        parser.add_argument("--cpus_per_chr",
                            type=int,
                            default=1,
                            help="Requests X CPUs (cores) for each chromosome.")
        parser.add_argument("--mem_per_chr",
                            type=int,
                            default=1,
                            help="Requests X GB RAM total for each chromosome.")
        parser.add_argument("--time_merge",
                            type=str,
                            default="short",
                            choices=["short", "medium", "long"],
                            help="Sets the work allocation time a.k.a. "
                                 "walltime for the merge job.")
        parser.add_argument("--cpus_merge",
                            type=int,
                            default=4,
                            help="Requests X CPUs (cores) for the merge job.")
        parser.add_argument("--mem_merge",
                            type=int,
                            default=16,
                            help="Requests X GB RAM total for the merge job.")
        parser.add_argument("--naive_merge",
                            action='store_false',
                            help="Concatenate VCF or BCF files without "
                                 "recompression. This is very fast but "
                                 "requires that all files are of the same "
                                 "type (all VCF or all BCF) and have the same "
                                 "headers. This is because all tags and "
                                 "chromosome names in the BCF body rely on "
                                 "the order of the contig and tag definitions "
                                 "in the header. A header check compatibility "
                                 "is performed and the program throws an "
                                 "error if it is not safe to use the option.")
        parser.add_argument("--outdir",
                            type=str,
                            required=True,
                            help="The path for the output directory.")
        parser.add_argument("--outfolder",
                            type=str,
                            required=True,
                            help="The name of the output directory.")
        parser.add_argument("--dryrun",
                            action='store_true',
                            help="Add this flag to disable submitting the"
                                 "job files.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Look for available files")
        files = []
        if "CHR" in self.vcf_file_format:
            for chr in self.chromosomes:
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
        chr_vcf_files = []
        for (chr, file) in files:
            joblog_path, jobfile_path, chr_vcf_file = self.create_job_file(
                chr=chr,
                infile=file
            )
            job_files.append((joblog_path, jobfile_path))
            chr_vcf_files.append(chr_vcf_file)

        if not self.dryrun:
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

        print("Generating merge and convert script")
        self.create_merge_and_convert_job_file(
            infiles=chr_vcf_files
        )

    def create_job_file(self, chr, infile):
        job_name = "VCFFilter_CHR{}".format(chr)
        joblog_path = os.path.join(self.jobs_outdir, job_name + ".out")

        norm_outfile = os.path.join(self.norm_outdir, os.path.basename(infile).replace(".vcf.gz", "_norm.vcf.gz"))
        filter_outfile = os.path.join(self.filtered_outdir, os.path.basename(norm_outfile).replace(".vcf.gz", "_vcffilter"))
        bgzip_outfile = os.path.join(self.bgzip_outdir, os.path.basename(filter_outfile) + ".vcf.gz")
        annotated_file = os.path.join(self.annotated_outdir, os.path.basename(bgzip_outfile).replace(".vcf.gz", "_annotated.vcf.gz"))

        lines = ["#!/bin/bash",
                 "#SBATCH --job-name={}".format(job_name),
                 "#SBATCH --output={}".format(joblog_path),
                 "#SBATCH --error={}".format(joblog_path),
                 "#SBATCH --time={}".format(self.time_per_chr),
                 "#SBATCH --cpus-per-task={}".format(self.cpus_per_chr),
                 "#SBATCH --mem={}gb".format(self.mem_per_chr),
                 "#SBATCH --nodes=1",
                 "#SBATCH --open-mode=append",
                 "#SBATCH --export=NONE",
                 "#SBATCH --get-user-env=L",
                 "",
                 "module load BCFtools/1.16-GCCcore-7.3.0",
                 "module load Python/3.7.4-GCCcore-7.3.0-bare",
                 "",
                 "bcftools norm -m -any {} -o {}".format(infile, norm_outfile),
                 "",
                 "python3 {} {} {}".format(self.filter_script, norm_outfile, filter_outfile),
                 "",
                 "bcftools view {}-filtered.vcf.gz -Oz -o {}".format(filter_outfile, bgzip_outfile),
                 "bcftools index {}".format(bgzip_outfile),
                 # "bgzip -c {}-filtered.vcf.gz > {}".format(filter_outfile, bgzip_outfile),
                 # "tabix -p vcf {}".format(bgzip_outfile),
                 "",
                 "bcftools annotate -a {} -c ID {} -o {}".format(self.annotations, bgzip_outfile, annotated_file),
                 ]

        jobfile_path = os.path.join(self.jobdir, job_name + ".sh")
        self.write_lines_to_file(
            lines=lines,
            filepath=jobfile_path
        )
        return joblog_path, jobfile_path, annotated_file
    def create_merge_and_convert_job_file(self, infiles):
        job_name = "MERGE_CHR_AND_TO_PGEN"
        joblog_path = os.path.join(self.jobs_outdir, job_name + ".out")

        merged_outfile = os.path.join(self.filtered_outdir, os.path.basename(infiles[0]).replace(".vcf.gz", "_allchromosomes.vcf.gz"))
        if "CHR" in self.vcf_file_format and set([chr_file.replace("CHR", "") for chr_file in infiles]) == 0:
            merged_outfile = os.path.join(self.merged_outdir, os.path.basename(infiles[0]).replace("CHR", "all"))
        plink_outfile = os.path.join(self.plink_outdir, os.path.basename(merged_outfile).replace(".vcf.gz", ""))

        lines = ["#!/bin/bash",
                 "#SBATCH --job-name={}".format(job_name),
                 "#SBATCH --output={}".format(joblog_path),
                 "#SBATCH --error={}".format(joblog_path),
                 "#SBATCH --time={}".format(self.time_merge),
                 "#SBATCH --cpus-per-task={}".format(self.cpus_merge),
                 "#SBATCH --mem={}gb".format(self.mem_merge),
                 "#SBATCH --nodes=1",
                 "#SBATCH --open-mode=append",
                 "#SBATCH --export=NONE",
                 "#SBATCH --get-user-env=L",
                 "",
                 "module load BCFtools/1.16-GCCcore-7.3.0",
                 "bcftools concat {}{} -o {}".format(" ".join(infiles), " --naive" if self.naive_merge else "", merged_outfile),
                 "",
                 "module load PLINK/2.0-alpha2-20191006",
                 "plink2 --vcf {} --make-pgen --out {}".format(merged_outfile, plink_outfile),
                 ""
                 ]

        self.write_lines_to_file(
            lines=lines,
            filepath=os.path.join(self.jobdir, job_name + ".sh")
        )

    @staticmethod
    def write_lines_to_file(lines, filepath):
        with open(filepath, "w") as f:
            for line in lines:
                f.write(line + "\n")
        f.close()
        print("\tSaved jobfile: {}".format(os.path.basename(filepath)))

    @staticmethod
    def run_command(command):
        print("\t" + " ".join(command))
        subprocess.call(command)

    def print_arguments(self):
        print("Arguments:")
        print("  > filter_script: {}".format(self.filter_script))
        print("  > vcf_indir: {}".format(self.vcf_indir))
        print("  > vcf_file_format: {}".format(self.vcf_file_format))
        print("  > exclude: {}".format(self.exclude))
        print("  > annotations: {}".format(self.annotations))
        print("  > per chromosome jobs:")
        print("    > Time: {}".format(self.time_per_chr))
        print("    > CPU-per-task: {}Gb".format(self.cpus_per_chr))
        print("    > Memory: {}Gb".format(self.mem_per_chr))
        print("  > merge job:")
        print("    > Time: {}".format(self.time_merge))
        print("    > CPU-per-task: {}Gb".format(self.cpus_merge))
        print("    > Memory: {}Gb".format(self.mem_merge))
        print("    > Naive: {}".format(self.naive_merge))
        print("  > outdir: {}".format(self.outdir))
        print("  > dryrun: {}".format(self.dryrun))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()