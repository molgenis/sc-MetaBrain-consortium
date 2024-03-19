#!/usr/bin/env python3

"""
File:         filter_WGS_VCF_files.py
Created:      2022/10/19
Last Changed: 2023/01/13
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
from datetime import datetime
import subprocess
import argparse
import os
import re
import time

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
./filter_WGS_VCF_files.py -h
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
        self.exclude_chr = getattr(arguments, 'exclude_chr')
        self.exclude_snp = getattr(arguments, 'exclude_snp')
        self.sex_path = getattr(arguments, 'sex')
        self.annotations = getattr(arguments, 'annotations')
        outdir = getattr(arguments, 'outdir')
        outfolder = getattr(arguments, 'outfolder')
        self.dryrun = getattr(arguments, 'dryrun')

        if not self.vcf_file_format.endswith(".vcf.gz"):
            print("Error: --vcf_file_format must end with '.vcf.gz'.")
            exit()
        if self.vcf_file_format.count('CHR') > 1:
            print("Error: --vcf_file_format must only contain 'CHR' once.")
            exit()

        self.chromosomes = [CHR for CHR in CHROMOSOMES if (self.exclude_chr is None) or (CHR not in self.exclude_chr)]
        self.default_id = "%CHROM:%POS:%REF{10}:%ALT{10}"

        # Pre-process the dataset output directory.
        date_str = datetime.now().strftime("%Y-%m-%d")
        if not re.match("^[0-9]{4}-(0[1-9]|1[0-2])-(0[1-9]|[1-2][0-9]|3[0-1])", outfolder):
            outfolder = "{}-{}".format(date_str, outfolder)
        self.outdir = os.path.join(outdir, outfolder)
        self.norm_outdir = os.path.join(self.outdir, "1-bcftools_norm")
        self.filter_outdir = os.path.join(self.outdir, "2-custom_vcffilter")
        self.merged_outdir = os.path.join(self.outdir, "3-bcftools_concat")
        self.annotate_outdir = os.path.join(self.outdir, "4-bcftools_annotate")
        self.plink_outdir = os.path.join(self.outdir, "5-plink2_makepgen")
        self.jobdir = os.path.join(self.outdir, "jobs")
        self.jobs_outdir = os.path.join(self.jobdir, "output")
        for dir in [self.outdir, self.norm_outdir, self.filter_outdir,
                    self.merged_outdir, self.annotate_outdir,
                    self.plink_outdir, self.jobdir, self.jobs_outdir]:
            if not os.path.exists(dir):
                os.makedirs(dir)

        if ("X" in self.chromosomes or "Y" in self.chromosomes) and (self.sex_path is None):
            print("Error, use --sex argument when analyzing sex chromosomes")
            exit()

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
        parser.add_argument("--exclude_chr",
                            nargs="*",
                            type=str,
                            choices=CHROMOSOMES,
                            default=None,
                            help="Exclude certain chromosomes.")
        parser.add_argument("--exclude_snp",
                            nargs="*",
                            type=str,
                            default=None,
                            help="Exclude certain variants.")
        parser.add_argument("--sex",
                            type=str,
                            required=False,
                            help="The sample-sex data file.")
        parser.add_argument("--annotations",
                            type=str,
                            required=True,
                            help="Bgzip-compressed and tabix-indexed file "
                                 "with annotations..")
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

        print("\nSubmitting job scripts")
        filter_vcf_outfiles = []
        filter_jobids = []
        for (chr, file) in files:
            print("")

            ####################################################################

            norm_vcf_outfile = os.path.join(self.norm_outdir, os.path.basename(file).replace(".vcf.gz", "_norm.vcf.gz"))
            norm_jobfile, norm_logfile = self.create_jobfile(
                job_name="BCFTOOLS_NORM_CHR{}".format(chr),
                module_load=["BCFtools/1.16-GCCcore-7.3.0"],
                commands=["bcftools norm -m -any {} -o {}".format(file, norm_vcf_outfile)]
            )

            norm_jobid = self.submit_job(
                jobfile=norm_jobfile,
                logfile=norm_logfile
            )

            ####################################################################

            filter_vcf_outfile = os.path.join(self.filter_outdir, os.path.basename(norm_vcf_outfile).replace(".vcf.gz", "_vcffilter"))
            filter_jobfile, filter_logfile = self.create_jobfile(
                job_name="CUSTOMVCFFILTER_CHR{}".format(chr),
                module_load=["Python/3.7.4-GCCcore-7.3.0-bare"],
                commands=["python3 {} -i {}{} -o {}".format(self.filter_script, norm_vcf_outfile, " --sex {}".format(self.sex_path) if chr in ["X", "Y"] else "", filter_vcf_outfile)],
                time="medium"
            )
            filter_vcf_outfiles.append(filter_vcf_outfile)

            filter_jobid = self.submit_job(
                jobfile=filter_jobfile,
                logfile=filter_logfile,
                depend=norm_jobid
            )
            filter_jobids.append(filter_jobid)

            ####################################################################

        if len(files) < 2:
            return

        ####################################################################

        print("\nSubmitting merge job script")
        concat_outfile = os.path.join(self.merged_outdir, self.vcf_file_format.replace(".vcf.gz", "_norm_vcffilter_concat.vcf.gz"))
        concat_jobfile, concat_logfile = self.create_jobfile(
            job_name="BCFTOOLS_CONCAT_AND_INDEX",
            module_load=["BCFtools/1.16-GCCcore-7.3.0"],
            commands=["bcftools concat {} -o {}".format(" ".join(["{}-filtered.vcf.gz".format(filter_vcf_outfile) for filter_vcf_outfile in filter_vcf_outfiles]), concat_outfile),
                      "bcftools index {}".format(concat_outfile)]
        )

        concat_jobid = self.submit_job(
            jobfile=concat_jobfile,
            logfile=concat_logfile,
            depend=filter_jobids
        )

        ####################################################################

        print("\nSubmitting annotate script")
        annotate_outfile = os.path.join(self.annotate_outdir, os.path.basename(concat_outfile).replace(".vcf.gz", "_annotate.vcf.gz"))
        annotate_jobfile, annotate_logfile = self.create_jobfile(
            job_name="BCFTOOLS_ANNOTATE",
            module_load=["BCFtools/1.16-GCCcore-7.3.0"],
            commands=["bcftools annotate -a {} -c ID {} --set-id +'{}' -o {}".format(self.annotations, concat_outfile, self.default_id, annotate_outfile)],
            time="medium",
            cpu=1,
            mem=2,
        )

        annotate_jobid = self.submit_job(
            jobfile=annotate_jobfile,
            logfile=annotate_logfile,
            depend=concat_jobid
        )

        ####################################################################

        print("\nSubmitting Plink job script")
        exclude_str = ""
        if self.exclude_snp:
            exclude_filepath = os.path.join(self.outdir, "exclude_variants.txt")
            exclude_str = " --exclude {}".format(exclude_filepath)
            self.write_lines_to_file(lines=[self.exclude_snp], filepath=exclude_filepath)
        plink_outfile = os.path.join(self.plink_outdir, os.path.basename(annotate_outfile).replace(".vcf.gz", ""))
        plink_cpu = 4
        plink_mem = 6
        plink_jobfile, plink_logfile = self.create_jobfile(
            job_name="PLINK2_MAKEPGEN",
            module_load=["PLINK/2.0-alpha2-20191006"],
            commands=["plink2 --vcf {} --make-pgen{} --threads {} "
                      "--memory {} --out {}".format(annotate_outfile,
                                                    exclude_str,
                                                    plink_cpu * 2,
                                                    plink_mem * 1000,
                                                    plink_outfile)],
            cpu=plink_cpu,
            mem=plink_mem
        )

        _ = self.submit_job(
            jobfile=plink_jobfile,
            logfile=plink_logfile,
            depend=annotate_jobid
        )

    def create_jobfile(self, job_name, commands, time="short", cpu=1,
                       mem=1, module_load=None):
        joblog_path = os.path.join(self.jobs_outdir, job_name + ".out")

        lines = ["#!/bin/bash",
                 "#SBATCH --job-name={}".format(job_name),
                 "#SBATCH --output={}".format(joblog_path),
                 "#SBATCH --error={}".format(joblog_path),
                 "#SBATCH --time={}".format(TIME_DICT[time]),
                 "#SBATCH --cpus-per-task={}".format(cpu),
                 "#SBATCH --mem={}gb".format(mem),
                 "#SBATCH --nodes=1",
                 "#SBATCH --open-mode=append",
                 "#SBATCH --export=NONE",
                 "#SBATCH --get-user-env=L",
                 ""]

        if module_load is not None and isinstance(module_load, list):
            for module in module_load:
                lines.append("module load {}".format(module))

        for command in commands:
            lines.extend(["", command])

        lines.extend(["", "echo 'Job finished'"])

        jobfile_path = os.path.join(self.jobdir, job_name + ".sh")
        self.write_lines_to_file(
            lines=lines,
            filepath=jobfile_path
        )
        return jobfile_path, joblog_path

    @staticmethod
    def write_lines_to_file(lines, filepath):
        with open(filepath, "w") as f:
            for line in lines:
                f.write(line + "\n")
        f.close()
        print("\tSaved jobfile: {}".format(os.path.basename(filepath)))

    def submit_job(self, jobfile, logfile, depend=None):
        if self.dryrun:
            return None

        if os.path.exists(logfile):
            for line in open(logfile, 'r'):
                if line == "Job finished":
                    print("\t\tSkipping '{}'".format(os.path.basename(jobfile)))
                    return None

        command = ['sbatch', jobfile]
        if depend is not None and (isinstance(depend, str) or isinstance(depend, list)):
            if isinstance(depend, str):
                depend = [depend]
            command.insert(1, '--depend=afterok:{}'.format(":".join(depend)))
            time.sleep(1)

        print("\t\t" + " ".join(command))
        sucess, output = self.run_command(command=command)
        print("\t\t{}".format(output))
        if not sucess:
            print("\t\tError, subprocess failed.")
            exit()

        return output.replace("Submitted batch job ", "").replace("\n", "")

    @staticmethod
    def run_command(command):
        success = False
        try:
            output = subprocess.check_output(command, stderr=subprocess.STDOUT).decode()
            success = True
        except subprocess.CalledProcessError as e:
            output = e.output.decode()
        except Exception as e:
            output = str(e)
        return success, output

    def print_arguments(self):
        print("Arguments:")
        print("  > filter_script: {}".format(self.filter_script))
        print("  > vcf_indir: {}".format(self.vcf_indir))
        print("  > vcf_file_format: {}".format(self.vcf_file_format))
        print("  > exclude chr: {}".format(self.exclude_chr))
        print("  > exclude SNP: {}".format(self.exclude_snp))
        print("  > Sex input: {}".format(self.sex_path))
        print("  > outdir: {}".format(self.outdir))
        print("  > dryrun: {}".format(self.dryrun))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()