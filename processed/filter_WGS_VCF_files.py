#!/usr/bin/env python3

"""
File:         filter_WGS_VCF_files.py
Created:      2022/10/19
Last Changed: 2022/11/03
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

### AMP-AD ###
./filter_WGS_VCF_files.py \
    --filter_script /groups/umcg-biogen/tmp01/input/rawdata/2017-GTExV8Genotypes/customVCFFilterV4.py \
    --vcf_indir /groups/umcg-biogen/tmp01/input/AMP-AD/2017-12-08-joint-WGS \
    --vcf_file_format NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_CHR.recalibrated_variants.vcf.gz \
    --exclude X Y others \
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
        self.filter_outdir = os.path.join(self.outdir, "2-custom_vcffilter")
        self.merged_outdir = os.path.join(self.outdir, "3-bcftools_concat")
        self.plink_outdir = os.path.join(self.outdir, "4-plink2_makepgen")
        self.jobdir = os.path.join(self.outdir, "jobs")
        self.jobs_outdir = os.path.join(self.jobdir, "output")
        for dir in [self.outdir, self.norm_outdir, self.filter_outdir,
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

        print("\nSubmitting job scripts")
        filter_vcf_outfiles = []
        filter_jobids = []
        for (chr, file) in files:
            ####################################################################

            norm_vcf_outfile = os.path.join(self.norm_outdir, os.path.basename(file).replace(".vcf.gz", "_norm.vcf.gz"))
            norm_jobfile, norm_logfile = self.create_jobfile(
                job_name="BCFTOOLS_NORM_CHR{}".format(chr),
                module_load=["BCFtools/1.16-GCCcore-7.3.0"],
                commands=["bcftools norm -m -any {} -o {}".format(file, norm_vcf_outfile)]
            )

            completed = False
            if os.path.exists(norm_logfile):
                    for line in open(norm_logfile, 'r'):
                        if line.startswith("Lines   total/split/realigned/skipped:"):
                            completed = True

            norm_jobid = None
            if completed:
                print("\tSkipping '{}'".format(os.path.basename(norm_jobfile)))
            else:
                norm_jobid = self.submit_job(jobfile=norm_jobfile)
            time.sleep(1)

            ####################################################################

            filter_vcf_outfile = os.path.join(self.filter_outdir, os.path.basename(norm_vcf_outfile).replace(".vcf.gz", "_vcffilter"))
            filter_jobfile, filter_logfile = self.create_jobfile(
                job_name="CUSTOMVCFFILTER_CHR{}".format(chr),
                module_load=["Python/3.7.4-GCCcore-7.3.0-bare"],
                commands=["python3 {} {} {}".format(self.filter_script, norm_vcf_outfile, filter_vcf_outfile)],
                time="medium"
            )
            filter_vcf_outfiles.append(filter_vcf_outfile)

            completed = False
            if os.path.exists(filter_logfile):
                    for line in open(filter_logfile, 'r'):
                        if "Done. How about that!" in line:
                            completed = True

            filter_jobid = None
            if completed:
                print("\tSkipping '{}'".format(os.path.basename(filter_jobfile)))
            else:
                filter_jobid = self.submit_job(jobfile=filter_jobfile, depend=norm_jobid)
            filter_jobids.append(filter_jobid)

            ####################################################################

        if len(files) < 2:
            return

        print("\nGenerating merge and convert script")
        merged_outfile = os.path.join(self.merged_outdir,
                                      os.path.basename(filter_vcf_outfiles[0]).replace(".vcf.gz", "_allchromosomes.vcf.gz"))
        if "CHR" in self.vcf_file_format and set([chr_file.replace("CHR", "") for chr_file in filter_vcf_outfiles]) == 0:
            merged_outfile = os.path.join(self.merged_outdir,
                                          os.path.basename(filter_vcf_outfiles[0]).replace("CHR", "all"))

        plink_outfile = os.path.join(self.plink_outdir, os.path.basename(merged_outfile).replace(".vcf.gz", ""))
        merge_jobfile, _ = self.create_jobfile(
            job_name="CONCAT_AND_MAKEPGEN",
            module_load=["BCFtools/1.16-GCCcore-7.3.0",
                         "PLINK/2.0-alpha2-20191006"],
            commands=["bcftools concat {}{} -o {}".format(" ".join(filter_vcf_outfiles),
                                                          " --naive" if self.naive_merge else "",
                                                          merged_outfile),
                      "plink2 --vcf {} --make-pgen --out {}".format(merged_outfile,
                                                                    plink_outfile)],
            time="medium",
            cpu=4,
            mem=32
        )

        # _ = self.submit_job(jobfile=merge_jobfile, depend=filter_jobids)

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

    def submit_job(self, jobfile, depend=None):
        if self.dryrun:
            return None

        command = ['sbatch', jobfile]
        if depend is not None and (isinstance(depend, str) or isinstance(depend, list)):
            if isinstance(depend, str):
                depend = [depend]
            command.insert(1, '--depend=afterok:{}'.format(":".join(depend)))

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
        print("  > exclude: {}".format(self.exclude))
        print("  > Naive merge: {}".format(self.naive_merge))
        print("  > outdir: {}".format(self.outdir))
        print("  > dryrun: {}".format(self.dryrun))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()