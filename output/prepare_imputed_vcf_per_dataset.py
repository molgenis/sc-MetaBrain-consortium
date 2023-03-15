#!/usr/bin/env python3

"""
File:         prepare_imputed_vcf_per_dataset.py
Created:      2023/02/23
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
import subprocess
import argparse
import os
import time

# Third party imports.

# Local application imports.

# Metadata
__program__ = "Prepare Imputed VCF Per Dataset"
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

### Mathys ###
./prepare_imputed_vcf_per_dataset.py \
    --work_dir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2023-02-02-AMP_AD/Step1-Imputation \
    --singularity /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2022-10-07-Imputation/WG1-pipeline-QC_imputation.sif \
    --dataset Mathys \
    --samples /groups/umcg-biogen/tmp01/input/processeddata/single-cell/Mathys2019/Mathys_genotype_samples.txt \
    --dryrun
    
"""

TIME_DICT = {
    "short": "05:59:59",
    "medium": "23:59:00",
    "long": "6-23:59:00"
}


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.work_dir = getattr(arguments, 'work_dir')
        self.bind_paths = ",".join(getattr(arguments, 'bind_paths'))
        self.singularity = getattr(arguments, 'singularity')
        self.dataset = getattr(arguments, 'dataset')
        self.samples = getattr(arguments, 'samples')
        self.dryrun = getattr(arguments, 'dryrun')

        self.vcf_all_merged_dataset_dir = os.path.join(self.work_dir, "vcf_all_merged", self.dataset)
        self.vcf_4_demultiplex_dataset_dir = os.path.join(self.work_dir, "vcf_4_demultiplex", self.dataset)
        self.metrics_dataset_dir = os.path.join(self.work_dir, "metrics", self.dataset)
        self.jobdir = os.path.join(self.work_dir, "jobs", self.dataset)
        self.jobs_outdir = os.path.join(self.jobdir, "output")
        for dir in [self.vcf_all_merged_dataset_dir,
                    self.vcf_4_demultiplex_dataset_dir,
                    self.metrics_dataset_dir,
                    self.jobdir,
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
        parser.add_argument("--work_dir",
                            type=str,
                            required=True,
                            default=None,
                            help="The working directory.")
        parser.add_argument("--bind_paths",
                            nargs="*",
                            type=str,
                            default=["/groups/umcg-biogen/tmp01/",
                                     "/groups/umcg-biogen/tmp01/umcg-mvochteloo/simulated_home:/home/umcg-mvochteloo"],
                            help="List of paths to bind to Singularity. "
                                 "Default: ['/groups/umcg-biogen/tmp01/']")
        parser.add_argument("--singularity",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--dataset",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--samples",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--dryrun",
                            action='store_true',
                            help="Add this flag to disable submitting the"
                                 "job files.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        ####################################################################

        print("\nSubmitting extract dataset job script")
        imputed_vcf = os.path.join(self.work_dir, "vcf_all_merged", "imputed_hg38.vcf.gz")
        dataset_imputed_vcf = os.path.join(self.vcf_all_merged_dataset_dir, "{}_imputed_hg38.vcf.gz".format(self.dataset))
        extract_jobfile, extract_logfile = self.create_jobfile(
            job_name="EXTRACT_DATASET",
            module_load=["Anaconda3"],
            commands=["source activate wg1_snakemake",
                      "singularity exec --bind {} {} bcftools view -S {} {} -Oz -o {}".format(self.bind_paths, self.singularity, self.samples, imputed_vcf, dataset_imputed_vcf),
                      "singularity exec --bind {} {} bcftools index {}".format(self.bind_paths, self.singularity, dataset_imputed_vcf)
                      ],
            cpu=1,
            mem=8,
            time="medium"
        )

        extract_jobid = self.submit_job(
            jobfile=extract_jobfile,
            logfile=extract_logfile
        )

        del imputed_vcf, extract_jobfile, extract_logfile

        ####################################################################

        print("\nSubmitting filter 4 demultiplexing script")
        dataset_info_filled_vcf = os.path.join(self.vcf_all_merged_dataset_dir, "{}_imputed_hg38_info_filled.vcf.gz".format(self.dataset))
        dataset_qc_filtered_vcf = os.path.join(self.vcf_all_merged_dataset_dir, "{}_imputed_hg38_R2_0.3_MAF0.05.vcf.gz".format(self.dataset))
        dataset_location_filtered = os.path.join(self.vcf_all_merged_dataset_dir, "{}_imputed_hg38_R2_0.3_MAF0.05_exons".format(self.dataset))
        dataset_location_filtered_complete = os.path.join(self.vcf_all_merged_dataset_dir, "{}_imputed_hg38_R2_0.3_MAF0.05_exons_complete_cases".format(self.dataset))
        filter4demultiplexing_jobfile, filter4demultiplexing_logfile = self.create_jobfile(
            job_name="FILTER4DEMULTIPLEXING",
            module_load=["Anaconda3"],
            commands=["source activate wg1_snakemake",
                      "singularity exec --bind {} {} bcftools +fill-tags -Oz --output {} {}".format(self.bind_paths, self.singularity, dataset_info_filled_vcf, dataset_imputed_vcf),
                      "singularity exec --bind {} {} bcftools filter --include 'MAF>=0.05 & R2>=0.3' -Oz --output {} {}".format(self.bind_paths, self.singularity, dataset_qc_filtered_vcf, dataset_info_filled_vcf),
                      "singularity exec --bind {} {} vcftools --gzvcf {} --max-alleles 2 --remove-indels --bed /opt/hg38exonsUCSC.bed --recode --recode-INFO-all --out {}".format(self.bind_paths, self.singularity, dataset_qc_filtered_vcf, dataset_location_filtered),
                      "singularity exec --bind {} {} vcftools --recode --recode-INFO-all --vcf {}.recode.vcf --max-missing 1 --out {}".format(
                          self.bind_paths, self.singularity, dataset_location_filtered, dataset_location_filtered_complete)
                      ],
            cpu=4,
            mem=20

        )

        filter4demultiplexing_jobid = self.submit_job(
            jobfile=filter4demultiplexing_jobfile,
            logfile=filter4demultiplexing_logfile,
            depend=extract_jobid
        )

        del dataset_info_filled_vcf, dataset_qc_filtered_vcf, dataset_location_filtered, filter4demultiplexing_jobfile, filter4demultiplexing_logfile

        ####################################################################

        print("\nSubmitting sort 4 demultiplexing script")
        dataset_complete_cases_sorted_vcf = os.path.join(self.vcf_4_demultiplex_dataset_dir, "{}_imputed_hg38_R2_0.3_MAF0.05_exons_sorted.vcf".format(self.dataset))
        sort4demultiplexing_mem = 25
        sort4demultiplexing_jobfile, sort4demultiplexing_logfile = self.create_jobfile(
            job_name="SORT4DEMULTIPLEXING",
            module_load=["Anaconda3"],
            commands=["source activate wg1_snakemake",
                      "singularity exec --bind {} {} java -Xmx{}g -Xms{}g -jar /opt/picard/build/libs/picard.jar SortVcf I={}.recode.vcf O={}".format(self.bind_paths, self.singularity, sort4demultiplexing_mem - 2, sort4demultiplexing_mem - 2, dataset_location_filtered_complete, dataset_complete_cases_sorted_vcf)],
            cpu=5,
            mem=sort4demultiplexing_mem
        )

        _ = self.submit_job(
            jobfile=sort4demultiplexing_jobfile,
            logfile=sort4demultiplexing_logfile,
            depend=filter4demultiplexing_jobid
        )

        ####################################################################

        print("\nSubmitting count SNPs script")
        count_snps_jobfile, count_snps_logfile = self.create_jobfile(
            job_name="COUNT_SNPS",
            module_load=["Anaconda3"],
            commands=["source activate wg1_snakemake",
                      "singularity exec --bind {} {} Rscript /opt/WG1-pipeline-QC/Imputation/scripts/SNP_numbers.R {} {}".format(self.bind_paths, self.singularity, self.work_dir, self.metrics_dataset_dir)]
        )

        _ = self.submit_job(
            jobfile=count_snps_jobfile,
            logfile=count_snps_logfile,
            depend=filter4demultiplexing_jobid
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
                 "#SBATCH --tmp={}gb".format(mem),
                 "#SBATCH --nodes=1",
                 "#SBATCH --export=ALL",
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
        print("  > work_dir: {}".format(self.work_dir))
        print("  > bind_paths: {}".format(self.bind_paths))
        print("  > singularity: {}".format(self.singularity))
        print("  > dataset: {}".format(self.dataset))
        print("  > samples: {}".format(self.samples))
        print("  > dryrun: {}".format(self.dryrun))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()