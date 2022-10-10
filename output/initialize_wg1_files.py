#!/usr/bin/env python3

"""
File:         initialize_wg1_files.py
Created:      2022/10/07
Last Changed: 2022/10/10
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
import argparse
import os

# Third party imports.

# Local application imports.

# Metadata
__program__ = "Initialize Workgroup1 Files"
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

### ImputationTestDataset ###
./initialize_wg1_files.py \
    --work_dir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC \
    --imputation_subdir 2022-10-07-Imputation \
    --demultiplexing_subdir 2022-10-10-DemultiplexingAndDoubletRemoval \
    --ref_dir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/hg38 \
    --plink_dir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2022-10-07-Imputation/ImputationTestDataset_plink \
    --samplesheet_filepath /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2022-10-10-DemultiplexingAndDoubletRemoval/TestData4PipelineSmall/samplesheet.txt \
    --scRNAseq_dir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2022-10-10-DemultiplexingAndDoubletRemoval/TestData4PipelineSmall \
    --snp_genotypes_filepath /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2022-10-10-DemultiplexingAndDoubletRemoval/TestData4PipelineSmall/test_dataset.vcf \
    --individual_list_dir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2022-10-10-DemultiplexingAndDoubletRemoval/TestData4PipelineSmall/individuals_list_dir \
    --dataset_outdir ImputationTestDataset
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        work_dir = getattr(arguments, 'work_dir')
        imputation_subdir = getattr(arguments, 'imputation_subdir')
        demultiplexing_subdir = getattr(arguments, 'demultiplexing_subdir')
        self.ref_dir = getattr(arguments, 'ref_dir')
        self.plink_dir = getattr(arguments, 'plink_dir')
        self.samplesheet_filepath = getattr(arguments, 'samplesheet_filepath')
        self.scRNAseq_dir = getattr(arguments, 'scRNAseq_dir')
        self.snp_genotypes_filepath = getattr(arguments, 'snp_genotypes_filepath')
        self.individual_list_dir = getattr(arguments, 'individual_list_dir')
        self.individual_list_dir = getattr(arguments, 'individual_list_dir')
        self.bind_paths = ", ".join(getattr(arguments, 'bind_paths'))
        dataset_outdir = getattr(arguments, 'dataset_outdir')

        # Define general variables.
        self.imputation_dir = os.path.join(work_dir, imputation_subdir)
        self.imputation_snakefile = os.path.join(self.imputation_dir, "Snakefile")
        self.demultiplexing_dir = os.path.join(work_dir, demultiplexing_subdir)
        self.demultiplexing_snakefile = os.path.join(self.demultiplexing_dir, "Snakefile")

        date_str = datetime.now().strftime("%Y-%m-%d")
        dataset_work_dir = os.path.join(work_dir, "{}-{}".format(date_str, dataset_outdir))

        # Define step 1 variables.
        self.dataset_imputation_output_dir = os.path.join(dataset_work_dir, "Step1-Imputation")
        self.dataset_imputation_logdir = os.path.join(self.dataset_imputation_output_dir, "log")
        self.dataset_imputation_configfile = os.path.join(self.dataset_imputation_output_dir, "PreImputation.yaml")

        # Define step 2 variables.
        self.dataset_demultiplexing_output_dir = os.path.join(dataset_work_dir, "Step2-DemultiplexingAndDoubletRemoval")
        self.dataset_demultiplexing_logdir = os.path.join(self.dataset_demultiplexing_output_dir, "log")
        self.dataset_demultiplexing_configfile = os.path.join(self.dataset_demultiplexing_output_dir, "sceQTL-Gen_Demultiplex.yaml")

        # Creating output directories
        for outdir in [dataset_work_dir,
                       self.dataset_imputation_output_dir,
                       self.dataset_imputation_logdir,
                       self.dataset_demultiplexing_output_dir,
                       self.dataset_demultiplexing_logdir,
                       ]:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

        self.time_dict = {
            "short": "05:59:59",
            "medium": "23:59:00",
            "long": "6-23:59:00"
        }

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
                            help="The working directory.")
        parser.add_argument("--imputation_subdir",
                            type=str,
                            required=True,
                            help="The imputation step subdirectory.")
        parser.add_argument("--demultiplexing_subdir",
                            type=str,
                            required=True,
                            help="The demultiplexing and doublet removal step "
                                 "subdirectory.")
        parser.add_argument("--ref_dir",
                            type=str,
                            required=True,
                            help="This is the path to the directory containing "
                                 "the imputation references provided for the "
                                 "SNP processing and imputation")
        parser.add_argument("--plink_dir",
                            type=str,
                            required=True,
                            help="Absolute path to directory that contains "
                                 "plink2 pgen, pvar and psam files for "
                                 "pipeline on hg19.")
        parser.add_argument("--samplesheet_filepath",
                            type=str,
                            required=True,
                            help="Tab separated file that has a header. Each "
                                 "line has a pool name used for the scRNA-seq "
                                 "directories and the number of individuals in "
                                 "each pool.")
        parser.add_argument("--scRNAseq_dir",
                            type=str,
                            required=True,
                            help="The parent directory that has directories "
                                 "for each pool and the scRNA-seq output "
                                 "below it.")
        parser.add_argument("--snp_genotypes_filepath",
                            type=str,
                            required=True,
                            help="The path to the genotype file that has just "
                                 "SNPs that are within exons and with a minor"
                                 " allele frequency > 0.05 (5%) in this "
                                 "sample.")
        parser.add_argument("--individual_list_dir",
                            type=str,
                            required=True,
                            help="Directory that has a different file for "
                                 "each pool. Each file contains a list of "
                                 "the individuals IDs that are in the pool "
                                 "separated by line (match the genotype "
                                 "individual IDs).")
        parser.add_argument("--bind_paths",
                            nargs="*",
                            type=str,
                            required=False,
                            default=["/groups/umcg-biogen/tmp01/"],
                            help="List of paths to bind to Singularity. "
                                 "Default: ['/groups/umcg-biogen/tmp01/']")
        parser.add_argument("--dataset_outdir",
                            type=str,
                            required=True,
                            help="The name of the output directory where you "
                                 "would like all outputs/results saved.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Generating step 1 files")
        self.generate_step1_files(
            output_dir=self.dataset_imputation_output_dir,
            log_dir=self.dataset_imputation_logdir,
            snakefile=self.imputation_snakefile,
            configfile=self.dataset_imputation_configfile
        )

        print("Generating step 2 files")
        self.generate_step2_files(
            output_dir=self.dataset_demultiplexing_output_dir,
            log_dir=self.dataset_demultiplexing_logdir,
            snakefile=self.demultiplexing_snakefile,
            configfile=self.dataset_demultiplexing_configfile
        )

    def generate_step1_files(self, output_dir, log_dir, snakefile, configfile):
        self.write_imputation_configfile(
            output_dir=output_dir,
            outpath=configfile
        )

        self.write_dry_run_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir
        )

        self.write_run_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir,
            log_dir=log_dir
        )

        self.write_run_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir,
            log_dir=log_dir,
            jobs=20,
            outfile="run_20jobs"
        )

        self.write_report_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir
        )


    def write_imputation_configfile(self, output_dir, outpath):
        yaml_lines = []
        for line in open(os.path.join(self.imputation_dir, "PreImputation.yaml"), 'r'):
            for label, argument in (("ref_dir", self.ref_dir),
                                    ("singularity_image", os.path.join(self.imputation_dir, "WG1-pipeline-QC_imputation.sif")),
                                    ("plink_dir", self.plink_dir),
                                    ("bind_paths", self.bind_paths),
                                    ("output_dir", output_dir)):
                if label in line:
                    line = "  {}: {}\n".format(label, argument)
                line = line.replace("\n", "")
            yaml_lines.append(line)

        self.write_lines_to_file(
            lines=yaml_lines,
            path=outpath
        )

    def generate_step2_files(self, output_dir, log_dir, snakefile, configfile):
        self.write_demultiplexing_configfile(
            output_dir=output_dir,
            outpath=configfile
        )

        self.write_dry_run_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir
        )

        self.write_run_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir,
            log_dir=log_dir,
            jobs=4
        )

    def write_demultiplexing_configfile(self, output_dir, outpath):
        yaml_lines = []
        for line in open(os.path.join(self.demultiplexing_dir, "sceQTL-Gen_Demultiplex.yaml"), 'r'):
            for label, argument in (("ref_dir", self.ref_dir),
                                    ("singularity_image", os.path.join(self.demultiplexing_dir, "WG1-pipeline-QC_wgpipeline.sif")),
                                    ("bind_path", self.bind_paths),
                                    ("samplesheet_filepath", self.samplesheet_filepath),
                                    ("scRNAseq_dir", self.scRNAseq_dir),
                                    ("snp_genotypes_filepath", self.snp_genotypes_filepath),
                                    ("individual_list_dir", self.individual_list_dir),
                                    ("output_dir", output_dir)):
                if label in line:
                    line = "  {}: {}\n".format(label, argument)
                line = line.replace("\n", "")
            yaml_lines.append(line)

        self.write_lines_to_file(
            lines=yaml_lines,
            path=outpath
        )

    def write_dry_run_script(self, snakefile, configfile, output_dir, cores=1):
        lines = [
            '#!/bin/bash',
            '',
            'snakemake \\',
            '    --snakefile {} \\'.format(snakefile),
            '    --configfile {} \\'.format(configfile),
            '    --dryrun \\',
            '    --cores {} \\'.format(cores),
            '    --reason'
        ]
        self.write_lines_to_file(
            lines=lines,
            path=os.path.join(output_dir, "dry_run.sh")
        )

    def write_run_script(self, snakefile, configfile, output_dir, log_dir,
                         jobs=1, restart_times=2, latency_wait=180, nodes=1,
                         time="short", outfile="run"):
        if time not in self.time_dict.keys():
            print("Error, time not recognized.")
            exit()

        lines = [
            '#!/bin/bash',
            '',
            'nohup \\',
            '  snakemake \\',
            '    --snakefile {} \\'.format(snakefile),
            '    --configfile {} \\'.format(configfile),
            '    --rerun-incomplete \\',
            '    --jobs {} \\'.format(jobs),
            '    --use-singularity \\',
            '    --restart-times {} \\'.format(restart_times),
            '    --keep-going \\',
            '    --latency-wait {} \\'.format(latency_wait),
            '    --cluster \\',
            '       "sbatch \\',
            '       --qos regular \\',
            '       -N {threads} \\',
            '       --nodes={} \\'.format(nodes),
            '       --mem={resources.mem_per_thread_gb}G \\',
            '       --tmp={resources.disk_per_thread_gb}G \\',
            '       -o {} \\'.format(os.path.join(log_dir, '{rule}.out')),
            '       --export ALL \\',
            '       --time={}" \\'.format(self.time_dict[time]),
            '    > {}/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &'.format(log_dir)
        ]
        self.write_lines_to_file(
            lines=lines,
            path=os.path.join(output_dir, "{}.sh".format(outfile))
        )

    def write_report_script(self, snakefile, configfile, output_dir):
        lines = [
            '#!/bin/bash',
            '',
            'nohup \\',
            '  snakemake \\',
            '    --snakefile {} \\'.format(snakefile),
            '    --configfile {} \\'.format(configfile),
            '    --report {} \\'.format(output_dir, 'imputation_report.html')
        ]
        self.write_lines_to_file(
            lines=lines,
            path=os.path.join(output_dir, "report.sh")
        )

    @staticmethod
    def write_lines_to_file(lines, path):
        with open(path, "w") as f:
            for line in lines:
                f.write(line + "\n")
        f.close()
        print("\tSaved file: {}".format(os.path.basename(path)))

    def print_arguments(self):
        print("General arguments")
        print("  > Imputation directory: {}".format(self.imputation_dir))
        print("  > Imputation snakefile: {}".format(self.imputation_snakefile))
        print("  > Demultiplexing directory: {}".format(self.demultiplexing_dir))
        print("  > Demultiplexing snakefile: {}".format(self.demultiplexing_snakefile))
        print("  > Reference directory: {}".format(self.ref_dir))
        print("")

        print("[Step 1] Imputation arguments:")
        print("  > Output directory: {}".format(self.dataset_imputation_output_dir))
        print("  > Log directory: {}".format(self.dataset_imputation_logdir))
        print("  > Plink directory: {}".format(self.plink_dir))
        print("  > Configuration file: {}".format(self.dataset_imputation_configfile))
        print("")

        print("[Step2] Demultiplexing and doublet removal arguments:")
        print("  > Output directory: {}".format(self.dataset_demultiplexing_output_dir))
        print("  > Log directory: {}".format(self.dataset_demultiplexing_logdir))
        print("  > Samplesheet file: {}".format(self.samplesheet_filepath))
        print("  > scRNA-seq directory: {}".format(self.scRNAseq_dir))
        print("  > SNP-genotypes file: {}".format(self.snp_genotypes_filepath))
        print("  > Individual list directory: {}".format(self.individual_list_dir))
        print("  > Configuration file: {}".format(self.dataset_demultiplexing_configfile))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()