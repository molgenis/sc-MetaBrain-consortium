#!/usr/bin/env python3

"""
File:         initialize_wg1_files.py
Created:      2022/10/07
Last Changed: 2023/02/02
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
import argparse
import os
import re

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
    --ref_dir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/hg38 \
    --dataset_outdir 2022-10-13-ImputationTestDataset \
    --imputation_subdir 2022-10-07-Imputation \
    --plink_dir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2022-10-07-Imputation/ImputationTestDataset_plink \
    --demultiplexing_subdir 2022-10-10-DemultiplexingAndDoubletRemoval \
    --samplesheet_filepath /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2022-10-10-DemultiplexingAndDoubletRemoval/TestData4PipelineSmall/samplesheet.txt \
    --scRNAseq_dir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2022-10-10-DemultiplexingAndDoubletRemoval/TestData4PipelineSmall \
    --snp_genotypes_filepath /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2022-10-10-DemultiplexingAndDoubletRemoval/TestData4PipelineSmall/test_dataset.vcf \
    --individual_list_dir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2022-10-10-DemultiplexingAndDoubletRemoval/TestData4PipelineSmall/individuals_list_dir
    
### AMP-AD ###
./initialize_wg1_files.py \
    --step 1 \
    --work_dir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC \
    --ref_dir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/hg38 \
    --dataset_outdir 2023-02-02-AMP_AD \
    --imputation_subdir 2022-10-07-Imputation \
    --plink_dir /groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/2022-11-03-FilteredGenotypes/5-plink2_makepgen_after_fill_all_eur
    
./initialize_wg1_files.py \
    --step 2 \
    --work_dir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC \
    --ref_dir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/hg38 \
    --dataset_outdir 2023-02-02-AMP_AD \
    --demultiplexing_subdir 2022-10-10-DemultiplexingAndDoubletRemoval \
    --demultiplexing_singularity /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2022-10-10-DemultiplexingAndDoubletRemoval/WG1-pipeline-QC_wgpipeline.simg \
    --samplesheet_filepath /groups/umcg-biogen/tmp01/input/processeddata/single-cell/Mathys2019/samplesheet.txt \
    --scRNAseq_dir /groups/umcg-biogen/tmp01/input/processeddata/single-cell/Mathys2019 \
    --snp_genotypes_filepath /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2023-02-02-AMP_AD/Step1-Imputation/vcf_all_merged/Mathys/Mathys_imputed_hg38_R2_0.3_MAF0.05.vcf \
    --individual_list_dir /groups/umcg-biogen/tmp01/input/processeddata/single-cell/Mathys2019/individual_list_dir
    
### Roche ###
./initialize_wg1_files.py \
    --step 1 \
    --work_dir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-08-21-WorkGroup1QC \
    --ref_dir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-08-21-WorkGroup1QC/hg38 \
    --dataset_outdir 2023-08-23-Roche \
    --imputation_subdir 2023-08-21-Imputation \
    --plink_dir 
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.step = getattr(arguments, 'step')
        self.work_dir = getattr(arguments, 'work_dir')
        self.ref_dir = getattr(arguments, 'ref_dir')
        self.bind_paths = ",".join(getattr(arguments, 'bind_paths'))
        dataset_outdir = getattr(arguments, 'dataset_outdir')

        # Pre-process the dataset output directory.
        date_str = datetime.now().strftime("%Y-%m-%d")
        if not re.match("^[0-9]{4}-(0[1-9]|1[0-2])-(0[1-9]|[1-2][0-9]|3[0-1])", dataset_outdir):
            dataset_outdir = "{}-{}".format(date_str, dataset_outdir)
        self.dataset_outdir = dataset_outdir

        # Step 1 arguments.
        self.imputation_subdir = getattr(arguments, 'imputation_subdir')
        self.imputation_singularity = getattr(arguments, 'imputation_singularity')
        self.imputation_config = getattr(arguments, 'imputation_config')
        self.plink_dir = getattr(arguments, 'plink_dir')

        if self.step is None or self.step == 1:
            for label, value in [("--imputation_subdir", self.imputation_subdir),
                                 ("--plink_dir", self.plink_dir)]:
                if value is None:
                    print("Argument {} is required when --step equals {}.".format(label, self.step))
                    exit()

        # Step 2 arguments.
        self.demultiplexing_subdir = getattr(arguments, 'demultiplexing_subdir')
        self.demultiplexing_singularity = getattr(arguments, 'demultiplexing_singularity')
        self.demultiplexing_config = getattr(arguments, 'demultiplexing_config')
        self.samplesheet_filepath = getattr(arguments, 'samplesheet_filepath')
        self.scRNAseq_dir = getattr(arguments, 'scRNAseq_dir')
        self.snp_genotypes_filepath = getattr(arguments, 'snp_genotypes_filepath')
        self.individual_list_dir = getattr(arguments, 'individual_list_dir')

        if self.step is None or self.step == 2:
            for label, value in [("--demultiplexing_subdir", self.demultiplexing_subdir),
                                 ("--samplesheet_filepath", self.samplesheet_filepath),
                                 ("--scRNAseq_dir", self.scRNAseq_dir),
                                 ("--snp_genotypes_filepath", self.snp_genotypes_filepath),
                                 ("--individual_list_dir", self.individual_list_dir)]:
                if value is None:
                    print("Argument {} is required when --step equals {}.".format(label, self.step))
                    exit()

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # general arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit")
        parser.add_argument("--step",
                            type=int,
                            choices=[1, 2],
                            default=None,
                            help="Which step to prepare. Default: None.")
        parser.add_argument("--work_dir",
                            type=str,
                            required=True,
                            default=None,
                            help="The working directory.")
        parser.add_argument("--ref_dir",
                            type=str,
                            required=True,
                            default=None,
                            help="This is the path to the directory containing "
                                 "the imputation references provided for the "
                                 "SNP processing and imputation")
        parser.add_argument("--bind_paths",
                            nargs="*",
                            type=str,
                            default=["/groups/umcg-biogen/tmp01/",
                                     "/groups/umcg-biogen/tmp01/umcg-mvochteloo/simulated_home:/home/umcg-mvochteloo"],
                            help="List of paths to bind to Singularity. "
                                 "Default: ['/groups/umcg-biogen/tmp01/']")
        parser.add_argument("--dataset_outdir",
                            type=str,
                            required=True,
                            default=None,
                            help="The name of the output directory where you "
                                 "would like all outputs/results saved.")

        # arguments for step 1.
        parser.add_argument("--imputation_subdir",
                            type=str,
                            default=None,
                            help="The imputation step subdirectory.")
        parser.add_argument("--imputation_singularity",
                            type=str,
                            default="WG1-pipeline-QC_imputation.sif",
                            help="")
        parser.add_argument("--imputation_config",
                            type=str,
                            default="PreImputation.yaml",
                            help="")
        parser.add_argument("--plink_dir",
                            type=str,
                            default=None,
                            help="Absolute path to directory that contains "
                                 "plink2 pgen, pvar and psam files for "
                                 "pipeline on hg19.")

        # argument for step 2.
        parser.add_argument("--demultiplexing_subdir",
                            type=str,
                            default=None,
                            help="The demultiplexing and doublet removal step "
                                 "subdirectory.")
        parser.add_argument("--demultiplexing_singularity",
                            type=str,
                            default="WG1-pipeline-QC_wgpipeline.simg",
                            help="The complete path to the singularity image "
                                 "that has all the softwares")
        parser.add_argument("--demultiplexing_config",
                            type=str,
                            default="sceQTL-Gen_Demultiplex.yaml",
                            help="")
        parser.add_argument("--samplesheet_filepath",
                            type=str,
                            default=None,
                            help="Tab separated file that has a header. Each "
                                 "line has a pool name used for the scRNA-seq "
                                 "directories and the number of individuals in "
                                 "each pool.")
        parser.add_argument("--scRNAseq_dir",
                            type=str,
                            default=None,
                            help="The parent directory that has directories "
                                 "for each pool and the scRNA-seq output "
                                 "below it.")
        parser.add_argument("--snp_genotypes_filepath",
                            type=str,
                            default=None,
                            help="The path to the genotype file that has just "
                                 "SNPs that are within exons and with a minor"
                                 " allele frequency > 0.05 (5%) in this "
                                 "sample.")
        parser.add_argument("--individual_list_dir",
                            type=str,
                            default=None,
                            help="Directory that has a different file for "
                                 "each pool. Each file contains a list of "
                                 "the individuals IDs that are in the pool "
                                 "separated by line (match the genotype "
                                 "individual IDs).")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        dataset_work_dir = os.path.join(self.work_dir, self.dataset_outdir)
        if not os.path.exists(dataset_work_dir):
            os.makedirs(dataset_work_dir)

        if self.step is None or self.step == 1:
            print("Generating step 1 files")

            output_dir = os.path.join(dataset_work_dir, "Step1-Imputation")
            log_dir = os.path.join(output_dir, "log")
            for outdir in [output_dir, log_dir]:
                if not os.path.exists(outdir):
                    os.makedirs(outdir)

            self.generate_step1_files(
                output_dir=output_dir,
                log_dir=log_dir,
                singularity=os.path.join(self.work_dir, self.imputation_subdir, self.imputation_singularity),
                snakefile=os.path.join(self.work_dir, self.imputation_subdir, "Snakefile"),
                configtemplate=os.path.join(self.work_dir, self.imputation_subdir, self.imputation_config),
                configfile=os.path.join(output_dir, self.imputation_config)
            )

            print("")

        if self.step is None or self.step == 2:
            print("Generating step 2 files")

            output_dir = os.path.join(dataset_work_dir, "Step2-DemultiplexingAndDoubletRemoval")
            log_dir = os.path.join(output_dir, "log")
            for outdir in [output_dir, log_dir]:
                if not os.path.exists(outdir):
                    os.makedirs(outdir)

            self.generate_step2_files(
                output_dir=output_dir,
                log_dir=log_dir,
                singularity=os.path.join(self.work_dir, self.demultiplexing_subdir, self.demultiplexing_singularity),
                snakefile=os.path.join(self.work_dir, self.demultiplexing_subdir, "Snakefile"),
                configtemplate=os.path.join(self.work_dir, self.demultiplexing_subdir, self.demultiplexing_config),
                configfile=os.path.join(output_dir, self.demultiplexing_config)
            )

            print("")

    def generate_step1_files(self, output_dir, log_dir, singularity, snakefile,
                             configtemplate, configfile):
        """
        Step 1) run 'dry_run.sh'
        Step 2) run 'build_dag1.sh'
        Step 3) run 'run_20jobs_short.sh'
        Step 4) check 'pca_sex_checks/check_sex_update_remove.tsv'
        Step 5) check 'pca_sex_checks/ancestry_update_remove.tsv'
        Step 6) run 'dry_run.sh'
        Step 7) run 'build_dag2.sh'
        Step 8) run 'run_22jobs_until_prehasing.sh'
        Step 9) run 'run_22jobs_until_combine.sh'
        Step 10) run 'run_22jobs_short.sh'
        Step 11) run 'report.sh'
        """
        config_arguments = (
            ("ref_dir", self.ref_dir),
            ("singularity_image", singularity),
            ("plink_dir", self.plink_dir),
            ("bind_paths", self.bind_paths),
            ("output_dir", output_dir)
        )
        self.write_configfile(
            template=configtemplate,
            arguments=config_arguments,
            outpath=configfile
        )

        self.write_dry_run_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir
        )

        self.write_unlock_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir
        )

        for i in range(1, 3):
            self.write_build_dag_script(
                snakefile=snakefile,
                configfile=configfile,
                output_dir=output_dir,
                outfile="dag{}".format(i)
            )

        for time in ["short", "medium", "long"]:
            self.write_run_script(
                snakefile=snakefile,
                configfile=configfile,
                output_dir=output_dir,
                log_dir=log_dir,
                jobs=22,
                time=time,
                outfile="run_22jobs_{}".format(time)
            )

        self.write_run_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir,
            log_dir=log_dir,
            until="eagle_prephasing",
            jobs=22,
            time="short",
            outfile="run_22jobs_until_prehasing"
        )

        self.write_run_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir,
            log_dir=log_dir,
            until="combine_vcfs_ancestry",
            jobs=22,
            time="long",
            outfile="run_22jobs_until_combine"
        )

        self.write_report_script(
            snakefile=snakefile,
            configfile=configfile,
            report=os.path.join(output_dir, "imputation_report.html"),
            log_dir=log_dir,
            output_dir=output_dir
        )

    def generate_step2_files(self, output_dir, log_dir, singularity, snakefile,
                             configtemplate, configfile):
        """

        Step 1) run 'dry_run.sh'
        Step 2) run 'build_dag1.sh'
        Step 3) run 'run_22jobs_short.sh'
        Step 4) run 'dry_run.sh'
        Step 5) run 'build_dag2.sh'
        Step 6) run 'run_22jobs_short.sh'
        Step 7) check 'manual_selections/DoubletDetection/DoubletDetection_manual_selection.tsv'
        Step 8) check 'manual_selections/scrublet/scrublet_percentile_manual_selection.tsv'
        Step 9) run 'dry_run.sh'
        Step 10) run 'build_dag3.sh'
        Step 11) run 'run_22jobs_short.sh'
        Step 12) run 'report.sh'
        """
        config_arguments = (
            ("ref_dir", self.ref_dir),
            ("singularity_image", singularity),
            ("bind_path", self.bind_paths),
            ("samplesheet_filepath", self.samplesheet_filepath),
            ("scRNAseq_dir", self.scRNAseq_dir),
            ("snp_genotypes_filepath", self.snp_genotypes_filepath),
            ("individual_list_dir", self.individual_list_dir),
            ("output_dir", output_dir)
        )
        self.write_configfile(
            template=configtemplate,
            arguments=config_arguments,
            outpath=configfile
        )

        self.write_dry_run_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir
        )

        self.write_unlock_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir
        )

        for i in range(1, 4):
            self.write_build_dag_script(
                snakefile=snakefile,
                configfile=configfile,
                output_dir=output_dir,
                outfile="dag{}".format(i)
            )

        for time in ["short", "medium", "long"]:
            self.write_run_script(
                snakefile=snakefile,
                configfile=configfile,
                output_dir=output_dir,
                log_dir=log_dir,
                jobs=22,
                time=time,
                outfile="run_20jobs_{}".format(time)
            )

        self.write_report_script(
            snakefile=snakefile,
            configfile=configfile,
            report=os.path.join(output_dir, "demultiplexing_report.html"),
            log_dir=log_dir,
            output_dir=output_dir
        )

    def write_configfile(self, template, arguments, outpath):
        yaml_lines = []
        for line in open(template, 'r'):
            line = line.replace("\n", "")
            for label, argument in arguments:
                if line.statswith("  {}".format(label)):
                    line = "  {}: {}".format(label, argument)
            yaml_lines.append(line)

        self.write_lines_to_file(
            lines=yaml_lines,
            path=outpath
        )

    def write_build_dag_script(self, snakefile, configfile, output_dir,
                               outfile="dag"):
        lines = [
            '#!/bin/bash',
            '',
            'snakemake \\',
            '    --snakefile {} \\'.format(snakefile),
            '    --configfile {} \\'.format(configfile),
            '    --dag | \\',
            '    dot -Tsvg \\',
            '        > {}.svg'.format(outfile)
        ]
        self.write_lines_to_file(
            lines=lines,
            path=os.path.join(output_dir, "build_{}.sh".format(outfile))
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

    def write_unlock_script(self, snakefile, configfile, output_dir):
        lines = [
            '#!/bin/bash',
            '',
            'snakemake \\',
            '    --snakefile {} \\'.format(snakefile),
            '    --configfile {} \\'.format(configfile),
            '    --unlock'
        ]
        self.write_lines_to_file(
            lines=lines,
            path=os.path.join(output_dir, "unlock.sh")
        )

    def write_run_script(self, snakefile, configfile, output_dir, log_dir,
                         jobs=1, restart_times=2, latency_wait=30, nodes=1,
                         time="short", until=None, outfile="run"):
        time_dict = {
            "short": "05:59:00",
            "medium": "23:59:00",
            "long": "6-23:59:00"
        }

        if time not in time_dict.keys():
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
            '       --nodes={} \\'.format(nodes),
            '       --cpus-per-task={threads} \\',
            '       --mem=\$(({resources.mem_per_thread_gb} * {threads}))G \\',
            '       --tmp=\$(({resources.disk_per_thread_gb} * {threads}))G \\',
            '       --time={} \\'.format(time_dict[time]),
            '       --output={} \\'.format(os.path.join(log_dir, '{rule}.out')),
            '       --error={} \\'.format(os.path.join(log_dir, '{rule}.out')),
            '       --export=ALL \\',
            '       --qos=regular" \\',
            '    > {}/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &'.format(log_dir),
            '',
            'echo "Check status of command with:" ps -p $! -u'
        ]

        if until is not None:
            lines.insert(6, '    --until {} \\'.format(until))

        self.write_lines_to_file(
            lines=lines,
            path=os.path.join(output_dir, "{}.sh".format(outfile))
        )

    def write_report_script(self, snakefile, configfile, report, log_dir,
                            output_dir):
        lines = [
            '#!/bin/bash',
            '',
            'nohup \\',
            '  snakemake \\',
            '    --snakefile {} \\'.format(snakefile),
            '    --configfile {} \\'.format(configfile),
            '    --report {} \\'.format(report),
            '    > {}/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &'.format(log_dir),
            '',
            'echo "Check status of command with:" ps -p $! -u'
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
        print("  > Step:                         {}".format("1 & 2" if self.step is None else str(self.step)))
        print("  > Working directory:            {}".format(self.work_dir))
        print("  > Reference directory:          {}".format(self.ref_dir))
        print("  > Bind paths:                   {}".format(self.bind_paths))
        print("  > Dataset output directory:     {}".format(self.dataset_outdir))
        print("")

        if self.step is None or self.step == 1:
            print("[Step 1] Imputation arguments:")
            print("  > Imputation subdirectory:      {}".format(self.imputation_subdir))
            print("  > Imputation singularity:       {}".format(self.imputation_singularity))
            print("  > Imputation configuration:     {}".format(self.imputation_config))
            print("  > Plink directory:              {}".format(self.plink_dir))
            print("")

        if self.step is None or self.step == 2:
            print("[Step2] Demultiplexing and doublet removal arguments:")
            print("  > Demultiplexing subdirectory:  {}".format(self.demultiplexing_subdir))
            print("  > Demultiplexing singularity:   {}".format(self.demultiplexing_singularity))
            print("  > Demultiplexing configuration: {}".format(self.demultiplexing_config))
            print("  > Samplesheet file:             {}".format(self.samplesheet_filepath))
            print("  > scRNA-seq directory:          {}".format(self.scRNAseq_dir))
            print("  > SNP-genotypes file:           {}".format(self.snp_genotypes_filepath))
            print("  > Individual list directory:    {}".format(self.individual_list_dir))
            print("")


if __name__ == '__main__':
    m = main()
    m.start()