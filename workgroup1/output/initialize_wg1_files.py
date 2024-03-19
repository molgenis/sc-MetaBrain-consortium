#!/usr/bin/env python3

"""
File:         initialize_wg1_files.py
Created:      2022/10/07
Last Changed: 2023/10/05
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
./initialize_wg1_files.py -h
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.step = getattr(arguments, 'step')
        self.work_dir = getattr(arguments, 'work_dir')
        self.ref_dir = getattr(arguments, 'ref_dir')
        self.bind_path = getattr(arguments, 'bind_path')
        self.sif_path = getattr(arguments, 'sif_path')
        dataset_outdir = getattr(arguments, 'dataset_outdir')
        self.repo_dir = getattr(arguments, 'repo_dir')
        self.exclude_temp_in_sbatch = getattr(arguments, 'exclude_temp_in_sbatch')

        # Pre-process the dataset output directory.
        date_str = datetime.now().strftime("%Y-%m-%d")
        if not re.match("^[0-9]{4}-(0[1-9]|1[0-2])-(0[1-9]|[1-2][0-9]|3[0-1])", dataset_outdir):
            dataset_outdir = "{}-{}".format(date_str, dataset_outdir)
        self.dataset_outdir = dataset_outdir

        # Step 1 arguments.
        self.imputation_config = getattr(arguments, 'imputation_config')
        self.genotype_path = getattr(arguments, 'genotype_path')
        self.psam = getattr(arguments, 'psam')
        self.genome_build = getattr(arguments, 'genome_build')
        self.dataset_samples = getattr(arguments, 'dataset_samples')
        if self.dataset_samples != "{}":
            self.dataset_samples = "{" + ", ".join([x.replace(":", ": ") for x in self.dataset_samples]) + "}"
        self.is_wgs = getattr(arguments, 'is_wgs')
        self.force_ancestry = getattr(arguments, 'force_ancestry')

        if self.step is None or self.step == 1:
            for label, value in [("--imputation_config", self.imputation_config),
                                 ("--genotype_path", self.genotype_path),
                                 ("--psam", self.psam),
                                 ("--genome_build", self.genome_build),
                                 ("--dataset_samples", self.dataset_samples),
                                 ("--is_wgs", self.is_wgs),
                                 ("--force_ancestry", self.force_ancestry)]:
                if value is None:
                    print("Argument {} is required when --step equals {}.".format(label, self.step))
                    exit()

        # Step 2 arguments.
        self.demultiplexing_config = getattr(arguments, 'demultiplexing_config')
        self.poolsheet_path = getattr(arguments, 'poolsheet_path')
        self.individual_list_dir = getattr(arguments, 'individual_list_dir')
        self.individual_coupling = getattr(arguments, 'individual_coupling')
        if self.individual_list_dir is None:
            self.individual_list_dir = ""
        self.vcf = getattr(arguments, 'vcf')
        if self.vcf is None:
            self.vcf = "[]"
        self.is_multiplexed = getattr(arguments, 'is_multiplexed')
        self.sc_data_type = getattr(arguments, 'sc_data_type')
        self.rename_bam = getattr(arguments, 'rename_bam')
        self.check_sample_swaps = getattr(arguments, 'no_sample_swap_check')

        if self.step is None or self.step == 2:
            for label, value in [("--poolsheet_path", self.poolsheet_path),
                                 ("--is_multiplexed", self.is_multiplexed),
                                 ("--sc_data_type", self.sc_data_type),
                                 ("--rename_bam", self.rename_bam),
                                 ("--check_sample_swaps", self.check_sample_swaps)]:
                if value is None:
                    print("Argument {} is required when --step equals {}.".format(label, self.step))
                    exit()
            if self.is_multiplexed:
                for label, value in [("--individual_list_dir", self.individual_list_dir),
                                     ("--vcf", self.vcf)]:
                    if value is None:
                        print("Argument {} is required when --is_multiplexed equals True.".format(label))
                        exit()
            if not self.is_multiplexed:
                for label, value in [("--individual_coupling", self.individual_coupling)]:
                    if value is None:
                        print("Argument {} is required when --is_multiplexed equals False.".format(label))
                        # exit()

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
        parser.add_argument("--bind_path",
                            type=str,
                            required=True,
                            default=None,
                            help="List of paths to bind to Singularity.")
        parser.add_argument("--sif_path",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--dataset_outdir",
                            type=str,
                            required=True,
                            default=None,
                            help="The name of the output directory where you "
                                 "would like all outputs/results saved.")
        parser.add_argument("--repo_dir",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--exclude_temp_in_sbatch",
                            action='store_true',
                            help="")

        # arguments for step 1.
        parser.add_argument("--imputation_config",
                            type=str,
                            default="sceQTL-Gen_Imputation.yaml",
                            help="")
        parser.add_argument("--genotype_path",
                            type=str,
                            default=None,
                            help="")
        parser.add_argument("--genome_build",
                            type=str,
                            default=None,
                            choices=['hg18', 'b36', 'hg19', 'b37', 'hg38', 'b38'],
                            help="")
        parser.add_argument("--psam",
                            type=str,
                            default="",
                            help="")
        parser.add_argument("--dataset_samples",
                            nargs="*",
                            type=str,
                            default="{}",
                            help="")
        parser.add_argument("--is_wgs",
                            action='store_true',
                            help="")
        parser.add_argument("--force_ancestry",
                            type=str,
                            default="null",
                            help="")
        parser.add_argument("--exclude_y_chr_in_sex_check",
                            action='store_true',
                            help="")

        # argument for step 2.
        parser.add_argument("--demultiplexing_config",
                            type=str,
                            default="sceQTL-Gen_Demultiplex.yaml",
                            help="")
        parser.add_argument("--poolsheet_path",
                            type=str,
                            default=None,
                            help="Tab separated file that has a header. Each "
                                 "line has a pool name used for the scRNA-seq "
                                 "directories and the number of individuals in "
                                 "each pool.")
        parser.add_argument("--individual_list_dir",
                            type=str,
                            default=None,
                            help="Directory that has a different file for "
                                 "each pool. Each file contains a list of "
                                 "the individuals IDs that are in the pool "
                                 "separated by line (match the genotype "
                                 "individual IDs).")
        parser.add_argument("--individual_coupling",
                            type=str,
                            default="",
                            help="")
        parser.add_argument("--vcf",
                            nargs="*",
                            type=str,
                            default=None,
                            help="")
        parser.add_argument("--is_multiplexed",
                            action='store_true',
                            default=False,
                            help="")
        parser.add_argument("--sc_data_type",
                            type=str,
                            choices=["single-cell", "single-nucleus"],
                            default="",
                            help="")
        parser.add_argument("--rename_bam",
                            action='store_true',
                            default=False,
                            help="")
        parser.add_argument("--no_sample_swap_check",
                            action='store_false',
                            default=True,
                            help="")

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
            slurm_log_dir = os.path.join(output_dir, "slurm_log")
            for outdir in [output_dir, log_dir, slurm_log_dir]:
                if not os.path.exists(outdir):
                    os.makedirs(outdir)

            self.generate_step1_files(
                repo_dir=os.path.join(self.repo_dir),
                output_dir=output_dir,
                log_dir=log_dir,
                slurm_log_dir=slurm_log_dir,
                snakefile=os.path.join(self.repo_dir, "Imputation", "Snakefile"),
                configtemplate=os.path.join(self.repo_dir, "Imputation", self.imputation_config),
                configfile=os.path.join(output_dir, self.imputation_config)
            )

            print("")

        if self.step is None or self.step == 2:
            print("Generating step 2 files")

            output_dir = os.path.join(dataset_work_dir, "Step2-DemultiplexingAndDoubletRemoval")
            log_dir = os.path.join(output_dir, "log")
            slurm_log_dir = os.path.join(output_dir, "slurm_log")
            for outdir in [output_dir, log_dir, slurm_log_dir]:
                if not os.path.exists(outdir):
                    os.makedirs(outdir)

            self.generate_step2_files(
                repo_dir=os.path.join(self.repo_dir),
                output_dir=output_dir,
                log_dir=log_dir,
                slurm_log_dir=slurm_log_dir,
                snakefile=os.path.join(self.repo_dir, "Demultiplexing", "Snakefile"),
                configtemplate=os.path.join(self.repo_dir, "Demultiplexing", self.demultiplexing_config),
                configfile=os.path.join(output_dir, self.demultiplexing_config)
            )

            print("")

    def generate_step1_files(self, repo_dir, output_dir, log_dir, slurm_log_dir, snakefile, configtemplate, configfile):
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
            ("bind_path", self.bind_path),
            ("singularity_image", self.sif_path),
            ("repo_dir", repo_dir),
            ("genotype_path", self.genotype_path),
            ("psam", self.psam),
            ("genome_build", self.genome_build),
            ("dataset_samples", self.dataset_samples),
            ("ref_dir", self.ref_dir),
            ("is_wgs", self.is_wgs),
            ("force_ancestry", self.force_ancestry),
            ("output_dir", output_dir)
        )
        self.write_configfile(
            template=configtemplate,
            arguments=config_arguments,
            outpath=configfile
        )

        cluster_status_script = self.write_cluster_status_script(
            output_dir=output_dir
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

        for jobs in [1, 3, 24, 100, 1000, 10000]:
            self.write_run_script(
                snakefile=snakefile,
                configfile=configfile,
                output_dir=output_dir,
                log_dir=log_dir,
                slurm_log_dir=slurm_log_dir,
                cluster_status_script=cluster_status_script,
                jobs=jobs,
                outfile="run_{}jobs".format(jobs)
            )

        self.write_run_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir,
            log_dir=log_dir,
            slurm_log_dir=slurm_log_dir,
            cluster_status_script=cluster_status_script,
            jobs=3,
            restart_times=0,
            qos="priority",
            outfile="run_dev"
        )

        self.write_run_local_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir
        )

        self.write_report_script(
            snakefile=snakefile,
            configfile=configfile,
            report=os.path.join(output_dir, "imputation_report.html"),
            output_dir=output_dir
        )

    def generate_step2_files(self, repo_dir, output_dir, log_dir, slurm_log_dir, snakefile, configtemplate, configfile):
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
            ("bind_path", self.bind_path),
            ("singularity_image", self.sif_path),
            ("repo_dir", repo_dir),
            ("ref_dir", self.ref_dir),
            ("poolsheet_path", self.poolsheet_path),
            ("individual_list_dir", self.individual_list_dir),
            ("individual_coupling", self.individual_coupling),
            ("vcf", self.vcf),
            ("output_dir", output_dir),
            ("is_multiplexed", self.is_multiplexed),
            ("sc_data_type", self.sc_data_type),
            ("rename_bam", self.rename_bam),
            ("check_sample_swaps", self.check_sample_swaps)
        )
        self.write_configfile(
            template=configtemplate,
            arguments=config_arguments,
            outpath=configfile
        )

        self.write_check_bam_script(
            output_dir=output_dir
        )

        cluster_status_script = self.write_cluster_status_script(
            output_dir=output_dir
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

        for jobs in [1, 2, 22, 100, 1000, 10000]:
            self.write_run_script(
                snakefile=snakefile,
                configfile=configfile,
                output_dir=output_dir,
                log_dir=log_dir,
                slurm_log_dir=slurm_log_dir,
                cluster_status_script=cluster_status_script,
                jobs=jobs,
                outfile="run_{}jobs".format(jobs)
            )

        self.write_run_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir,
            log_dir=log_dir,
            slurm_log_dir=slurm_log_dir,
            cluster_status_script=cluster_status_script,
            jobs=3,
            restart_times=0,
            qos="priority",
            outfile="run_dev"
        )

        self.write_run_local_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir
        )

        self.write_report_script(
            snakefile=snakefile,
            configfile=configfile,
            report=os.path.join(output_dir, "demultiplexing_report.html"),
            output_dir=output_dir
        )

    def write_configfile(self, template, arguments, outpath):
        yaml_lines = []
        for line in open(template, 'r'):
            line = line.replace("\n", "")
            for label, argument in arguments:
                if line.startswith("  {}:".format(label)):
                    line = "  {}: {}".format(label, argument)
            yaml_lines.append(line)

        self.write_lines_to_file(
            lines=yaml_lines,
            path=outpath
        )

    def write_cluster_status_script(self, output_dir):
        """
        https://www.embl.org/groups/bioinformatics-rome/blog/2022/05/snakemake-profile-5-handling-memory-and-timeout-errors/
        """
        outfile = os.path.join(output_dir, "status-sacct.sh")

        lines = [
            '#!/usr/bin/env bash',
            '# Check status of Slurm job',
            'jobid="$1"',
            'if [[ "$jobid" == Submitted ]]',
            'then',
            '  echo smk-simple-slurm: Invalid job ID: "$jobid" >&2',
            '  echo smk-simple-slurm: Did you remember to add the flag --parsable to your sbatch call? >&2',
            '  exit 1',
            'fi',
            'output=`sacct -j "$jobid" --format State --noheader | head -n 1 | awk \'{print $1}\'`',
            'if [[ $output =~ ^(COMPLETED).* ]]',
            'then',
            '  echo success',
            'elif [[ $output =~ ^(RUNNING|PENDING|COMPLETING|CONFIGURING|SUSPENDED).* ]]',
            'then',
            '  echo running',
            'else',
            '  echo failed',
            'fi']

        self.write_lines_to_file(
            lines=lines,
            path=outfile
        )

        return outfile

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
            '        > {}.svg'.format(os.path.join(output_dir, outfile))
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
                         slurm_log_dir, cluster_status_script, nodes=1,
                         jobs=1, restart_times=0, latency_wait=60,
                         qos="regular", until=None, outfile="run"):
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
            '       --job-name=snakemake_{rule}_{wildcards} \\',
            '       --nodes={} \\'.format(nodes),
            '       --cpus-per-task={threads} \\',
            '       --mem=\$(({resources.mem_per_thread_gb} * {threads}))G \\',
            '       --tmp=\$(({resources.disk_per_thread_gb} * {threads}))G \\',
            '       --time={resources.time} \\',
            '       --output={}'.format(slurm_log_dir) + '/{rule}_{wildcards}.out \\',
            '       --error={}'.format(slurm_log_dir) + '/{rule}_{wildcards}.out \\',
            '       --export=NONE \\',
            '       --qos={} \\'.format(qos),
            '       --parsable" \\',
            '    --cluster-status {} \\'.format(cluster_status_script),
            '    > {}/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &'.format(log_dir),
            '',
            'echo "Check status of command with:" ps -p $! -u'
        ]

        if self.exclude_temp_in_sbatch:
            lines.remove('       --tmp=\$(({resources.disk_per_thread_gb} * {threads}))G \\')

        if until is not None:
            lines.insert(6, '    --until {} \\'.format(until))

        self.write_lines_to_file(
            lines=lines,
            path=os.path.join(output_dir, "{}.sh".format(outfile))
        )

    def write_run_local_script(self, snakefile, configfile, output_dir):
        lines = [
            '#!/bin/bash',
            '',
            'snakemake \\',
            '  --snakefile {} \\'.format(snakefile),
            '  --configfile {} \\'.format(configfile),
            '  --cores 1'
        ]
        self.write_lines_to_file(
            lines=lines,
            path=os.path.join(output_dir, "run_local.sh")
        )

    def write_report_script(self, snakefile, configfile, report, output_dir):
        lines = [
            '#!/bin/bash',
            '',
            'snakemake \\',
            '  --snakefile {} \\'.format(snakefile),
            '  --configfile {} \\'.format(configfile),
            '  --report {}'.format(report)
        ]
        self.write_lines_to_file(
            lines=lines,
            path=os.path.join(output_dir, "report.sh")
        )


    def write_check_bam_script(self, output_dir):
        bam_index = None
        bam_files = []
        with open(self.poolsheet_path, 'r') as f:
            for line in f:
                values = line.strip("\n").split("\t")

                if bam_index is None:
                    if "Bam" in values:
                        bam_index = values.index("Bam")
                    continue

                bam_files.append(values[bam_index])
        f.close()

        if len(bam_files) == 0:
            return

        lines = [
            '#!/bin/bash',
            '',
            'declare -a arr=({})'.format(" ".join(bam_files)),
            '',
            'echo "Checking for \'chr\' prefix in BAM files:"',
            '',
            'sum=0',
            'for bam in "${arr[@]}"',
            'do',
            '   n=$(singularity exec --bind {bind} {sif} samtools view -H $bam | grep SN:chr[0-9XY]* | wc -l)'.format(bind=self.bind_path, sif=self.sif_path),
            '   if [[ $n -eq 0 ]]; then echo "    $bam"; else echo "    $bam has prefix"; fi',
            '   sum=$(($sum + $n))',
            'done',
            'if [[ $sum -gt 0 ]]; then echo "Please set \'rename_bam\' to True."; else echo "Please set \'rename_bam\' to False."; fi',
            '',
        ]
        self.write_lines_to_file(
            lines=lines,
            path=os.path.join(output_dir, "check_bam.sh")
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
        print("  > Bind path:                    {}".format(self.bind_path))
        print("  > Singularity path:             {}".format(self.sif_path))
        print("  > Dataset output directory:     {}".format(self.dataset_outdir))
        print("  > Repository directory:         {}".format(self.repo_dir))
        print("  > Exclude TEMP in SBATCH:       {}".format(self.exclude_temp_in_sbatch))
        print("")

        if self.step is None or self.step == 1:
            print("[Step 1] Imputation arguments:")
            print("  > Imputation configuration:     {}".format(self.imputation_config))
            print("  > Genotype path:                {}".format(self.genotype_path))
            print("  > PSAM file:                    {}".format(self.psam))
            print("  > Genome build:                 {}".format(self.genome_build))
            print("  > Dataset samples:              {}".format(self.dataset_samples))
            print("  > Is WGS:                       {}".format(self.is_wgs))
            print("  > Force ancestry:               {}".format(self.force_ancestry))
            print("")

        if self.step is None or self.step == 2:
            print("[Step2] Demultiplexing and doublet removal arguments:")
            print("  > Demultiplexing configuration: {}".format(self.demultiplexing_config))
            print("  > Poolsheet file:               {}".format(self.poolsheet_path))
            print("  > VCF:                          {}".format(", ".join(self.vcf) if isinstance(self.vcf, list) else ""))
            print("  > Individual list directory:    {}".format(self.individual_list_dir))
            print("  > Individual coupling file:     {}".format(self.individual_coupling))
            print("  > Is multiplexed:               {}".format(self.is_multiplexed))
            print("  > Single-cell data type:        {}".format(self.sc_data_type))
            print("")


if __name__ == '__main__':
    m = main()
    m.start()