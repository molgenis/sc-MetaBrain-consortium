#!/usr/bin/env python3

"""
File:         initialize_cd_files.py
Created:      2024/03/25
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
from datetime import datetime
import argparse
import os
import re

# Third party imports.

# Local application imports.

# Metadata
__program__ = "Initialize Combine Datasets Files"
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
./initialize_cd_files.py -h
"""

class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.work_dir = getattr(arguments, 'work_dir')
        dataset_outdir = getattr(arguments, 'dataset_outdir')
        self.comb_datasets_dir = getattr(arguments, 'comb_datasets_dir')
        self.comb_datasets_config = getattr(arguments, 'comb_datasets_config')
        self.bind_path = getattr(arguments, 'bind_path')
        self.sif_path = getattr(arguments, 'sif_path')
        self.repo_dir = getattr(arguments, 'repo_dir')
        self.datasheet_path = getattr(arguments, 'datasheet_path')
        self.variants_path = getattr(arguments, 'variants_path')
        self.exclude_temp_in_sbatch = getattr(arguments, 'exclude_temp_in_sbatch')

        # Pre-process the dataset output directory.
        date_str = datetime.now().strftime("%Y-%m-%d")
        if not re.match("^[0-9]{4}-(0[1-9]|1[0-2])-(0[1-9]|[1-2][0-9]|3[0-1])", dataset_outdir):
            dataset_outdir = "{}-{}".format(date_str, dataset_outdir)
        self.dataset_outdir = dataset_outdir

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
        parser.add_argument("--work_dir",
                            type=str,
                            required=True,
                            default=None,
                            help="The working directory.")
        parser.add_argument("--dataset_outdir",
                            type=str,
                            required=True,
                            default=None,
                            help="The name of the output directory where you "
                                 "would like all outputs/results saved.")
        parser.add_argument("--comb_datasets_dir",
                            type=str,
                            required=True,
                            default=None,
                            help=".")
        parser.add_argument("--comb_datasets_config",
                            type=str,
                            default="scMetaBrain_CombineDatasets.yaml",
                            help="")
        parser.add_argument("--repo_dir",
                            type=str,
                            required=True,
                            default=None,
                            help="")
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
        parser.add_argument("--datasheet_path",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--variants_path",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--exclude_temp_in_sbatch",
                            action='store_true',
                            help="")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Generating files")

        output_dir = os.path.join(self.work_dir, self.dataset_outdir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        log_dir = os.path.join(output_dir, "log")
        slurm_log_dir = os.path.join(output_dir, "slurm_log")
        for outdir in [output_dir, log_dir, slurm_log_dir]:
            if not os.path.exists(outdir):
                os.makedirs(outdir)
        snakefile = os.path.join(self.comb_datasets_dir, "Snakefile")
        configfile = os.path.join(output_dir, self.comb_datasets_config)

        config_arguments = (
            ("bind_path", self.bind_path),
            ("singularity_image", self.sif_path),
            ("repo_dir", os.path.join(self.repo_dir)),
            ("singularity_image", self.sif_path),
            ("datasheet_path", os.path.join(self.datasheet_path)),
            ("variants_path", self.variants_path),
            ("output_dir", os.path.join(output_dir))
        )
        self.write_configfile(
            template=os.path.join(self.comb_datasets_dir, self.comb_datasets_config),
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

        self.write_build_dag_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir,
            outfile="dag"
        )

        for jobs in [10000]:
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

        # self.write_report_script(
        #     snakefile=snakefile,
        #     configfile=configfile,
        #     report=os.path.join(output_dir, "combine_datasets_report.html"),
        #     output_dir=output_dir
        # )

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

    @staticmethod
    def write_lines_to_file(lines, path):
        with open(path, "w") as f:
            for line in lines:
                f.write(line + "\n")
        f.close()
        print("\tSaved file: {}".format(os.path.basename(path)))

    def print_arguments(self):
        print("General arguments")
        print("  > Working directory:              {}".format(self.work_dir))
        print("  > Dataset output directory:       {}".format(self.dataset_outdir))
        print("  > Comb. Datasets directory:       {}".format(self.comb_datasets_dir))
        print("  > Comb. Datasets configuration:   {}".format(self.comb_datasets_config))
        print("  > Bind path:                      {}".format(self.bind_path))
        print("  > Singularity path:               {}".format(self.sif_path))
        print("  > Repository folder:              {}".format(self.repo_dir))
        print("  > Datasheet path:                 {}".format(self.datasheet_path))
        print("  > Variants path:                  {}".format(self.variants_path))
        print("  > Exclude TEMP in SBATCH:         {}".format(self.exclude_temp_in_sbatch))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()