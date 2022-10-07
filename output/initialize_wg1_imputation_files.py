#!/usr/bin/env python3

"""
File:         initialize_wg1_imputation_files.py
Created:      2022/10/07
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
import os

# Third party imports.

# Local application imports.

# Metadata
__program__ = "Initialize Workgroup1 Imputation Files"
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
./initialize_wg1_imputation_files.py \
    --workdir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2022-10-07-Imputation/ \
    --plink_dir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2022-10-07-Imputation/ImputationTestDataset_plink \
    --output_dir 2022-10-07-ImputationTestDataset
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.workdir = getattr(arguments, 'workdir')
        self.yaml_template = getattr(arguments, 'yaml_template')
        self.ref_dir = getattr(arguments, 'ref_dir')
        self.singularity_image = getattr(arguments, 'singularity_image')
        self.plink_dir = getattr(arguments, 'plink_dir')
        self.bind_paths = ", ".join(getattr(arguments, 'bind_paths'))
        self.output_dir = os.path.join(self.workdir, getattr(arguments, 'output_dir'))

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
                            help="The working directory.")
        parser.add_argument("--yaml_template",
                            type=str,
                            required=False,
                            default="PreImputation.yaml",
                            help="The yaml template filename. Default: "
                                 "'PreImputation.yaml'.")
        parser.add_argument("--ref_dir",
                            type=str,
                            required=False,
                            default="hg38",
                            help="The directory the imputation references "
                                 "provided for the SNP processing and "
                                 "imputation. Default: 'hg38'.")
        parser.add_argument("--singularity_image",
                            type=str,
                            required=False,
                            default="WG1-pipeline-QC_imputation.sif",
                            help="The singularity image that has all the "
                                 "softwares. Default: "
                                 "'WG1-pipeline-QC_imputation.sif'")
        parser.add_argument("--plink_dir",
                            type=str,
                            required=True,
                            help="Absolute path to directory that contains "
                                 "plink2 pgen, pvar and psam files for "
                                 "pipeline on hg19.")
        parser.add_argument("--bind_paths",
                            nargs="*",
                            type=str,
                            required=False,
                            default=["/groups/umcg-biogen/tmp01/"],
                            help="List of paths to bind to Singularity. "
                                 "Default: ['/groups/umcg-biogen/tmp01/']")
        parser.add_argument("--output_dir",
                            type=str,
                            required=True,
                            help="The name of the output directory where you "
                                 "would like all outputs/results saved.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        # Define directory / file variables.
        log_dir = os.path.join(self.output_dir, "log")
        imputation_snakefile = os.path.join(self.workdir, "Snakefile")

        print("Making output directories")
        for outdir in [self.output_dir, log_dir]:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

        print("Generating files")
        ########################################################################

        # Reading template and update lines.
        yaml_lines = []
        for line in open(os.path.join(self.workdir, self.yaml_template), 'r'):
            for label, argument in (("ref_dir", os.path.join(self.workdir, self.ref_dir)),
                                    ("singularity_image", os.path.join(self.workdir, self.singularity_image)),
                                    ("plink_dir", self.plink_dir),
                                    ("bind_paths", self.bind_paths),
                                    ("output_dir", self.output_dir)):
                if label in line:
                    line = "  {}: {}\n".format(label, argument)
                line = line.replace("\n", "")
            yaml_lines.append(line)

        imputation_config = os.path.join(self.output_dir, self.yaml_template)
        self.write_lines_to_file(
            lines=yaml_lines,
            path=imputation_config
        )

        ########################################################################

        dry_lines = [
            '#!/bin/bash',
            '',
            'snakemake \\',
            '    --snakefile {} \\'.format(imputation_snakefile),
            '    --configfile {} \\'.format(imputation_config),
            '    --dryrun \\',
            '    --cores 1 \\',
            '    --reason'
        ]
        self.write_lines_to_file(
            lines=dry_lines,
            path=os.path.join(self.output_dir, "dry_run.sh")
        )

        ########################################################################

        run_phase1_lines = [
            '#!/bin/bash',
            '',
            'nohup \\',
            '  snakemake \\',
            '    --snakefile {} \\'.format(imputation_snakefile),
            '    --configfile {} \\'.format(imputation_config),
            '    --rerun-incomplete \\',
            '    --jobs 48 \\',
            '    --use-singularity \\',
            '    --restart-times 2 \\',
            '    --keep-going \\',
            '    --cluster \\',
            '       "sbatch \\',
            '       --qos debug \\',
            '       -N 1 \\',
            '       --ntasks 1 \\',
            '       --cpus-per-task 48 \\',
            '       -o {} \\'.format(os.path.join(log_dir, '%{rule}.out')),
            '       --export ALL" \\',
            '    > {}/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &'.format(log_dir)
        ]
        self.write_lines_to_file(
            lines=run_phase1_lines,
            path=os.path.join(self.output_dir, "run_phase1.sh")
        )

        ########################################################################

        run_phase1_with_latency_lines = [
            '#!/bin/bash',
            '',
            'nohup \\',
            '  snakemake \\',
            '    --snakefile {} \\'.format(imputation_snakefile),
            '    --configfile {} \\'.format(imputation_config),
            '    --rerun-incomplete \\',
            '    --jobs 1 \\',
            '    --use-singularity \\',
            '    --restart-times 2 \\',
            '    --keep-going \\',
            '    --latency-wait 30 \\',
            '    --cluster \\',
            '       "sbatch \\',
            '       --qos regular \\',
            '       -N {threads} \\',
            '        --mem={resources.mem_per_thread_gb}G \\',
            '        --tmp={resources.disk_per_thread_gb}G \\',
            '       -o {} \\'.format(os.path.join(log_dir, '%{rule}.out')),
            '        --export ALL \\',
            '        --time=05:59:59" \\',
            '    > {}/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &'.format(log_dir)
        ]
        self.write_lines_to_file(
            lines=run_phase1_with_latency_lines,
            path=os.path.join(self.output_dir, "run_phase1_with_latency_wait.sh")
        )

        ########################################################################

        run_phase2_lines = [
            '#!/bin/bash',
            '',
            'nohup \\',
            '  snakemake \\',
            '    --snakefile {} \\'.format(imputation_snakefile),
            '    --configfile {} \\'.format(imputation_config),
            '    --rerun-incomplete \\',
            '    --jobs 20 \\',
            '    --use-singularity \\',
            '    --restart-times 2 \\',
            '    --keep-going \\',
            '    --latency-wait 30 \\',
            '    --cluster \\',
            '        "qsub -S /bin/bash \\',
            '            -q short.q \\',
            '            -r yes \\',
            '            -pe smp {threads} \\',
            '            -l tmp_requested={resources.disk_per_thread_gb}G \\',
            '            -l mem_requested={resources.mem_per_thread_gb}G \\',
            '            -e {} \\'.format(log_dir),
            '            -o {} \\'.format(log_dir),
            '            -j y \\',
            '            -V" \\',
            '    > {}/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &'.format(log_dir)
        ]
        self.write_lines_to_file(
            lines=run_phase2_lines,
            path=os.path.join(self.output_dir, "run_phase2.sh")
        )

        ########################################################################

        report_lines = [
            '#!/bin/bash',
            '',
            'nohup \\',
            '  snakemake \\',
            '    --snakefile {} \\'.format(imputation_snakefile),
            '    --configfile {} \\'.format(imputation_config),
            '    --report {} \\'.format(self.output_dir, 'imputation_report.html')
        ]
        self.write_lines_to_file(
            lines=report_lines,
            path=os.path.join(self.output_dir, "report.sh")
        )

    @staticmethod
    def write_lines_to_file(lines, path):
        with open(path, "w") as f:
            for line in lines:
                f.write(line + "\n")
        f.close()
        print("\tSaved file: {}".format(os.path.basename(path)))

    def print_arguments(self):
        print("Arguments:")
        print("  > yaml_template: {}".format(self.yaml_template))
        print("  > ref_dir: {}".format(self.ref_dir))
        print("  > singularity_image: {}".format(self.singularity_image))
        print("  > plink_dir: {}".format(self.plink_dir))
        print("  > bind_paths: {}".format(self.bind_paths))
        print("  > output_dir: {}".format(self.output_dir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()