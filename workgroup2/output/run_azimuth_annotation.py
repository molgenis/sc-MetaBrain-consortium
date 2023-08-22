#!/usr/bin/env python3

"""
File:         run_azimuth_annotation.py
Created:      2023/04/07
Last Changed: 2023/04/10
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
import subprocess
import argparse
import glob
import os

# Third party imports.
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Run Azimuth Annotation"
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

TIME_DICT = {
    "short": "05:59:00",
    "medium": "23:59:00",
    "long": "6-23:59:00"
}

"""
Syntax: 

### Mathys2019 ###
./run_azimuth_annotation.py \
    --rscript /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-02-02-WorkGroup2CellType/2023-04-05-Mathys2019/azimuth_annotation.R \
    --workdir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-02-02-WorkGroup2CellType/2023-04-05-Mathys2019 \
    --query /groups/umcg-biogen/tmp01/input/processeddata/single-cell/datasets/Mathys2019/2023-04-05-SeuratObject/Mathys2019.rds \
    --refdir /groups/umcg-biogen/tmp01/input/processeddata/single-cell/screference/Bakken2020/reference \
    --dry_run
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.rscript = getattr(arguments, 'rscript')
        self.workdir = getattr(arguments, 'workdir')
        self.query = getattr(arguments, 'query')
        self.refdir = getattr(arguments, 'refdir')
        self.time = getattr(arguments, 'time')
        self.cpus_per_task = getattr(arguments, 'cpus_per_task')
        self.mem = getattr(arguments, 'mem')
        self.dry_run = getattr(arguments, 'dry_run')

        # Set variables.
        self.jobdir = os.path.join(self.workdir, "jobs")
        self.jobs_outdir = os.path.join(self.jobdir, "output")
        for dir in [self.workdir, self.jobdir, self.jobs_outdir]:
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
        parser.add_argument("--rscript",
                            type=str,
                            required=True,
                            help="The run_azimuth r script.")
        parser.add_argument("--workdir",
                            type=str,
                            required=True,
                            help="The working directory.")
        parser.add_argument("--query",
                            type=str,
                            required=True,
                            help="The query seurat file.")
        parser.add_argument("--refdir",
                            type=str,
                            required=True,
                            help="The reference seurat file directory.")
        parser.add_argument("--time",
                            type=str,
                            required=False,
                            default="short",
                            help="The maximum walltime for each job.")
        parser.add_argument("--cpus_per_task",
                            type=int,
                            required=False,
                            default=4,
                            help="Restricts script to use specified amount "
                                 "of cores. (default: 14)")
        parser.add_argument("--mem",
                            type=int,
                            required=False,
                            default=32,
                            help="Restricts script to use specified amount "
                                 "of memory (in GB). (default: 32)")
        parser.add_argument("--dry_run",
                            action='store_true',
                            help="Only create the job files, don't submit them.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Starting job files.")
        for reference in glob.glob(os.path.join(self.refdir, "*")):
            if not reference.lower().endswith(".rds"):
                continue

            folder = "all"
            if "_" in reference:
                folder = reference.split("_")[1].split(".")[0]
            print("  Running azimuth with reference '{}'".format(folder))
            outpath = os.path.join(self.workdir, folder)
            if not os.path.exists(outpath):
                os.makedirs(outpath)

            jobfile, logfile = self.create_jobfile(
                job_name="AZIMUTH_ANNOTATE_{}".format(folder),
                commands=["Rscript {} {} {} {} {}".format(self.rscript, self.workdir, self.query, reference, folder)],
                time=self.time,
                cpu=self.cpus_per_task,
                mem=self.mem,
                module_load=["R", "HDF5/1.12.2-gompi-2022a"],
            )

            skip = False
            if os.path.exists(logfile):
                for line in open(logfile, 'r'):
                    if line == "Job finished":
                        print("\t\tSkipping '{}'".format(os.path.basename(jobfile)))
                        skip = True
            if skip:
                continue

            self.submit_job(jobfile=jobfile)


    @staticmethod
    def load_file(inpath, header=0, index_col=None, sep=","):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    def filter_arguments(self):
        arguments = {}
        for name, value, default in self.defaults:
            if value != default and default != "required":
                arguments[name] = value
        return arguments

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

    def submit_job(self, jobfile):
        if self.dry_run:
            return None

        command = ['sbatch', jobfile]
        print("\t" + " ".join(command))
        subprocess.call(command)

    def print_arguments(self):
        print("Arguments:")
        print("  > Rscript: {}".format(self.rscript))
        print("  > Working directory: {}".format(self.workdir))
        print("  > Query object: {}".format(self.query))
        print("  > Reference directory: {}".format(self.refdir))
        print("  > Time: {}".format(self.time))
        print("  > CPUs per task: {}".format(self.cpus_per_task))
        print("  > Memory: {}".format(self.mem))
        print("  > dryrun: {}".format(self.dry_run))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()