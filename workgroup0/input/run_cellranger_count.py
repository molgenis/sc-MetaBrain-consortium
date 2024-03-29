#!/usr/bin/env python3

"""
File:         run_cellranger_count.py
Created:      2022/09/22
Last Changed: 2023/07/04
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
import os

# Third party imports.
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Run Cell Ranger Count"
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
./run_cellranger_count.py -h
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.env = getattr(arguments, 'env')
        self.workdir = getattr(arguments, 'workdir')
        self.sample_id_table = getattr(arguments, 'sample_id_table')
        self.sample_col = getattr(arguments, 'sample_col')
        self.id_col = getattr(arguments, 'id_col')
        self.cellranger = getattr(arguments, 'cellranger')
        time = getattr(arguments, 'time')
        if time in TIME_DICT:
            time = TIME_DICT[time]
        self.time = time
        self.break_lock = getattr(arguments, 'break_lock')
        self.dry_run = getattr(arguments, 'dry_run')

        # Safe the Cell Ranger arguments.
        self.fastqs = getattr(arguments, 'fastqs')
        self.libraries = getattr(arguments, 'libraries')
        self.transcriptome = getattr(arguments, 'transcriptome')
        self.feature_ref = getattr(arguments, 'feature_ref')
        self.target_panel = getattr(arguments, 'target_panel')
        self.no_target_umi_filter = getattr(arguments, 'no_target_umi_filter')
        self.expect_cells = getattr(arguments, 'expect_cells')
        self.force_cells = getattr(arguments, 'force_cells')
        self.include_introns = getattr(arguments, 'include_introns')
        self.nosecondary = getattr(arguments, 'nosecondary')
        self.no_bam = getattr(arguments, 'no_bam')
        self.no_libraries = getattr(arguments, 'no_libraries')
        self.chemistry = getattr(arguments, 'chemistry')
        self.r1_length = getattr(arguments, 'r1_length')
        self.r2_length = getattr(arguments, 'r2_length')
        self.lanes = getattr(arguments, 'lanes')
        self.localcores = getattr(arguments, 'localcores')
        self.localmem = getattr(arguments, 'localmem')
        self.check_library_compatibility = getattr(arguments, 'check_library_compatibility')

        # Set variables.
        self.jobdir = os.path.join(self.workdir, "jobs")
        self.jobs_outdir = os.path.join(self.jobdir, "output")
        for dir in [self.workdir, self.jobdir, self.jobs_outdir]:
            if not os.path.exists(dir):
                os.makedirs(dir)

        self.defaults = [
            ("env", self.env, "required"),
            ("workdir", self.workdir, "required"),
            ("sample_id_table", self.sample_id_table, "required"),
            ("sample_col", self.sample_col, "required"),
            ("id_col", self.id_col, "required"),
            ("cellranger", self.cellranger, "required"),
            ("time", self.time, "required"),
            ("fastqs", self.fastqs, "required"),
            ("break_lock", self.dry_run, "required"),
            ("dry_run", self.dry_run, "required"),
            ("libraries", self.libraries, None),
            ("transcriptome", self.transcriptome, "required"),
            ("feature-ref", self.feature_ref, None),
            ("target-panel", self.target_panel, None),
            ("no-target-umi-filter", self.no_target_umi_filter, False),
            ("expect-cells", self.expect_cells, None),
            ("force-cells", self.force_cells, None),
            ("include-introns", self.include_introns, True),
            ("nosecondary", self.nosecondary, False),
            ("no-bam", self.no_bam, False),
            ("no-libraries", self.no_libraries, False),
            ("chemistry", self.chemistry, "auto"),
            ("r1-length", self.r1_length, None),
            ("r2-length", self.r2_length, None),
            ("lanes", self.lanes, None),
            ("localcores", self.localcores, None),
            ("localmem", self.localmem, None),
            ("check-library-compatibility", self.check_library_compatibility, True)
        ]

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
        parser.add_argument("--env",
                            type=str,
                            required=True,
                            help="The python environment.")
        parser.add_argument("--workdir",
                            type=str,
                            required=True,
                            help="The working directory.")
        parser.add_argument("--sample_id_table",
                            type=str,
                            required=True,
                            help="The sample-ID link matrix.")
        parser.add_argument("--sample_col",
                            type=str,
                            required=True,
                            help="The sample ID column in --sample-id-table.")
        parser.add_argument("--id_col",
                            type=str,
                            required=True,
                            help="The ID column in --sample-id-table.")
        parser.add_argument("--cellranger",
                            type=str,
                            required=True,
                            help="The path to the cellranger software.")
        parser.add_argument("--time",
                            type=str,
                            required=False,
                            default="short",
                            help="The time for each job.")
        parser.add_argument("--fastqs",
                            type=str,
                            required=True,
                            help="Path of the fastq_path folder generated by "
                                 "cellranger mkfastq.")
        parser.add_argument("--break_lock",
                            action='store_true',
                            help="Remove locks from locked directories.")
        parser.add_argument("--dry_run",
                            action='store_true',
                            help="Only create the job files, don't submit them.")

        parser.add_argument("--libraries",
                            type=str,
                            default=None,
                            help="Path to a libraries.csv file declaring "
                                 "FASTQ paths and library types of input "
                                 "libraries. Required for gene expression + "
                                 "Feature Barcode analysis. When using this "
                                 "argument, --fastqs and --sample must not be "
                                 "used. This argument should not be used when "
                                 "performing gene expression-only analysis; "
                                 "use --fastqs instead.")
        parser.add_argument("--transcriptome",
                            type=str,
                            required=True,
                            help="Path to the Cell Ranger-compatible "
                                 "transcriptome reference.")
        parser.add_argument("--feature_ref",
                            type=str,
                            default=None,
                            help="Path to a Feature Reference CSV file "
                                 "declaring the Feature Barcode reagents used "
                                 "in the experiment. Required for Feature "
                                 "Barcode analysis.")
        parser.add_argument("--target_panel",
                            type=str,
                            default=None,
                            help="Path to a Target Panel CSV file declaring "
                                 "the target panel used, if any. Required for "
                                 "Targeted Gene Expression analysis. Default "
                                 "analysis will exclude intronic mapped "
                                 "reads, which is the recommended mode for "
                                 "targeted assays.")
        parser.add_argument("--no_target_umi_filter",
                            action='store_true',
                            help="Add this flag to disable targeted UMI "
                                 "filtering. See Targeted algorithms for "
                                 "details.")
        parser.add_argument("--expect_cells",
                            type=int,
                            default=None,
                            help="Override the pipeline’s "
                                 "auto-estimate. See cell calling algorithm "
                                 "overview for details on how this parameter "
                                 "is used. If used, enter the expected number "
                                 "of recovered cells.")
        parser.add_argument("--force_cells",
                            type=int,
                            default=None,
                            help="Force pipeline to use this number of cells, "
                                 "bypassing the cell detection algorithm. Use "
                                 "this if the number of cells estimated by "
                                 "Cell Ranger is not consistent with the "
                                 "barcode rank plot.")
        parser.add_argument("--include_introns",
                            action='store_false',
                            help="Set to false to exclude intronic reads in "
                                 "count. Including introns in analysis is "
                                 "recommended to maximize sensitivity, "
                                 "except when --target-panel is used. "
                                 "Default: true.")
        parser.add_argument("--nosecondary",
                            action='store_true',
                            help="Add this flag to skip secondary analysis of "
                                 "the feature-barcode matrix (dimensionality "
                                 "reduction, clustering, and visualization). "
                                 "Set this if you plan to use cellranger "
                                 "reanalyze or your own custom analysis.")
        parser.add_argument("--no_bam",
                            action='store_true',
                            help="Add this flag to skip BAM file generation. "
                                 "This will reduce the total computation time "
                                 "for the pipestance and the size of the "
                                 "output directory. If unsure, we recommend "
                                 "not using this option, as BAM files can be "
                                 "useful for troubleshooting and downstream "
                                 "analysis. Default: false.")
        parser.add_argument("--no_libraries",
                            action='store_true',
                            help="Add this flag to proceed with processing "
                                 "data using a --feature-ref but no Feature "
                                 "Barcode data specified with the --libraries "
                                 "flag.")
        parser.add_argument("--chemistry",
                            type=str,
                            default="auto",
                            choices=["auto",
                                     "threeprime",
                                     "fiveprime",
                                     "SC3Pv2",
                                     "SC3Pv3",
                                     "SC3Pv3LT",
                                     "SC3Pv3HT",
                                     "SC5P-PE",
                                     "SC5P-R2",
                                     "SC3Pv1"],
                            help="Assay configuration. NOTE: by default the "
                                 "assay configuration is detected "
                                 "automatically, which is the recommended "
                                 "mode. You should only specify chemistry if "
                                 "there is an error in automatic detection.")
        parser.add_argument("--r1_length",
                            type=int,
                            default=None,
                            help="Limit the length of the input Read 1 "
                                 "sequence of Gene Expression (and any Feature "
                                 "Barcode) library to the first N bases, where "
                                 "N is a user-supplied value. Note that the "
                                 "length includes the 10x Barcode and UMI "
                                 "sequences so do not set this below 26 for "
                                 "Single Cell 3′ v2 or Single Cell 5′. This "
                                 "and --r2-length are useful options for "
                                 "determining the optimal read length for "
                                 "sequencing.")
        parser.add_argument("--r2_length",
                            type=int,
                            default=None,
                            help="Limit the length of the input R2 sequence "
                                 "to the first N bases, where N is a "
                                 "user-supplied value. Trimming occurs before "
                                 "sequencing metrics are computed and "
                                 "therefore, limiting R2 read length may "
                                 "affect Q30 scores.")
        parser.add_argument("--lanes",
                            type=str,
                            default=None,
                            help="Lanes associated with this sample")
        parser.add_argument("--localcores",
                            type=int,
                            default=None,
                            help="Restricts cellranger to use specified number "
                                 "of cores to execute pipeline stages. By "
                                 "default, cellranger will use all of the "
                                 "cores available on your system.")
        parser.add_argument("--localmem",
                            type=int,
                            default=None,
                            help="Restricts cellranger to use specified amount "
                                 "of memory (in GB) to execute pipeline stages."
                                 " By default, cellranger will use 90% of the "
                                 "memory available on your system.")
        parser.add_argument("--check_library_compatibility",
                            action='store_false',
                            help="This option allows users to disable the "
                                 "check that evaluates 10x Barcode overlap "
                                 "between libraries when multiple libraries "
                                 "are specified (e.g., Gene Expression + "
                                 "Antibody Capture). Setting this option to "
                                 "false will disable the check across all "
                                 "library combinations. We recommend running "
                                 "this check (default), however if the "
                                 "pipeline errors out, users can bypass the "
                                 "check to generate outputs for "
                                 "troubleshooting. Default: true")


        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading sample-ID table.")
        sid_df = self.load_file(self.sample_id_table)
        print("\t  > Found {} samples".format(sid_df.shape[0]))

        print("Creating job files.")
        arguments = self.filter_arguments()
        info = []
        for _, row in sid_df.iterrows():
            sample = row[self.sample_col]
            id_value = row[self.id_col]

            jobfile_path, logfile_path = self.create_job_file(sample=sample,
                                                              id_value=id_value,
                                                              arguments=arguments)
            info.append((jobfile_path, logfile_path, id_value))

        print("Create aggregation CSV")
        aggregation_data = []
        for jobfile_path, _, id_value in info:
            aggregation_data.append([id_value, os.path.join(self.workdir, str(id_value), "outs", "molecule_info.h5")])
        aggregation_df = pd.DataFrame(aggregation_data, columns=["sample_id", "molecule_h5"])
        self.save_file(df=aggregation_df,
                       outpath=os.path.join(self.workdir, "aggregation.csv"))

        print("Starting job files.")
        for jobfile_path, logfile_path, id_value in info:
            success, log_exists, error_reason = self.parse_logfile(logfile_path)

            if success:
                print("\tSample '{}' already completed successfully!".format(id_value))
            else:
                if log_exists:
                    print("\tSample '{}' failed due to {}".format(id_value, error_reason))

                lock_file = os.path.join(self.workdir, str(id_value), "_lock")
                if os.path.exists(lock_file):
                    if self.break_lock:
                        print("\tWarning, folder is locked. Breaking lock and continue sbatch.")
                        os.remove(lock_file)
                        self.run_command(command=['sbatch', jobfile_path])
                    else:
                        print("\tWarning, folder is locked. Sbatch not possible.")
                else:
                    self.run_command(command=['sbatch', jobfile_path])


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

    def create_job_file(self, sample, id_value, arguments):
        job_name = "cellranger_count_{}".format(sample)
        logfile_path = os.path.join(self.jobs_outdir, job_name + ".out")

        lines = ["#!/bin/bash",
                 "#SBATCH --job-name={}".format(job_name),
                 "#SBATCH --output={}".format(logfile_path),
                 "#SBATCH --error={}".format(logfile_path),
                 "#SBATCH --time={}".format(self.time),
                 "#SBATCH --cpus-per-task={}".format(1 if self.localcores is None else self.localcores),
                 "#SBATCH --mem={}gb".format(8 if self.localmem is None else self.localmem),
                 "#SBATCH --nodes=1",
                 "#SBATCH --open-mode=append",
                 "#SBATCH --export=NONE",
                 "#SBATCH --get-user-env=L",
                 "",
                 "export PATH={}:$PATH".format(self.cellranger),
                 "source {}".format(self.env),
                 "",
                 "cd {} || exit".format(self.workdir),
                 "",
                 "cellranger count \\",
                 "  --id={} \\".format(id_value),
                 "  --fastqs={} \\".format(self.fastqs),
                 "  --sample={} \\".format(sample),
                 "  --transcriptome={} \\".format(self.transcriptome)]

        for name, value in arguments.items():
            if isinstance(value, bool):
                lines.append("  --{} \\".format(name))
            else:
                lines.append("  --{}={} \\".format(name, value))

        if lines[-1].endswith("\\"):
            lines[-1] = lines[-1].strip("\\")

        lines.extend(["", "deactivate", ""])

        jobfile_path = os.path.join(self.jobdir, job_name + ".sh")
        with open(jobfile_path, "w") as f:
            for line in lines:
                f.write(line + "\n")
        f.close()
        print("\tSaved jobfile: {}".format(os.path.basename(jobfile_path)))
        return jobfile_path, logfile_path

    @staticmethod
    def parse_logfile(logfile_path):
        success = False
        log_exists = os.path.exists(logfile_path)
        error_reason = "UNKNOWN"
        if log_exists:
            with open(logfile_path, 'r') as f:
                for line in f:
                    if 'Pipestance completed successfully!' in line:
                        success = True
                    if 'Pipestance failed' in line:
                        error_reason = "Pipestance error: UNKNOWN"
                    if 'slurmstepd: error:' in line:
                        if "TIME LIMIT" in line:
                            error_reason = "SLURM error: TIME LIMIT"
                        elif "MEMORY LIMIT" in line:
                            error_reason = "SLURM error: MEMORY LIMIT"
            f.close()

        return success, log_exists, error_reason

    def run_command(self, command):
        print("\t" + " ".join(command))
        if self.dry_run:
            print("\t\tNot executed due to dry run")
            return

        subprocess.call(command)

    @staticmethod
    def save_file(df, outpath, header=True, index=False, sep=","):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    def print_arguments(self):
        print("Arguments:")
        for name, value, default in self.defaults:
            if value != default:
                print("  > {}: {}".format(name, value))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()