#!/usr/bin/env python3

"""
File:         run_cellbender.py
Created:      2023/03/23
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
import glob
import os

# Third party imports.
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Run CellBender"
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
./run_cellbender.py \
    --workdir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-03-23-CellBender/2023-03-23-Mathys2019 \
    --inputdir /groups/umcg-biogen/tmp01/input/processeddata/single-cell/Mathys2019 \
    --dry_run \
    --epochs 250

"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.workdir = getattr(arguments, 'workdir')
        self.inputdir = getattr(arguments, 'inputdir')
        self.time = TIME_DICT[getattr(arguments, 'time')]
        self.cpus_per_task = getattr(arguments, 'cpus_per_task')
        self.mem = getattr(arguments, 'mem')
        self.dry_run = getattr(arguments, 'dry_run')

        # Safe the CellBender arguments.
        self.expected_cells = getattr(arguments, 'expected_cells')
        self.total_droplets_included = getattr(arguments, 'total_droplets_included')
        self.model = getattr(arguments, 'model')
        self.epochs = getattr(arguments, 'epochs')
        self.cuda = getattr(arguments, 'cuda')
        self.low_count_threshold = getattr(arguments, 'low_count_threshold')
        self.z_dim = getattr(arguments, 'z_dim')
        self.z_layers = getattr(arguments, 'z_layers')
        self.training_fraction = getattr(arguments, 'training_fraction')
        self.empty_drop_training_fraction = getattr(arguments, 'empty_drop_training_fraction')
        self.blacklist_genes = getattr(arguments, 'blacklist_genes')
        self.fpr = getattr(arguments, 'fpr')
        self.exclude_antibody_capture = getattr(arguments, 'exclude_antibody_capture')
        self.learning_rate = getattr(arguments, 'learning_rate')
        self.final_elbo_fail_fraction = getattr(arguments, 'final_elbo_fail_fraction')
        self.epoch_elbo_fail_fraction = getattr(arguments, 'epoch_elbo_fail_fraction')
        self.num_training_tries = getattr(arguments, 'num_training_tries')
        self.learning_rate_retry_mult = getattr(arguments, 'learning_rate_retry_mult')
        self.posterior_batch_size = getattr(arguments, 'posterior_batch_size')
        self.cells_posterior_reg_calc = getattr(arguments, 'cells_posterior_reg_calc')

        # Set variables.
        self.jobdir = os.path.join(self.workdir, "jobs")
        self.jobs_outdir = os.path.join(self.jobdir, "output")
        for dir in [self.workdir, self.jobdir, self.jobs_outdir]:
            if not os.path.exists(dir):
                os.makedirs(dir)

        self.defaults = [
            ("workdir", self.workdir, "required"),
            ("inputdir", self.inputdir, "required"),
            ("time", self.time, "required"),
            ("cpus-per-task", self.cpus_per_task, "required"),
            ("mem", self.mem, "required"),
            ("dry_run", self.dry_run, "required"),
            ("expected-cells", self.expected_cells, None),
            ("total-droplets-included", self.total_droplets_included, 25000),
            ("model", self.model, "full"),
            ("epochs", self.epochs, 150),
            ("cuda", self.cuda, False),
            ("low-count-threshold", self.low_count_threshold, 15),
            ("z-dim", self.z_dim, 100),
            ("z-layers", self.z_layers, [500]),
            ("training-fraction", self.training_fraction, 0.9),
            ("empty-drop-training-fraction", self.empty_drop_training_fraction, 0.5),
            ("blacklist-genes", self.blacklist_genes, []),
            ("fpr", self.fpr, [0.01]),
            ("exclude-antibody-capture", self.exclude_antibody_capture, False),
            ("learning-rate", self.learning_rate, 1e-4),
            ("final-elbo-fail-fraction", self.final_elbo_fail_fraction, None),
            ("epoch-elbo-fail-fraction", self.epoch_elbo_fail_fraction, None),
            ("num-training-tries", self.num_training_tries, 1),
            ("learning-rate-retry-mult", self.learning_rate_retry_mult, 0.5),
            ("posterior-batch-size", self.posterior_batch_size, 20),
            ("cells-posterior-reg-calc", self.cells_posterior_reg_calc, 100)
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
        parser.add_argument("--workdir",
                            type=str,
                            required=True,
                            help="The working directory.")
        parser.add_argument("--inputdir",
                            type=str,
                            required=True,
                            help="The directory path containing the folders "
                                 "that yield the data on which to run tool, as "
                                 "either _raw_gene_barcode_matrices_h5.h5 file,"
                                 " or as the path containing the raw "
                                 "matrix.mtx, barcodes.tsv, and genes.tsv "
                                 "files. Supported for outputs of CellRanger "
                                 "v2 and v3.")
        parser.add_argument("--time",
                            type=str,
                            required=False,
                            choices=["short", "medium", "long"],
                            default="short",
                            help="The maximum walltime for each job.")
        parser.add_argument("--cpus_per_task",
                            type=int,
                            required=False,
                            default=1,
                            help="Restricts CellBender to use specified amount "
                                 "of cores. (default: 1)")
        parser.add_argument("--mem",
                            type=int,
                            required=False,
                            default=8,
                            help="Restricts CellBender to use specified amount "
                                 "of memory (in GB). (default: 8)")
        parser.add_argument("--dry_run",
                            action='store_true',
                            help="Only create the job files, don't submit them.")

        parser.add_argument("--expected_cells",
                            type=int,
                            required=False,
                            default=None,
                            help="Number of cells expected in the dataset "
                                 "(a rough estimate within a factor of 2 is "
                                 "sufficient).")
        parser.add_argument("--total_droplets_included",
                            type=int,
                            required=False,
                            default=25000,
                            help="The number of droplets from the rank-ordered "
                                 "UMI plot that will be analyzed. The largest "
                                 "'total_droplets' droplets will have their "
                                 "cell probabilities inferred as an output. "
                                 "(default: 25000).")
        parser.add_argument("--model",
                            type=str,
                            required=False,
                            choices=["simple", "ambient", "swapping", "full"],
                            default="full",
                            help="Which model is being used for count data. "
                                 "'simple' does not model either ambient RNA "
                                 "or random barcode swapping (for debugging "
                                 "purposes -- not recommended). 'ambient' "
                                 "assumes background RNA is incorporated "
                                 "into droplets. 'swapping' assumes "
                                 "background RNA comes from random barcode "
                                 "swapping. 'full' uses a combined ambient "
                                 "and swapping model. (default: full).")
        parser.add_argument("--epochs",
                            type=int,
                            required=False,
                            default=150,
                            help="Number of epochs to train. (default: 150).")
        parser.add_argument("--cuda",
                            action='store_true',
                            help="Including the flag --cuda will run the "
                                 "inference on a GPU.")
        parser.add_argument("--low_count_threshold",
                            type=int,
                            default=15,
                            help="Droplets with UMI counts below this number "
                                 "are completely excluded from the analysis. "
                                 "This can help identify the correct prior "
                                 "for empty droplet counts in the rare case "
                                 "where empty counts are extremely high "
                                 "(over 200). (default: 15)")
        parser.add_argument("--z_dim",
                            type=int,
                            required=False,
                            default=100,
                            help="Dimension of latent variable z. (default: "
                                 "100)")
        parser.add_argument("--z_layers",
                            nargs="+",
                            type=int,
                            required=False,
                            default=[500],
                            help="Dimension of hidden layers in the encoder "
                                 "for z. (default: [500])")
        parser.add_argument("--training_fraction",
                            type=float,
                            required=False,
                            default=0.9,
                            help="Training detail: the fraction of the data "
                                 "used for training. The rest is never seen "
                                 "by the inference algorithm. Speeds up "
                                 "learning. (default: 0.9)")
        parser.add_argument("--empty_drop_training_fraction",
                            type=float,
                            required=False,
                            default=0.5,
                            help="Training detail: the fraction of the "
                                 "training data each epoch that is drawn "
                                 "(randomly sampled) from surely empty "
                                 "droplets. (default: 0.5)")
        parser.add_argument("--blacklist_genes",
                            nargs="+",
                            type=int,
                            required=False,
                            default=[],
                            help="Integer indices of genes to ignore entirely. "
                                 "In the output count matrix, the counts for "
                                 "these genes will be set to zero.")
        parser.add_argument("--fpr",
                            nargs="+",
                            type=float,
                            required=False,
                            default=[0.01],
                            help="Target false positive rate in (0, 1). A "
                                 "false positive is a true signal count that "
                                 "is erroneously removed. More background "
                                 "removal is accompanied by more signal "
                                 "removal at high values of FPR. You can "
                                 "specify multiple values, which will create "
                                 "multiple output files. (default: [0.01])")
        parser.add_argument("--exclude_antibody_capture",
                            action='store_true',
                            help="Including the flag "
                                 "--exclude-antibody-capture will cause "
                                 "remove-background to operate on gene "
                                 "counts only, ignoring other features.")
        parser.add_argument("--learning_rate",
                            type=float,
                            required=False,
                            default=1e-4,
                            help="Training detail: lower learning rate for "
                                 "inference. A OneCycle learning rate schedule "
                                 "is used, where the upper learning rate is "
                                 "ten times this value. (For this value, "
                                 "probably do not exceed 1e-3). (default: "
                                 "0.0001)")
        parser.add_argument("--final_elbo_fail_fraction",
                            type=float,
                            required=False,
                            default=None,
                            help="Training is considered to have failed if "
                                 "(best_test_ELBO - final_test_ELBO)/"
                                 "(best_test_DLBO - initial_train_ELBO) > "
                                 "FINAL_ELBO_FAIL_FRACTION.(default: do not "
                                 "fail training based on final_training_ELBO)")
        parser.add_argument("--epoch_elbo_fail_fraction",
                            type=float,
                            required=False,
                            default=None,
                            help="Training is considered to have failed if "
                                 "(previous_epoch_test_ELBO - "
                                 "current_epoch_test_ELBO)/"
                                 "(previous_epoch_test_ELBO - "
                                 "initial_train_ELBO) > "
                                 "EPOCH_ELBO_FAIL_FRACTION.(default: do not "
                                 "fail training based on epoch_training_ELBO)")
        parser.add_argument("--num_training_tries",
                            type=int,
                            required=False,
                            default=1,
                            help=" Number of times to attempt to train the "
                                 "model. Each subsequent attempt the learning "
                                 "rate is multiplied by "
                                 "LEARNING_RATE_RETRY_MULT. (default: 1)")
        parser.add_argument("--learning_rate_retry_mult",
                            type=float,
                            required=False,
                            default=0.5,
                            help="Learning rate is multiplied by this amount "
                                 "each time (default: 0.5)")
        parser.add_argument("--posterior_batch_size",
                            type=int,
                            required=False,
                            default=20,
                            help="Size of batches when creating the posterior. "
                                 "Reduce this to avoid running out of GPU "
                                 "memory creating the posterior (will be "
                                 "slower). (default: 20)")
        parser.add_argument("--cells_posterior_reg_calc",
                            type=int,
                            required=False,
                            default=100,
                            help="Number of cells used to estimate posterior "
                                 "regularization lambda. (default: 100)")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Starting job files.")
        arguments = self.filter_arguments()
        for path in glob.glob(os.path.join(self.inputdir, "*")):
            folder = os.path.basename(path)
            input = self.find_input(folder)
            if input is None:
                continue

            outpath = os.path.join(self.workdir, folder)
            if not os.path.exists(outpath):
                os.makedirs(outpath)
            output = os.path.join(self.workdir, folder, "cellbender_remove_background_output.h5")
            if os.path.exists(output):
                print("\tSkipping sample '{}'".format(folder))

            jobfile_path = self.create_job_file(sample=folder,
                                                input=input,
                                                output=output,
                                                arguments=arguments)

            command = ['sbatch', jobfile_path]
            self.run_command(command)


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

    def find_input(self, folder):
        raw_feature_bc_matrix_path = os.path.join(self.inputdir, folder, "outs", "raw_feature_bc_matrix.h5")
        if os.path.exists(raw_feature_bc_matrix_path):
            return raw_feature_bc_matrix_path

        raw_gene_bc_matrices_h5_path = os.path.join(self.inputdir, folder, "outs", "raw_gene_bc_matrices_h5.h5")
        if os.path.exists(raw_gene_bc_matrices_h5_path):
            return raw_gene_bc_matrices_h5_path

        raw_feature_bc_matrix_folder = os.path.join(self.inputdir, folder, "outs", "raw_feature_bc_matrix")
        if os.path.exists(raw_feature_bc_matrix_folder):
            for subdir, _, files in os.walk(raw_feature_bc_matrix_folder):
                if 'matrix.mtx' in files and 'barcodes.tsv' in files and 'genes.tsv' in files:
                    return subdir
                elif 'matrix.mtx.gz' in files and 'barcodes.tsv.gz' in files and 'genes.tsv.gz' in files:
                    return subdir

        raw_gene_bc_matrices_folder = os.path.join(self.inputdir, folder, "outs", "raw_gene_bc_matrices")
        if os.path.exists(raw_gene_bc_matrices_folder):
            for subdir, _, files in os.walk(raw_gene_bc_matrices_folder):
                if 'matrix.mtx' in files and 'barcodes.tsv' in files and 'genes.tsv' in files:
                    return subdir
                elif 'matrix.mtx.gz' in files and 'barcodes.tsv.gz' in files and 'genes.tsv.gz' in files:
                    return subdir

        return None

    def create_job_file(self, sample, input, output, arguments):
        job_name = "cellbender_remove_background_{}".format(sample)

        lines = ["#!/bin/bash",
                 "#SBATCH --job-name={}".format(job_name),
                 "#SBATCH --output={}".format(os.path.join(self.jobs_outdir, job_name + ".out")),
                 "#SBATCH --error={}".format(os.path.join(self.jobs_outdir, job_name + ".out")),
                 "#SBATCH --time={}".format(self.time),
                 "#SBATCH --cpus-per-task={}".format(self.cpus_per_task),
                 "#SBATCH --mem={}gb".format(self.mem),
                 "#SBATCH --nodes=1",
                 "#SBATCH --open-mode=append",
                 "#SBATCH --export=NONE",
                 "#SBATCH --get-user-env=L",
                 "",
                 "module load Anaconda3",
                 "source activate cellbender",
                 "",
                 "cd {} || exit".format(self.workdir),
                 "",
                 "cellbender remove-background \\",
                 "  --input={} \\".format(input),
                 "  --output={} \\".format(output)
                 ]

        for name, value in arguments.items():
            if isinstance(value, bool):
                lines.append("  --{} \\".format(name))
            else:
                lines.append("  --{}={} \\".format(name, value))

        if lines[-1].endswith("\\"):
            lines[-1] = lines[-1].strip("\\")

        lines.extend(["", "conda deactivate", ""])

        jobfile_path = os.path.join(self.jobdir, job_name + ".sh")
        with open(jobfile_path, "w") as f:
            for line in lines:
                f.write(line + "\n")
        f.close()
        print("\tSaved jobfile: {}".format(os.path.basename(jobfile_path)))
        return jobfile_path

    def run_command(self, command):
        if self.dry_run:
            return

        print("\t" + " ".join(command))
        # subprocess.call(command)

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