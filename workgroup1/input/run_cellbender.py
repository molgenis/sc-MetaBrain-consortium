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
import math
import os
import re

# Third party imports.
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
    --time medium \
    --dry_run

./run_cellbender.py \
    --workdir /scratch/p301710/2023-03-28-scMetaBrainConsortium/output/2023-03-28-CellBender/2023-03-28-Mathys2019 \
    --inputdir /scratch/p301710/2023-03-28-scMetaBrainConsortium/input/processeddata/Mathys2019/ \
    --cuda \
    --num_training_tries 3 \
    --partition gpu \
    --time 01:59:00 \
    --dry_run
    
./run_cellbender.py \
    --workdir /scratch/p301710/2023-03-28-scMetaBrainConsortium/output/2023-03-28-CellBender/2023-03-28-Mathys2019 \
    --inputdir /scratch/p301710/2023-03-28-scMetaBrainConsortium/input/processeddata/Mathys2019/ \
    --cuda \
    --num_training_tries 3 \
    --partition gpu \
    --time 01:59:00 \
    --rerun 10101327 11409232 11336574 20249897 11302830 \
    --dry_run 

"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.workdir = getattr(arguments, 'workdir')
        self.inputdir = getattr(arguments, 'inputdir')
        self.partition = getattr(arguments, 'partition')
        self.gpus_per_node = getattr(arguments, 'gpus_per_node')
        time = getattr(arguments, 'time')
        self.time = TIME_DICT[time] if time in TIME_DICT else time
        self.cpus_per_task = getattr(arguments, 'cpus_per_task')
        self.mem = getattr(arguments, 'mem')
        self.dry_run = getattr(arguments, 'dry_run')
        self.extensions = getattr(arguments, 'extension')
        self.rerun = getattr(arguments, 'rerun')

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
            ("partition", self.partition, "required"),
            ("gpus-per-node", self.gpus_per_node, "required"),
            ("time", self.time, "required"),
            ("cpus-per-task", self.cpus_per_task, "required"),
            ("mem", self.mem, "required"),
            ("dry_run", self.dry_run, "required"),
            ("extensions", self.extensions, "required"),
            ("rerun", self.rerun, "required"),
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
        parser.add_argument("--partition",
                            type=str,
                            required=False,
                            choices=["regular", "parallel", "gpu", "himem", "gelifes"],
                            default=None,
                            help="The partition to submit to.")
        parser.add_argument("--gpus_per_node",
                            type=str,
                            required=False,
                            choices=["a100.20gb:1", "a100:1", "a100:2"],
                            default="a100.20gb:1",
                            help="Request a specific GPU type when using "
                                 "--partition=gpu.")
        parser.add_argument("--time",
                            type=str,
                            required=False,
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
                            default=4,
                            help="Restricts CellBender to use specified amount "
                                 "of memory (in GB). (default: 4)")
        parser.add_argument("--dry_run",
                            action='store_true',
                            help="Only create the job files, don't submit them.")
        parser.add_argument("--extension",
                            nargs="+",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")
        parser.add_argument("--rerun",
                            nargs="+",
                            type=str,
                            required=False,
                            default=None,
                            help="The name of sample to run again.")

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

        stats_data = []
        training_data = []
        test_data = []

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

            sample_rerun = False
            if os.path.exists(output) and ((self.rerun is not None) and (folder in self.rerun)):
                sample_rerun = True
                arguments["learning-rate"] = self.learning_rate * self.learning_rate_retry_mult

            jobfile_path = self.create_job_file(sample=folder,
                                                input=input,
                                                output=output,
                                                arguments=arguments)

            if sample_rerun:
                print("\tRerunning sample '{}'".format(folder))
                command = ['sbatch', jobfile_path]
                self.run_command(command)
            elif os.path.exists(output):
                print("\tSample '{}' already finished".format(folder))

                resultfile = os.path.join(self.workdir, folder, "cellbender_remove_background_output.h5")
                if not os.path.exists(resultfile):
                    continue
                stats_s, training_s, test_s = self.parse_results(resultfile, folder)

                logfile = os.path.join(self.workdir, folder, "cellbender_remove_background_output.log")
                if os.path.exists(logfile):
                    stats_s = self.add_stats_from_log(stats_s, logfile)

                stats_data.append(stats_s)
                training_data.append(training_s)
                test_data.append(test_s)

                if folder in ['20254740', '20282398']:
                    print(stats_s)
            else:
                command = ['sbatch', jobfile_path]
                self.run_command(command)

        if len(stats_data) == 0:
            exit()

        print("Summarizing results.")
        stats_df = pd.concat(stats_data, axis=1)
        self.save_file(df=stats_df, outpath=os.path.join(self.workdir, "cellbender_remove_background_stats.txt.gz"))

        training_df = pd.concat(training_data, axis=1)
        self.save_file(df=training_df, outpath=os.path.join(self.workdir, "cellbender_remove_background_training_loss.txt.gz"))

        test_df = pd.concat(test_data, axis=1)
        self.save_file(df=test_df, outpath=os.path.join(self.workdir, "cellbender_remove_background_test_loss.txt.gz"))

        print("Visualising results.")
        self.visualise_stats(stats_df)

        plot_df = self.build_plot_df(training_df, test_df)
        validation_df = self.validate_training_procedure(training_df, plot_df)
        self.save_file(df=validation_df, outpath=os.path.join(self.workdir, "cellbender_remove_background_training_validation.txt.gz"))

        self.visualise_training_procedure(plot_df)

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
        # raw_feature_bc_matrix_path = os.path.join(self.inputdir, folder, "outs", "raw_feature_bc_matrix.h5")
        raw_feature_bc_matrix_path = os.path.join(self.inputdir, folder, "raw_feature_bc_matrix.h5")
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
                 "source activate CellBender",
                 "",
                 "cd {} || exit".format(self.workdir),
                 "",
                 "cellbender remove-background \\",
                 "  --input={} \\".format(input),
                 "  --output={} \\".format(output)
                 ]

        if self.partition is not None:
            lines.insert(5, "#SBATCH --partition={}".format(self.partition))
        if self.cuda is not None:
            lines.insert(6, "#SBATCH --gpus-per-node={}".format(self.gpus_per_node))

        for name, value in arguments.items():
            if isinstance(value, bool):
                lines.append("  --{} \\".format(name))
            else:
                lines.append("  --{}={} \\".format(name, value))

        if lines[-1].endswith("\\"):
            lines[-1] = lines[-1].strip("\\")

        lines.extend(["", "conda deactivate", "", "echo 'Job finished'"])

        jobfile_path = os.path.join(self.jobdir, job_name + ".sh")
        with open(jobfile_path, "w") as f:
            for line in lines:
                f.write(line + "\n")
        f.close()
        print("\tSaved jobfile: {}".format(os.path.basename(jobfile_path)))
        return jobfile_path

    @staticmethod
    def parse_results(resultsfile, folder):
        hf = h5py.File(resultsfile, 'r')

        # stats
        stats_data = {}
        for key in hf.get('matrix').keys():
            if key == "overdispersion_mean_and_scale":
                overdispersion_mean_and_scale = np.array(hf.get('matrix/overdispersion_mean_and_scale'))
                stats_data["overdispersion_mean"] = float(overdispersion_mean_and_scale[0])
                stats_data["overdispersion_scale"] = float(overdispersion_mean_and_scale[1])
            elif key == "contamination_fraction_params":
                contamination_fraction_params = np.array(hf.get('matrix/contamination_fraction_params'))
                stats_data["contamination_fraction_rho_alpha"] = float(contamination_fraction_params[0])
                stats_data["contamination_fraction_rho_beta"] = float(contamination_fraction_params[1])
            elif key == "lambda_multiplier":
                stats_data["lambda_multiplier"] = float(np.array(hf.get('matrix/lambda_multiplier')))
            elif key == "target_false_positive_rate":
                stats_data["target_false_positive_rate"] = float(np.array(hf.get('matrix/target_false_positive_rate')))
            elif key == "fraction_data_used_for_testing":
                stats_data["fraction_data_used_for_testing"] = float(np.array(hf.get('matrix/fraction_data_used_for_testing')))
            else:
                pass
        stats_s = pd.Series(stats_data)
        stats_s.name= folder

        # training loss
        training_s = pd.Series(np.array(hf.get('matrix/training_elbo_per_epoch')) * -1)
        training_s.index = np.arange(1, training_s.shape[0] + 1, 1)
        training_s.name = folder

        # test loss
        test_s = pd.Series(np.array(hf.get('matrix/test_elbo')) * -1)
        test_s.index = pd.Series(np.array(hf.get('matrix/test_epoch')))
        test_s.name = folder

        hf.close()

        return stats_s, training_s, test_s

    @staticmethod
    def add_stats_from_log(stats_s, logfile):
        with open(logfile, 'r') as f:
            for line in f:
                if re.search("Including ([0-9]+) genes that have nonzero counts.", line):
                    stats_s["nonzero_genes"] = int(re.search("Including ([0-9]+) genes that have nonzero counts.", line).group(1))
                elif re.search("Prior on counts in empty droplets is ([0-9]+)", line):
                    stats_s["prior_counts_empty"] = int(re.search("Prior on counts in empty droplets is ([0-9]+)", line).group(1))
                elif re.search("Prior on counts for cells is ([0-9]+)", line):
                    stats_s["prior_counts_cell"] = int(re.search("Prior on counts for cells is ([0-9]+)", line).group(1))
                elif re.search("Excluding barcodes with counts below ([0-9]+)", line):
                    stats_s["barcodes_threshold"] = int(re.search("Excluding barcodes with counts below ([0-9]+)", line).group(1))
                elif re.search("Using ([0-9]+) probable cell barcodes, plus an additional ([0-9]+) barcodes, and ([0-9]+) empty droplets.", line):
                    match = re.search("Using ([0-9]+) probable cell barcodes, plus an additional ([0-9]+) barcodes, and ([0-9]+) empty droplets.", line)
                    stats_s["probable_cell_barcodes"] = int(match.group(1))
                    stats_s["additional_barcodes"] = int(match.group(2))
                    stats_s["empty_droplets"] = int(match.group(3))
                elif re.search("Largest surely-empty droplet has ([0-9]+) UMI counts.", line):
                    stats_s["max_empty_droplet_count"] = int(re.search("Largest surely-empty droplet has ([0-9]+) UMI counts.", line).group(1))
                # elif re.search("\[epoch ([0-9]+)]  average training loss: ([0-9]+.[0-9]+)", line):
                #     match = re.search("\[epoch ([0-9]+)]  average training loss: ([0-9]+.[0-9]+)", line)
                #     epoch = int(match.group(1))
                #     training_loss = float(match.group(2))
                #     training_data[epoch] = training_loss
                # elif re.search("\[epoch ([0-9]+)] average test loss: ([0-9]+.[0-9]+)", line):
                #     match = re.search("\[epoch ([0-9]+)] average test loss: ([0-9]+.[0-9]+)", line)
                #     epoch = int(match.group(1))
                #     test_loss = float(match.group(2))
                #     test_data[epoch] = test_loss
                elif re.search("Optimal posterior regularization factor = ([0-9]+.[0-9]+)", line):
                    stats_s["optimal_regular_factor"] = float(re.search("Optimal posterior regularization factor = ([0-9]+.[0-9]+)", line).group(1))
                elif "Learning failed.  Retrying with learning-rate " in line:
                    stats_s.drop(labels=[label for label in ["nonzero_genes",
                                                             "prior_counts_empty",
                                                             "prior_counts_cell",
                                                             "barcodes_threshold",
                                                             "probable_cell_barcodes",
                                                             "additional_barcodes",
                                                             "empty_droplets",
                                                             "max_empty_droplet_count",
                                                             "optimal_regular_factor"
                                                             ] if
                                         label in stats_s.index], inplace=True)
                else:
                    pass
        f.close()

        return stats_s

    def run_command(self, command):
        if self.dry_run:
            return

        print("\t" + " ".join(command))
        subprocess.call(command)

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    def visualise_stats(self, stats_df):
        df = stats_df.copy()
        df = df.transpose()
        df["total_barcodes"] = df["probable_cell_barcodes"] + df["additional_barcodes"]

        df.reset_index(drop=False, inplace=True)
        df = df.melt(id_vars=["index"])
        df.columns = ["sample", "variable", "value"]

        df["variable"] = df["variable"].map({
            "overdispersion_mean": "overdispersion (mean)",
            "overdispersion_scale": "overdispersion (scale)",
            "contamination_fraction_rho_alpha": "contamination fraction distribution rho alpha",
            "contamination_fraction_rho_beta": "contamination fraction distribution rho beta",
            "target_false_positive_rate": "target false positive rate",
            "lambda_multiplier": "lambda multiplier",
            "fraction_data_used_for_testing": "fraction data used for testing",
            "nonzero_genes": "genes with nonzero counts",
            "prior_counts_empty": "prior counts in empty droplets",
            "prior_counts_cell": "prior counts for cells",
            "barcodes_threshold": "excluding barcodes <N counts",
            "probable_cell_barcodes": "N probable cell barcodes",
            "additional_barcodes": "N additional barcodes",
            "total_barcodes": "N barcodes",
            "empty_droplets": "N empty droplets",
            "max_empty_droplet_count": "largest surely-empty droplet UMI",
            "optimal_regular_factor": "optimal posterior regularization factor"
        })

        self.barplot(
            df=df,
            panels=[panel for panel in ["overdispersion (mean)",
                                        "overdispersion (scale)",
                                        "contamination fraction distribution rho alpha",
                                        "contamination fraction distribution rho beta",
                                        "target false positive rate",
                                        "lambda multiplier",
                                        "fraction data used for testing",
                                        "genes with nonzero counts",
                                        "prior counts in empty droplets",
                                        "prior counts for cells",
                                        "excluding barcodes <N counts",
                                        "N probable cell barcodes",
                                        "N additional barcodes",
                                        "N barcodes",
                                        "N empty droplets",
                                        "largest surely-empty droplet UMI",
                                        "optimal posterior regularization factor"]
                    if panel in df["variable"].unique()],
            x="sample",
            y="value",
            panel="variable",
            filename="cellbender_remove_background_stats"
        )


    def barplot(self, df, panels, x="x", y="y", panel=None, palette=None,
                ylabel="", title="", filename="plot"):
        nplots = len(panels)
        ncols = math.ceil(np.sqrt(nplots))
        nrows = math.ceil(nplots / ncols)

        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharex='all',
                                 sharey='none',
                                 figsize=(6 * ncols, 6 * nrows))
        sns.set(color_codes=True)

        row_index = 0
        col_index = 0
        for i in range(ncols * nrows):
            if nrows == 1 and ncols == 1:
                ax = axes
            elif nrows == 1 and ncols > 1:
                ax = axes[col_index]
            elif nrows > 1 and ncols == 1:
                ax = axes[row_index]
            else:
                ax = axes[row_index, col_index]

            if i < nplots:
                data = df.loc[df[panel] == panels[i], :]
                if data.shape[0] == 0:
                    continue

                sns.barplot(
                    data=data,
                    x=x,
                    y=y,
                    color="black",
                    palette=palette,
                    ax=ax
                )

                ax.set_xticklabels(ax.get_xmajorticklabels(), rotation=90)
                ax.set_ylim(int(data[y].min() - (data[y].max() * 0.05)), ax.get_ylim()[1])

                ax.set_xlabel("",
                              fontsize=10,
                              fontweight='bold')
                ax.set_ylabel("",
                              fontsize=10,
                              fontweight='bold')
                ax.set_title(panels[i],
                             fontsize=14,
                             fontweight='bold')

            else:
                ax.set_xticklabels(ax.get_xmajorticklabels(), rotation=90)

            col_index += 1
            if col_index > (ncols - 1):
                col_index = 0
                row_index += 1

        fig.suptitle(title,
                     fontsize=40,
                     fontweight='bold')

        plt.tight_layout()
        for extension in self.extensions:
            outpath = os.path.join(self.workdir, "{}.{}".format(filename, extension))
            fig.savefig(outpath)
            print("\tSaved figure: {}".format(os.path.basename(outpath)))
        plt.close()


    @staticmethod
    def build_plot_df(training_df, test_df):
        training_df_dfm = training_df.copy()
        training_df_dfm.reset_index(drop=False, inplace=True)
        training_dfm = training_df_dfm.melt(id_vars=["index"])
        training_dfm["hue"] = "Train"

        test_dfm = test_df.copy()
        test_dfm.reset_index(drop=False, inplace=True)
        test_dfm = test_dfm.melt(id_vars=["index"])
        test_dfm["hue"] = "Test"

        dfm = pd.concat([training_dfm, test_dfm], axis=0)
        dfm.columns = ["epoch", "sample", "loss", "group"]
        dfm["sample"] = dfm["sample"].astype(str)

        return dfm

    def validate_training_procedure(self, training_df, plot_df):
        train_pct_change_df = training_df.pct_change() * 100
        train_pct_change_df[train_pct_change_df < 0] = np.nan
        validation_df = train_pct_change_df.describe().transpose()
        validation_df.sort_values(by=["mean", "std"], inplace=True)
        print(validation_df)

        samples = []
        subtitles = []
        for index, row in validation_df.iterrows():
            samples.append(index)
            subtitles.append(" [{:.2f}% ±{:.2f}]".format(row["mean"], row["std"]))

        self.lineplot_per_sample(
            df=plot_df,
            samples=samples,
            subtitles=subtitles,
            x="epoch",
            y="loss",
            hue="group",
            xlabel="Epoch",
            ylabel="ELBO",
            filename="cellbender_remove_background_training_procedure_per_sample"
        )

        return validation_df

    def lineplot_per_sample(self, df, samples, subtitles=None, x="x", y="y",
                            sample_column="sample", hue=None, palette=None,
                            xlabel="", ylabel="", filename="plot"):

        if subtitles is not None and len(subtitles) != len(samples):
            print("Error, subtitles are not the same length as the samples")
            return

        nplots = len(samples)
        ncols = math.ceil(np.sqrt(nplots))
        nrows = math.ceil(nplots / ncols)

        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharex='none',
                                 sharey='none',
                                 figsize=(6 * ncols, 6 * nrows))
        sns.set(color_codes=True)

        row_index = 0
        col_index = 0
        for i in range(ncols * nrows):
            if nrows == 1 and ncols == 1:
                ax = axes
            elif nrows == 1 and ncols > 1:
                ax = axes[col_index]
            elif nrows > 1 and ncols == 1:
                ax = axes[row_index]
            else:
                ax = axes[row_index, col_index]

            if i < nplots:
                data = df.loc[df[sample_column] == samples[i], :]
                if data.shape[0] == 0:
                    continue

                subtitle = None
                if subtitles is not None:
                    subtitle = subtitles[i]

                self.lineplot(
                    fig=fig,
                    ax=ax,
                    data=data,
                    x=x,
                    y=y,
                    hue=hue,
                    palette=palette,
                    xlabel=xlabel,
                    ylabel=ylabel,
                    title="{}{}".format(samples[i], subtitle)
                )

            else:
                ax.set_axis_off()

            col_index += 1
            if col_index > (ncols - 1):
                col_index = 0
                row_index += 1

        plt.tight_layout()
        for extension in self.extensions:
            outpath = os.path.join(self.workdir, "{}.{}".format(filename, extension))
            fig.savefig(outpath)
            print("\tSaved figure: {}".format(os.path.basename(outpath)))
        plt.close()

    def lineplot(self, fig, ax, data, x, y, hue, palette=None,
                 xlabel="", ylabel="", title=""):

        sns.despine(fig=fig, ax=ax)

        sns.lineplot(data=data,
                     x=x,
                     y=y,
                     hue=hue,
                     style=hue,
                     markers=True,
                     palette=palette,
                     ax=ax)

        ax.set_xlabel(xlabel,
                       fontsize=10,
                       fontweight='bold')
        ax.set_ylabel(ylabel,
                       fontsize=10,
                       fontweight='bold')
        ax.set_title(title,
                      fontsize=14,
                      fontweight='bold')

        ax.invert_yaxis()
        xaxis_labels = np.arange(min(data[x]), max(data[x]), 20)
        ax.xaxis.set_ticks(xaxis_labels)
        ax.xaxis.set_ticklabels(xaxis_labels - 1, rotation=0)

    def visualise_training_procedure(self, plot_df):
        sns.set()
        sns.set_style("ticks")
        fig, ax = plt.subplots(figsize=(6, 6))

        self.lineplot(
            fig=fig,
            ax=ax,
            data=plot_df,
            x="epoch",
            y="loss",
            hue="group",
            xlabel="Epoch",
            ylabel="ELBO",
            title="Progress of the training procedure"
        )

        plt.tight_layout()
        for extension in self.extensions:
            outpath = os.path.join(self.workdir, "cellbender_remove_background_training_procedure.{}".format(extension))
            fig.savefig(outpath)
            print("\tSaved figure: {}".format(os.path.basename(outpath)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        for name, value, default in self.defaults:
            if value != default:
                print("  > {}: {}".format(name, value))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()