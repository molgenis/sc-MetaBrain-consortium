#!/usr/bin/env python3

"""
File:         run_cellbender.py
Created:      2023/03/23
Last Changed: 2023/07/07
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
import math
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
### Zhou 2020 ###
./run_cellbender.py \
    --workdir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-08-28-CellBender-v0.3.0/2023-09-10-Zhou2020-test \
    --inputdir /groups/umcg-biogen/tmp02/input/processeddata/single-cell/Zhou2020/ \
    --gres gpu:a40:1 \
    --cuda \
    --dry_run
    
./run_cellbender.py \
    --workdir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-08-28-CellBender-v0.3.0/2023-09-10-Zhou2020-ExtendedEpochLowLearn \
    --inputdir /groups/umcg-biogen/tmp02/input/processeddata/single-cell/Zhou2020/ \
    --gres gpu:a40:1 \
    --epochs 300 \
    --cuda \
    --learning_rate 1e-5 \
    --dry_run

./run_cellbender.py \
    --workdir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-08-28-CellBender-v0.3.0/2023-09-10-Zhou2020-ExtendedEpochLowLearn-CellRangerExpectedCells \
    --inputdir /groups/umcg-biogen/tmp02/input/processeddata/single-cell/Zhou2020/ \
    --gres gpu:a40:1 \
    --expected_cells -1 \
    --epochs 300 \
    --cuda \
    --learning_rate 1e-5 \
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
        self.gres = getattr(arguments, 'gres')
        time = getattr(arguments, 'time')
        self.time = TIME_DICT[time] if time in TIME_DICT else time
        self.cpus_per_task = getattr(arguments, 'cpus_per_task')
        self.mem = getattr(arguments, 'mem')
        preflights = getattr(arguments, 'preflights')
        self.dry_run = getattr(arguments, 'dry_run')
        self.rerun = getattr(arguments, 'rerun')

        # Format preflights.
        self.preflights = []
        if preflights:
            command = []
            for word in preflights:
                if word == ",":
                    self.preflights.append(" ".join(command))
                    command = []
                else:
                    command.append(word)
            self.preflights.append(" ".join(command))

        # Safe the CellBender arguments.
        self.cuda = getattr(arguments, 'cuda')
        self.expected_cells = getattr(arguments, 'expected_cells')
        self.total_droplets_included = getattr(arguments, 'total_droplets_included')
        self.force_cell_umi_prior = getattr(arguments, 'force_cell_umi_prior')
        self.force_empty_umi_prior = getattr(arguments, 'force_empty_umi_prior')
        self.model = getattr(arguments, 'model')
        self.epochs = getattr(arguments, 'epochs')
        self.low_count_threshold = getattr(arguments, 'low_count_threshold')
        self.z_dim = getattr(arguments, 'z_dim')
        self.z_layers = getattr(arguments, 'z_layers')
        self.training_fraction = getattr(arguments, 'training_fraction')
        self.empty_drop_training_fraction = getattr(arguments, 'empty_drop_training_fraction')
        self.ignore_features = getattr(arguments, 'ignore_features')
        self.fpr = getattr(arguments, 'fpr')
        self.exclude_feature_types = getattr(arguments, 'exclude_feature_types')
        self.projected_ambient_count_threshold = getattr(arguments, 'projected_ambient_count_threshold')
        self.learning_rate = getattr(arguments, 'learning_rate')
        self.checkpoint_mins = getattr(arguments, 'checkpoint_mins')
        self.final_elbo_fail_fraction = getattr(arguments, 'final_elbo_fail_fraction')
        self.epoch_elbo_fail_fraction = getattr(arguments, 'epoch_elbo_fail_fraction')
        self.num_training_tries = getattr(arguments, 'num_training_tries')
        self.learning_rate_retry_mult = getattr(arguments, 'learning_rate_retry_mult')
        self.posterior_batch_size = getattr(arguments, 'posterior_batch_size')
        self.posterior_regularization = getattr(arguments, 'posterior_regularization')
        self.alpha = getattr(arguments, 'alpha')
        self.q = getattr(arguments, 'q')
        self.estimator = getattr(arguments, "estimator")
        self.estimator_multiple_cpu = getattr(arguments, 'estimator_multiple_cpu')
        self.constant_learning_rate = getattr(arguments, "constant_learning_rate")
        self.cpu_threads = getattr(arguments, 'cpu_threads')

        if self.partition == "gpu" and self.gpus_per_node is None:
            self.gpus_per_node = "a100.20gb:1"

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
            ("gres", self.gres, "required"),
            ("time", self.time, "required"),
            ("cpus-per-task", self.cpus_per_task, "required"),
            ("preflights", self.preflights, "required"),
            ("mem", self.mem, "required"),
            ("dry_run", self.dry_run, "required"),
            ("rerun", self.rerun, "required"),
            ("cuda", self.cuda, False),
            ("expected-cells", self.expected_cells, None),
            ("total-droplets-included", self.total_droplets_included, None),
            ("force-cell-umi-prior", self.force_cell_umi_prior, None),
            ("force-empty-umi-prior", self.force_empty_umi_prior, None),
            ("model", self.model, "full"),
            ("epochs", self.epochs, 150),
            ("low-count-threshold", self.low_count_threshold, 5),
            ("z-dim", self.z_dim, 64),
            ("z-layers", self.z_layers, [512]),
            ("training-fraction", self.training_fraction, 0.9),
            ("empty-drop-training-fraction", self.empty_drop_training_fraction, 0.2),
            ("ignore-features", self.ignore_features, []),
            ("fpr", self.fpr, [0.01]),
            ("exclude-feature-types", self.exclude_feature_types, []),
            ("projected-ambient-count-threshold", self.projected_ambient_count_threshold, 0.1),
            ("learning-rate", self.learning_rate, 1e-4),
            ("checkpoint-mins", self.checkpoint_mins, 7),
            ("final-elbo-fail-fraction", self.final_elbo_fail_fraction, None),
            ("epoch-elbo-fail-fraction", self.epoch_elbo_fail_fraction, None),
            ("num-training-tries", self.num_training_tries, 1),
            ("learning-rate-retry-mult", self.learning_rate_retry_mult, 0.2),
            ("posterior-batch-size", self.posterior_batch_size, 128),
            ("posterior-regularization", self.posterior_regularization, None),
            ("alpha", self.alpha, None),
            ("q", self.q, None),
            ("estimator", self.estimator, "mckp"),
            ("estimator-multiple-cpu", self.estimator_multiple_cpu, False),
            ("constant-learning-rate", self.constant_learning_rate, False),
            ("cpu-threads", self.cpu_threads, None)
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
        parser.add_argument("--bind_paths",
                            nargs="*",
                            type=str,
                            default=["/groups/umcg-biogen/tmp02/",
                                     "/groups/umcg-biogen/tmp02/users/umcg-mvochteloo/simulated_home:/home/umcg-mvochteloo"],
                            help="List of paths to bind to Singularity. "
                                 "Default: ['/groups/umcg-biogen/tmp01/']")
        parser.add_argument("--singularity",
                            type=str,
                            required=False,
                            default="/groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-08-28-CellBender-v0.3.0/cellbender-v0.3.0.simg",
                            help="")
        parser.add_argument("--partition",
                            type=str,
                            required=False,
                            choices=["regular", "parallel", "gpu", "himem",
                                     "gelifes", "gpu_a40"],
                            default=None,
                            help="The partition to submit to.")
        parser.add_argument("--gpus_per_node",
                            type=str,
                            required=False,
                            choices=["a100.20gb:1", "a100:1", "a100:2"],
                            default=None,
                            help="Request a specific GPU type on Habrok when"
                                 " using --partition=gpu.")
        parser.add_argument("--gres",
                            type=str,
                            required=False,
                            choices=["gpu:a40:1"],
                            default=None,
                            help="Request a specific GPU type on Nibbler when "
                                 "using --partition=gpu.")
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
                            default=-1,
                            help="Restricts CellBender to use specified amount "
                                 "of memory (in GB). (default: -1)")
        parser.add_argument("--preflights",
                            nargs="+",
                            type=str,
                            default=[],
                            help="")
        parser.add_argument("--dry_run",
                            action='store_true',
                            help="Only create the job files, don't submit them.")
        parser.add_argument("--rerun",
                            nargs="+",
                            type=str,
                            required=False,
                            default=None,
                            help="The name of sample to run again.")

        parser.add_argument("--cuda",
                            action="store_true",
                            help="Including the flag --cuda will run the "
                                 "inference on a GPU.")
        parser.add_argument("--checkpoint", # new in v0.3.0
                            type=str,
                            default='ckpt.tar.gz',
                            help="Checkpoint tarball produced by the same version "
                                 "of CellBender remove-background.  If present, "
                                 "and the workflow hashes match, training will "
                                 "restart from this checkpoint.")
        parser.add_argument("--expected_cells",
                            type=int,
                            default=None,
                            help="Number of cells expected in the dataset "
                                 "(a rough estimate within a factor of 2 "
                                 "is sufficient).")
        parser.add_argument("--total_droplets_included", # no default 25k in v0.3.0
                            type=int,
                            default=None,
                            help="The number of droplets from the "
                                 "rank-ordered UMI plot that will have their "
                                 "cell probabilities inferred as an output. "
                                 "Include the droplets which might contain "
                                 "cells. Droplets beyond TOTAL_DROPLETS_INCLUDED "
                                 "should be 'surely empty' droplets.")
        parser.add_argument("--force_cell_umi_prior", # new in v0.3.0
                            type=float,
                            default=None,
                            help="Ignore CellBender's heuristic prior estimation, "
                                 "and use this prior for UMI counts in cells.")
        parser.add_argument("--force_empty_umi_prior", # new in v0.3.0
                            type=float,
                            default=None,
                            help="Ignore CellBender's heuristic prior estimation, "
                                 "and use this prior for UMI counts in empty droplets.")
        parser.add_argument("--model", # naive added in v0.3.0
                            type=str,
                            default="full",
                            choices=["naive", "simple", "ambient", "swapping", "full"],
                            help="Which model is being used for count data. "
                                 "'naive' subtracts the estimated ambient profile."
                                 " 'simple' does not model either ambient "
                                 "RNA or random barcode swapping (for "
                                 "debugging purposes -- not recommended).  "
                                 "'ambient' assumes background RNA is "
                                 "incorporated into droplets.  'swapping' "
                                 "assumes background RNA comes from random "
                                 "barcode swapping (via PCR chimeras).  'full' "
                                 "uses a combined ambient and swapping model.")
        parser.add_argument("--epochs",
                            type=int,
                            default=150,
                            help="Number of epochs to train.")
        parser.add_argument("--low_count_threshold", # default decreased from 15 to 5 in v0.3.0
                            type=int,
                            default=5,
                            help="Droplets with UMI counts below this "
                                 "number are completely excluded from the "
                                 "analysis.  This can help identify the "
                                 "correct prior for empty droplet counts "
                                 "in the rare case where empty counts "
                                 "are extremely high (over 200)."
                                 "(default: 5)")
        parser.add_argument("--z_dim", # default decreased from 100 to 64 in v0.3.0
                            type=int,
                            default=64,
                            help="Dimension of latent variable z.")
        parser.add_argument("--z_layers", # default increased from 500 to 512 in v0.3.0
                            nargs="+",
                            type=int,
                            default=[512],
                            help="Dimension of hidden layers in the encoder for z.")
        parser.add_argument("--training_fraction",
                            type=float,
                            default=0.9,
                            help="Training detail: the fraction of the "
                                 "data used for training. The rest is never "
                                 "seen by the inference algorithm. Speeds up "
                                 "learning. (default: 0.9)")
        parser.add_argument("--empty_drop_training_fraction",
                            type=float,
                            default=0.2,
                            help="Training detail: the fraction of the "
                                 "training data each epoch "
                                 "that is drawn (randomly sampled) from "
                                 "surely empty droplets. (default: 0.2)")
        parser.add_argument("--ignore_features", # renamed from blaklist-genes in v0.3.0
                            nargs="+",
                            type=int,
                            default=[],
                            help="Integer indices of features to ignore "
                                 "entirely.  In the output count matrix, the "
                                 "counts for these features will be unchanged.")
        parser.add_argument("--fpr",
                            nargs="+",
                            default=[0.01],
                            help="Target 'delta' false positive rate in [0, 1). "
                                 "Use 0 for a cohort of samples which will be "
                                 "jointly analyzed for differential expression. "
                                 "A false positive is a true signal count that is "
                                 "erroneously removed.  More background removal "
                                 "is accompanied by more signal removal "
                                 "at high values of FPR.  You can specify "
                                 "multiple values, which will create multiple "
                                 "output files.")
        parser.add_argument("--exclude_feature_types",
                            type=str,
                            nargs="+",
                            default=[],
                            help="Feature types to ignore during the analysis.  "
                                 "These features will be left unchanged in the "
                                 "output file.")
        parser.add_argument("--projected_ambient_count_threshold", # new in v0.3.0
                            type=float,
                            default=0.1,
                            help="Controls how many features are included in the "
                                 "analysis, which can lead to a large speedup. "
                                 "If a feature is expected to have less than "
                                 "PROJECTED_AMBIENT_COUNT_THRESHOLD counts total "
                                 "in all cells (summed), then that gene is "
                                 "excluded, and it will be unchanged in the "
                                 "output count matrix.  For example, "
                                 "PROJECTED_AMBIENT_COUNT_THRESHOLD = 0 will "
                                 "include all features which have even a single "
                                 "count in any empty droplet.")
        parser.add_argument("--learning_rate",
                            type=float,
                            default=1e-4,
                            help="Training detail: lower learning rate for "
                                 "inference. A OneCycle learning rate schedule "
                                 "is used, where the upper learning rate is ten "
                                 "times this value. (For this value, probably "
                                 "do not exceed 1e-3).")
        parser.add_argument("--checkpoint_mins", # new in v0.3.0
                            type=float,
                            default=7.,
                            help="Checkpoint file will be saved periodically, "
                                 "with this many minutes between each checkpoint.")
        parser.add_argument("--final_elbo_fail_fraction",
                            type=float,
                            help="Training is considered to have failed if "
                                 "(best_test_ELBO - final_test_ELBO)/(best_test_ELBO "
                                 "- initial_test_ELBO) > FINAL_ELBO_FAIL_FRACTION.  Training will "
                                 "automatically re-run if --num-training-tries > "
                                 "1.  By default, will not fail training based "
                                 "on final_training_ELBO.")
        parser.add_argument("--epoch_elbo_fail_fraction",
                            type=float,
                            help="Training is considered to have failed if "
                                 "(previous_epoch_test_ELBO - current_epoch_test_ELBO)"
                                 "/(previous_epoch_test_ELBO - initial_train_ELBO) "
                                 "> EPOCH_ELBO_FAIL_FRACTION.  Training will "
                                 "automatically re-run if --num-training-tries > "
                                 "1.  By default, will not fail training based "
                                 "on epoch_training_ELBO.")
        parser.add_argument("--num_training_tries",
                            type=int,
                            default=1,
                            help="Number of times to attempt to train the model.  "
                                 "At each subsequent attempt, the learning rate is "
                                 "multiplied by LEARNING_RATE_RETRY_MULT.")
        parser.add_argument("--learning_rate_retry_mult", # default decreased from 0.5 to 0.2 in v0.3.0
                            type=float,
                            default=0.2,
                            help="Learning rate is multiplied by this amount each "
                                 "time a new training attempt is made.  (This "
                                 "parameter is only used if training fails based "
                                 "on EPOCH_ELBO_FAIL_FRACTION or "
                                 "FINAL_ELBO_FAIL_FRACTION and NUM_TRAINING_TRIES"
                                 " is > 1.)")
        parser.add_argument("--posterior_batch_size", # default increased from 20 to 128
                            type=int,
                            default=128,
                            help="Training detail: size of batches when creating "
                                 "the posterior.  Reduce this to avoid running "
                                 "out of GPU memory creating the posterior "
                                 "(will be slower).")
        parser.add_argument("--posterior_regularization", # new in v0.3.0
                            type=str,
                            default=None,
                            choices=["PRq", "PRmu", "PRmu_gene"],
                            help="Posterior regularization method. (For experts: "
                                 "not required for normal usage, see "
                                 "documentation). PRq is approximate quantile-"
                                 "targeting. PRmu is approximate mean-targeting "
                                 "aggregated over genes (behavior of v0.2.0). "
                                 "PRmu_gene is approximate mean-targeting per "
                                 "gene.")
        parser.add_argument("--alpha", # new in v0.3.0
                            type=float,
                            default=None,
                            help="Tunable parameter alpha for the PRq posterior "
                                 "regularization method (not normally used: see "
                                 "documentation).")
        parser.add_argument("--q", # new in v0.3.0
                            type=float,
                            default=None,
                            help="Tunable parameter q for the CDF threshold "
                                 "estimation method (not normally used: see "
                                 "documentation).")
        parser.add_argument("--estimator", # new in v0.3.0
                            type=str,
                            default="mckp",
                            choices=["map", "mean", "cdf", "sample", "mckp"],
                            help="Output denoised count estimation method. (For "
                                 "experts: not required for normal usage, see "
                                 "documentation).")
        parser.add_argument("--estimator_multiple_cpu", # new in v0.3.0
                            action="store_true",
                            default=False,
                            help="Including the flag --estimator-multiple-cpu will "
                                 "use more than one CPU to compute the MCKP "
                                 "output count estimator in parallel (does "
                                 "nothing for other estimators).")
        parser.add_argument("--constant_learning_rate", # new in v0.3.0
                            action="store_true",
                            default=False,
                            help="Including the flag --constant-learning-rate will "
                                 "use the ClippedAdam optimizer instead of the "
                                 "OneCycleLR learning rate schedule, which is "
                                 "the default.  Learning is faster with the "
                                 "OneCycleLR schedule.  However, training can "
                                 "easily be continued from a checkpoint for more "
                                 "epochs than the initial command specified when "
                                 "using ClippedAdam.  On the other hand, if using "
                                 "the OneCycleLR schedule with 150 epochs "
                                 "specified, it is not possible to pick up from "
                                 "that final checkpoint and continue training "
                                 "until 250 epochs.")
        parser.add_argument("--cpu_threads", # new in v0.3.0
                            type=int,
                            default=None,
                            help="Number of threads to use when pytorch is run "
                                 "on CPU. Defaults to the number of logical cores.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Creating job files.")
        metrics_data = []
        arguments = self.filter_arguments()
        info = []
        aggregation_data = []
        for path in glob.glob(os.path.join(self.inputdir, "*")):
            folder = os.path.basename(path)
            if not os.path.isdir(path):
                continue

            input = self.find_input(folder)
            if input is None:
                print("\t\tWarning, no input file found for '{}'.".format(folder))
                continue

            metrics_path = self.find_metrics(folder)
            if metrics_path is None:
                print("\t\tWarning, no metrics file found for '{}'.".format(folder))
                continue

            # Loading CellRanger metrics summary file.
            metrics_df = self.load_file(metrics_path, header=0, index_col=None)
            metrics_df = metrics_df.replace(',', '', regex=True).replace('%', '', regex=True).astype(float)
            metrics_df.index = [folder]
            metrics_data.append(metrics_df)

            mem = self.mem
            if self.mem == -1:
                # Using a predicted memory usage based on the CellRanger number of cells.
                cellranger_expected_cells = int(metrics_df.loc[folder, "Estimated Number of Cells"])
                mem = int(math.ceil(3.5 + cellranger_expected_cells * 0.0015))

            if self.expected_cells == -1:
                # Using CellRanger expected number of cells.
                cellranger_expected_cells = int(metrics_df.loc[folder, "Estimated Number of Cells"])
                arguments["expected-cells"] = cellranger_expected_cells

            # Create output.
            outpath = os.path.join(self.workdir, folder)
            if not os.path.exists(outpath):
                os.makedirs(outpath)
            output = os.path.join(self.workdir, folder, "cellbender_remove_background_output.h5")
            aggregation_data.append([folder, output])

            # Check if sample is being rerun.
            sample_rerun = False
            if os.path.exists(output) and ((self.rerun is not None) and (folder in self.rerun)):
                sample_rerun = True
                arguments["learning-rate"] = self.learning_rate * self.learning_rate_retry_mult

            # Create job files.
            jobfile_path, logfile_path = self.create_job_file(sample=folder,
                                                              mem=mem,
                                                              input=input,
                                                              output=output,
                                                              arguments=arguments)
            info.append((jobfile_path, logfile_path, folder, sample_rerun))

        print("Saving aggregated files")
        aggregation_df = pd.DataFrame(aggregation_data, columns=["sample_id", "molecule_h5"])
        print(aggregation_df)
        self.save_file(df=aggregation_df,
                       outpath=os.path.join(self.workdir, "aggregation.csv"),
                       sep=",",
                       index=False)

        metrics_df = pd.concat(metrics_data, axis=0)
        print(metrics_df)
        self.save_file(df=metrics_df, outpath=os.path.join(self.workdir, "metrics_summary.txt.gz"))

        print("Starting job files.")
        for jobfile_path, logfile_path, folder, sample_rerun in info:
            if os.path.exists(logfile_path):
                success, missing = self.parse_logfile(logfile_path)
                if success:
                    if sample_rerun:
                        print("\tRerunning sample '{}' with adjusted learning rate".format(folder))
                    else:
                        print("\tSample '{}' already completed successfully!".format(folder))
                        continue
                else:
                    if len(missing) != 0:
                        print("\tSample '{}' failed due to missing '{}'.".format(folder, ", ".join(missing)))
                        continue
                    else:
                        print("\tSample '{}' failed. Rerunning sample.".format(folder))

            command = ['sbatch', jobfile_path]
            self.run_command(command)

    def filter_arguments(self):
        arguments = {}
        for name, value, default in self.defaults:
            if value != default and default != "required":
                arguments[name] = value
        return arguments

    def find_input(self, folder):
        raw_feature_bc_matrix_path = os.path.join(self.inputdir, folder, "raw_feature_bc_matrix.h5")
        if os.path.exists(raw_feature_bc_matrix_path):
            return raw_feature_bc_matrix_path

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

    def find_metrics(self, folder):
        metrics_outs_path = os.path.join(self.inputdir, folder, "outs", "metrics_summary.csv")
        if os.path.exists(metrics_outs_path):
            return metrics_outs_path

        return None

    @staticmethod
    def load_file(inpath, header=0, index_col=None, sep=","):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    def create_job_file(self, sample, mem, input, output, arguments):
        job_name = "cellbender_remove_background_{}".format(sample)
        logfile_path = os.path.join(self.jobs_outdir, job_name + ".out")

        lines = ["#!/bin/bash",
                 "#SBATCH --job-name={}".format(job_name),
                 "#SBATCH --output={}".format(logfile_path),
                 "#SBATCH --error={}".format(logfile_path),
                 "#SBATCH --time={}".format(self.time),
                 "#SBATCH --cpus-per-task={}".format(self.cpus_per_task),
                 "#SBATCH --mem={}gb".format(mem),
                 "#SBATCH --nodes=1",
                 "#SBATCH --open-mode=append",
                 "#SBATCH --export=NONE",
                 "#SBATCH --get-user-env=L",
                 "",
                 "cd {} || exit".format(os.path.join(self.workdir, sample)),
                 "",
                 "singularity exec \\",
                 "  --nv \\",
                 "  --bind /groups/umcg-biogen/tmp02/,/groups/umcg-biogen/tmp02/users/umcg-mvochteloo/simulated_home:/home/umcg-mvochteloo \\",
                 "  /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-08-28-CellBender-v0.3.0/cellbender-v0.3.0.simg \\",
                 "  cellbender remove-background \\",
                 "    --input={} \\".format(input),
                 "    --output={} \\".format(output)
                 ]

        insert_index = 5
        if self.partition is not None:
            lines.insert(insert_index, "#SBATCH --partition={}".format(self.partition))
            insert_index += 1
        if self.cuda is not None and self.gpus_per_node is not None:
            lines.insert(insert_index, "#SBATCH --gpus-per-node={}".format(self.gpus_per_node))
            insert_index += 1
        if self.cuda is not None and self.gres is not None:
            lines.insert(insert_index, "#SBATCH --gres={}".format(self.gres))
            insert_index += 1

        for command in self.preflights:
            lines.insert(insert_index + 7, command)
            insert_index += 1

        for name, value in arguments.items():
            if isinstance(value, bool):
                lines.append("    --{} \\".format(name))
            else:
                lines.append("    --{}={} \\".format(name, value))

        if lines[-1].endswith("\\"):
            lines[-1] = lines[-1].strip("\\")

        lines.extend(["", "echo 'Job finished'"])

        for command in self.preflights:
            if "conda" in command and "activate" in command:
                lines.extend(["", "conda deactivate"])
            elif "activate" in command:
                lines.extend(["", "deactivate"])

        jobfile_path = os.path.join(self.jobdir, job_name + ".sh")
        with open(jobfile_path, "w") as f:
            for line in lines:
                f.write(line + "\n")
        f.close()
        print("\tSaved jobfile: {}".format(os.path.basename(jobfile_path)))
        return jobfile_path, logfile_path

    @staticmethod
    def parse_logfile(logfile_path):
        inference_completed = False
        posterior_added = False
        saved_summary_plots = False
        saved_barcodes = False
        saved_output = 0
        saved_metrics = False
        saved_report = False
        completed = False
        with open(logfile_path, 'r') as f:
            for line in f:
                if 'Running remove-background' in line:
                    inference_completed = False
                    posterior_added = False
                    saved_summary_plots = False
                    saved_barcodes = False
                    saved_output = 0
                    saved_metrics = False
                    saved_report = False
                    completed = False
                if 'Inference procedure complete.' in line:
                    inference_completed = True
                if 'Added posterior object to checkpoint file.' in line:
                    posterior_added = True
                if 'Saved summary plots as' in line:
                    saved_summary_plots = True
                if 'Saved cell barcodes in' in line:
                    saved_barcodes = True
                if 'Succeeded in writing CellRanger format output to file' in line:
                    saved_output += 1
                if 'Saved output metrics as' in line:
                    saved_metrics = True
                if 'Succeeded in writing report to' in line:
                    saved_report = True
                if 'Completed remove-background.' in line:
                    completed = True
        f.close()

        missing = [label for variable, label in [(inference_completed, "inference"),
                                                 (posterior_added, "posterior"),
                                                 (saved_summary_plots, "summary"),
                                                 (saved_barcodes, "barcodes"),
                                                 (saved_output == 2, "output"),
                                                 (saved_metrics, "metrics"),
                                                 (saved_report, "report")] if not variable]

        return inference_completed and posterior_added and saved_summary_plots and saved_barcodes and saved_output == 2 and saved_metrics and saved_report and completed, missing

    def run_command(self, command):
        print("\t" + " ".join(command))
        if self.dry_run:
            print("\t\tNot executed due to dry run")
            return

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

    def print_arguments(self):
        print("Arguments:")
        for name, value, default in self.defaults:
            if value != default:
                print("  > {}: {}".format(name, value))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()