#!/usr/bin/env python3

"""
File:         visualise_cellbender_results.py
Created:      2023/07/07
Last Changed: 2023/07/11
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
import glob
import math
import os
import re

# Third party imports.
import numpy as np
import pandas as pd
import h5py
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Visualise CellBender Results"
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

### Mathys2019 ###
./visualise_cellbender_results.py \
    --workdir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-03-23-CellBender/2023-03-23-Mathys2019
    
./visualise_cellbender_results.py \
    --workdir /scratch/p301710/2023-03-28-scMetaBrainConsortium/output/2023-03-28-CellBender/2023-03-28-Mathys2019
    
./visualise_cellbender_results.py \
    --workdir /scratch/p301710/2023-03-28-scMetaBrainConsortium/output/2023-03-28-CellBender/2023-07-07-Mathys2019


## Zhou 2020 ###
./visualise_cellbender_results.py \
    --workdir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-03-28-CellBender/2023-07-07-Zhou2020

./visualise_cellbender_results.py \
    --workdir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-03-28-CellBender/2023-07-10-Zhou2020-CellRangerExpectedCells

"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.workdir = getattr(arguments, 'workdir')
        self.extensions = getattr(arguments, 'extension')

        # Set variables.
        self.plot_outdir = os.path.join(self.workdir, "plot")
        if not os.path.exists(self.plot_outdir):
            os.makedirs(self.plot_outdir)

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
        parser.add_argument("--extension",
                            nargs="+",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        stats_data = []
        training_data = []
        test_data = []
        cell_prob_data = []
        latent_gene_encoding_data = []

        print("Loading results.")
        for path in glob.glob(os.path.join(self.workdir, "*")):
            folder = os.path.basename(path)
            output = os.path.join(self.workdir, folder, "cellbender_remove_background_output.h5")
            if not os.path.exists(output):
                continue
            print("\tProcessing '{}'".format(folder))

            resultfile = os.path.join(self.workdir, folder, "cellbender_remove_background_output.h5")
            if not os.path.exists(resultfile):
                continue
            stats_s, training_s, test_s, cell_prob_s, latent_gene_encoding_df = self.parse_results(resultfile, folder)
            stats_s = pd.Series({"passed": True})
            stats_s.name = folder

            logfile = os.path.join(self.workdir, folder, "cellbender_remove_background_output.log")
            if os.path.exists(logfile):
                stats_s = self.add_stats_from_log(stats_s, logfile)

            barcodesfile = os.path.join(self.workdir, folder, "cellbender_remove_background_output_cell_barcodes.csv")
            if os.path.exists(barcodesfile):
                stats_s["barcodes"] = self.count_barcodes(barcodesfile)

            if not stats_s["passed"]:
                print("\t  Error in log file! Please rerun.")
                continue

            stats_data.append(stats_s)
            training_data.append(training_s)
            test_data.append(test_s)
            cell_prob_data.append(cell_prob_s)
            latent_gene_encoding_data.append(latent_gene_encoding_df)

        if len(stats_data) == 0:
            exit()

        print("Summarizing results.")
        stats_df = pd.concat(stats_data, axis=1)
        self.save_file(df=stats_df, outpath=os.path.join(self.workdir, "cellbender_remove_background_stats.txt.gz"))

        training_df = pd.concat(training_data, axis=1)
        self.save_file(df=training_df, outpath=os.path.join(self.workdir, "cellbender_remove_background_training_loss.txt.gz"))

        test_df = pd.concat(test_data, axis=1)
        self.save_file(df=test_df, outpath=os.path.join(self.workdir, "cellbender_remove_background_test_loss.txt.gz"))

        cell_prob_df = pd.concat(cell_prob_data, axis=1)
        self.save_file(df=cell_prob_df, outpath=os.path.join(self.workdir, "cellbender_remove_background_latent_cell_probability.txt.gz"))

        latent_gene_encoding_df = pd.concat(latent_gene_encoding_data, axis=0)
        self.save_file(df=latent_gene_encoding_df, outpath=os.path.join(self.workdir, "cellbender_remove_background_latent_gene_encoding.txt.gz"))

        print("Visualising results.")
        validation_df, sample_order = self.visualise_training_procedure(training_df, test_df)
        self.save_file(df=validation_df, outpath=os.path.join(self.workdir, "cellbender_remove_background_training_validation.txt.gz"))

        self.visualise_stats(stats_df, sample_order=sample_order)
        self.visualise_cell_determination(cell_prob_df=cell_prob_df, sample_order=sample_order)
        self.visualise_latent_gene_encoding(latent_gene_encoding_df=latent_gene_encoding_df, sample_order=sample_order)

    @staticmethod
    def parse_results(resultsfile, folder):
        hf = h5py.File(resultsfile, 'r')

        # stats
        stats_s = pd.Series({"passed": False,
                             "overdispersion_mean": float(np.array(hf.get('matrix/overdispersion_mean_and_scale'))[0]),
                             "overdispersion_scale": float(np.array(hf.get('matrix/overdispersion_mean_and_scale'))[1]),
                             "contamination_fraction_rho_alpha": float(np.array(hf.get('matrix/contamination_fraction_params'))[0]),
                             "contamination_fraction_rho_beta": float(np.array(hf.get('matrix/contamination_fraction_params'))[1]),
                             "lambda_multiplier": float(np.array(hf.get('matrix/lambda_multiplier'))),
                             "target_false_positive_rate": float(np.array(hf.get('matrix/target_false_positive_rate'))),
                             "fraction_data_used_for_testing": float(np.array(hf.get('matrix/fraction_data_used_for_testing')))
                             })
        stats_s.name= folder

        # training loss
        training_s = pd.Series(np.array(hf.get('matrix/training_elbo_per_epoch')))
        training_s.index = np.arange(0, training_s.shape[0], 1)
        training_s.name = folder

        # test loss
        test_s = pd.Series(np.array(hf.get('matrix/test_elbo')))
        test_s.index = pd.Series(np.array(hf.get('matrix/test_epoch')))
        test_s.name = folder

        # cel probabilities.
        p = np.array(hf.get('matrix/latent_cell_probability'))
        cell_prob_s = pd.Series(p)
        cell_prob_s.name = folder

        # latent encoding.
        z = np.array(hf.get('matrix/latent_gene_encoding'))
        d = np.array(hf.get('matrix/latent_scale'))
        pca = PCA(n_components=2)
        if p is None:
            p = np.ones_like(d)
        z_pca = pca.fit_transform(z[p >= 0.5])
        latent_gene_encoding_df = pd.DataFrame(z_pca[:, :2], columns=["PC0", "PC1"])
        latent_gene_encoding_df["sample"] = folder

        hf.close()

        return stats_s, training_s, test_s, cell_prob_s, latent_gene_encoding_df


    @staticmethod
    def add_stats_from_log(stats_s, logfile):
        input_keys = stats_s.keys()
        start_datetime = None
        end_datetime = None
        with open(logfile, 'r') as f:
            for line in f:
                if line.startswith("cellbender remove-background"):
                    for part in line.split("--")[1:]:
                        if "=" in part:
                            key, value = part.split("=")
                            try:
                                stats_s["parameter:" + key] = float(value)
                            except ValueError:
                                pass
                        else:
                            stats_s[part] = 1
                elif "Training succeeded" in line:
                    stats_s["passed"] = True
                elif re.search("(\d{4})-(\d{2})-(\d{2}) (\d{2}):(\d{2}):(\d{2})", line):
                    tmp_datetime = datetime.strptime(re.search("(\d{4})-(\d{2})-(\d{2}) (\d{2}):(\d{2}):(\d{2})", line).group(0), "%Y-%m-%d %H:%M:%S")
                    if start_datetime is None:
                        start_datetime = tmp_datetime
                    else:
                        end_datetime = tmp_datetime
                elif re.search("Including ([0-9]+) genes that have nonzero counts.", line):
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
                    stats_s.drop(labels=[label for label in stats_s.index if label not in input_keys], inplace=True)
                    start_datetime = end_datetime
                    end_datetime = None
                else:
                    pass
        f.close()

        if start_datetime is not None and end_datetime is not None:
            stats_s["time"] = (end_datetime - start_datetime).total_seconds() / 60.0

        return stats_s

    @staticmethod
    def count_barcodes(barcodesfiles):
        count = 0
        with open(barcodesfiles, 'r') as f:
            for _ in f:
                count += 1
        f.close()

        return count

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

    def visualise_training_procedure(self, training_df, test_df):
        train_pct_change_df = training_df.pct_change() * 100
        train_pct_change_df[train_pct_change_df < 0] = np.nan
        validation_df = train_pct_change_df.describe().transpose()
        validation_df.sort_values(by=["mean", "std"], inplace=True)

        samples = []
        subtitles = []
        for index, row in validation_df.iterrows():
            samples.append(index)
            subtitles.append(" [{:.2f}% ±{:.2f}]".format(row["mean"], row["std"]))

        plot_df = self.build_plot_df(training_df, test_df)
        self.plot_per_sample(
            df=plot_df,
            samples=samples,
            plottype="line",
            subtitles=subtitles,
            x="epoch",
            y="loss",
            hue="group",
            xlabel="Epoch",
            ylabel="ELBO",
            title="Progress of the training procedure",
            filename="cellbender_remove_background_training_procedure"
        )

        return validation_df, samples

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

    def visualise_stats(self, stats_df, sample_order):
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
            "total_barcodes": "N total barcodes",
            "empty_droplets": "N empty droplets",
            "max_empty_droplet_count": "largest surely-empty droplet UMI",
            "optimal_regular_factor": "optimal posterior regularization factor",
            "time": "time (minutes)",
            "barcodes": "N barcodes"
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
                                        "N total barcodes",
                                        "N empty droplets",
                                        "largest surely-empty droplet UMI",
                                        "optimal posterior regularization factor",
                                        "time (minutes)",
                                        "N barcodes"]
                    if panel in df["variable"].unique()],
            x="sample",
            y="value",
            panel="variable",
            order=sample_order,
            filename="cellbender_remove_background_stats"
        )


    def barplot(self, df, panels, x="x", y="y", panel=None, palette=None, order=None,
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
                    order=order,
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

        for extension in self.extensions:
            outpath = os.path.join(self.plot_outdir, "{}.{}".format(filename, extension))
            fig.savefig(outpath)
            print("\tSaved figure: {}".format(os.path.basename(outpath)))
        plt.close()

    def visualise_cell_determination(self, cell_prob_df, sample_order):
        plot_df = cell_prob_df.copy()
        plot_df.reset_index(drop=False, inplace=True)
        plot_df = plot_df.melt(id_vars=["index"])
        plot_df.columns = ["index", "sample", "value"]

        self.plot_per_sample(
            df=plot_df,
            samples=sample_order,
            plottype="scatter",
            x="index",
            y="value",
            color="red",
            alpha=0.3,
            xlabel="Barcode index",
            ylabel="Cell probability",
            title="Determination of which barcodes contain cells",
            filename="cellbender_remove_background_latent_cell_probability"
        )

    def visualise_latent_gene_encoding(self, latent_gene_encoding_df, sample_order):
        self.plot_per_sample(
            df=latent_gene_encoding_df,
            samples=sample_order,
            plottype="scatter",
            x="PC0",
            y="PC1",
            alpha=0.3,
            xlabel="PC 0",
            ylabel="PC 1",
            title="PCA of latent encoding of cell gene expression",
            filename="cellbender_remove_background_latent_gene_encoding"
        )


    def plot_per_sample(self, df, samples, plottype, subtitles=None, x="x", y="y", alpha=1.0,
                            sample_column="sample", hue=None, color="black", palette=None,
                            xlabel="", ylabel="", title="", filename="plot"):
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
                data = df.loc[df[sample_column] == samples[i], :].dropna()
                if data.shape[0] == 0:
                    continue

                subtitle = ""
                if subtitles is not None:
                    subtitle = subtitles[i]

                sns.despine(fig=fig, ax=ax)

                if plottype == "line":
                    sns.lineplot(data=data,
                                 x=x,
                                 y=y,
                                 hue=hue,
                                 style=hue,
                                 color=color,
                                 alpha=alpha,
                                 markers=["o"] * len(data[hue].unique()) if hue in data else "o",
                                 markeredgewidth=0.0,
                                 palette=palette,
                                 ax=ax)
                elif plottype == "scatter":
                    sns.scatterplot(data=data,
                                    x=x,
                                    y=y,
                                    hue=hue,
                                    color=color,
                                    alpha=alpha,
                                    palette=palette,
                                    ax=ax)

                ax.set_xlabel(xlabel,
                              fontsize=10,
                              fontweight='bold')
                ax.set_ylabel(ylabel,
                              fontsize=10,
                              fontweight='bold')
                ax.set_title("{}{}".format(samples[i], subtitle),
                             fontsize=14,
                             fontweight='bold')


            else:
                ax.set_axis_off()

            col_index += 1
            if col_index > (ncols - 1):
                col_index = 0
                row_index += 1

        fig.suptitle(title,
                     fontsize=16,
                     color="#000000",
                     weight='bold')

        for extension in self.extensions:
            outpath = os.path.join(self.plot_outdir, "{}.{}".format(filename, extension))
            fig.savefig(outpath)
            print("\tSaved figure: {}".format(os.path.basename(outpath)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Working directory: {}".format(self.workdir))
        print("  > Extensions: {}".format(", ".join(self.extensions)))
        print("  > Plot output directory: {}".format(self.plot_outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()