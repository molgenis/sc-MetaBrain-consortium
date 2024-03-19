#!/usr/bin/env python3

"""
File:         compare_vcf.py
Created:      2024/01/25
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
import argparse
import gzip
import math
import random
import os

# Third party imports.
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Compare VCF"
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
./compare_vcf.py -h
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.vcf1 = getattr(arguments, 'vcf1')
        self.vcf2 = getattr(arguments, 'vcf2')
        self.n = getattr(arguments, 'n')
        self.title = getattr(arguments, 'title')
        self.extensions = getattr(arguments, 'extension')

        # Creating output directory
        self.outdir = os.path.join(str(os.path.dirname(os.path.abspath(__file__))), 'plots')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

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
        parser.add_argument("--vcf1",
                            type=str,
                            required=True,
                            help="")
        parser.add_argument("--vcf2",
                            type=str,
                            required=True,
                            help="")
        parser.add_argument("--n",
                            type=int,
                            required=False,
                            default=None,
                            help="")
        parser.add_argument("--title",
                            type=str,
                            required=False,
                            default="",
                            help="")
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

        print("Loading VCF 1 indices")
        vcf1_indices = self.load_indices(self.vcf1)
        print("\tN = {:,}".format(len(vcf1_indices)))

        # vcf1_indices_dict = {}
        # for vcf1_index in vcf1_indices:
        #     chr, pos, ref, alt = vcf1_index.split(":")
        #     if pos in vcf1_indices_dict:
        #         vcf1_indices_dict[pos].append(vcf1_index)
        #     else:
        #         vcf1_indices_dict[pos] = [vcf1_index]

        print("\nLoading VCF 2 indices")
        vcf2_indices = self.load_indices(self.vcf2)
        print("\tN = {:,}".format(len(vcf2_indices)))

        print("\nFinding overlap")
        overlap = vcf1_indices & vcf2_indices
        print("\tN = {:,}".format(len(overlap)))
        if len(overlap) == 0:
            return
        if self.n is not None and len(overlap) > self.n:
            overlap = random.sample(overlap, self.n)
            print("\tLimiting to N = {:,}".format(len(overlap)))
        del vcf1_indices, vcf2_indices

        print("\nLoading VCF 1")
        df1 = self.load_data(self.vcf1, query_indices=overlap)
        print(df1)
        df1.to_csv(os.path.join(self.outdir, "2024-01-26-df1.txt.gz"), sep="\t", header=True, index=True)
        # df1 = pd.read_csv(os.path.join(self.outdir, "df1.txt.gz"), sep="\t", header=0, index_col=0)

        print("\nLoading VCF 2")
        df2 = self.load_data(self.vcf2, query_indices=overlap)
        print(df2)
        df2.to_csv(os.path.join(self.outdir, "2024-01-26-df2.txt.gz"), sep="\t", header=True, index=True)
        # df2 = pd.read_csv(os.path.join(self.outdir, "df2.txt.gz"), sep="\t", header=0, index_col=0)

        print("\nComparing:")
        print("  Fully identical:", df1.equals(df2))
        print("")

        if df1.equals(df2):
            return

        if df1.index.tolist() != df2.index.tolist():
            print("  Indices differ, filtering on overlapping")
            overlap = [index for index in df1.index if index in df2.index]
            df1 = df1.loc[overlap, :]
            df2 = df2.loc[overlap, :]
            print("  Data frames identical after row filter:", df1.equals(df2))
        else:
            print("  Indices are identical")
        print("")

        if df1.columns.tolist() != df2.columns.tolist():
            print("  Columns differ, filtering on overlapping")
            overlap = [column for column in df1.columns if column in df2.columns]
            df1 = df1.loc[:, overlap]
            df2 = df2.loc[:, overlap]
            print("  Data frames identical after col filter:", df1.equals(df2))
        else:
            print("  Columns are identical")
        print("")

        print(df1)
        print(df2)

        if not df1.equals(df2):
            plot_data = {}
            comparison_data = []
            comparison_data2 = []
            for col in df1.columns:
                if col in df2.columns:
                    df = df1[[col]].merge(df2[[col]], left_index=True, right_index=True)
                    df.columns = ["x", "y"]
                    df.replace(-1, np.nan, inplace=True)
                    df.dropna(inplace=True)
                    df[col] = df["x"] == df["y"]
                    n_same = df[col].sum()
                    n_tested = df.shape[0]
                    pcnt_same = (100 / n_tested) * n_same
                    print("  {:,} / {:,} [{:.2f}%] identical".format(n_same, n_tested, pcnt_same))

                    if "GT"in col:
                        comparison_data.append([col, n_same, n_tested, pcnt_same])
                        comparison_data2.append(df[col])
                        # if len(comparison_data) > 5:
                        #     break

                    # if df1[col].equals(df2[col]):
                    #     print("  ", col, "identical")
                    # elif df1[col].values.tolist() == df2[col].values.tolist():
                    #     print("  ", col, "identical values, different indices")
                    # elif pd.api.types.is_numeric_dtype(df1[col]) and pd.api.types.is_numeric_dtype(df2[col]):
                    #     corr_df = df1[[col]].merge(df2[[col]], left_index=True, right_index=True)
                    #     corr_df.columns = ["x", "y"]
                    #     corr_df.dropna(inplace=True)
                    #     pearson_coef, pearson_p = stats.pearsonr(corr_df["x"], corr_df["y"])
                    #     print("  ", col, "different", pearson_coef, pearson_p)
                    #
                    #     plot_data[col] = corr_df
                    # else:
                    #     print("  ", col, "different values, different indices")

            comparison_df = pd.DataFrame(comparison_data, columns=["Column", "Same", "Tested", "Percentage"])
            print(comparison_df)
            print(comparison_df.mean(axis=0))
            print(comparison_df["Percentage"].min())
            print(comparison_df["Percentage"].max())
            print(comparison_df["Percentage"].mean())

            comparison2_df = pd.concat(comparison_data2, axis=1)
            comparison2_df.fillna(True, inplace=True)
            print(comparison2_df)
            print(comparison2_df.sum(axis=0).min())
            print(comparison2_df.sum(axis=0).max())
            print(comparison2_df.sum(axis=0).mean())
            print(comparison2_df.sum(axis=0).value_counts())
            print("")
            print(comparison2_df.sum(axis=1).min())
            print(comparison2_df.sum(axis=1).max())
            print(comparison2_df.sum(axis=1).mean())
            print(comparison2_df.sum(axis=1).value_counts())
            print(comparison2_df.sum(axis=1).value_counts()[0])
            print(comparison2_df.loc[comparison2_df.sum(axis=1) == 0, :])

            if len(plot_data.keys()) > 0:
                self.visualise_data(plot_data=plot_data, title=self.title, filename=self.title.lower())

        print("\n----------------------------------------\n")

    @staticmethod
    def load_indices(inpath, nrows=None, query_indices=None):
        chrom_index = None
        pos_index = None
        ref_index = None
        alt_index = None
        indices = set()
        with gzip.open(inpath, 'rt') as f:
            for i, line in enumerate(f):
                if (i == 0) or (i % 100000 == 0):
                    print("  parsed {:,} lines".format(i))
                if line.startswith("##"):
                    continue

                if nrows is not None and len(indices) >= nrows:
                    break

                values = line.strip("\n").split("\t")
                if line.startswith("#CHROM"):
                    chrom_index = values.index("#CHROM")
                    pos_index = values.index("POS")
                    ref_index = values.index("REF")
                    alt_index = values.index("ALT")
                    continue

                alleles = [values[ref_index], values[alt_index]]
                alleles.sort()
                index = "{}:{}:{}:{}".format(values[chrom_index], values[pos_index], alleles[0], alleles[1])
                indices.add(index)

                # if query_indices is not None and values[pos_index] in query_indices and index not in query_indices[values[pos_index]]:
                #     print(index)
                #     print(query_indices[values[pos_index]])
                #     print("")
        f.close()

        return indices

    @staticmethod
    def load_data(inpath, nrows=None, query_indices=None):
        data = {}
        n_skipped = 0
        n_saved = 0

        chrom_index = None
        pos_index = None
        id_index = None
        ref_index = None
        alt_index = None
        qual_index = None
        filter_index = None
        info_index = None
        format_index = None
        samples_index = None
        samples = []
        with gzip.open(inpath, 'rt') as f:
            for i, line in enumerate(f):
                if (i == 0) or (i % 100000 == 0):
                    print("  parsed {:,} lines, skipped {:,} variants, saved {:,} variants".format(i, n_skipped, n_saved))

                # if nrows is None and n_saved > 0:
                #     break

                if nrows is not None and n_saved >= nrows:
                    break
                if line.startswith("##"):
                    continue

                values = line.strip("\n").split("\t")
                if line.startswith("#CHROM"):
                    chrom_index = values.index("#CHROM")
                    pos_index = values.index("POS")
                    id_index = values.index("ID")
                    ref_index = values.index("REF")
                    alt_index = values.index("ALT")
                    qual_index = values.index("QUAL")
                    filter_index = values.index("FILTER")
                    info_index = values.index("INFO")
                    format_index = values.index("FORMAT")
                    for samples_index, column in enumerate(values):
                        if column not in ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]:
                            break
                    samples = values[samples_index:]
                    continue

                alleles = [values[ref_index], values[alt_index]]
                alleles.sort()
                index = "{}:{}:{}:{}".format(values[chrom_index], values[pos_index], alleles[0], alleles[1])
                if query_indices is not None and index not in query_indices:
                    n_skipped += 1
                    continue

                data[index] = {
                    "#CHROM": values[chrom_index],
                    "POS": values[pos_index],
                    "ID": values[id_index],
                    "REF": values[ref_index],
                    "ALT": values[alt_index],
                    "QUAL": values[qual_index],
                    "FILTER": values[filter_index]
                }

                for field in values[info_index].split(";"):
                    split_field = field.split("=")
                    if len(split_field) == 1:
                        key = split_field[0]
                        value = 1
                    else:
                        key = split_field[0]
                        value = split_field[1]
                    data[index][key] = value

                gt_index = values[format_index].split(":").index("GT")
                for sample, value in zip(samples, values[samples_index:]):
                    gt_value = value.split(":")[gt_index]
                    if gt_value == "./.":
                        data[index][sample + "-GT"] = -1.0
                    elif "/" in gt_value:
                        data[index][sample + "-GT"] = sum([float(dosage) for dosage in gt_value.split("/")])
                    elif "|" in gt_value:
                        data[index][sample + "-GT"] = sum([float(dosage) for dosage in gt_value.split("|")])
                    else:
                        print(gt_value)
                        exit()

                n_saved += 1

        f.close()

        print("  parsed {:,} lines, skipped {:,} variants, saved {:,} variants".format(i, n_skipped, n_saved))

        return pd.DataFrame(data).T

    @staticmethod
    def visualise_data(plot_data, title="", filename="plot"):
        keys = list(plot_data.keys())
        keys.sort()

        nplots = len(keys)
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
                data = plot_data[keys[i]]
                if data.shape[0] == 0:
                    continue

                sns.despine(fig=fig, ax=ax)

                sns.regplot(x="x", y="y", data=data, ci=95,
                            scatter_kws={'facecolors': "#808080",
                                         'edgecolors': "#808080"},
                            line_kws={"color": "#b22222"},
                            ax=ax
                            )

                ax.axline((0, 0), slope=1, ls='--', color="#000000", alpha=0.15, zorder=-1)

                # Set annotation.
                pearson_coef, _ = stats.pearsonr(data["y"], data["x"])
                ax.annotate(
                    'total N = {:,}'.format(data.shape[0]),
                    xy=(0.03, 0.94),
                    xycoords=ax.transAxes,
                    color="#000000",
                    fontsize=14,
                    fontweight='bold')
                ax.annotate(
                    'total r = {:.4f}'.format(pearson_coef),
                    xy=(0.03, 0.88),
                    xycoords=ax.transAxes,
                    color="#000000",
                    fontsize=14,
                    fontweight='bold')

                ax.set_xlabel("Drew",
                              fontsize=10,
                              fontweight='bold')
                ax.set_ylabel("Martijn",
                              fontsize=10,
                              fontweight='bold')
                ax.set_title(keys[i],
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

        if not os.path.exists("plot"):
            os.makedirs("plot")

        plt.tight_layout()
        outpath = os.path.join("plot/{}.png".format(filename))
        fig.savefig(outpath)
        print("\tSaved figure: {}".format(os.path.basename(outpath)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > VCF1: {}".format(self.vcf1))
        print("  > VCF2: {}".format(self.vcf2))
        print("  > N {:,}".format(self.n))
        print("  > Title: {}".format(self.title))
        print("  > Extensions: {}".format(", ".join(self.extensions)))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()