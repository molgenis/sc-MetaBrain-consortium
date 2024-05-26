#!/usr/bin/env python3

"""
File:         create_wg1_metadata.py
Created:      2024/04/03
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
import glob
import os

# Third party imports.
import numpy as np
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Create WG1 Metadata"
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
./create_wg1_metadata.py -h
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.wg1_workdir = getattr(arguments, 'wg1_workdir')
        self.metadata_workdir = getattr(arguments, 'metadata_workdir')
        self.datasets = getattr(arguments, 'datasets')
        self.rochedir = getattr(arguments, 'rochedir')
        self.outdir = getattr(arguments, 'outdir')
        self.extensions = getattr(arguments, 'extension')

        # Creating output directory
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

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
        parser.add_argument("--wg1_workdir",
                            type=str,
                            required=True,
                            help="The WG1 working directory.")
        parser.add_argument("--metadata_workdir",
                            type=str,
                            required=True,
                            help="The metadata directory.")
        parser.add_argument("--datasets",
                            nargs="*",
                            type=str,
                            required=True,
                            default=None,
                            help="The WG1 datasets to load.")
        parser.add_argument("--rochedir",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--outdir",
                            type=str,
                            required=True,
                            help="The output directory.")
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

        print("Loading sample mapping files")
        roche_df = self.load_file(os.path.join(self.rochedir, "samples.individual_ids.eqtl.txt"))
        roche_columbia_df = self.load_file(os.path.join(self.rochedir, "key_columbia_rosmap.txt"))
        roche_columbia_ind_mapping = dict(zip(roche_columbia_df["individual_id_anon"], roche_columbia_df["individual_id"]))
        roche_df["individual_id"] = roche_df["individual_id"].replace(roche_columbia_ind_mapping)
        print(roche_df)
        print(len(roche_df["individual_id"].unique()))

        scmetabrain_df = self.load_sample_mapping(
            workdir=self.metadata_workdir,
            datasets=self.datasets
        )
        print(scmetabrain_df)

        sample_mapping_df = roche_df.merge(scmetabrain_df, how="left")
        sample_mapping_df["sc_metabrain_sample_id"] = sample_mapping_df["sc_metabrain_sample_id"].astype(str)
        print(sample_mapping_df.loc[sample_mapping_df["dataset"] == "Mathys2019", :])
        print(len(sample_mapping_df["sc_metabrain_sample_id"].unique()))
        print(len(sample_mapping_df["individual_id"].unique()))
        self.save_file(sample_mapping_df[["sc_metabrain_sample_id", "individual_id", "dataset"]], outpath=os.path.join(self.outdir,  "GTE.txt"))

        sample_mapping = dict(zip(sample_mapping_df["sc_metabrain_sample_id"], sample_mapping_df["individual_id"]))

        print("Loading metadata")
        metadata = self.load_metadata(
            workdir=self.wg1_workdir,
            datasets=self.datasets,
            sample_mapping=sample_mapping
        )
        print(metadata)

        self.save_file(metadata, outpath=os.path.join(self.outdir,  "Final_Assignments_demultiplexing_doublets.tsv.gz"))

    def load_sample_mapping(self, workdir, datasets):
        df_list = []
        for dataset in datasets:
            print("  loading '{}' results".format(dataset))
            full_link_path = os.path.join(workdir, dataset, dataset + "_full_link_table.tsv")
            if not os.path.exists(full_link_path):
                print("Error, could not find *_full_link_table.tsv input file.")
                continue
            df = self.load_file(full_link_path)
            if dataset == "Mathys2019":
                df["sample_id"] = [specimen_id.split("_")[0] for specimen_id in df["specimenID"]]
                df = df[["sample_id", "projid"]]
            elif dataset == "Zhou2020":
                df = df[["specimenID", "individualID"]]
            elif (dataset == "RocheAD2022") or (dataset == "RocheMS2022"):
                df["sample_id"] = [id_value.split("_")[0] for id_value in df["id"]]
                df = df[["sample_id", "sampleID"]]
            elif dataset == "RocheColumbia2022":
                df = df[["sample_id_anon", "SampleID"]]
            else:
                print("Not implemented")
                exit()

            df.columns = ["sample_id", "sc_metabrain_sample_id"]
            df["dataset"] = dataset
            df_list.append(df)

        return pd.concat(df_list, axis=0)

    def load_metadata(self, workdir, datasets, sample_mapping = None):
        df_list = []
        for dataset in datasets:
            dataset_list = []

            print("  loading '{}' results".format(dataset))
            subdirs = glob.glob(os.path.join(workdir, "*"))
            subdirs = [os.path.basename(subdir) for subdir in subdirs if dataset in os.path.basename(subdir)]
            subdirs.sort()
            dataset_folder = None
            if len(subdirs) == 1:
                dataset_folder = subdirs[0]
            elif len(subdirs) > 1:
                dataset_folder = subdirs[-1]
            else:
                print("Error, could not find dataset name")
                exit()
            print(dataset_folder)

            assignments = set()
            for path in glob.glob(os.path.join(workdir, dataset_folder, "Step2-DemultiplexingAndDoubletRemoval", "*")):
                pool = str(os.path.basename(path))
                if pool in ["CombinedResults", "manual_selection", "QC_figures", "log", "slurm_log", "genotypes"] or pool.endswith(".sh") or pool.endswith(".yaml") or pool.endswith(".yaml~"):
                    continue

                doubletfinder_path = os.path.join(workdir, dataset_folder, "Step2-DemultiplexingAndDoubletRemoval", pool, "DoubletFinderRun1", "DoubletFinder_doublets_singlets.tsv.gz")
                if not os.path.exists(doubletfinder_path):
                    print("Error, could not find DoubletFinder input file: {}.".format(doubletfinder_path))
                    continue
                doubletfinder_df = self.load_file(doubletfinder_path)

                scdblfinder_path = os.path.join(workdir, dataset_folder, "Step2-DemultiplexingAndDoubletRemoval", pool, "scDblFinderRun1", "scDblFinder_doublets_singlets.tsv.gz")
                if not os.path.exists(scdblfinder_path):
                    print("Error, could not find scDblFinder input file: {}.".format(scdblfinder_path))
                    continue
                scdblfinder_df = self.load_file(scdblfinder_path)


                df = doubletfinder_df.merge(scdblfinder_df, on="Barcode")
                df["Pool"] = pool
                df["Barcode"] = [barcode.split("-")[0] + "_" + pool for barcode in df["Barcode"]]
                df["Assignment"] = np.nan
                if sample_mapping is not None and pool in sample_mapping:
                    df["Assignment"] = sample_mapping[pool]
                    assignments.add(sample_mapping[pool])
                df["DropletType"] = "doublet"
                df.loc[(df["DoubletFinder_DropletType"] == "singlet") & (df["scDblFinder_DropletType"] == "singlet"), "DropletType"] = "singlet"
                df["Dataset"] = dataset

                dataset_list.append(df[["Pool", "Barcode", "Assignment", "DropletType"]])
                df_list.append(df[["Pool", "Barcode", "Assignment", "DropletType", "Dataset"]])


            dataset_df = pd.concat(dataset_list, axis=0)
            print(dataset, len(assignments))
            self.save_file(dataset_df, outpath=os.path.join(self.outdir, dataset + "_Final_Assignments_demultiplexing_doublets.tsv.gz"))

        return pd.concat(df_list, axis=0)
    @staticmethod
    def load_file(inpath, header=0, index_col=None, sep="\t", skiprows=None,
                  nrows=None, dtype=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         skiprows=skiprows, nrows=nrows, dtype=dtype)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def save_file(df, outpath, header=True, index=False, sep="\t"):
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
        print("  > WG1 working directory: {}".format(self.wg1_workdir))
        print("  > Metadata working directory: {}".format(self.metadata_workdir))
        print("  > Datasets: {}".format(", ".join(self.datasets)))
        print("  > Roche directory: {}".format(self.rochedir))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()