#!/usr/bin/env python3

"""
File:         combine_wg3_eqtl_files.py
Created:      2023/04/24
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
import h5py

# Local application imports.

"""
Syntax:
./combine_wg3_eqtl_files.py -h
"""

# Metadata
__program__ = "Combine Workgroup 3 eQTL files"
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


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.work_dir = getattr(arguments, 'work_dir')
        self.verbose = getattr(arguments, 'verbose')

        self.outdir = os.path.join(self.work_dir, "merged")
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
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
        parser.add_argument("--verbose",
                            action='store_true',
                            help="Print all info.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("  Loading Workgroup 3 Data")
        for _, cell_type in self.cell_type_dict.items():
            print("  Processing {}".format(cell_type))

            results_file_prefix = "qtl_results_"

            h5_files = glob.glob(os.path.join(self.wg3_folder, "output", self.annotation_level, cell_type, "qtl", "{}*.h5".format(results_file_prefix)))
            h5_files = [(int(os.path.basename(h5_file).split("_")[2]), h5_file) for h5_file in h5_files]
            h5_files.sort(key=lambda x: x[0])
            if len(h5_files) == 0:
                print("  Warning, no h5 files found for this cell type")
                return None

            dfs = []
            prev_chr = None
            for chr, h5_filepath in h5_files:
                if prev_chr is None or chr != prev_chr:
                    print("\tLoading chromosome {}".format(os.path.basename(h5_filepath).split("_")[2]))
                h5_df = self.load_h5_file(filepath=h5_filepath,
                                          cell_type=cell_type,
                                          results_file_prefix=results_file_prefix)
                if h5_df is None:
                    continue

                dfs.append(h5_df)
                prev_chr = chr
                del h5_df

            if len(dfs) == 0:
                print("  Warning, no data loaded for this cell type")
                return None

            df = pd.concat(dfs, axis=0)

            if df.shape[0] == 0:
                print("  Warning, no data loaded for this cell type")
                return None

            print("  Saving file")
            outpath = os.path.join(self.outdir, cell_type)
            self.save_file(df=df, outpath=outpath + ".txt.gz")
            self.save_file(df=df, outpath=outpath + ".pkl")
            # self.save_file(df=df, outpath=outpath + ".xlsx")

    def load_h5_file(self, filepath, cell_type, results_file_prefix):
        analysis_subset = os.path.basename(filepath).replace(results_file_prefix, "").replace(".h5", "")
        snp_metadata_file = os.path.join(self.wg3_folder, "output", self.annotation_level, cell_type, "qtl", "snp_metadata_{}.txt".format(analysis_subset))
        feature_metadata_file = os.path.join(self.wg3_folder, "output", self.annotation_level, cell_type, "qtl", "feature_metadata_{}.txt".format(analysis_subset))

        if os.path.exists(feature_metadata_file):
            ffea_df = self.load_file(feature_metadata_file, header=0, index_col=None)
        elif os.path.exists(feature_metadata_file + ".gz"):
            ffea_df = self.load_file(feature_metadata_file + ".gz", header=0, index_col=None)
        else:
            print("  Warning, skipping: '{}' missing feature metadata file.".format(analysis_subset))
            return None

        if os.path.exists(snp_metadata_file):
            fsnp_df = self.load_file(snp_metadata_file, header=0, index_col=None)
        elif os.path.exists(snp_metadata_file + ".gz"):
            fsnp_df = self.load_file(snp_metadata_file + ".gz", header=0, index_col=None)
        else:
            print("  Warning, skipping: '{}' missing SNP metadata file.".format(analysis_subset))
            return None

        ffea_df = ffea_df.rename(index=str,
                                 columns={"chromosome": "feature_chromosome",
                                          "start": "feature_start",
                                          "end": "feature_end"})
        fsnp_df = fsnp_df.rename(index=str,
                                 columns={"chromosome": "snp_chromosome",
                                          "position": "snp_position"})

        ########################################################################

        snp_id_df = fsnp_df["snp_id"].str.split(":", n=None, expand=True)
        snp_id_df.columns = ["snp_chromosome", "snp_position", "alleleA", "alleleB"]
        snp_id_df.drop(["snp_chromosome", "snp_position"], axis=1, inplace=True)
        fsnp_df = pd.concat([fsnp_df, snp_id_df], axis=1)
        fsnp_df['allele1'] = np.minimum(fsnp_df['alleleA'], fsnp_df['alleleB'])
        fsnp_df['allele2'] = np.maximum(fsnp_df['alleleA'], fsnp_df['alleleB'])
        fsnp_df['alleles'] = fsnp_df['allele1'] + fsnp_df['allele2']
        fsnp_df.drop(['alleleA', 'alleleB', 'allele1', 'allele2'], axis=1, inplace=True)

        ########################################################################

        frez = h5py.File(filepath, 'r')
        frez_keys = [k.replace('_i_', '') for k in list(frez.keys())]
        if len(frez_keys) == 0:
            return None

        dfs = []
        for frez_key in frez_keys:
            frez_df = pd.DataFrame(np.array(frez[frez_key]))
            frez_df['feature_id'] = frez_key
            dfs.append(frez_df)
        df = pd.concat(dfs, axis=0)
        df['snp_id'] = df['snp_id'].astype(str)

        if self.verbose:
            print("\tLoaded h5 file: {} with shape: {}".format(os.path.basename(filepath), df.shape))

        df = pd.merge(df, ffea_df, on='feature_id', how='left')

        ########################################################################

        if(len(glob.glob(os.path.join(self.wg3_folder, "output", self.annotation_level, cell_type, "qtl", "snp_qc_metrics_naContaining_feature_*.txt"))) > 0):
            print("  Error, code not implemented yet")
            exit()
            # tmp_df = pd.DataFrame(columns=df.columns)
            # for key in frez_keys:
            #     qc_metrics_inpath = os.path.join(self.wg3_folder, "output", self.annotation_level, cell_type, "qtl", "snp_qc_metrics_naContaining_feature_{}.txt".format(key))
            #     if os.path.isfile(qc_metrics_inpath):
            #         fsnp_rel = self.load_file(qc_metrics_inpath, header=0, index_col=None)
            #         tmp_t = df.loc[df["feature_id"] == key]
            #         fsnp_t = fsnp_df.loc[:,["snp_id", "snp_chromosome", "snp_position", "assessed_allele"]]
            #         fsnp_t = pd.merge(fsnp_t, fsnp_rel, on='snp_id', how='right')
            #         tmp_t = pd.merge(tmp_t, fsnp_t, on='snp_id', how='left')
            #         tmp_df = tmp_df.append(tmp_t, sort=False)
            #     else:
            #         tmp_t = df.loc[df["feature_id"] == key]
            #         tmp_t = pd.merge(tmp_t, fsnp_df, on='snp_id', how='left')
            #         tmp_df = tmp_df.append(tmp_t,sort=False)
            #     data[key]=np.zeros(len(np.unique(list(frez_keys))),dtype='object')+np.nan
            # df = tmp_df
            # del tmp_df
        else:
            df = pd.merge(df, fsnp_df, on='snp_id', how='left')

        df['empirical_feature_p_value'] = df['empirical_feature_p_value'].astype(float)
        df['p_value'] = df['p_value'].astype(float)

        del ffea_df, fsnp_df

        return df

    def save_file(self, df, outpath, header=True, index=True, sep="\t", na_rep="NA",
                  sheet_name="Sheet1"):
        if outpath.endswith('pkl'):
            df.to_pickle(outpath)
        elif outpath.endswith('xlsx'):
            df.to_excel(outpath,
                        sheet_name=sheet_name,
                        na_rep=na_rep,
                        header=header,
                        index=index)
        else:
            compression = 'infer'
            if outpath.endswith('.gz'):
                compression = 'gzip'

            df.to_csv(outpath, sep=sep, index=index, header=header,
                      compression=compression)

        if self.verbose:
            print("\tSaved dataframe: {} "
                  "with shape: {}".format(os.path.basename(outpath),
                                          df.shape))

    def print_arguments(self):
        print("Arguments:")
        print("  > Working directory:    {}".format(self.work_dir))
        print("  > Verbose:              {}".format(self.verbose))
        print("  > Output directory:     {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
