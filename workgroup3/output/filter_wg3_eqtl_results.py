#!/usr/bin/env python3

"""
File:         filter_wg3_eqtl_results.py
Created:      2023/04/27
Last Changed:
Author:       M.Vochteloo

Copyright (C) 2020 M.Vochteloo
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
./filter_wg3_eqtl_results.py \
    --input_dir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-11-WorkGroup3eQTLAndDEA/2023-04-12-Mathys2019/output/L1/AST \
    --snp_gene /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-25-mbQTLBryoisRepl/combos/ASTcombos.txt
    
./filter_wg3_eqtl_results.py \
    --input_dir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-11-WorkGroup3eQTLAndDEA/2023-04-12-Mathys2019/output/L1/EX \
    --snp_gene /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-25-mbQTLBryoisRepl/combos/EXcombos.txt

./filter_wg3_eqtl_results.py \
    --input_dir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-11-WorkGroup3eQTLAndDEA/2023-04-12-Mathys2019/output/L1/IN \
    --snp_gene /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-25-mbQTLBryoisRepl/combos/INcombos.txt

./filter_wg3_eqtl_results.py \
    --input_dir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-11-WorkGroup3eQTLAndDEA/2023-04-12-Mathys2019/output/L1/MIC \
    --snp_gene /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-25-mbQTLBryoisRepl/combos/MICcombos.txt
    
./filter_wg3_eqtl_results.py \
    --input_dir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-11-WorkGroup3eQTLAndDEA/2023-04-12-Mathys2019/output/L1/OLI \
    --snp_gene /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-25-mbQTLBryoisRepl/combos/OLIcombos.txt
    
./filter_wg3_eqtl_results.py \
    --input_dir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-11-WorkGroup3eQTLAndDEA/2023-04-12-Mathys2019/output/L1/OPC \
    --snp_gene /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-25-mbQTLBryoisRepl/combos/OPCcombos.txt
    
./filter_wg3_eqtl_results.py \
    --input_dir /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-11-WorkGroup3eQTLAndDEA/2023-04-12-Mathys2019/output/L1/PER \
    --snp_gene /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-04-25-mbQTLBryoisRepl/combos/PERcombos.txt
"""

# Metadata
__program__ = "Filter Workgroup 3 eQTL Results"
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
        self.input_dir = getattr(arguments, 'input_dir')
        self.snp_gene = getattr(arguments, 'snp_gene')
        self.verbose = getattr(arguments, 'verbose')

        self.results_file_prefix = "qtl_results_"

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
        parser.add_argument("--input_dir",
                            type=str,
                            required=True,
                            help="")
        parser.add_argument("--snp_gene",
                            type=str,
                            required=True,
                            help=".")
        parser.add_argument("--verbose",
                            action='store_true',
                            help="Print all info.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading filter list")
        snp_gene_df = self.load_file(self.snp_gene, header=None, index_col=None)
        snp_gene_df.index = snp_gene_df.iloc[:, 0] + "_" + snp_gene_df.iloc[:, 1]
        snps = set(snp_gene_df.iloc[:, 0].values)
        feature_ids = set(snp_gene_df.iloc[:, 1].values)
        feature_to_snp_dict = dict(zip(snp_gene_df.iloc[:, 1], snp_gene_df.iloc[:, 0]))
        print("\t  filtering for {:,.0f} eQTLs".format(snp_gene_df.shape[0]))
        print("")

        print("Loading workgroup 3 eQTL results data")
        h5_files = glob.glob(os.path.join(self.input_dir, "qtl", "{}*.h5".format(self.results_file_prefix)))
        h5_files = [(int(os.path.basename(h5_file).split("_")[2]), h5_file) for h5_file in h5_files]
        h5_files.sort(key=lambda x: x[0])
        if len(h5_files) == 0:
            print("  Warning, no h5 files found for this cell type")
            return None

        prev_chr = None
        chr_rows = 0
        dfs = []
        for chr, h5_filepath in h5_files:
            if prev_chr is None:
                print("\tLoading chromosome {}".format(os.path.basename(h5_filepath).split("_")[2]))
            elif chr != prev_chr:
                print("\t  {:,.0f} rows stored".format(chr_rows))
                print("\tLoading chromosome {}".format(os.path.basename(h5_filepath).split("_")[2]))
                chr_rows = 0
            prev_chr = chr

            h5_df = self.load_h5_file(filepath=h5_filepath,
                                      snps=snps,
                                      feature_ids=feature_ids,
                                      feature_to_snp_dict=feature_to_snp_dict)
            if h5_df is None:
                continue

            chr_rows += h5_df.shape[0]

            dfs.append(h5_df)
            del h5_df

        if len(dfs) == 0:
            print("Warning, no data loaded for this cell type")
            exit()

        df = pd.concat(dfs, axis=0)
        print("\t  found {:,.0f}/{:,.0f} eQTLs".format(df.shape[0], snp_gene_df.shape[0]))

        if df.shape[0] == 0:
            print("Warning, no data loaded for this cell type")
            exit()
        print("")

        print("Saving file")
        self.save_file(df=df, outpath=os.path.join(self.input_dir, "{}BryoiseQTLs.txt.gz".format(self.results_file_prefix)))
        print("")

    def load_h5_file(self, filepath, snps, feature_ids, feature_to_snp_dict):
        analysis_subset = os.path.basename(filepath).replace(self.results_file_prefix, "").replace(".h5", "")
        snp_metadata_file = os.path.join(self.input_dir, "qtl", "snp_metadata_{}.txt".format(analysis_subset))
        feature_metadata_file = os.path.join(self.input_dir, "qtl", "feature_metadata_{}.txt".format(analysis_subset))

        if os.path.exists(feature_metadata_file):
            ffea_df = self.load_file(feature_metadata_file, header=0, index_col=None)
        elif os.path.exists(feature_metadata_file + ".gz"):
            ffea_df = self.load_file(feature_metadata_file + ".gz", header=0, index_col=None)
        else:
            print("  Warning, skipping: '{}' missing feature metadata file.".format(analysis_subset))
            return None

        overlapping_features = set(ffea_df["feature_id"]).intersection(feature_ids)
        if len(overlapping_features) == 0:
            return None

        if os.path.exists(snp_metadata_file):
            fsnp_df = self.load_file(snp_metadata_file, header=0, index_col=None)
        elif os.path.exists(snp_metadata_file + ".gz"):
            fsnp_df = self.load_file(snp_metadata_file + ".gz", header=0, index_col=None)
        else:
            print("  Warning, skipping: '{}' missing SNP metadata file.".format(analysis_subset))
            return None

        overlapping_snps = set(fsnp_df["snp_id"]).intersection(snps)
        if len(overlapping_snps) == 0:
            return None

        ffea_df = ffea_df.rename(index=str,
                                 columns={"chromosome": "feature_chromosome",
                                          "start": "feature_start",
                                          "end": "feature_end"})
        fsnp_df = fsnp_df.rename(index=str,
                                 columns={"chromosome": "snp_chromosome",
                                          "position": "snp_position"})

        frez = h5py.File(filepath, 'r')
        frez_keys = [k.replace('_i_', '') for k in list(frez.keys())]

        dfs = []
        for frez_key in frez_keys:
            if frez_key not in overlapping_features or frez_key not in feature_to_snp_dict:
                continue
            frez_df = pd.DataFrame(np.array(frez[frez_key]))
            frez_df['snp_id'] = frez_df['snp_id'].astype(str)

            eqtl_snp = feature_to_snp_dict[frez_key]
            if eqtl_snp not in overlapping_snps:
                continue

            eqtl_df = frez_df.loc[frez_df["snp_id"] == eqtl_snp, :].copy()
            eqtl_df['feature_id'] = frez_key
            dfs.append(eqtl_df)
            del frez_df

        if len(dfs) == 0:
            return None

        df = pd.concat(dfs, axis=0)

        if self.verbose:
            print("\tLoaded h5 file: {} with shape: {}".format(os.path.basename(filepath), df.shape))

        df = pd.merge(df, ffea_df, on='feature_id', how='left')

        ########################################################################

        if(len(glob.glob(os.path.join(self.input_dir, "qtl", "snp_qc_metrics_naContaining_feature_*.txt"))) > 0):
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
            df = pd.merge(df, fsnp_df, on='snp_id', how='inner')

        df['empirical_feature_p_value'] = df['empirical_feature_p_value'].astype(float)
        df['p_value'] = df['p_value'].astype(float)

        del ffea_df, fsnp_df

        return df

    def load_file(self, inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)

        if self.verbose:
            print("\tLoaded dataframe: {} "
                  "with shape: {}".format(os.path.basename(inpath),
                                          df.shape))
        return df

    def save_file(self, df, outpath, header=True, index=False, sep="\t"):
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
        print("  > Input directory:    {}".format(self.input_dir))
        print("  > SNP-gene file:      {}".format(self.snp_gene))
        print("  > Verbose:            {}".format(self.verbose))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
