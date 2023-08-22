#!/usr/bin/env python3

"""
File:         replace_PSAM_file.py
Created:      2022/11/15
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
import os

# Third party imports.
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Replace PSAM file"
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

### AMP-AD ###
./replace_PSAM_file.py \
    --original /groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/2022-11-03-FilteredGenotypes/5-plink2_makepgen/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_CHR.recalibrated_variants_norm_vcffilter_concat_annotate.psam \
    --new /groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/AMP_AD_annotdata.psam
    
./replace_PSAM_file.py \
    --original /groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/2022-11-03-FilteredGenotypes/5-plink2_makepgen_after_fill/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_CHR.recalibrated_variants_norm_vcffilter_concat_annotate.psam \
    --new /groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/AMP_AD_annotdata.psam
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.original_path = getattr(arguments, 'original')
        self.new_path = getattr(arguments, 'new')

        for path in [self.original_path, self.new_path]:
            if not path.endswith(".psam"):
                print("Error, input files need to be '.psam'.")
                exit()

        self.required_columns = ['#FID', 'IID', 'PAT', 'MAT', 'SEX', 'Provided_Ancestry', 'genotyping_platform', 'array_available', 'wgs_available', 'wes_available', 'age', 'age_range', 'Study', 'smoking_status', 'hormonal_contraception_use_currently', 'menopause', 'pregnancy_status']

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
        parser.add_argument("--original",
                            type=str,
                            required=True,
                            help="The path to the original PSAM file"
                                 "outputted by PLINK.")
        parser.add_argument("--new",
                            type=str,
                            required=True,
                            help="The path to the new PSAM file.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data.")
        original_df = self.load_file(self.original_path)
        print(original_df)

        new_df = self.load_file(self.new_path)
        new_df.fillna("NA", inplace=True)
        print(new_df)

        print("Validating data.")
        if original_df["#IID"].tolist() != new_df["IID"].tolist():
            print("Error, IID columns do not match.")
            exit()
        if new_df.columns.tolist() != self.required_columns:
            print("Error, new PSAM file columns do not match.")
            exit()

        print("Moving old PSAM file.")
        original_dir, original_filename = os.path.split(self.original_path)
        movedir = os.path.join(original_dir, 'original_psam')
        if not os.path.exists(movedir):
            os.makedirs(movedir)
        replace_outpath = os.path.join(movedir, original_filename)

        self.save_file(df=original_df, outpath=replace_outpath)

        print("Overwriting original file.")
        self.save_file(df=new_df, outpath=self.original_path)

    @staticmethod
    def load_file(inpath, header=0, index_col=None, sep="\t"):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col)
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
        print("  > Original: {}".format(self.original_path))
        print("  > New: {}".format(self.new_path))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()