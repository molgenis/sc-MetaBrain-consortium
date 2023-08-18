#!/usr/bin/env python3

"""
File:         download_from_synapse.py
Created:      2023/01/19
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
import argparse
import os
import getpass

# Third party imports.
import synapseclient
import synapseutils

# Local application imports.

# Metadata
__program__ = "Download From Synapse"
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
  
./download_from_synapse.py 

### Mathys 2019 ###
./download_from_synapse.py \
    --synapse_ids syn18681734 syn18638475 syn18642926 \
    --folders processed fastq metadata \
    --outdir /groups/umcg-biogen/tmp01/input/rawdata/single-cell/Mathys2019
    
### Cain 2020 ###
./download_from_synapse.py \
    --synapse_ids syn21589957 syn17055069 \
    --folders processed fastq \
    --outdir /groups/umcg-biogen/tmp01/input/rawdata/single-cell/Cain2020

### Zhou 2020 ### 
./download_from_synapse.py \
    --synapse_ids syn21682218 syn21126450 \
    --folders metadata fastq \
    --outdir /groups/umcg-biogen/tmp02/input/rawdata/single-cell/Zhou2020
    
### Fujita 2022 ### 
./download_from_synapse.py \
    --synapse_ids syn31512863 \
    --folders Fastq \
    --outdir /groups/umcg-biogen/tmp01/input/rawdata/single-cell/Fujita2022
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.synapse_ids = getattr(arguments, 'synapse_ids')
        self.folders = getattr(arguments, 'folders')
        self.out_dir = getattr(arguments, 'outdir')

        if len(self.synapse_ids) != len(self.folders):
            print("Error, -synapse_ids and --folders must have equal lengths.")
            exit()

        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

        os.chdir(self.out_dir)

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # general arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit")
        parser.add_argument("--synapse_ids",
                            nargs="*",
                            type=str,
                            required=True,
                            help="The synapse ID to download.")
        parser.add_argument("--folders",
                            nargs="*",
                            type=str,
                            required=True,
                            help="The subfolders to download.")
        parser.add_argument("--outdir",
                            type=str,
                            required=True,
                            help="The path to the output directory.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        # Log in on Synapse.
        user = input("Synapse username:")
        password = getpass.getpass('Synapse password:')

        syn = synapseclient.Synapse()
        syn.login(user, password)

        # Download.
        for synapse_id, folder in zip(self.synapse_ids, self.folders):
            _ = synapseutils.syncFromSynapse(syn,
                                             entity=synapse_id,
                                             path=folder + "/")

    def print_arguments(self):
        print("General arguments")
        for synapse_id, folder in zip(self.synapse_ids, self.folders):
            print("  > entity: {}\tpath: {}/".format(synapse_id, folder))
        print("  > Output directory: {}".format(self.out_dir))
        print("")

if __name__ == '__main__':
    m = main()
    m.start()