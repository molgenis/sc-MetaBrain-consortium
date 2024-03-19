#!/usr/bin/env python3

"""
File:         download_from_synapse.py
Created:      2023/01/19
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
./download_from_synapse.py -h
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.entities = getattr(arguments, 'entities')
        self.folders = getattr(arguments, 'folders')
        self.out_dir = getattr(arguments, 'outdir')

        if len(self.entities) != len(self.folders):
            print("Error, -input and --folders must have equal lengths.")
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
        parser.add_argument("--entities",
                            nargs="*",
                            type=str,
                            required=False,
                            help="The synapse ID to download or a file containing synapse IDs.")
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
        authToken = getpass.getpass('Synapse authentication token:')

        syn = synapseclient.Synapse()
        syn.login(authToken=authToken)

        # Download.
        for entity, folder in zip(self.entities, self.folders):
            if os.path.exists(entity):
                with open(entity, 'r') as f:
                    for line in f:
                        entity_code, entity_file = line.strip("\n").split("  ")
                        if os.path.exists(folder + "/" + entity_file):
                            print("Skipping '{}' - '{}'".format(entity_code, entity_file))
                            continue
                        self.download_file(syn=syn, entity=entity_code, path=folder + "/")
                f.close()
            else:
                self.download_file(syn=syn, entity=entity, path=folder + "/")

    @staticmethod
    def download_file(syn, entity, path):
        try:
            _ = synapseutils.syncFromSynapse(syn,
                                             entity=entity,
                                             ifcollision="keep.local",
                                             path=path)
        except synapseclient.core.exceptions.SynapseUnmetAccessRestrictions as e:
            print(e)
        except synapseclient.core.exceptions.SynapseMd5MismatchError as e:
            print(e)

    def print_arguments(self):
        print("General arguments")
        for entity, folder in zip(self.entities, self.folders):
            print("  > entity: {}\tpath: {}/".format(entity, folder))
        print("  > Output directory: {}".format(self.out_dir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()