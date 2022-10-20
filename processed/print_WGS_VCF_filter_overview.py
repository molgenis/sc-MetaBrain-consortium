#!/usr/bin/env python3

"""
File:         print_WGS_VCF_filter_overview.py
Created:      2022/10/22
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
import re

# Third party imports.
import pandas as pd
import numpy as np

# Local application imports.

# Metadata
__program__ = "Print WGS VCF Filter Overview"
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
./print_WGS_VCF_filter_overview.py \
    --workdir /groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/2022-10-20-FilteredGenotypes
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.workdir = getattr(arguments, 'workdir')

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
                            help="The path to the working directory")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Printing overview:")
        data = []
        for chr in [str(chr) for chr in range(1, 23)] + ["X", "Y", "others"]:
            logfile_path = os.path.join(self.workdir, "jobs", "output", "VCFFilter_CHR{}.out".format(chr))
            if os.path.exists(logfile_path):
                data.append([chr, True] + self.read_log_file(filepath=logfile_path))
            else:
                data.append([chr, False, 0, 0, np.nan, False, "00:00:00"])
        data.append([])
        df = pd.DataFrame(data, columns=["chromosome", "started", "parsed", "written", "kept (%)", "finished", "elapsed"])
        df.iloc[df.shape[0] - 1, :] = ["total", df["started"].all(), df["parsed"].sum(), df["written"].sum(), df["kept (%)"].mean(), df["finished"].all(), ""]
        print(df)

    @staticmethod
    def read_log_file(filepath):
        with open(filepath, 'r') as f:
            parsed = 0
            written = 0
            pct_kept = np.nan
            finished = False
            elapsed = ""
            for line in f:
                if re.match("([0-9]+) lines parsed, ([0-9]+) written", line):
                    match = re.search("([0-9]+) lines parsed, ([0-9]+) written", line)
                    parsed = int(match.group(1))
                    written = int(match.group(2))
                    pct_kept = np.round((written / parsed) * 100, 1)
                elif re.match("Done. How about that!", line):
                    finished = True
                elif re.match("[0-9]{7}.batch. ([0-9]{2}:[0-9]{2}:[0-9]{2})", line):
                    match = re.search("[0-9]{7}.batch. ([0-9]{2}:[0-9]{2}:[0-9]{2})", line)
                    elapsed = match.group(1)
                else:
                    pass
        f.close()

        return [parsed, written, pct_kept, finished, elapsed]

    def print_arguments(self):
        print("Arguments:")
        print("  > Working directory: {}".format(self.workdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()