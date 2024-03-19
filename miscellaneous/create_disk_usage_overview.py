#!/usr/bin/env python3

"""
File:         create_disk_usage_overview.py
Created:      2024/03/15
Last Changed:
Author:       M.Vochteloo

Copyright (C) 2024 University Medical Center Groningen.
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
import argparse
import os
import subprocess
from datetime import datetime

# Third party imports.

# Local application imports.

# Metadata
__program__ = "Create Disk Usage Overview"
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
./create_disk_usage_overview.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--basedir", type=str, required=False, nargs="*", default=["input", "output"], help="")
parser.add_argument("--skip_subfolder", type=str, required=False, nargs="*", default=None, help="")
parser.add_argument("--max_depth", type=int, required=False, default=None, help="")
parser.add_argument("--dry_run", action='store_true', help="")
args = parser.parse_args()

if args.skip_subfolder is None:
    args.skip_subfolder = []

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

def split_path(path):
    items = path.split("/")
    return items[2], items[3], items[4], "./" + "/".join(items[5:])

def run_du_command(command):
    try:
        output = subprocess.check_output(command, stderr=subprocess.STDOUT).decode()
    except subprocess.CalledProcessError as e:
        print("  Warning, cannot read some directories.")
        output = e.output.decode()
    except Exception as e:
        print("  Error, subprocess failed.")
        print(str(e))
        exit()

    return output

def convert_size_to_gb(size):
    multiplication = 1
    unit = None
    if size.endswith("K"):
        unit = "K"
        multiplication = (1 / 1e6)
    elif size.endswith("M"):
        unit = "M"
        multiplication = (1 / 1e3)
    elif size.endswith("G"):
        unit = "G"
        multiplication = 1
    elif size.endswith("T"):
        unit = "T"
        multiplication = 1e3

    if unit is None:
        print("Error, unexpected unit for size value '{}'".format(size))

    size_as_number = float(size.replace(unit, ""))
    return round(size_as_number * multiplication, 6), size_as_number, unit

def get_depth_of_path(path):
    return len([part for part in path.split("/") if part != ""])

def parse_du_output(output):
    data = []
    lines = output.split("\n")
    for line in lines:
        if line == "" or line.startswith("du: cannot read directory"):
            continue
        size, path = line.split("\t")
        size_in_gb, size, unit = convert_size_to_gb(size=size)
        depth = get_depth_of_path(path=path)
        data.append([size_in_gb, size, unit, depth, path])
    return data

date_str = datetime.now().strftime("%Y-%m-%d")
outpath = os.path.join(str(os.path.dirname(os.path.abspath(__file__))), date_str + '-du.tsv')
fhout = open(os.path.join(outpath), 'w')
outstring = "SizeInGB\tSize\tUnit\tDepth\tGroup\tDisk\tFolder\tQueryPath\tFullPath"
print(outstring)
fhout.write(outstring + "\n")

for basedir in args.basedir:
    print("Processing basedir '{}'".format(basedir))
    if not os.path.exists(basedir):
        print("Error, directory does not exist.")
        exit()

    group, disk, base_folder, remaining_path = split_path(path=basedir)

    command = ["du", "-h", basedir]
    if args.max_depth is not None:
        command.extend(["--max-depth", str(args.max_depth)])
    print("  Command: " + " ".join(command))
    if args.dry_run:
        print("  Not executed due to dry-run")
        continue
    output = run_du_command(command=command)
    data = parse_du_output(output=output)

    for line in data:
        values = line[:-1] + [group, disk, base_folder, remaining_path] + [line[-1]]
        outstring = "\t".join([str(value) for value in values])
        print(outstring)
        fhout.write(outstring + "\n")

    print("")

fhout.close()
