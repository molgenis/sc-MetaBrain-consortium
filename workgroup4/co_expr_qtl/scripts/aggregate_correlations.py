#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import gzip
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--indir", required=True, type=str,  help="")
parser.add_argument("--chr", required=False, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

os.makedirs(os.path.dirname(args.out), exist_ok=True)


def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


print("Loading file ...")

print("Done.")
