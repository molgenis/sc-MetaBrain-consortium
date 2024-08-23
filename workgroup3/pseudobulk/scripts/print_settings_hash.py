#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import hashlib

parser = argparse.ArgumentParser(description="")
parser.add_argument("--ncount_rna", required=False, type=int, default=500, help="")
parser.add_argument("--nfeature_rna", required=False, type=int, default=0, help="")
parser.add_argument("--percent_rb", required=False, type=int, default=0, help="")
parser.add_argument("--percent_mt", required=False, type=int, default=5, help="")
parser.add_argument("--malat1", required=False, type=int, default=0, help="")
parser.add_argument("--min_cells", required=False, type=int, default=5, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

def get_md5_hash(s):
    return hashlib.md5(s.encode()).hexdigest()

print("Settings hash:")
settings = {
    "nCount_RNA": args.ncount_rna,
    "nFeature_RNA": args.nfeature_rna,
    "percent.rb": args.percent_rb,
    "percent.mt": args.percent_mt,
    "MALAT1": args.malat1,
    "minCells": args.min_cells
}
md5_hash = get_md5_hash(str(settings))
print("\t{}: {}".format(md5_hash, settings))

print("\nDone")