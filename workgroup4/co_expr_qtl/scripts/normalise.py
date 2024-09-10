#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import gzip
import json
import numpy as np
import scipy as sp
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--counts", required=True, type=str, help="")
parser.add_argument("--avg_rd", required=True, type=str, help="")
parser.add_argument("--log1p", action="store_true", default=False, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

dir = os.path.dirname(args.outfile)
if dir != "":
    os.makedirs(dir, exist_ok=True)


def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)

def do_pf(mtx, sf):
    pf = mtx.sum(axis=1).A.ravel()
    pf = sp.sparse.diags(sf/pf) @ mtx
    return pf

# def norm_pf_log_pf(mtx):
#     pf_log_pf = do_pf(np.log1p(do_pf(mtx)))
#     return pf_log_pf

print("Loading count ...")
mtx = None

print("Loading average read count ...")
fh = open(args.avg_rd)
sf = json.load(fh)["avg_read_count"]
fh.close()
print("\taverage read count: {:.4f}".format(sf))

print("Proportional fit Normalise ...")
mtx = do_pf(mtx=mtx, sf=sf)

if args.log1p:
    print("log1p transform ...")
    mtx = np.log1p(mtx)

print("Saving count ...")

print("Done.")
