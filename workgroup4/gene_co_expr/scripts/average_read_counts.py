#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import gzip
import json
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--read_counts", required=True, type=str, nargs="*", help="")
parser.add_argument("--outfile", required=True, type=str, help="")
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

def parse_file(fpath):
    n = 0
    sum = 0
    with gzopen(fpath, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            values = line.rstrip("\n").split("\t")
            sum += float(values[1])
            n += 1
    f.close()
    return sum, n

print("Parsing files ...")
total_sum = 0
total_n = 0
for rc_fpath in args.read_counts:
    sum, n = parse_file(fpath=rc_fpath)
    total_sum += sum
    total_n += n
del sum, n

print("Calculating average read counts ...")
if total_n == 0:
    raise ValueError("Devision by 0.")
avg_read_count = total_sum / total_n

print("Saving file ...")

with gzopen(args.outfile, 'w') as f:
    json.dump({"avg_read_count": avg_read_count}, f, indent=4)
f.close()

print("Done.")
