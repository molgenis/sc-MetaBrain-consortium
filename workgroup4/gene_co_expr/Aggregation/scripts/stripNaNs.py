#!/usr/bin/env python
# Author: A. Kooijmans and M. Vochteloo

import argparse
import gzip
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--infile", required=True, type=str, help="Input file")
parser.add_argument("--outfile", required=True, type=str, help="Output file")
parser.add_argument("--minobs", required=False, type=int, default=15, help="Min nr non-nan observations [default: 15]")
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
        if mode == "w":
            return gzip.open(file, mode + 't', 4)
        else:
            return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


fh = gzopen(args.infile,'r')
fho = gzopen(args.outfile,'w')
header = fh.readline()
fho.write(header)

nrsamples = len(header.rstrip("\n").split("\t")) - 1
if args.minobs > nrsamples:
    print(f"Warning: minobs {args.minobs} smaller than number of samples: {nrsamples}. Replacing minimum with {nrsamples}")
    args.minobs = nrsamples

print(f"Removing features with < {args.minobs} observations")
fho.write(fh.readline()) # TODO: is this correct?
lctr = 0
wctr = 0
for line in fh:
    elems = line.rstrip("\n").split("\t")
    # sample = elems[0]
    nrobs = 0
    for i in range(1, len(elems)):
        if elems[i][0] != "N":
            nrobs += 1
    if nrobs >= args.minobs:
        fho.write(line)
        wctr += 1
    lctr += 1
    if lctr % 10000 == 0:
        perc = (wctr / lctr) * 100
        print(f"{lctr} read, {wctr} written {perc:2f} %", end='\r')
perc = (wctr / lctr) * 100
print(f"{lctr} read, {wctr} written {perc:2f} %", end='\n')

fh.close()
fho.close()

print("Done")
