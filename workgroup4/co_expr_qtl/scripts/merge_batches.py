#!/usr/bin/env python
# Author: A. Kooijmans, H.J. Westra, and M. Vochteloo

import argparse
import gzip

parser = argparse.ArgumentParser(description="")
parser.add_argument("--input", required=True, nargs="*", type=str,  help="")
parser.add_argument("--output", required=True, type=str,  help="")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        if mode == "w":
            return gzip.open(file, mode + 't', 4)
        else:
            return gzip.open(file, mode + 't')
    else:
        return open(file, mode)
    
files = args.input.sort()
print(f"Merging {len(files)} files")

fout = gzopen(args.output,'w')
header = None
all_lines = 0
for file in files:
    fin = gzopen(file,'r')
    line = 0
    if header == None:
        print("Writing header")
        header = fin.readline()
        fout.write(header)
        line += 1
    else:
        # Skip header
        fin.readline()
        line += 1
    for line in fin:
        fout.write(line)
        line += 1
        all_lines += 1
    fin.close()
    print(f"{line} lines written for {file}", end = "\r")

fin.flush()
fin.close()
print(f"\n{all_lines} lines written")
print("\nDone.\n")