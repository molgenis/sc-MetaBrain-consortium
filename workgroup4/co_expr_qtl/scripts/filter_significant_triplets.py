#!/usr/bin/env python
# Author: A. Kooijmans, H.J. Westra, and M. Vochteloo

import glob
import argparse
import gzip

parser = argparse.ArgumentParser(description="")
parser.add_argument("--significant_egenes", required=True, type=str,  help="")
parser.add_argument("--egene_column", required=True, type=str,  help="")
parser.add_argument("--signif_column", required=True, type=str,  help="")
parser.add_argument("--indir", required=True, type=str,  help="")
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

print(f"Extracting significance using {args.signif_column}")
fin = gzopen(args.significant_egenes,"rt")
header = fin.readline().strip().split("\t")
egene_pval = {}
for line in fin:
	values = line.strip().split("\t")
	egene = values[header.index(args.egene_column)]
	pval = values[header.index(args.signif_column)]
	egene_pval[egene] = float(pval)

print(f"{len(egene_pval)} significant top effects found")

print(f"Looking for files in: {args.indir}")
pattern = f"{args.indir}*-dump-AllEffects.txt.gz"
files = glob.glob(pattern)
print(f"{len(files)} file(s) detected")
print(files)

fout = gzopen(args.output,'w')
ctr = 0
for file in files:
	fin = gzopen(file,'r')
	all_lines = 0
	lctr = 0
	if ctr == 0:
		line = fin.readline()
		fout.write(line)
		header = line.strip().split("\t")
		all_lines += 1
		for line in fin:
			all_lines += 1
			values = line.strip().split("\t")
			egene = values[header.index("Group")]
			pval = float(values[header.index("MetaP")])
			# pval = float(line.strip().split("\t")[13])
			threshold = float(egene_pval[egene])
			if pval < threshold:
				print(f"{egene}: {pval} < {threshold}")
				fout.write(line)
				lctr += 1
	else:
		fin.readline()
		all_lines += 1
		for line in fin:
			values = line.strip().split("\t")
			all_lines += 1
			egene = values[header.index("Group")]
			pval = float(values[header.index("MetaP")])
			threshold = float(egene_pval[egene])
			if pval < threshold:
				fout.write(line)
				lctr += 1

	print(f"Found {lctr}/{all_lines} significant effects for {file}")
	fin.close()
	ctr += 1

fout.flush()
fout.close()
print("Done.")