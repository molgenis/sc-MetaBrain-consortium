#!/usr/bin/env python
# Author: A. Kooijmans, H.J. Westra, and M. Vochteloo

import glob
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="")
parser.add_argument("--input", required=True, type=str,  help="")
parser.add_argument("--out", required=True, type=str,  help="")
parser.add_argument("--batches", required=True, type=str,  help="")
parser.add_argument("--signif_column", default = "qval", type=str,  help="")
parser.add_argument("--alpha", default = 0.05, type=float,  help="")
parser.add_argument("--snp_genepair_triplets", required=True, nargs="*", type=str,  help="")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)

# Read in coeqtl results
print(f"Filtering significant effects using a threshold of {args.signif_column} < {args.alpha}")
df = pd.read_csv(args.input,sep="\t")
df = df[df[args.signif_column] < args.alpha]
df = df.sort_values(by = args.signif_column, ascending=True)
df.to_csv(f"{args.out}-TopEffects-significant.txt",index=False,sep="\t")

print("Getting significant groups")
fin = gzopen(f"{args.out}-TopEffects-significant.txt",'r')
fin.readline()
significant = set()
for line in fin:
	elems = line.strip().split("\t")
	grp = elems[0]
	significant.add(grp)
fin.close()
print(f"{len(significant)} significant groups found")

print("\nExtracting gene pairs from previously created batches")
batches = glob.glob(f"{args.batches}chr*/batches/*.txt")
print(f"{len(batches)} batches found")

fout = gzopen(f"{args.out}-genepairs_to_dump.txt",'w')
counter = 0
pairs = set()
for batch in batches:
	counter += 1
	fin = gzopen(batch,'r')
	for line in fin:
		elems = line.split("_")
		g1 = elems[0]
		if g1 in significant:
			fout.write(line)
			pairs.add(line.strip())
	if counter % 100 == 0:
		print(f"{counter}/{len(batches)} files processed ", end='\r')
	fin.close()
fout.close()
print(f"{counter}/{len(batches)} files processed - "+batch, end='\n')

print("\nWriting gene pair SNP triplets")
snp_genepair = args.snp_genepair_triplets
print(f"{len(snp_genepair)} files found")
fout = gzopen(f"{args.out}-snp_genepairs_to_dump.txt",'w')
for file in snp_genepair:
	fin = gzopen(file,'r')
	for line in fin:
		elems = line.strip().split("\t")
		id = elems[1]
		if id in pairs:
			fout.write(line)
	fin.close()
fout.close()

print("\nDone.\n")