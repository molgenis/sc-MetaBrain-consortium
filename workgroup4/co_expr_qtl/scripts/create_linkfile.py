#!/usr/bin/env python
# Author: A. Kooijmans, H.J. Westra, and M. Vochteloo

import gzip
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--corr", required=True, type=str,  help="")
parser.add_argument("--wg0", required=True, type=str,  help="")
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


print("Exctracting sample IDs from correlation file header")
fin = gzopen(args.corr,"r")
samples = fin.readline().strip().split("\t")[1:]
fin.close()

print(f"{len(samples)} sample IDs found")

print("\nLinking sample IDs to dataset using the WG0 poolsheet")
df = pd.read_csv(args.wg0,sep="\t")
pool_dataset = df.groupby('Pool')['Dataset'].apply(list).to_dict()
print(f"{len(pool_dataset)} pool-dataset links found")

print("\nWriting linkfile")
fho = gzopen(args.output,'w')
seen = set()
wctr = 0
for sample in samples:
    pool,id = sample.split(".")
    if id not in seen:
        seen.add(id)
        dataset = pool_dataset[pool][0]
        fho.write(id+"\t"+sample+"\t"+dataset+"\n")
        wctr += 1
fho.close()

print(f"{wctr} unique samples written")
print("\nDone.\n")