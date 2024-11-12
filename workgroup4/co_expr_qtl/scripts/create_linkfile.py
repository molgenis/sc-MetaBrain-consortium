import gzip
import pandas as pd
import argparse
import sys
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--sample_ids", required=True, type=str,  help="")
parser.add_argument("--wg0_file_directories", required=True, type=str,  help="")
parser.add_argument("--output", required=True, type=str,  help="")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")


fh = open(args.sample_ids,"rt")
elems = []
for line in fh:
    elems.append(line.strip())
fh.close()

print(f"{len(elems)} sample ids")

fho = open(args.output,'w')
# header = fh.readline()
seen = set()

df = pd.read_csv(args.wg0_file_directories,sep="\t")
pool_dataset = df.groupby('Pool')['Dataset'].apply(list).to_dict()

print(f"{len(pool_dataset)} pool-dataset links")

wctr = 0
for elem in elems:
    pool,ind = elem.split(".")
    if ind not in seen:
        seen.add(ind)
        dataset = pool_dataset[pool][0]
        fho.write(ind+"\t"+elem+"\t"+dataset+"\n")
        wctr += 1
fho.close()

print(f"{wctr} unique samples written")
print("Done.")