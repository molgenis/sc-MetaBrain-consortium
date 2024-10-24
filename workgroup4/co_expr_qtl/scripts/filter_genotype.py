import pandas as pd
import gzip
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("--top_effects", required=True, type=str,  help="eQTL top effects file")
parser.add_argument("--genotype", required=True, type=str,  help="genotype file")
parser.add_argument("--output", required=True, type=str,  help="output file")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

# Load genotype file
top_effects = pd.read_csv(args.top_effects,sep="\t",index_col=None)
top_effects = top_effects[top_effects["qval"] < 0.05] # Only load keep significant effects
top_snps = set(top_effects["SNP"])
del top_effects

# Load genotype file
fin = gzip.open(args.genotype,"rt")
fout = gzip.open(args.output,"wt")

print("Filtering genotype file...")
counter = 0
for line in fin:
    counter += 1
    if line.startswith("#"):
        fout.write(line)
    else:
        snp = line.split("\t")[2]
        if snp in top_snps:
            fout.write(line)
    if counter % 100000 == 0:
        print(f"{counter} lines procesed.",end="\r")
        sys.stdout.flush()

print("\nDone.\n")
fin.close()
fout.close()