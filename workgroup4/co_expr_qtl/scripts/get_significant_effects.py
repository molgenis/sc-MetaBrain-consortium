import glob
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="")
parser.add_argument("--input", required=True, type=str,  help="")
parser.add_argument("--out", required=True, type=str,  help="")
parser.add_argument("--batches", required=True, type=str,  help="")
parser.add_argument("--snp_genepair_triplets", required=True, nargs="*", type=str,  help="")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

# Read in coeqtl results
print("Filtering significant effects")
df = pd.read_csv(args.input,sep="\t")

print("Extracting significant genes")
df = df[df["qval"]<0.05]
df = df.sort_values(by='qval', ascending=True)
df.to_csv(f"{args.out}EX-TopEffectswithQval-significant.txt",index=False,sep="\t")


fh = open(f"{args.out}EX-TopEffectswithQval-significant.txt",'r')
fh.readline()

significant = set()

for line in fh:
	elems = line.strip().split("\t")
	grp = elems[0]
	significant.add(grp)

fh.close()

print(f"Number of significant eGenes: {len(significant)}")

batches = glob.glob(f"{args.batches}chr*/batches/*.txt")
print()
fho = open(f"{args.out}genepairs_to_dump.txt",'w')
bctr = 0

pairs = set()

for batch in batches:
	bctr += 1
	fh = open(batch,'r')
	for line in fh:
		elems = line.split("_")
		g1 = elems[0]
		if g1 in significant:
			fho.write(line)
			pairs.add(line.strip())
	if bctr % 100 == 0:
		print(f"{bctr}/{len(batches)} files processed ", end='\r')
	fh.close()

print(f"{bctr}/{len(batches)} files processed - "+batch, end='\n')
fho.close()

sglfiles = args.snp_genepair_triplets
fho = open(f"{args.out}snp_genepairs_to_dump.txt",'w')
for sglfile in sglfiles:
	print(f"{sglfile} linkfile!")
	fh = open(sglfile,'r')
	for line in fh:
		elems = line.strip().split("\t")
		id = elems[1]
		if id in pairs:
			fho.write(line)
	fh.close()
fho.close()
