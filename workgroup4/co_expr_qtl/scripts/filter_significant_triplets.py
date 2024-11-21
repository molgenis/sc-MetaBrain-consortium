import glob
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--significant_egenes", required=True, nargs="*", type=str,  help="")
parser.add_argument("--egene_column", required=True, type=str,  help="")
parser.add_argument("--signif_column", required=True, type=str,  help="")
parser.add_argument("--indir", required=True, type=str,  help="")
parser.add_argument("--output", required=True, type=str,  help="")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

print(f"Extracting significance using {args.signif_column}")
fin = open(args.significant_egenes,"rt")
header = fin.readline().strip().split("\t")
fin.readline()
egene_pval = {}
for line in fin:
	egene = line.strip().split("\t")[header.index(args.egene_col)]
	pval = line.strip().split("\t")[header.index(args.signif_column)]
	egene_pval[egene] = pval
fin.close()

print(f"Looking for files in: {args.indir}")
pattern = f"{args.indir}*-dump-TopEffects.txt"
files = glob.glob(pattern)
print(f"{len(files)} files detected")

fho = open(args.out,'w')
print("open outputfile")
ctr = 0
wctr = 0
for file in files:
	fh = open(file,'r')
	lctr = 0
	if ctr == 0:
		fho.write(fh.readline())
		for line in fh:
			egene = line.strip().split("\t")[0]
			pval = float(line.strip().split("\t")[13])
			threshold = float(egene_pval[egene])
			if pval < threshold:
				fho.write(fh.readline())
				lctr += 1
	else:
		fh.readline()
		for line in fh:
			egene = line.strip().split("\t")[0]
			pval = float(line.strip().split("\t")[13])
			threshold = float(egene_pval[egene])
			if pval < threshold:
				fho.write(fh.readline())
				lctr += 1

	print(f"{file}, {lctr} lines")
	fh.close()
	ctr += 1

fho.flush()
fho.close()
print("Done.")