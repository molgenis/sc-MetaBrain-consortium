import glob
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--significant_egenes", required=True, nargs="*", type=str,  help="")
parser.add_argument("--input", required=True, nargs="*", type=str,  help="")
parser.add_argument("--output", required=True, type=str,  help="")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

fin = open(args.significant_egenes,open="rt")
fin.readline()
egene_pval = {}
for line in fin:
	egene = line.strip().split("\t")[0]
	pval = line.strip().split("\t")[-1]
	egene_pval[egene] = pval
fin.close()

files = args.input
out = args.output

print(f"{len(files)} files detected")

fho = open(out,'w')
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
			pval = line.strip().split("\t")[13]
			threshold = egene_pval[egene]
			if pval < threshold:
				fho.write(fh.readline())
				lctr += 1
	else:
		fh.readline()
		for line in fh:
			egene = line.strip().split("\t")[0]
			pval = line.strip().split("\t")[13]
			threshold = egene_pval[egene]
			if pval < threshold:
				fho.write(fh.readline())
				lctr += 1

	print(f"{file}, {lctr} lines")
	fh.close()
	ctr += 1

fho.flush()
fho.close()
print("Done.")