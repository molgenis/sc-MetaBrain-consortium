import glob
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--input", required=True, nargs="*", type=str,  help="")
parser.add_argument("--output", required=True, type=str,  help="")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

files = args.input
out = args.output

# print(f"Looking for files ending with /chr*/output/chr*-batch-*-TopEffects.txt in {dir}")
# files = glob.glob(f"{dir}/chr*/output/chr*-batch-*-TopEffects.txt")
# files.sort()
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
		lctr += 1
	else:
		fh.readline()
	#	print("skip header")
		lctr += 1

	for line in fh:
		fho.write(line)
		wctr+=1
	#	print("wrote kn")
		lctr += 1
	print(f"{file}, {lctr} lines")

	fh.close()
	ctr += 1
fho.flush()
fho.close()
print(f"{wctr} lines written?")
print("Done.")