import glob
import sys

if len(sys.argv) < 3:
	print("Usage: indir output.txt")
	sys.close()

files = sys.argv[1]
out = sys.argv[2]

# print(f"Looking for files ending with /chr*/output/chr*-batch-*-TopEffects.txt in {dir}")
# files = glob.glob(f"{dir}/chr*/output/chr*-batch-*-TopEffects.txt")
# files.sort()
print(f"{len(files)} files detected")

fho = open(out,'w')
print("open outputfule")
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