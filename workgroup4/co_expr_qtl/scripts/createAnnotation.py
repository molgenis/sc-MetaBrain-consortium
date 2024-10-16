import gzip
import sys


def fopen(file,mode):
	if file.endswith(".gz"):
		if mode == 'r':
			return gzip.open(file,'rt')
		else:
			return gzip.open(file,'wt',4)
	return open(file,mode)


infile=sys.argv[1]
geneannotation=sys.argv[2]
outfile=sys.argv[3]


print("Reading gene annotation: "+geneannotation)
fh = fopen(geneannotation,'r')
data = {}
header = fh.readline().lower().strip().split("\t")
gcol = -1
chrcol = -1
stacol = -1
stocol = -1
strandcol = -1

for i in range(len(header)):
	if header[i] == "genesymbol":
		gcol = i
	elif header[i] == "chr":
		chrcol = i
	elif header[i] == "chrstart":
		stacol = i
	elif header[i] == "chrend":
		stocol = i
	elif header[i] == "strand":
		strandcol = i

print(f"gene -> {gcol}")
print(f"chr -> {chrcol}")
print(f"start -> {stacol}")
print(f"stop -> {stocol}")
print(f"strand -> {strandcol}")

for line in fh:
	elems = line.strip().split("\t")
	gene = elems[gcol]
	chr = elems[chrcol]
	sta = elems[stacol]
	sto = elems[stocol]
	strand = elems[strandcol]
	data[gene] = [chr,sta,sto,strand]
fh.close()

print(f"{len(geneannotation)} annotations read.")

print("Reading genepairs: "+infile)
fh = fopen(infile,'r')
print("Writing: "+outfile)
fho = fopen(outfile,'w')
geneswitherrors = set()
fho.write("Gene\tGeneSymbol\tChr\tChrStart\tChrEnd\tStrand\n")
wctr = 0
lctr = 0
for line in fh:
	p = line.strip().split("\t",2)[0] # also allow coexpression matrix as input
	pelems = p.split("_")
	g1 = pelems[0]
	annot = data.get(g1)
	if annot is None and g1 not in geneswitherrors:
		print(f"Warning: no annoation for {g1}")
		geneswitherrors.add(g1)
	elif annot is not None:
		fho.write(p+"\t"+p+"\t"+ "\t".join(annot) +"\n" )
		wctr += 1
	lctr += 1
	if lctr % 10000 == 0:
		print(f'{lctr} lines read, {wctr} lines written.',end='\r')
print(f'{lctr} lines read, {wctr} lines written.',end='\n')
fho.close()
fh.close()
print()