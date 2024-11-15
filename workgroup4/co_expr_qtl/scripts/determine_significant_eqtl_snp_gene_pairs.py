import sys

eqtlfile = sys.argv[1]
outfile = sys.argv[2]

fh = open(eqtlfile,'r')
fho = open(outfile,'wt')
fho.write("snp\tgene\tchr\n")
fh.readline()
for line in fh:
  elems = line.strip().split("\t")
  gen = elems[0]
  chr = elems[1]
  snp = elems[5]
  qval = float(elems[-2])
  if qval < 0.05:
     fho.write(snp+"\t"+gen+"\t"+chr+"\n")
fho.close()
fh.close()