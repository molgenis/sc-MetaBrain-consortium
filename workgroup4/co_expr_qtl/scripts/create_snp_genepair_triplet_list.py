import sys


eqtlfile = sys.argv[1]
querychr = int(sys.argv[2])
genepairfile=sys.argv[3]
tripletout=sys.argv[4]
tripletgroupout=sys.argv[5]


print("Making co-eQTL snp - genepair triplets")
for arg in sys.argv:
	print(f"{arg}")

sys.stdout.flush()

snpmap = {}
fh = open(eqtlfile,'r')
fh.readline()
for line in fh:
    elems = line.strip().split("\t")
    snp = elems[0]
    chr = int(elems[2])
    gene = elems[1]
    if int(chr) == querychr:
        snpmap[gene] = snp
fh.close()
genes = snpmap.keys()

print(f"{len(snpmap)} snp-gene pairs loaded")

fh = open(genepairfile,'r')
fh.readline()
fho = open(tripletout,'w')
fho2 = open(tripletgroupout,'w')

wctr = 0
for line in fh:
    id = line.strip()
    idpt = id.split("_")
    gene1 = idpt[0]
    gene2 = idpt[1]

    if gene1 in genes:
        snp = snpmap.get(gene1)
        if snp is not None:
            fho.write(snp+"\t"+id+"\n")
            fho2.write(id+"\t"+gene1+"\n")
            wctr += 1
#    if gene2 in genes:
#        snp = snpmap.get(gene2)
#	if snp is not None:
#	        fho.write(snp+"\t"+id+"\n")
#      		fho2.write(id+"\t"+gene2+"\n")

fho.close()
fho2.close()
fh.close()
print(f"{wctr} lines written - done.")
print()
sys.stdout.flush()
print("Done.")