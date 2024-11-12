import gzip
import sys
import os
import glob

debug = False

def gzopen(file):
	if file.endswith(".gz"):
		return gzip.open(file,'rt')
	else:
		return open(file,'r')

if len(sys.argv) < 9:
    print("Usage: createbatches.py batchnameprefix expfile.txt.gz gte.txt genotype.vcf.gz genelist.txt.gz annotation.txt.gz " +
          "snp_genepairs.txt group_annotation.txt template.sh nrmaxgenesperbatch outdir")
    actr = 1
    for arg in sys.argv:
        print(str(actr)+"\t"+arg)
        actr += 1
    sys.exit()

batchnameprefix = sys.argv[1]
expfile = sys.argv[2]
gte = sys.argv[3]
genotype = sys.argv[4]
genelist = sys.argv[5]
annotation = sys.argv[6]
snpgenepairs=sys.argv[7]
groupsfile = sys.argv[8]
template = sys.argv[9]
nrgenes = int(sys.argv[10])
out = sys.argv[11]

print()
print("-----------------------")
print("Batch creat0r")
print("-----------------------")

vctr = 1
for var in sys.argv:
	print(str(vctr)+" - "+var)
	vctr+=1
print()
if not out.endswith("/"):
    out = out + "/"


def writeJob(exp, gte, gt, template, batchfile, jobfile, outprefix, logprefix, chr, annotation, snpgenepairs, groups, jobname):
    if debug:
        print("Writing job: "+jobfile)
    fh = open(template,'r')
    lines = fh.readlines()
    fh.close()
    fho = open(jobfile,'w')
    for line in lines:
        line = line.replace("JOBNAME",jobname)
        line = line.replace("GENOTYPE",gt)
        line = line.replace("GTE",gte)
        line = line.replace("EXPRESSION",exp)
        line = line.replace("CHROM",str(chr))
        line = line.replace("BATCHFILE",batchfile)
        line = line.replace("OUTPREFIX",outprefix)
        line = line.replace("LOGPREFIX",logprefix)
        line = line.replace("ANNOTATION", annotation)
        line = line.replace("SNPGENELIMIT", snpgenepairs)
        line = line.replace("GROUPS", groups)
        fho.write(line)
    fho.close()


def checkDir(path):
    if os.path.exists(path):
        # delete contents
        files = glob.glob(path+"*")
        for file in files:
            if debug:
                print("Removing: "+file)
            os.remove(file)
    else:
        if debug:
            print("Creating dir: "+path)
        os.mkdir(path)


abspath = os.path.abspath(out)
checkDir(abspath+"/batches/")
checkDir(abspath+"/output/")
# checkDir(abspath+"/jobs/")
checkDir(abspath+"/logs/")

# read the ids from the expression file
print("Reading: "+expfile)
fh = gzopen(expfile)
genesPresentInData = set()
line_counter = 0
fh.readline()
for line in fh:
	gene = line.split('\t', 1)[0]
	genesPresentInData.add(gene)
	line_counter += 1
	if line_counter % 10000 == 0:
		print(f"{line_counter} lines parsed", end="\r")
print(f"{line_counter} lines parsed", end="\n")
fh.close()
print("Genes in expfile: {}".format(len(genesPresentInData)))

# read gene set
print("Reading: "+genelist)
geneLimitSet = set()
fh = gzopen(genelist)
for line in fh:
	gene = line.strip()
	if gene in genesPresentInData:
		geneLimitSet.add(gene)
fh.close()
print("Genes in genelist: {}".format(len(geneLimitSet)))

# read annotation
print("Annotation: " + annotation)
fh = gzopen(annotation)
header=fh.readline().lower().strip().split("\t")
gcol = -1
chrcol = -1
poscol = -1

for i in range(len(header)):
	if header[i] == "gene":
		gcol = i
	elif header[i] == "chr":
		chrcol = i
	elif header[i] == "chrstart":
		poscol = i

genesPerChr = {}
annotread = 0
for line in fh:
    elems = line.strip().split("\t")
    gene = elems[gcol]
    if gene in geneLimitSet:
        chr = int(elems[chrcol])
        if chr < 23:
            pos = int(elems[poscol])
            chrgenes = genesPerChr.get(chr)
            if chrgenes is None:
                chrgenes = []
            chrgenes.append(gene)
            genesPerChr[chr] = chrgenes
            annotread = annotread + 1
fh.close()
print("Annotation read for {} splice events".format(annotread))

# read gene groupings
print("Reading: "+groupsfile)
groups = None
geneToGroup = None
if os.path.exists(groupsfile):
	groups = {}
	geneToGroup = {}
	fh = None
	if groupsfile.endswith(".gz"):
		fh = gzip.open(groupsfile,'rt')
	else:
		fh = open(groupsfile,'r')
	for line in fh:
		elems = line.strip().split("\t")
		id = elems[0]
		grp = elems[1]
		grpset = groups.get(grp)
		if grpset is None:
			grpset = set()
		grpset.add(id)
		groups[grp] = grpset
		geneToGroup[id] = grp
	fh.close()
	print("{} groups loaded from {}".format(len(groups), groupsfile))

chrs = []
for chr in genesPerChr.keys():
	chrs.append(chr)
chrs.sort()

#for chr in genesPerChr.keys():
for chr in chrs:
	print("Writing batches for chr: "+str(chr))
	if groups is not None:
		bctr = 1
		gctr = 0
		bgctr = 0
	
		chrgenes = genesPerChr.get(chr)

		# determine groups on this chromosome
		groupsOnChr = set()
		ok = True
		for gene in chrgenes:
			grp = geneToGroup.get(gene)
			if grp is None:
				print("Error: groups defined, but "+gene+" not in a group")
				ok = False
			else:
				groupsOnChr.add(grp)
		if not ok:
			sys.exit(0)
		print("{} groups for chr {}".format(len(groupsOnChr), chr))
		groupsOnChrArr = []
		for grp in groupsOnChr:
			groupsOnChrArr.append(grp)
		groupsOnChrArr.sort()

		# write the first job
		batchname = "chr"+str(chr)+"-batch-"+str(bctr)
		batchfile = abspath+"/batches/"+batchname+".txt"
		jobfile = abspath+"/jobs/"+batchname+".sh"
		outprefix = abspath+"/output/"+batchname
		logprefix = abspath+"/logs/"+batchname
		chrgenotype = genotype.replace("CHR", str(chr))
		jobname = batchnameprefix +"-"+batchname
		# writeJob(expfile, gte, chrgenotype, template, batchfile, jobfile, outprefix, logprefix, chr, annotation, snpgenepairs,groupsfile, jobname)
		bgout = open(batchfile, 'w')

		grpctr = 0
		while grpctr < len(groupsOnChrArr):
			currentGroup = groupsOnChrArr[grpctr]
			currentGroupGenes = groups.get(currentGroup)

			# write all genes of this group into batch
			for gene in currentGroupGenes:
				bgout.write(gene+"\n")
				bgctr += 1
			grpctr += 1
			# check if there is a next group
			if grpctr < len(groupsOnChrArr):
				nextgrp = groupsOnChrArr[grpctr]
				nextgrpGenes = groups.get(nextgrp)
				# check if this would overflow the batch
				if bgctr >= nrgenes or bgctr + len(nextgrpGenes) >= nrgenes:
					# if so, close thecurrent batch.. and 
					bgout.close()

					# start a new batch
					bctr += 1

					batchname = "chr"+str(chr)+"-batch-"+str(bctr)
					batchfile = abspath+"/batches/"+batchname+".txt"
					# print a warning if next batch is large!
					if len(nextgrpGenes) > nrgenes:
						print("Warning: group "+nextgrp+" will have a large batch: n="+str(len(nextgrpGenes))+ " - "+batchfile)

					jobfile = abspath+"/jobs/"+batchname+".sh"
					outprefix = abspath+"/output/"+batchname
					logprefix = abspath+"/logs/"+batchname
					chrgenotype = genotype.replace("CHR", str(chr))	
					jobname = batchnameprefix +"-"+batchname
					# writeJob(expfile, gte, chrgenotype, template, batchfile, jobfile, outprefix, logprefix, chr, annotation, snpgenepairs, groupsfile, jobname)
					bgout = open(batchfile, 'w')
					bgctr = 0
			print("Writing group: {}".format(grpctr), end='\r')
		# close current batch, if any
		if not bgout.closed:
			bgout.close()
		print()
		print("Done with chr {}".format(chr))
	else:
		# no groups defined	
		print("No groups defined. Writing jobs")
		bctr = 1
		gctr = 0
		bgctr = 0
	
		chrgenes = genesPerChr.get(chr)

		# write the first job
		batchname = "chr"+str(chr)+"-batch-"+str(bctr)
		batchfile = abspath+"/batches/"+batchname+".txt"
		jobfile = abspath+"/jobs/"+batchname+".sh"
		outprefix = abspath+"/output/"+batchname
		logprefix = abspath+"/logs/"+batchname
		chrgenotype = genotype.replace("CHR", str(chr))
		jobname = batchnameprefix +"-"+batchname
		# writeJob(expfile, gte, chrgenotype, template, batchfile, jobfile, outprefix, logprefix, chr, annotation, snpgenepairs, groupsfile, jobname)
		bgout = open(batchfile, 'w')
		while gctr < len(chrgenes):
			gene = chrgenes[gctr]
			bgout.write(gene+"\n")
			bgctr = bgctr + 1
			# batch is overflowing
			if bgctr > nrgenes:
				# close old batch
				bgout.close()

				# create new batch
				bctr += 1
				batchname = "chr"+str(chr)+"-batch-"+str(bctr)
				batchfile = abspath+"/batches/"+batchname+".txt"
				jobfile = abspath+"/jobs/"+batchname+".sh"
				outprefix = abspath+"/output/"+batchname
				logprefix = abspath+"/logs/"+batchname
				chrgenotype = genotype.replace("CHR", str(chr))
				jobname = batchnameprefix +"-"+batchname
				# writeJob(expfile, gte, chrgenotype, template, batchfile, jobfile, outprefix, logprefix, chr, annotation, snpgenepairs, groupsfile, jobname)
				bgout = open(batchfile, 'w')
				bgctr = 0
			gctr += 1

		# close current batch, if any
		if not bgout.closed:
			bgout.close()
print()
print("Done creating batches")
print()
sys.exit()