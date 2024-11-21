#!/usr/bin/env python
# Author: A. Kooijmans, H.J. Westra, and M. Vochteloo

import argparse
import os
import gzip
import glob

parser = argparse.ArgumentParser(description="")
parser.add_argument("--genes", required=True, type=str, help="")
parser.add_argument("--nrgenes", required=True, type=int, help="")
parser.add_argument("--annotation", required=True, type=str, help="")
parser.add_argument("--expgroups", required=False, type=str, help="")
parser.add_argument("--chr", required=False, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)
    
def checkDir(path):
    if os.path.exists(path):
        files = glob.glob(path+"*")
        for file in files:
            os.remove(file)
    else:
        os.mkdir(path)

# Check directories
os.makedirs(os.path.dirname(args.out), exist_ok=True)
checkDir(args.out +"batches/")

# Load genes
print("Loading genes")
fin = gzopen(args.genes,"r")
n_duplicated = 0
genes = set()
for line in fin:
    gene = line.strip()
    if gene in genes:
            print(f"\tWarning, gene '{gene}' is duplicated. Skipping gene.")
            n_duplicated += 1
            continue
    genes.add(gene)
fin.close()
print(f"Loaded {len(genes)} genes, {n_duplicated} genes duplicates found")
del n_duplicated

# Loading groups
groups = None
geneToGroup = None
if args.expgroups is not None:
    print("\nLoading gene groups")
    groups = {}
    geneToGroup = {}
    fin = gzopen(args.expgroups, "r")
    for line in fin:
        values = line.strip().split("\t")
        gene_pair = values[0]
        grp = values[1]
        grp_set = groups.get(grp)
        if grp_set == None:
            grp_set = set()
        grp_set.add(gene_pair)
        groups[grp] = grp_set
        geneToGroup[gene_pair] = grp
    fin.close()
    print(f"{len(groups)} groups loaded")

print("\nLoading annotation")
fin = gzopen(args.annotation, mode = "r")
header = fin.readline().strip().split("\t")
genesPerChr = {}
for line in fin:
    values = line.strip().split("\t")
    gene_pair = values[header.index("Gene")]
    chr = values[header.index("Chr")]
    if args.chr is not None:
        if chr == args.chr:
            if chr not in genesPerChr.keys():
                genesPerChr[chr] = []
            genesPerChr[chr].append(gene_pair)
fin.close()

chromosomes = []
for chr in genesPerChr.keys():
	chromosomes.append(chr)
chromosomes.sort()

print("\nWriting batches")
for chr in chromosomes:
    if groups is not None:
        batch_ctr = 1
        batch_gene_ctr = 0
        
        chrgenes = genesPerChr.get(chr)
        groupsOnChr = set()
        for gene in chrgenes:
            group = geneToGroup.get(gene)
            if group is None:
                print(f"Error: groups defined, but {gene} not in a group")
                exit()
            else:
                groupsOnChr.add(group)
        groupsOnChrArr = sorted(groupsOnChr)

        print(f"{len(groupsOnChrArr)} groups for chr {chr}")

        batchname = f"chr{chr}-batch-{batch_ctr}"
        batchfile = f"{args.out}batches/{batchname}.txt"

        batch_out = open(batchfile,"wt")

        group_ctr = 0
        while group_ctr < len(groupsOnChrArr):
            currentGroup = groupsOnChrArr[group_ctr]
            currentGroupGenes = groups.get(currentGroup)
            # write all genes of this group into batch
            for gene in currentGroupGenes:
                batch_out.write(f"{gene}\n")
                batch_gene_ctr += 1
            
            group_ctr += 1

            # check if there is a next group
            if group_ctr < len(groupsOnChrArr):
                nextgroup = groupsOnChrArr[group_ctr]
                nextgroupGenes = groups.get(nextgroup)
                # check if this would overflow the batch
                if batch_gene_ctr >= args.nrgenes or batch_gene_ctr + len(nextgroupGenes) >= args.nrgenes:
                    batch_out.close()
                    # Start a new batch
                    batch_ctr += 1
                    batchname = f"chr{chr}-batch-{batch_ctr}"
                    batchfile = f"{args.out}batches/{batchname}.txt"
                    # print a warning if next batch is large!
                    if len(nextgroupGenes) > args.nrgenes:
                        print(f"Warning: group {nextgroup} will  have a large batch: n = {len(nextgroupGenes)} - {batchname}")
                    batch_out = open(batchfile,"wt")
                    batch_gene_ctr = 0

            print("Writing group: {}".format(group_ctr), end='\r')

        # Close current batch if any
        if not batch_out.closed:
            batch_out.close()

        print(f"Finished writing batches for chr {chr}")

    else:
        print("\nNo groups defined, writing batches")
        batch_ctr = 0
        gene_ctr = 0
        batch_gene_ctr = 0

        chrgenes = genesPerChr.get(chr)
        batchname = f"chr{chr}-batch-{batch_ctr}"
        batchfile = f"{args.out}batches/{batchname}.txt"
        batch_out = open(batchfile,"wt")
        while gene_ctr < len(chrgenes):
            gene = chrgenes [gene_ctr]
            batch_out.write(f"{gene}\n")
            batch_gene_ctr += 1
            # Check if batch if overflowing
            if batch_gene_ctr > args.nrgenes:
                batch_out.close()
                # Create new batch
                batch_ctr += 1
                batchname = f"chr{chr}-batch-{batch_ctr}"
                batchfile = f"{args.out}batches/{batchname}.txt"
                batch_out = open(batchfile,"wt")
                batch_gene_ctr = 0
            gene_ctr += 1
        if not batch_out.closed:
            batch_out.close()
        print(f"{batch_ctr} batches written")

print("\nDone.\n")