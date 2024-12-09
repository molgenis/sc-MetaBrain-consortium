#!/usr/bin/env python
# Author: A. Kooijmans, H.J. Westra, and M. Vochteloo

import argparse
import sys
import gzip

parser = argparse.ArgumentParser(description="")
parser.add_argument("--gene_pairs", required=True, type=str,  help="")
parser.add_argument("--gene_annotation", required=True, type=str,  help="")
parser.add_argument("--output", required=True, type=str,  help="")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

sys.stdout.flush()

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        if mode == "w":
            return gzip.open(file, mode + 't', 4)
        else:
            return gzip.open(file, mode + 't')
    else:
        return open(file, mode)
    
print("Reading gene annotation file")
fin = gzopen(args.gene_annotation,"r")
header = fin.readline().strip().split("\t")
expected_header = ["Gene","GeneSymbol","Chr","ChrStart","ChrEnd","Strand"]

if header != expected_header:
    print(f"Error: Header does not match the expected format.\n"
          f"Found: {header}\n"
          f"Expected: {expected_header}")
    exit()

annotation = {}
counter = 0
for line in fin:
    values = line.strip().split("\t")
    gene = values[header.index("GeneSymbol")]
    chr = values[header.index("Chr")]
    start = values[header.index("ChrStart")]
    stop = values[header.index("ChrEnd")]
    strand = values[header.index("Strand")]
    annotation[gene] = [chr,start,stop,strand]
fin.close()
print(f"{len(annotation.keys())} genes loaded")

print("\nWriting annotation file for the provided gene pairs")
fin = gzopen(args.gene_pairs,"r")
fout = gzopen(args.output,"w")

fout.write("Gene\tGeneSymbol\tChr\tChrStart\tChrEnd\tStrand\n")
missing_annotation = {}
counter = 0
for line in fin:
    gene_pair = line.strip()
    gene1,gene2 = gene_pair.split("_")
    gene_annot = annotation[gene1]
    if gene_annot is None and gene1 not in missing_annotation:
        print(f"Warning, no annotation available for {gene1}")
        missing_annotation.add(gene1)
    elif gene_annot is not None:
        fout.write(gene_pair + "\t" + gene_pair + "\t" + "\t".join(gene_annot) + "\n")
        counter += 1
    if counter % 100000 == 0:
        print(f"{counter} lines written.",end='\r')
fin.close()
fout.close()
print(f"{counter} lines written.")
print("\nDone.\n")