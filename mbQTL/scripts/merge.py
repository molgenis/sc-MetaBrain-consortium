#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import gzip
import glob
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--input_dir", required=True, type=str, help="")
parser.add_argument("--outfile", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

os.makedirs(os.path.dirname(args.outfile), exist_ok=True)


def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


print("Loading top effects files...")
files = glob.glob(args.input_dir + "*-TopEffects.txt")
print("\tfound {:,} files".format(len(files)))
if len(files) > 0:
    n_elems = None
    header = None
    pos = {}
    gene_line = {}
    gene_pval = []
    for file in files:
        with gzopen(file, 'r') as fh:
            for fctr, line in enumerate(fh):
                line = line.rstrip("\n")
                elems = line.split("\t")
                if n_elems is None:
                    n_elems = len(elems)
                if len(elems) != n_elems:
                    print("Error, unexpected number of elements.")
                    print(line)
                    exit()

                if fctr == 0:
                    if header is None:
                        header = line
                        pos = {column: index for index, column in enumerate(elems)}
                    if line != header:
                        print("Error, header differs between files.")
                        print(line)
                        exit()
                    continue

                gene = elems[pos["Gene"]]
                if gene in gene_line:
                    if line != gene_line[gene]:
                        print("Error, gene '{}' is tested multiple times and results differ.".format(gene))
                        exit()
                    print("Warning, gene '{}' is tested multiple times. Results were identical so it is save to skip.".format(gene))
                    continue
                gene_line[gene] = line

                pval = float(elems[pos["BetaAdjustedMetaP"]])
                gene_pval.append([gene, pval])
        fh.close()

        # sort by qvalue
        gene_pval.sort(key=lambda p: p[1])

        # Save the output.
        fho = gzopen(args.outfile, 'w')
        fho.write(header + "\n")
        for (gene, _) in gene_pval:
            fho.write(gene_line[gene] + "\n")
        fho.close()

print("Done")