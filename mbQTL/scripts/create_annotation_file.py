#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import gzip
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--data", required=True, type=str, help="")
parser.add_argument("--feature_name", required=False, type=str, default="gene_name", help="")
parser.add_argument("--autosomes_only", action="store_true", default=False, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

os.makedirs(os.path.dirname(args.out), exist_ok=True)

is_gtf = args.data.endswith(".gtf")
if args.feature_name is None and is_gtf:
    print("Error, --feature_name is required in --data is a GTF file.")
    exit()
if args.feature_name is not None and not is_gtf:
    print("Warning, --feature_name is ignored since --data is a GTF file.")


def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


print("Parsing annotation file...")

tmp_fpath = args.out + ".WithDuplicates.txt.gz"
fhtmp = gzopen(tmp_fpath, mode="w")
fhtmp.write("\t".join(["Gene", "GeneSymbol", "Chr", "ChrStart", "ChrEnd", "Strand"]) + "\n")

accepted_chr = [str(i) for i in range(1, 23)]
if not args.autosomes_only:
    accepted_chr.extend(["X", "MT", "M", "Y"])

features = set()
duplicated = set()
n_chr_skip = 0
n_unsorted_pos = 0
n_flipped = 0
n_duplicated = 0
n_features = 0
pos = {}
if is_gtf:
    pos = {"seqname": 0, "source": 1, "feature": 2, "start": 3, "end": 4, "score": 5, "strand": 6, "frame": 7, "attribute": 8}
with gzopen(args.data, mode="r") as f:
    for i, line in enumerate(f):
        # If gtf file, skip the header.
        if is_gtf and line.startswith("#"):
            continue

        # Split the value on tab.
        values = line.rstrip("\n").split("\t")

        # If not GTF file, parse the header.
        if not is_gtf and i == 0:
            pos = {value: index for index, value in enumerate(values)}

        if is_gtf:
            if values[pos["feature"]] != "gene":
                continue
            attributes = {key: val.replace('"', '') for key, val in (attr.split(' ') for attr in values[pos["attribute"]].split("; "))}
            if args.feature_name not in attributes:
                print(f"Error, could not find attribute '{args.feature_name}'.")
                exit()

            # Extract the data we need.
            gene = attributes[args.feature_name]
            gene_symbol = attributes["gene_name"]
            chr = values[pos["seqname"]]
            chr_start = values[pos["start"]]
            chr_end = values[pos["end"]]
            strand = values[pos["strand"]]
        else:
            # Extract the data we need.
            gene = values[pos["Gene"]]
            gene_symbol = values[pos["GeneSymbol"]]
            chr = values[pos["Chr"]]
            chr_start = values[pos["ChrStart"]]
            chr_end = values[pos["ChrEnd"]]
            strand = values[pos["Strand"]]

        if gene in features:
            n_duplicated += 1
            duplicated.add(gene)

        if chr.lower().startswith("chr"):
            chr = chr.lower().replace("chr", "")

        # Check if the chromosome is accepted.
        if chr not in accepted_chr:
            n_chr_skip += 1
            continue

        # Make sure that the start and end position are numerically sorted.
        if int(chr_start) > int(chr_end):
            n_unsorted_pos += 1
            tmp_chr_start = chr_start
            chr_start = chr_end
            chr_end = tmp_chr_start
            del tmp_chr_start

        # Check the strand of the gene.
        if strand == "+":
            # Positive strand, tss is at the gene start position.
            fhtmp.write("\t".join([gene, gene_symbol, chr, chr_start, chr_end, strand]) + "\n")
        elif strand == "-":
            n_flipped += 1
            # Negative strand, tss is at the gene end position.
            fhtmp.write("\t".join([gene, gene_symbol, chr, chr_end, chr_start, strand]) + "\n")
        else:
            print("Error, strand is not '-' or '+'.")
            exit()

        # Save gene name for next round.
        features.add(gene)
        n_features += 1
f.close()
fhtmp.close()

print("  Saved {:,} features with {:,} duplicates".format(n_features, n_duplicated))
print("    {:,} were skipped due to chr filter".format(n_chr_skip))
print("    {:,} entries had a chromosome start lower than the chromosome end".format(n_unsorted_pos))
print("    {:,} negative strand effects flipped".format(n_flipped))

fhout = gzopen(args.out + ".txt.gz", mode="w")
n_removed = 0
with gzopen(tmp_fpath, mode="r") as f:
    for line in f:
        if line.rstrip("\n").split("\t")[0] in duplicated:
            n_removed += 1
            continue
        fhout.write(line)
f.close()
fhout.close()

print("  Removed {:,} features with multiple entries".format(n_removed))

print("Done.")
