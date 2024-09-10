#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import gzip
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--in_gtf", required=True, type=str, help="")
parser.add_argument("--gene_pairs", required=True, type=str, help="")
parser.add_argument("--feature_name", required=False, type=str, default="gene_name", help="")
parser.add_argument("--autosomes_only", action="store_true", default=False, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

os.makedirs(os.path.dirname(args.out), exist_ok=True)


def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


print("Parsing annotation file...")

tmp_fpath = args.out + ".WithDuplicates.txt.gz"
fhtmp = gzopen(tmp_fpath, mode="w")
fhtmp.write("\t".join(["Platform", "ArrayAddress", "Symbol", "Chr", "ChrStart", "ChrEnd", "Probe", "Strand"]) + "\n")

accepted_chr = [str(i) for i in range(1, 23)]
if not args.autosomes_only:
    accepted_chr.extend(["X", "MT", "M", "Y"])

features = set()
duplicated = set()
with gzopen(args.in_gtf, mode="r") as f:
    for line in f:
        if line.startswith("#"):
            continue

        seqname, source, feature, start, end, score, strand, frame, attribute = line.rstrip("\n").split("\t")
        if feature != "gene":
            continue
        attributes = {key: val.replace('"', '') for key, val in (attr.split(' ') for attr in attribute.split("; "))}

        if attributes[args.feature_name] in features:
            duplicated.add(attributes[args.feature_name])

        if seqname.startswith("chr"):
            seqname = seqname.replace("chr", "")

        if seqname not in accepted_chr:
            continue

        fhtmp.write("\t".join([source, attributes[args.feature_name], attributes["gene_name"], seqname, start, end, attributes["gene_id"], strand]) + "\n")
        features.add(attributes[args.feature_name])
f.close()
fhtmp.close()

print("  Saved {:,} features with {:,} duplicates".format(len(features), len(duplicated)))

fhout = gzopen(args.out + ".txt.gz", mode="w")
n_removed = 0
with gzopen(tmp_fpath, mode="r") as f:
    for line in f:
        if line.rstrip("\n").split("\t")[1] in duplicated:
            n_removed += 1
            continue
        fhout.write(line)
f.close()
fhout.close()

print("  Removed {:,} features with multiple entries".format(n_removed))

print("Done.")
