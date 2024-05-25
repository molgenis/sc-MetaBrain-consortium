#!/usr/bin/env python
# Author: M. Vochteloo
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--wg1_metadata", required=True, type=str, nargs="*", help="")
parser.add_argument("--vcf", required=True, type=str, nargs="*", help="")
parser.add_argument("--outfile", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import gzip

def load_wg1_samples(inpath):
    assignments = set()
    assignment_col = None
    with gzip.open(inpath, 'rt') as f:
        for line in f:
            values = line.rstrip("\n").split("\t")
            if assignment_col is None:
                assignment_col = values.index("Assignment")
                continue

            if values[assignment_col] == "":
                continue

            assignments.add(values[assignment_col])
    return assignments


def load_vcf_samples(inpath):
    header = None
    with gzip.open(inpath, 'rt') as f:
        for line in f:
            if line.startswith("##"):
                continue
            header = line.strip("\n").split("\t")
            break
    f.close()
    return set([col for col in header if col not in ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]])

print("Loading WG1 metadata samples:")
wg1_metadata_samples = set()
for wg1_metadata in args.wg1_metadata:
    samples = load_wg1_samples(inpath=wg1_metadata)
    print("\t{} [N = {:,}]".format(wg1_metadata, len(samples)))
    wg1_metadata_samples.update(samples)
print("WG1 metadata has {:,} samples".format(len(wg1_metadata_samples)))

print("Loading VCF samples:")
vcf_samples = set()
for vcf in args.vcf:
    samples = load_vcf_samples(inpath=vcf)
    print("\t{} [N = {:,}]".format(vcf, len(samples)))
    vcf_samples.update(samples)
print("VCF file has {:,} samples".format(len(vcf_samples)))

overlap = wg1_metadata_samples.intersection(vcf_samples)
print("Overlap is {:,} samples".format(len(overlap)))
if len(overlap) == 0:
    exit()

overlap = list(overlap)
overlap.sort()

with gzip.open(args.outfile, "wt") as f:
    for sample in overlap:
        f.write(sample + "\n")
f.close()
