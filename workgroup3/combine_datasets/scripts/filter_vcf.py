#!/usr/bin/env python
# Author: M. Vochteloo
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--vcf", required=True, type=str, help="")
parser.add_argument("--variants", required=False, type=str, default=None, help="")
parser.add_argument("--samples", required=False, type=str, default=None, help="")
parser.add_argument("--outfile", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import gzip

variants = None
if args.variants is not None:
    print("Loading variants...")
    variants = set()
    with gzip.open(args.variants, 'rt') as f:
        for line in f:
            variants.add(line.rstrip("\n"))
    f.close()
    print("\tLoaded {} variants".format(len(variants)))

samples = None
if args.samples is not None:
    print("Loading samples...")
    samples = set()
    with gzip.open(args.samples, 'rt') as f:
        for line in f:
            samples.add(line.rstrip("\n"))
    f.close()
    print("\tLoaded {} samples".format(len(samples)))


print("Filtering VCF file...")
fhin = gzip.open(args.vcf, 'rt')
fhout = gzip.open(args.outfile, 'wt')

n_parsed = 0
n_written = 0
indices = []
vcf_header_columns = {"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"}
id_column = None
for line in fhin:
    if (n_parsed > 0) and (n_parsed % 500000 == 0):
        print("  Processed {:,} variants of which {:,} were kept for {:,} samples".format(n_parsed, n_written, len(indices)), end='\r', flush=True)

    if line.startswith("##"):
        fhout.write(line)
        continue

    values = line.rstrip("\n").split("\t")

    if line.startswith("#"):
        for index, value in enumerate(values):
            if (samples is None) or (value in vcf_header_columns) or (value in samples):
                indices.append(index)
                continue

        fhout.write("##filter_vcf_variants.py --vcf {} --variants {} --samples {} --outfile {}\n".format(args.vcf, args.variants, args.samples, args.outfile))
        id_column = values.index("ID")
        fhout.write("\t".join([values[index] for index in indices]) + "\n")
        continue

    variant = values[id_column]
    if (variants is None) or (variant in variants):
        fhout.write("\t".join([values[index] for index in indices]) + "\n")
        n_written += 1
    n_parsed += 1

fhin.close()
fhout.close()
print("  Processed {:,} variants of which {:,} were kept for {:,} samples".format(n_parsed, n_written, len(indices)), flush=True)
print("\nDone")
