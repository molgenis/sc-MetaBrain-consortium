#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import gzip
import os


parser = argparse.ArgumentParser(description="")
parser.add_argument("--vcf", required=True, type=str, help="")
parser.add_argument("--variant_stats", required=False, type=str, default=None, help="")
parser.add_argument("--chr", required=False, type=str, default=None, help="")
parser.add_argument("--maf", required=False, type=float, default=None, help="")
parser.add_argument("--cr", required=False, type=float, default=None, help="")
parser.add_argument("--hwep", required=False, type=float, default=None, help="")
parser.add_argument("--r2", required=False, type=float, default=None, help="")
parser.add_argument("--variants", required=False, type=str, default=None, help="")
parser.add_argument("--samples", required=False, type=str, default=None, help="")
parser.add_argument("--outfile", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

dir = os.path.dirname(args.outfile)
if dir != "":
    os.makedirs(dir, exist_ok=True)

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


pass_qc_variants = None
n_chr_exlude = 0
n_maf_exlude = 0
n_call_exlude = 0
n_hwep_exlude = 0
n_r2_exlude = 0
if args.variant_stats is not None and (args.maf is not None or args.cr is not None or args.r2 is not None):
    print("Loading variants stats...")
    pass_qc_variants = set()
    nvariants = 0
    with gzopen(args.variant_stats, 'r') as f:
        for i, line in enumerate(f):
            values = line.rstrip("\n").split("\t")
            if i == 0:
                pos = {value: index for index, value in enumerate(values)}
                if len({"CHR", "MAF", "CALL", "HWE", "MACH_R2"}.difference(set(values))) > 0:
                    print("Error, --variant_stats file does not have all required columns.")
                    exit()
                continue

            if args.chr is not None and values[pos["CHR"]] != args.chr:
                n_chr_exlude += 1
                continue

            if args.maf is not None and float(values[pos["MAF"]]) < args.maf:
                n_maf_exlude += 1
                continue

            if args.cr is not None and float(values[pos["CALL"]]) < args.cr:
                n_call_exlude += 1
                continue

            if args.hwep is not None and float(values[pos["HWE"]]) < args.hwep:
                n_hwep_exlude += 1
                continue

            if args.r2 is not None and float(values[pos["MACH_R2"]]) < args.r2:
                n_r2_exlude += 1
                continue

            pass_qc_variants.add(values[pos["ID"]])
            nvariants += 1
    f.close()
    print("\tParsed {:,} variants\tCHR exclude: {:,}\tMAF exclude: {:,}\tCALL exclude: {:,}\tHWE exclude: {:,}\tMACH_R2 exclude: {:,}\tN kept: {:,}".format(nvariants, n_chr_exlude, n_maf_exlude, n_call_exlude, n_hwep_exlude, n_r2_exlude, nvariants))

    if nvariants == 0:
        print("Error, no variants kept from --variant_stats file.")
        exit()

select_variants = None
if args.variants is not None:
    print("Loading variants...")
    select_variants = set()
    nvariants = 0
    with gzopen(args.variants, 'r') as f:
        for line in f:
            select_variants.add(line.rstrip("\n"))
            nvariants += 1
    f.close()
    print("\tLoaded {:,} variants".format(nvariants))

    if nvariants == 0:
        print("Error, no variants in --variants file.")
        exit()

select_samples = None
if args.samples is not None:
    print("Loading samples...")
    select_samples = set()
    nsamples = 0
    with gzopen(args.samples, 'r') as f:
        for line in f:
            select_samples.add(line.rstrip("\n"))
            nsamples += 1
    f.close()
    print("\tLoaded {:,} samples".format(nsamples))

    if nsamples == 0:
        print("Error, no samples in --samples file.")
        exit()


print("Filtering VCF file...")
fhin = gzopen(args.vcf, 'r')
fhout = gzopen(args.outfile, 'w')

i = 0
n_qc_exclude = 0
n_subset_exclude = 0
n_written = 0
indices = []
vcf_header_columns = {"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"}
pos = {}
for i, line in enumerate(fhin):
    if (i == 0) or ((i % 5e5) == 0):
        print("\tParsed {:,} lines    QC exclude: {:,}    Subset exclude: {:,}    Saved: {:,}".format(i, n_qc_exclude, n_subset_exclude, n_written), end='\r', flush=True)

    if line.startswith("##"):
        fhout.write(line)
        continue

    values = line.rstrip("\n").split("\t")

    if line.startswith("#"):
        for index, value in enumerate(values):
            pos[value] = index
            if (select_samples is None) or (value in vcf_header_columns) or (value in select_samples):
                indices.append(index)
                continue

        fhout.write("##filter_vcf_variants.py --vcf {} --variants {} --samples {} --outfile {}\n".format(args.vcf, args.variants, args.samples, args.outfile))
        fhout.write("\t".join([values[index] for index in indices]) + "\n")
        continue

    if pass_qc_variants is None and args.chr is not None:
        chr = values[pos["#CHROM"]]
        if chr != args.chr:
            n_qc_exclude += 1
            continue

    variant = values[pos["ID"]]
    if pass_qc_variants is not None and variant not in pass_qc_variants:
        n_qc_exclude += 1
        continue

    if select_variants is not None and variant not in select_variants:
        n_subset_exclude += 1
        continue

    fhout.write("\t".join([values[index] for index in indices]) + "\n")
    n_written += 1

fhin.close()
fhout.close()
print("\tParsed {:,} lines    QC exclude: {:,}    Subset exclude: {:,}    Saved: {:,}".format(i, n_qc_exclude, n_subset_exclude, n_written), flush=True)

if n_written == 0:
    print("Error, no variants written to output file.")
    exit()

print("\nDone")
