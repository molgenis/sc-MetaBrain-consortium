#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import gzip

parser = argparse.ArgumentParser(description="")
parser.add_argument("--gte", required=True, type=str, default=None, help="")
parser.add_argument("--vcf", required=True, type=str, default=None, help="")
parser.add_argument("--exp", required=True, type=str, default=None, help="")
parser.add_argument("--cov", required=False, type=str, default=None, help="")
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


def get_individuals_from_vcf(fpath):
    vcf_columns = {"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"}
    header = None
    fh = gzopen(fpath, mode='r')
    for line in fh:
        if line.startswith("#CHROM"):
            header = line
            break
    fh.close()
    columns = set(header.strip("\n").split("\t"))
    return columns.symmetric_difference(vcf_columns)

def get_individuals_from_matrix(fpath):
    fh = gzopen(fpath, mode='r')
    header = fh.readline()
    fh.close()
    return set(header.rstrip("\n").split("\t")[1:])

def load_gte(fpath):
    geno_samples = set()
    expr_samples = set()
    fh = gzopen(fpath, mode='r')
    for line in fh:
        values = line.rstrip("\n").split("\t")
        geno_samples.add(values[0])
        expr_samples.add(values[1])
    fh.close()
    return geno_samples, expr_samples

def load_cov(fpath, ref_expr_samples):
    # This can be GTE or matrix.
    _, cov_gte_expr_samples = load_gte(fpath=fpath)
    cov_matrix_expr_samples = get_individuals_from_matrix(fpath=fpath)

    # We assume it is the one with the highest overlap.
    cov_gte_expr_overlap = cov_gte_expr_samples.intersection(ref_expr_samples)
    cov_matrix_expr_overlap = cov_matrix_expr_samples.intersection(ref_expr_samples)
    if len(cov_gte_expr_overlap) > len(cov_matrix_expr_overlap):
        return cov_gte_expr_samples
    elif len(cov_gte_expr_overlap) < len(cov_matrix_expr_overlap):
        return cov_matrix_expr_overlap
    else:
        print("Error, unable to determine if cov is a GTE file or a matrix.")
        exit()


##################################################

print("Loading genotype-to-expression coupling (--gte):")
geno_samples, expr_samples = load_gte(fpath=args.gte)
print("\t{:,} genotype samples".format(len(geno_samples)))
print("\t{:,} expression samples".format(len(expr_samples)))

print("\nLoading genotype input (--vcf):")
vcf_geno_samples = get_individuals_from_vcf(fpath=args.vcf)
print("\t{:,} genotype samples".format(len(vcf_geno_samples)))

print("\nLoading expression matrix (--exp):")
exp_expr_samples = get_individuals_from_matrix(fpath=args.exp)
print("\t{:,} expression samples".format(len(exp_expr_samples)))

cov_expr_samples = expr_samples
if args.cov is not None:
    print("\nLoading covariates (--cov):")
    cov_expr_samples = load_cov(fpath=args.cov, ref_expr_samples=expr_samples)
    print("\t{:,} expression samples".format(len(cov_expr_samples)))

print("\nOverlap:")
geno_overlap = geno_samples.intersection(vcf_geno_samples)
expr_overlap = expr_samples.intersection(exp_expr_samples).intersection(cov_expr_samples)
print("\t{:,} genotype samples".format(len(geno_overlap)))
print("\t{:,} expression samples".format(len(expr_overlap)))

print("\nDone")