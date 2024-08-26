#!/usr/bin/env python
# Author: M. Vochteloo
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--data", required=True, type=str, help="")
parser.add_argument("--nom_pvalue_column", required=False, type=str, default="pvalue", help="")
parser.add_argument("--perm_pvalue_column", required=False, type=str, default="perm-pvalue", help="")
parser.add_argument("--qvalue_column", required=False, type=str, default="qvalue", help="")
parser.add_argument("--minimimal_reporting_p", required=False, type=float, default=0.05, help="")
parser.add_argument("--out", required=True, type=str, help="The output file where results will be saved.")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import numpy as np
import pandas as pd
import gzip
import glob
import os

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


print("Counting number of eQTLs\n")
data = []
for fpath in glob.glob(args.data):
    cov = fpath.split(os.sep)[-2]
    n = 0
    n_nom = 0
    n_perm = 0
    n_qval = 0
    indices = {}
    with gzopen(fpath, 'r') as f:
        for i, line in enumerate(f):
            values = line.rstrip("\n").split("\t")
            if i == 0:
                indices = {colname: index for index, colname in enumerate(values)}
                if args.nom_pvalue_column not in indices:
                    n_nom = np.nan
                if args.perm_pvalue_column not in indices:
                    n_perm = np.nan
                if args.qvalue_column not in indices:
                    n_qval = np.nan
                continue

            if args.nom_pvalue_column in indices and float(values[indices[args.nom_pvalue_column]]) < args.minimimal_reporting_p:
                n_nom += 1
            if args.perm_pvalue_column in indices and float(values[indices[args.perm_pvalue_column]]) < args.minimimal_reporting_p:
                n_perm += 1
            if args.qvalue_column in indices and float(values[indices[args.qvalue_column]]) < args.minimimal_reporting_p:
                n_qval += 1
            n += 1
    f.close()

    data.append([cov, n, n_nom, n_perm, n_qval])

colnames = ["Cov", "N-effects", args.nom_pvalue_column, args.perm_pvalue_column, args.qvalue_column]
df = pd.DataFrame(data, columns=colnames).sort_values(by=[args.nom_pvalue_column, "N-effects", "Cov"], ascending=[False, True, True])
print(df)

print("Writing output")
df.to_csv(args.out, sep="\t", header=True, index=False)

print("\nEND")