#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
from statsmodels.regression.linear_model import OLS
import numpy as np
import gzip
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--data", required=True, type=str, help="")
parser.add_argument("--cov", required=True, type=str, help="")
parser.add_argument("--header", required=False, default=0, help="")
parser.add_argument("--index_col", required=False, default=0, help="")
parser.add_argument("--no_index_name", dest="has_index_name", action="store_false", default=True, help="")
parser.add_argument("--allow_na", action="store_true", default=False, help="")
parser.add_argument("--ignore_categorical", action="store_true", default=False, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

os.makedirs(os.path.dirname(args.out), exist_ok=True)

sep = "\t"

if args.header is not None and args.header != 0:
    print("Error, header should be None or 0. Other values might show unintended behaviour.")
    exit()


def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


def load_file(inpath, header_only=False):
    fhin = gzopen(inpath, mode="r")

    header_index = args.header
    if header_only and header_index is None:
        header_index = 0

    colnames = []
    rownames = []
    data = []
    n_values = None
    for i, line in enumerate(fhin):
        values = line.rstrip("\n").split(sep)
        if n_values is None:
            n_values = len(values)
        if len(values) != n_values:
            print("  Error, unequal number of columns in the input file.")
            fhin.close()
            exit()

        if args.allow_na:
            values = ["nan" if value == '' else value for value in values]

        rowname = i
        if args.index_col is not None:
            rowname = values[args.index_col]
            del values[args.index_col]
        rownames.append(rowname)

        if i == header_index:
            colnames = values
            if header_only:
                break
            continue

        data.append(values)

    fhin.close()

    if not colnames:
        colnames = [i for i, _ in enumerate(data[0])]

    index_name = ""
    if args.has_index_name:
        index_name = rownames[0]
        del rownames[0]

    return colnames, index_name, rownames, np.array(data)


def calc_residuals(y, X):
    y_hat = np.empty(y.shape, dtype=np.float64)

    # Mask the NaN values.
    obs_mask = np.ones(y_hat.shape, dtype=bool)
    if args.allow_na:
        obs_mask = (~np.isnan(y)) & (np.isnan(X).sum(axis=1) == 0)
        y_hat[~obs_mask] = np.nan

    # Calculate residuals using OLS.
    y_hat[obs_mask] = OLS(y[obs_mask], X[obs_mask, :]).fit().resid

    # Cast back to list of strings.
    return y_hat.astype(str).tolist()


def regress_data(inpath, outpath, covariates, sample_mask):
    # TODO: this casting to array and back is a bit slow and not even required if sample_mask is all True.

    fhin = gzopen(inpath, mode="r")
    fhout = gzopen(outpath, mode="w")

    i = 0
    n_values = None
    for i, line in enumerate(fhin):
        if i % 1000 == 0:
            print("  Processed {:,} rows".format(i), end='\r')

        values = line.rstrip("\n").split(sep)
        if n_values is None:
            n_values = len(values)
        if len(values) != n_values:
            print("  Error, unequal number of columns in the input file.")
            fhin.close()
            exit()

        if args.allow_na:
            values = ["nan" if value == '' else value for value in values]

        index = i
        if args.index_col is not None:
            index = values[args.index_col]
            del values[args.index_col]

        if i == args.header:
            fhout.write(sep.join(["-"] + np.array(values, dtype=object)[sample_mask].tolist()) + "\n")
            continue

        adj_values = calc_residuals(y=np.array(values, dtype=np.float64)[sample_mask], X=covariates)
        fhout.write(sep.join([index] + adj_values) + "\n")

    fhin.close()
    fhout.close()

    print("  Processed {:,} rows\n".format(i))

# First load the covariates file completely and the header of the data file.
cov_colnames, _, cov_rownames, cov_m = load_file(args.cov)
print("Covariate file has {} rows and {} columns".format(len(cov_rownames), len(cov_colnames)))
if len(set(cov_colnames)) != len(cov_colnames):
    print("Error, covariate columns should have unique values")
    exit()
if len(set(cov_rownames)) != len(cov_rownames):
    print("Error, covariate indices should have unique values")
    exit()

data_colnames, _, _, _ = load_file(args.data, header_only=True)
if len(set(data_colnames)) != len(data_colnames):
    print("Error, data columns should have unique values")
    exit()

# Check if the samples overlap. If not, try transposing the covariates file. By default
# we expect the samples to be on the rows for the covariates and on the columns for the data.
samples = set(cov_rownames).intersection(set(data_colnames))
if len(samples) == 0:
    print("No matching samples detected between covariate file and dataset. Maybe your covariate file needs to be transposed? Will test that for you now:")

    samples = set(cov_colnames).intersection(set(data_colnames))
    print("Transposing the covariate file reveals: {} samples present.".format(len(samples)))
    if len(samples) == 0:
        exit()

    cov_colnames_tmp = cov_colnames.copy()
    cov_colnames = cov_rownames
    cov_rownames = cov_colnames_tmp
    cov_m = np.transpose(cov_m)
    del cov_colnames_tmp

print("Your covariate corrected dataset will have {} samples, after removing samples with missing covariate values.".format(len(samples)))

print("Casting all covariates to numeric")
cov_num_arrays = []
cov_num_colnames = []
for cov_index, cov_colname in enumerate(cov_colnames):
    try:
        cov_num_arrays.append(cov_m[:, cov_index].astype(np.float64))
        cov_num_colnames.append(cov_colname)
    except ValueError:
        if args.ignore_categorical:
            print("\t{} removed since ignore_categorical is True".format(cov_colname))
            continue

        n_unique = len(np.unique(cov_m[:, cov_index]))
        if n_unique == 1:
            print("\t{} removed due to no variance".format(cov_colname))
            continue
        if n_unique == cov_m.shape[0]:
            print("\t{} removed due to all unique values".format(cov_colname))
            continue

        print("\tCovariate {} is not numeric, one-hot encoding {:,} values ignoring most frequent value.".format(cov_colname, n_unique))

        # Count the unique categorical values.
        cat_cov_counts = list(zip(*np.unique(cov_m[:, cov_index], return_counts=True)))
        cat_cov_counts.sort(key=lambda x: -x[1])

        # We skip the most frequent one and on-hot encode the rest.
        for (cov_oh_colname, _) in cat_cov_counts[1:]:
            cov_oh_a = np.zeros(cov_m.shape[0], dtype=np.float64)
            mask = cov_m[:, cov_index] == cov_oh_colname
            cov_oh_a[mask] = 1

            # Save as normal.
            cov_num_arrays.append(cov_oh_a.astype(np.float64))
            cov_num_colnames.append(cov_colname + "_" + cov_oh_colname)
cov_m = np.transpose(np.vstack(cov_num_arrays))
cov_colnames = cov_num_colnames
print("\nUpdated covariates:\n{}\n".format(";".join(cov_colnames)))

print("Checking variance of covariates\n")
cov_colmask = np.std(cov_m, axis=0) != 0
cov_m = cov_m[:, cov_colmask]
cov_colnames = [colname for (colname, keep) in zip(cov_colnames, cov_colmask) if keep]
print("Remaining covariates:\n{}\n".format(";".join(cov_colnames)))

# Updating colorder to match data file. Better to sort this file so we
# do not need to sort the data file which is likely much bigger.
data_sample_mask = np.ones(len(data_colnames), dtype=bool)
if cov_rownames != data_colnames:
    print("Reordering covariate matrix based on data matrix columns.")
    cov_rownames_tmp = []
    cov_roworder = []
    for sample_index, sample in enumerate(data_colnames):
        if sample not in cov_rownames:
            data_sample_mask[sample_index] = False
            continue
        cov_rownames_tmp.append(sample)
        cov_roworder.append(cov_rownames.index(sample))
    cov_rownames = cov_rownames_tmp
    cov_m = cov_m[cov_roworder, :]
    del cov_rownames_tmp, cov_roworder

    if cov_rownames != list(np.array(data_colnames)[data_sample_mask]):
        print("Error, covariate samples do not match data samples.")
        exit()

# Adding the intercept as the first row.
cov_m = np.hstack([np.ones((cov_m.shape[0], 1)), cov_m])
cov_colnames = ["intercept"] + cov_colnames

# Define the output file.
out_fpath = args.out + ".CovariatesRemovedOLS.txt" + (".gz" if args.data.endswith(".gz") else "")

print("Calculating OLS residuals...")
regress_data(inpath=args.data, outpath=out_fpath, covariates=cov_m, sample_mask=data_sample_mask)

print("Saved: {}".format(out_fpath))
print("Done")
