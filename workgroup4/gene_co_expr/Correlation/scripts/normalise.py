#!/usr/bin/env python
# Author: M. Vochteloo and A. Kooijmans

import argparse
import gzip
import json
import numpy as np
import scipy as sp
import os
import scanpy
import h5py

parser = argparse.ArgumentParser(description="")
parser.add_argument("--counts", required=True, type=str, help="")
parser.add_argument("--avg_rd", required=True, type=str, help="")
parser.add_argument("--log1p", action="store_true", default=False, help="")
parser.add_argument("--weights_out", required=False, type=str, default=None, help="Path to the output weights file.")
parser.add_argument("--out", required=True, type=str, help="")

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

def do_pf(mtx, sf):
    pf = mtx.sum(axis=1).A.ravel()
    pf = sp.sparse.diags(sf/pf) @ mtx
    return pf

def save_filtered_counts_h5(fpath, adata):
    """
    Write hdf5 file from Cell Ranger v4 or later versions.
    Source: https://github.com/scverse/anndata/issues/595
    """
    if os.path.exists(fpath):
        raise FileExistsError(f"There already is a file `{fpath}`.")

    int_max = lambda x: int(max(np.floor(len(str(int(max(x)))) / 4), 1) * 4)
    str_max = lambda x: max([len(i) for i in x])

    # Save.
    w = h5py.File(fpath, 'w')
    grp = w.create_group("matrix")
    grp.create_dataset("barcodes", data=np.array(adata.obs_names, dtype=f'|S{str_max(adata.obs_names)}'))
    grp.create_dataset("data", data=np.array(adata.X.data, dtype=f'<i{int_max(adata.X.data)}'))
    ftrs = grp.create_group("features")
    # this group will lack the following keys:
    # '_all_tag_keys', 'feature_type', 'genome', 'id', 'name', 'pattern', 'read', 'sequence'
    ftrs.create_dataset("feature_type", data=np.array(adata.var.feature_types, dtype=f'|S{str_max(adata.var.feature_types)}'))
    ftrs.create_dataset("genome", data=np.array(adata.var.genome, dtype=f'|S{str_max(adata.var.genome)}'))
    ftrs.create_dataset("id", data=np.array(adata.var.gene_ids, dtype=f'|S{str_max(adata.var.gene_ids)}'))
    ftrs.create_dataset("name", data=np.array(adata.var.index, dtype=f'|S{str_max(adata.var.index)}'))
    grp.create_dataset("indices", data=np.array(adata.X.indices, dtype=f'<i{int_max(adata.X.indices)}'))
    grp.create_dataset("indptr", data=np.array(adata.X.indptr, dtype=f'<i{int_max(adata.X.indptr)}'))
    grp.create_dataset("shape", data=np.array(list(adata.X.shape)[::-1], dtype=f'<i{int_max(adata.X.shape)}'))


# def norm_pf_log_pf(mtx):
#     pf_log_pf = do_pf(np.log1p(do_pf(mtx)))
#     return pf_log_pf

print("Loading count ...")
adata = scanpy.read_10x_h5(args.counts)
mtx = adata.X
barcodes = adata.obs_names

print(f"Original matrix shape: {mtx.shape}")

print("Loading average read count ...")
fh = open(args.avg_rd, 'r')
sf = json.load(fh)["avg_read_count"]
fh.close()
print("\taverage read count: {:.4f}".format(sf))

print("Proportional fit Normalise ...")
mtx = do_pf(mtx=mtx, sf=sf)

if args.log1p:
    print("log1p transform ...")
    mtx = np.log1p(mtx)

if args.weights_out is not None:
    print("Calculating read counts ...")
    pf = mtx.A.sum(axis=1).ravel().astype(int)
    print(f"Matrix shape after log {mtx.shape}")

    with gzopen(args.weights_out, "w") as fh:
        for idx in range(0, len(barcodes)):
            fh.write(f"{barcodes[idx]}\t{pf[idx]}\n")
    fh.close()

print("Saving count ...")
adata.X = mtx
save_filtered_counts_h5(args.out, adata)

print("Done.")