#!/usr/bin/env python3

import scanpy as sc
from scipy import io
import gzip
import sys

adata = sc.read_h5ad(sys.argv[1])
out_dir = sys.argv[2]

with gzip.open(out_dir + '/barcodes.tsv.gz', 'wb') as f:
    for item in adata.obs_names:
        f.write('{}\n'.format(item).encode())

with gzip.open(out_dir + '/features.tsv.gz', 'wb') as f:
    for item in ['\t'.join([x, x, 'Gene Expression']) for x in adata.var_names]:
        f.write('{}\n'.format(item).encode())

adata.obs.to_csv(out_dir + '/metadata.csv.gz')

io.mmwrite(out_dir + '/matrix.mtx', adata.X.T)
