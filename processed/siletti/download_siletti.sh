#!/bin/bash

ml HDF5/1.12.2-gompi-2022a
ml Python/3.10.4-GCCcore-11.3.0
ml R
# Seurat_4.3.0
# SeuratObject_4.1.3
# Matrix_1.5-3

source /groups/umcg-biogen/tmp01/umcg-mvochteloo/sc_env_3.10.4/bin/activate
# anndata==0.8.0

cd /groups/umcg-biogen/tmp01/input/processeddata/single-cell/Siletti2022/ || exit

mkdir raw
curl -o raw/nonneurons.h5ad https://storage.googleapis.com/linnarsson-lab-human/Nonneurons.h5ad
curl -o raw/neurons.h5ad https://storage.googleapis.com/linnarsson-lab-human/Neurons.h5ad

mkdir -p cellranger/nonneurons
python3 export_h5ad_to_cellranger.py raw/nonneurons.h5ad cellranger/nonneurons
gzip cellranger/nonneurons/matrix.mtx

mkdir -p cellranger/neurons
python3 export_h5ad_to_cellranger.py raw/neurons.h5ad cellranger/neurons
gzip cellranger/neurons/matrix.mtx

mkdir -p seurat/nonneurons
Rscript cellranger_to_seurat.R cellranger/nonneurons seurat/nonneurons

mkdir -p seurat/neurons
Rscript cellranger_to_seurat.R cellranger/neurons seurat/neurons

Rscript combine_siletti_reference.R seurat



