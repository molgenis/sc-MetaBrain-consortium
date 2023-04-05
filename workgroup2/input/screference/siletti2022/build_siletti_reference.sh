#!/bin/bash

# ml HDF5/1.12.2-gompi-2022a
ml HDF5/1.12.1-gompi-2021b
# ml Python/3.10.4-GCCcore-11.3.0
ml Python/3.7.10-GCCcore-11.2.0
ml R
# Seurat_4.3.0
# Matrix_1.5-3

# source /groups/umcg-biogen/tmp01/umcg-mvochteloo/sc_env_3.10.4/bin/activate
source /mnt/depts/dept04/compbio/human_genetics/metabrain_sc/mvochtel/env/bin/activate
# anndata==0.8.0

# cd /groups/umcg-biogen/tmp01/input/processeddata/single-cell/screference/Siletti2022 || exit
cd /mnt/depts/dept04/compbio/human_genetics/metabrain_sc/input/processeddata/Siletti2022 || exit

mkdir raw
curl -o raw/nonneurons.h5ad https://storage.googleapis.com/linnarsson-lab-human/Nonneurons.h5ad
curl -o raw/neurons.h5ad https://storage.googleapis.com/linnarsson-lab-human/Neurons.h5ad

# Rscript single_cell_reference.R /groups/umcg-biogen/tmp01/input/processeddata/single-cell/screference/Siletti2022
Rscript single_cell_reference.R /mnt/depts/dept04/compbio/human_genetics/metabrain_sc/input/processeddata/Siletti2022



