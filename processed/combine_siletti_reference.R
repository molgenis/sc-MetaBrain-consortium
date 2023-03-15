# Title     : combine_siletti_reference.R
# Objective :
# Created by: mvochteloo
# Created on: 2023/03/15

# install.packages('Seurat')
# install.packages('SeuratDisk')
library(Seurat)
library(SeuratDisk)

Convert("h5ad/nonneurons.h5ad", "h5Seurat/nonneurons.h5Seurat")
nonneurons <- LoadH5Seurat("h5Seurat/nonneurons.h5Seurat")
Convert("h5ad/neurons.h5ad", "h5Seurat/neurons.h5Seurat")
neurons <- LoadH5Seurat("h5Seurat/neurons.h5Seurat")

reference <- merge(nonneurons, neurons)
reference <- SCTransform(reference)
set.seed(7777)
reference <- RunPCA(reference)
reference <- RunUMAP(reference, dims = 1:30, return.model = T)
reference <- FindNeighbors(reference, dims = 1:30, k.param = 20)
reference <- FindClusters(reference, resolution = 1.2)
saveRDS(reference, "siletti.rds")