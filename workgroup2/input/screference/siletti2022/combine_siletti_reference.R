#!/usr/bin/env Rscript

# Title     : combine_siletti_reference.R
# Objective :
# Created by: mvochteloo
# Created on: 2023/03/20

# install.packages('Seurat')
# install.packages('Matrix')
library(Seurat)
library(Matrix)

args <- commandArgs(trailingOnly=TRUE)
work_dir <- args[1]

nonneurons <- readRDS(paste0(work_dir, '/nonneurons/seurat_object.rds'))
nonneurons <- UpdateSeuratObject(nonneurons)
neurons <- readRDS(paste0(work_dir, '/neurons/seurat_object.rds'))
neurons <- UpdateSeuratObject(neurons)

reference <- merge(nonneurons, neurons)
reference <- SCTransform(reference)
set.seed(7777)
reference <- RunPCA(reference)
reference <- RunUMAP(reference, dims = 1:30, return.model = T)
reference <- FindNeighbors(reference, dims = 1:30, k.param = 20)
reference <- FindClusters(reference, resolution = 1.2)
saveRDS(reference, paste0(work_dir, '/siletti.rds'))