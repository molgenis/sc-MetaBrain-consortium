#!/usr/bin/env Rscript

# Title     : cellranger_to_seurat.R
# Objective :
# Created by: mvochteloo
# Created on: 2023/03/20

# install.packages('Seurat')
# install.packages('Matrix')
library(Seurat)
library(Matrix)

args <- commandArgs(trailingOnly=TRUE)
work_dir <- args[1]
out_dir <- args[2]

list.files(work_dir) # Should show barcodes.tsv.gz, genes.tsv.gz, and matrix.mtx
expression_matrix <- Read10X(data.dir = work_dir)

meta_data <- read.csv(paste0(work_dir, '/metadata.csv.gz'), header = T, stringsAsFactors = F)
rownames(meta_data) <- meta_data$CellID

seurat_object <- CreateSeuratObject(counts = expression_matrix, meta.data = meta_data)
saveRDS(seurat_object, paste0(out_dir, '/seurat_object.rds'))