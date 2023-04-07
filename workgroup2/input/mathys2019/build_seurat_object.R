#!/usr/bin/env Rscript

# Title     : build_seurat_object.R
# Objective :
# Created by: mvochteloo
# Created on: 2023/04/05

# install.packages('Seurat')
# install.packages('Matrix')
library(Seurat)
# library(Matrix)
# Attaching SeuratObject
library(Matrix)

work_dir <- "/groups/umcg-biogen/tmp01/input/processeddata/single-cell/datasets/Mathys2019/2022-09-26-CellRanger/Aggregate/outs"
expression_dir <- paste0(work_dir, "/count/filtered_feature_bc_matrix")
expression_matrix <- Read10X(data.dir = expression_dir)

barcodes <- read.table(paste0(expression_dir, "/barcodes.tsv.gz"), header = F, stringsAsFactors = F)
barcodes <- cbind(barcodes, data.frame(do.call('rbind', strsplit(as.character(barcodes$V1), '-', fixed = T))))
colnames(barcodes) <- c("TAG", "barcode", "sample.id")

aggregation <- read.table(paste0(work_dir, "/aggregation.csv"), sep=",", header = T, stringsAsFactors = F)
aggregation$sample.id <- 1:nrow(aggregation)
aggregation$molecule_h5 <- NULL
colnames(aggregation) <- c("projid", "sample.id")
barcodes <- merge(x = barcodes,
                  y = aggregation,
                  by = "sample.id",
                  all.x = TRUE,
                  sort = F)
barcodes$sample.id <- NULL

meta_data <- read.table("/groups/umcg-biogen/tmp01/input/ROSMAP-scRNAseq/filtered_column_metadata.txt.gz", header = T, stringsAsFactors = F)
meta_data$barcode <- data.frame(do.call('rbind', strsplit(as.character(meta_data$TAG), '.', fixed = T)))[, 1]
meta_data$TAG <- NULL
meta_data <- merge(x = barcodes,
                   y = meta_data,
                   by = c("barcode", "projid"),
                   all.x = TRUE,
                   sort = F)
rownames(meta_data) <- meta_data$TAG
meta_data$barcode <- NULL
meta_data$TAG <- NULL

seurat_object <- CreateSeuratObject(counts = expression_matrix, meta.data = meta_data)

seurat_object <- FindVariableFeatures(seurat_object)
# Calculating gene variances
# Calculating feature variances of standardized and clipped values
seurat_object <- SCTransform(seurat_object)
# Calculating cell attributes from input UMI matrix: log_umi
# Variance stabilizing transformation of count matrix of size 30542 by 70494
# Model formula is y ~ log_umi
# Get Negative Binomial regression parameters per gene
# Using 2000 genes, 5000 cells
# Found 73 outliers - those will be ignored in fitting/regularization step
# Second step: Get residuals using fitted parameters for 30542 genes
# Computing corrected count matrix for 30542 genes
# Calculating gene attributes
# Wall clock passed: Time difference of 15.0849 mins
# Determine variable features
# Place corrected count matrix in counts slot
# Centering data matrix
# Set default assay to SCT
set.seed(7777)
seurat_object <- RunPCA(seurat_object, verbose = FALSE)
seurat_object <- RunUMAP(seurat_object, reduction = "pca", dims = 1:50, return.model = T)
# Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
# To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
# This message will be shown once per session
# UMAP will return its model
# 14:12:31 UMAP embedding parameters a = 0.9922 b = 1.112
# 14:12:31 Read 70494 rows and found 50 numeric columns
# 14:12:31 Using Annoy for neighbor search, n_neighbors = 30
# 14:12:31 Building Annoy index with metric = cosine, n_trees = 50
# 14:12:41 Writing NN index file to temp file /local/1652931//RtmpXkWA3W/file609971d69f95
# 14:12:41 Searching Annoy index using 1 thread, search_k = 3000
# 14:13:06 Annoy recall = 100%
# 14:13:11 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
# 14:13:15 Initializing from normalized Laplacian + noise (using irlba)
# 14:13:31 Commencing optimization for 200 epochs, with 3148470 positive edges
# Using method 'umap'
# 14:15:07 Optimization finished
seurat_object <- FindNeighbors(seurat_object, dims = 1:50, k.param = 20)
# Computing nearest neighbor graph
# Computing SNN
seurat_object <- FindClusters(seurat_object, resolution = 1.2)
# Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#
# Number of nodes: 70494
# Number of edges: 2543550
#
# Running Louvain algorithm...
# Maximum modularity in 10 random starts: 0.9064
# Number of communities: 38
# Elapsed time: 34 seconds

saveRDS(seurat_object, paste0("/groups/umcg-biogen/tmp01/input/processeddata/single-cell/datasets/Mathys2019/2023-04-05-SeuratObject/Mathys2019.rds"))