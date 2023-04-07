#!/usr/bin/env Rscript

############################################################################################################################
# Authors: Martijn Vochteloo, based on code from Andrew Butler
# Source: https://github.com/satijalab/azimuth-references/blob/master/human_motorcortex/scripts/export.R
# Name: build_reference.R
# Function: export the human motor cortex reference as a Azimuth reference.
############################################################################################################################

library(Seurat)
library(Azimuth)
args <- commandArgs(trailingOnly = TRUE)

ref <- readRDS(file = args[1])

downsample <- args[2]
if (downsample != "all") {
  paste0("Downsample to ", args[2], " cells")
  ref <- subset(x = ref, downsample = strtoi(args[2]))
}

ref <- SCTransform(ref)
# Calculating cell attributes from input UMI matrix: log_umi
# Variance stabilizing transformation of count matrix of size 32359 by 76621
# Model formula is y ~ log_umi
# Get Negative Binomial regression parameters per gene
# Using 2000 genes, 5000 cells
# Found 66 outliers - those will be ignored in fitting/regularization step
#
# Second step: Get residuals using fitted parameters for 32359 genes
# Computing corrected count matrix for 32359 genes
# Calculating gene attributes
# Wall clock passed: Time difference of 12.21147 mins
# Determine variable features
# Place corrected count matrix in counts slot
# Centering data matrix
# Set default assay to SCT
set.seed(7777)
ref <- RunPCA(ref, verbose = FALSE)
ref <- RunUMAP(ref, reduction = "pca", dims = 1:50, return.model = T)
# Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
# To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
# This message will be shown once per session
# UMAP will return its model
# 11:23:13 UMAP embedding parameters a = 0.9922 b = 1.112
# 11:23:13 Read 76621 rows and found 50 numeric columns
# 11:23:13 Using Annoy for neighbor search, n_neighbors = 30
# 11:23:13 Building Annoy index with metric = cosine, n_trees = 50
# 11:23:18 Writing NN index file to temp file /tmp/Rtmp3H8V4Y/file3a41d29d3d71a
# 11:23:18 Searching Annoy index using 1 thread, search_k = 3000
# 11:23:38 Annoy recall = 100%
# 11:23:43 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
# 11:23:45 Initializing from normalized Laplacian + noise (using irlba)
# 11:24:03 Commencing optimization for 200 epochs, with 3195436 positive edges
# 11:25:29 Optimization finished
ref <- FindNeighbors(ref, dims = 1:50, k.param = 20)
# Computing nearest neighbor graph
# Computing SNN
ref <- FindClusters(ref, resolution = 1.2)
# Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#
# Number of nodes: 76621
# Number of edges: 2723456
#
# Running Louvain algorithm...
# Maximum modularity in 10 random starts: 0.9332
# Number of communities: 51
# Elapsed time: 17 seconds

# SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]], file = paste0(args[3], "/idx.annoy"))
saveRDS(object = ref, file = paste0(args[3], "/ref.Rds"))
