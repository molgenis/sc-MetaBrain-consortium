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
  ref <- subset(x = ref, downsample = strtoi(args[2]))
}

ref <- SCTransform(ref)
set.seed(7777)
ref <- RunPCA(ref, verbose = FALSE)
ref <- RunUMAP(ref, reduction = "pca", dims = 1:50, return.model = T)
# Only one graph name supplied, storing nearest-neighbor graph only
ref <- FindNeighbors(ref, dims = 1:50, k.param = 20)
ref <- FindClusters(ref, resolution = 1.2)

ref <- AzimuthReference(
  object = ref,
  refUMAP = "umap",
  refDR = "pca",
  refAssay = "integrated",
  metadata = c("class", "cluster", "major_subclass", "minor_subclass", "cross_species_cluster"),
  dims = 1:50,
  k.param = 31,
  reference.version = "1.0.0"
)

SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]], file = paste0(args[2], "idx.annoy"))
saveRDS(object = ref, file = paste0(args[2], "ref.Rds"))
