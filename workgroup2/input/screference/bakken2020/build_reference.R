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
appendix <- ""
if (downsample != "all") {
  ref <- subset(x = ref, downsample = strtoi(args[2]))
  appendix <- paste0("downsample", downsample)
}
ref <- RunUMAP(object = ref, reduction = "pca", dims = 1:50, return.model = TRUE)

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

SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]], file = paste0("reference/idx", appendix,".annoy"))
saveRDS(object = ref, file = paste0("reference/ref", appendix,".Rds"))