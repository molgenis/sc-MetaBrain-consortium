#!/usr/bin/env Rscript

############################################################################################################################
# Authors: Martijn Vochteloo, based on code from Andrew Butler
# Source: https://github.com/satijalab/azimuth-references/blob/master/human_motorcortex/scripts/integrate.R
# Source: Based on code from Andrew Butler: https://github.com/satijalab/azimuth-references/blob/master/human_motorcortex/scripts/export.R
# Name: integrate.R
# Function: integrate the human motor cortex reference files.
############################################################################################################################

.libPaths("/usr/local/lib/R/site-library")
library(Seurat)
library(feather)
args <- commandArgs(trailingOnly = TRUE)

dat <- readRDS(file = args[1])
metadata <- as.data.frame(x = read_feather(path = args[2]))
rownames(x = metadata) <- metadata$sample_id
metadata <- metadata[colnames(x = dat), ]

ob <- CreateSeuratObject(counts = dat, meta.data = metadata)
ob.list <- SplitObject(object = ob, split.by = "donor_id")
ob.list <- lapply(X = ob.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = ob.list, nfeatures = 3000)
ob.list <- PrepSCTIntegration(object.list = ob.list, anchor.features = features)
anchors <- FindIntegrationAnchors(
  object.list = ob.list,
  anchor.features = features,
  normalization.method = "SCT",
  dims = 1:30
)
mo.int <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  dims = 1:30
)
mo.int <- RunPCA(mo.int, verbose = FALSE)

saveRDS(object = mo.int, file = args[3])