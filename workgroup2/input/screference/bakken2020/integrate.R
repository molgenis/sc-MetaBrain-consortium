#!/usr/bin/env Rscript

############################################################################################################################
# Authors: Martijn Vochteloo, based on code from Andrew Butler
# Source: https://github.com/satijalab/azimuth-references/blob/master/human_motorcortex/scripts/integrate.R
# Source: Based on code from Andrew Butler: https://github.com/satijalab/azimuth-references/blob/master/human_motorcortex/scripts/export.R
# Name: integrate.R
# Function: integrate the human motor cortex reference files.
############################################################################################################################

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

get_minor_to_major_subclass_dict <- function(){
  cell_type_dict <- list()
  # major cell types
  cell_type_dict[["Astro"]] <- "AST"
  cell_type_dict[["Endo"]] <- "END"
  cell_type_dict[["L2/3 IT"]] <- "EX"
  cell_type_dict[["L5 ET"]] <- "EX"
  cell_type_dict[["L5 IT"]] <- "EX"
  cell_type_dict[["L5/6 NP"]] <- "EX"
  cell_type_dict[["L6 CT"]] <- "EX"
  cell_type_dict[["L6 IT"]] <- "EX"
  cell_type_dict[["L6 IT Car3"]] <- "EX"
  cell_type_dict[["L6b"]] <- "EX"
  cell_type_dict[["Lamp5"]] <- "IN"
  cell_type_dict[["Micro-PVM"]] <- "MIC"
  cell_type_dict[["Oligo"]] <- "OLI"
  cell_type_dict[["OPC"]] <- "OPC"
  cell_type_dict[["Pvalb"]] <- "IN"
  cell_type_dict[["Sncg"]] <- "IN"
  cell_type_dict[["Sst"]] <- "IN"
  cell_type_dict[["Sst Chodl"]] <- "IN"
  cell_type_dict[["Vip"]] <- "IN"
  cell_type_dict[["VLMC"]] <- "PER"
  return(cell_type_dict)
}
cell_type_dict <- get_minor_to_major_subclass_dict()

Idents(object = mo.int) <- "subclass_label"
mo.int$class <- as.factor(mo.int$class_label)
mo.int$cluster <- as.factor(mo.int$cluster_label)
major_subclass <- cell_type_dict[mo.int$subclass_label]
names(major_subclass) <- names(mo.int$subclass_label)
mo.int$major_subclass <- as.factor(unlist(major_subclass))
mo.int$minor_subclass <- as.factor(mo.int$subclass_label)
mo.int$cross_species_cluster <- as.factor(mo.int$cross_species_cluster_label)

saveRDS(object = mo.int, file = args[3])