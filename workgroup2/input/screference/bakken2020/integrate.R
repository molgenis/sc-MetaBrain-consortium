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
  cell_type_dict[["Astro"]] <- "astrocyte"
  cell_type_dict[["Endo"]] <- "endothelial cell"
  cell_type_dict[["L2/3 IT"]] <- "excitatory neuron"
  cell_type_dict[["L5 ET"]] <- "excitatory neuron"
  cell_type_dict[["L5 IT"]] <- "excitatory neuron"
  cell_type_dict[["L5/6 NP"]] <- "excitatory neuron"
  cell_type_dict[["L6 CT"]] <- "excitatory neuron"
  cell_type_dict[["L6 IT"]] <- "excitatory neuron"
  cell_type_dict[["L6 IT Car3"]] <- "excitatory neuron"
  cell_type_dict[["L6b"]] <- "excitatory neuron"
  cell_type_dict[["Lamp5"]] <- "inhibitory neuron"
  cell_type_dict[["Micro-PVM"]] <- "perivascular macrophage"
  cell_type_dict[["Oligo"]] <- "oligodendrocyte"
  cell_type_dict[["OPC"]] <- "oligodendrocyte precursor cell"
  cell_type_dict[["Pvalb"]] <- "inhibitory neuron"
  cell_type_dict[["Sncg"]] <- "inhibitory neuron"
  cell_type_dict[["Sst"]] <- "inhibitory neuron"
  cell_type_dict[["Sst Chodl"]] <- "inhibitory neuron"
  cell_type_dict[["Vip"]] <- "inhibitory neuron"
  cell_type_dict[["VLMC"]] <- "pericyte"
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