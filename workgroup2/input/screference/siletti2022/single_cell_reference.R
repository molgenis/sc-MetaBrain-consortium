#!/usr/bin/env Rscript

# Title     : single_cell_reference.R
# Objective :
# Created by: mvochteloo, adapted from scRNAsequest - scRef.R by Zhengyu Ouyang & Yu Sun: https://github.com/interactivereport/scRNAsequest/blob/bade84eb821617dad8a411bc2e43cf5a96910e27/src/scRef.R
# Created on: 2023/03/27

# install.packages('Seurat')
# install.packages("bit64")
# install.packages('rhdf5')
# install.packages('Matrix')
# if (!requireNamespace('remotes', quietly = TRUE) {
#   install.packages('remotes')
# }
# remotes::install_github('satijalab/azimuth', ref = 'master')
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("glmGamPoi")

library(Seurat)
library(bit64)
library(rhdf5)
library(Matrix)
library(Azimuth)
options(future.globals.maxsize=3145728000, stringsAsFactors=F)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]
# work_dir <- '/mnt/depts/dept04/compbio/human_genetics/metabrain_sc/input/processeddata/Siletti2022'

# X      H5I_GROUP
# layers H5I_GROUP  named list()
# obs    H5I_GROUP  BadCells / CellID / Paris10 / Paris21 / Paris9 / PunchcardClusters / ROIGroup / ROIGroupCoarse / ROIGroupFine / Split / TopSplit / assay_ontology_term_id / cell_cycle_score / cluster_id / development_stage_ontology_term_id / disease_ontology_term_id / dissection / donor_id / fraction_mitochondrial / fraction_unspliced / organism_ontology_term_id / roi / sample_id / self_reported_ethnicity_ontology_term_id / sex_ontology_term_id / supercluster_term / suspension_type / tissue_ontology_term_id / total_UMIs / total_genes
# obsm   H5I_GROUP  X_UMAP / X_tSNE
# obsp   H5I_GROUP  named list()
# uns    H5I_GROUP  batch_condition (sample_id) /  batch_condition (3.0.0) / title (Neurons)
# var    H5I_GROUP  Accession / Biotype / Chromosome / End / Gene / Start
# varm   H5I_GROUP  named list()
# varp   H5I_GROUP  named list()

getX <- function(strH5ad){
  message("\tobtaining X ...")
  X <- h5read(strH5ad, "X", bit64conversion='bit64')
  cID <- h5read(strH5ad, "obs/CellID")
  gID <- h5read(strH5ad, "var/Accession")
  M <- sparseMatrix(i=X$indices,
                    p=X$indptr,
                    x=as.numeric(X$data),
                    dims=c(length(gID), length(cID)),
                    dimnames=list(gID, cID))
  return(M)
}

getobs <- function(strH5ad){
  message("\tobtaining obs ...")
  obs <- h5read(strH5ad, "obs")
  cID <- h5read(strH5ad, "obs/CellID")
  names <- names(obs)[names(obs) != "CellID"]
  meta <- data.frame(matrix(ncol = length(names(obs)), nrow = length(cID)))
  dimnames(meta) <- list(cID, names(obs))
  for(one in names(obs)){
    if ("categories" %in% names(obs[[one]])) {
      meta[, one] <- obs[[one]][["categories"]][obs[[one]][["codes"]] + 1]
    } else {
      meta[, one] <- obs[[one]]
    }
  }
  return(meta)
}

nonneurons <- paste0(work_dir, '/raw/nonneurons.h5ad')
neurons <- paste0(work_dir, '/raw/neurons.h5ad')


message("Building seurat object ...")
nonneurons_object <- CreateSeuratObject(counts = getX(nonneurons),
                                        project = "Non-neurons",
                                        meta.data = getobs(nonneurons))
neurons_object <- CreateSeuratObject(counts = getX(neurons),
                                     project = "Neurons",
                                     meta.data = getobs(neurons))
seurat_object <- merge(nonneurons, y = neurons, add.cell.ids = c("Non-neurons", "Neurons"), project = "Siletti2022")
rm(nonneurons_object, neurons_object)
gc()

message("SCTransform ...")
# no batches so we can just do everything at once
seurat_object <- SCTransform(seurat_object,
                             method = 'glmGamPoi', # poisson
                             new.assay.name = "SCT",
                             return.only.var.genes = FALSE, # TRUE
                             verbose = FALSE) # NULL
DefaultAssay(seurat_object) <- "SCT"

message("Processing ...")
# set.seed(7777)
seurat_object <- RunPCA(seurat_object)

#https://github.com/satijalab/azimuth/wiki/Azimuth-Reference-Format
# seurat_object <- unifySCTmodel(seurat_object)
seurat_object <- FindNeighbors(seurat_object, dims = 1:50, reduction = "pca", verbose = FALSE) # dims = 1:30
seurat_object <- RunSPCA(seurat_object, npcs = ncol(seurat_object[["pca"]]), graph = 'SCT_snn')

seurat_object <- FindNeighbors(seurat_object, dims = 1:50, reduction = "spca", verbose = FALSE)
seurat_object <- RunUMAP(seurat_object, dims = 1:50, reduction = "spca", umap.method = "uwot", return.model = TRUE) # dims = 1:30

message("Creating Azimuth reference ...")
reference <- AzimuthReference(seurat_object,
                              refUMAP = "umap",
                              refDR = "spca",
                              refAssay = "SCT",
                              dims = 1:ncol(seurat_object[["pca"]]),
                              plotref = "umap",
                              metadata = c('BadCells', 'ClustersPremerge', 'Paris10', 'Paris21', 'Paris25', 'Paris9', 'PunchcardClusters', 'ROIGroup', 'ROIGroupCoarse', 'ROIGroupFine', 'Split', 'TopSplit', 'roi', 'organism_ontology_term_id', 'disease_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'assay_ontology_term_id', 'sex_ontology_term_id', 'development_stage_ontology_term_id', 'donor_id', 'suspension_type', 'dissection', 'fraction_mitochondrial', 'fraction_unspliced', 'cell_cycle_score', 'total_genes', 'total_UMIs', 'sample_id', 'cluster_id', 'supercluster_term', 'cell_type_ontology_term_id', 'tissue_ontology_term_id')
)
saveRDS(reference, paste0(work_dir, 'siletti2022.rds'))
