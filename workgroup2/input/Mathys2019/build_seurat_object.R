#!/usr/bin/env Rscript

# Title     : build_seurat_object.R
# Objective :
# Created by: mvochteloo
# Created on: 2023/04/05

# install.packages('Seurat')
# install.packages('Matrix')
library(Seurat)
library(Matrix)

work_dir <- "/groups/umcg-biogen/tmp01/input/processeddata/single-cell/datasets/Mathys2019/2022-09-26-CellRanger/Aggregate/outs"
expression_dir <- paste0(work_dir, "/count/filtered_feature_bc_matrix")
list.files(expression_dir) # Should show barcodes.tsv.gz, genes.tsv.gz, and matrix.mtx
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
saveRDS(seurat_object, paste0("/groups/umcg-biogen/tmp01/input/processeddata/single-cell/datasets/Mathys2019/2023-04-05-SeuratObject/Mathys2019.rds"))