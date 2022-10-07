# Title     : compare_Mathys2019_CellRanger_output.R
# Objective :
# Created by: mvochteloo
# Created on: 27/09/2022

# install.packages('Matrix')
# install.packages('Seurat')
library(Seurat)
library(Matrix)

### Mathys data download; CellRanger v2.0.0 ###

cellranger_output_expression_matrix <- Read10X(data.dir = '/groups/umcg-biogen/tmp01/input/rawdata/single-cell/Mathys2019/processed/CellRangerOutput')
cellranger_output_seurat_object <- Seurat::CreateSeuratObject(counts = cellranger_output_expression_matrix, project = "cellranger_output")
cellranger_output_seurat_object
# 32,738 genes across 35,389,440 cells

notfiltered_expression_matrix <- Read10X(data.dir = '/groups/umcg-biogen/tmp01/input/rawdata/single-cell/Mathys2019/processed/notfiltered_feature_bc_matrix')
notfiltered_seurat_object <- Seurat::CreateSeuratObject(counts = notfiltered_expression_matrix, project = "notfiltered")
notfiltered_seurat_object
# 18,192 genes across 80,660 cells

filtered_expression_matrix <- Read10X(data.dir = '/groups/umcg-biogen/tmp01/input/rawdata/single-cell/Mathys2019/processed/filtered_feature_bc_matrix')
filtered_seurat_object <- Seurat::CreateSeuratObject(counts = filtered_expression_matrix, project = "filtered")
filtered_seurat_object
# 17,926 genes across 70,634 cells

### Mathys data re-processing; CellRanger v7.0.1 ###

my_expression_matrix <- Read10X(data.dir = '/groups/umcg-biogen/tmp01/input/processeddata/single-cell/Mathys2019/2022-09-27-CellRanger_Aggregate/outs/count/filtered_feature_bc_matrix/')
my_seurat_object <- Seurat::CreateSeuratObject(counts = my_expression_matrix, project = "My")
my_seurat_object
# 36,601 genes across 70,494 cells

filtered_seurat_object <- RunPCA(filtered_seurat_object)
my_seurat_object <- RunPCA(my_seurat_object)

merged <- merge(filtered_seurat_object, my_seurat_object, add.cell.ids = c("MathysFiltered", "My"), project = "merged_seurat")