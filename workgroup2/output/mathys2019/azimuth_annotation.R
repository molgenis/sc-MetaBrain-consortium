#!/usr/bin/env Rscript

############################################################################################################################
# Authors: Martijn Vochteloo, based on code by Roy Oelen
# Name: azimuth_annotation.R
# Function: annotate Mathys 2019 dataset with Bakken 2020 as cell type reference
############################################################################################################################

library(Seurat)
library(Azimuth) # module load HDF5/1.12.2-gompi-2022a
library(SeuratData)
library(patchwork)
library(ggplot2)

# work_dir <- "/groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-02-02-WorkGroup2CellType/2023-04-05-Mathys2019"
# query_file <- "/groups/umcg-biogen/tmp01/input/processeddata/single-cell/datasets/Mathys2019/2023-04-05-SeuratObject/Mathys2019.rds"
# ref_dir <- "/groups/umcg-biogen/tmp01/input/processeddata/single-cell/screference/Bakken2020/reference"
# subset <- "all"

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]
query_file <- args[2]
ref_dir <- args[3]
subset <- args[4]

dir.create(file.path(work_dir, subset), showWarnings = FALSE)
dir.create(file.path(work_dir, subset, "data"), showWarnings = FALSE)
dir.create(file.path(work_dir, subset, "plots"), showWarnings = FALSE)
dir.create(file.path(work_dir, subset, "plots", "barplots"), showWarnings = FALSE)
dir.create(file.path(work_dir, subset, "plots", "dimplots"), showWarnings = FALSE)
dir.create(file.path(work_dir, subset, "plots", "prediction_scores"), showWarnings = FALSE)
dir.create(file.path(work_dir, subset, "plots", "prediction_scores", "feature_plots"), showWarnings = FALSE)
dir.create(file.path(work_dir, subset, "plots", "prediction_scores", "violin_plots"), showWarnings = FALSE)
dir.create(file.path(work_dir, subset, "plots", "tables"), showWarnings = FALSE)

query_broad_cell_type_colors <- c(
  "Ast" = "#D55E00",
  "End" = "#CC79A7",
  "Ex" = "#0072B2",
  "In" = "#56B4E9",
  "Mic" = "#E69F00",
  "Oli" = "#009E73",
  "Opc" = "#F0E442",
  "Per" = "#808080",
  "NA" = "#000000"
)
reference_major_subclass_colors <- c(
  "astrocyte" = "#D55E00",
  "endothelial cell" = "#CC79A7",
  "excitatory neuron" = "#0072B2",
  "inhibitory neuron" = "#56B4E9",
  "perivascular macrophage" = "#E69F00",
  "oligodendrocyte" = "#009E73",
  "oligodendrocyte precursor cell" = "#F0E442",
  "pericyte" = "#808080"
)

# read the query and the reference
query <- readRDS(query_file)
reference <- readRDS(paste0(ref_dir, "/", subset, "/ref.Rds"))

# find transfer anchors
anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:50
)
# Normalizing query using reference SCT model
# Projecting cell embeddings
# Finding neighborhoods
# Finding anchors
#         Found 14531 anchors

# do the reference mapping
query <- MapQuery(
  anchorset = anchors,
  query = query,
  reference = reference,
  refdata = list(
    cluster = "cluster",
    class = "class",
    major_subclass = "major_subclass",
    minor_subclass = "minor_subclass"
  ),
  reference.reduction = "pca",
  reduction.model = "umap"
)
# Finding integration vectors
# Finding integration vector weights
# Predicting cell labels
# Predicting cell labels
# Predicting cell labels
# Warning: Keys should be one or more alphanumeric characters followed by an underscore, setting key from predictionscoremajor_subclass_ to predictionscoremajorsubclass_
# Predicting cell labels
# Warning: Keys should be one or more alphanumeric characters followed by an underscore, setting key from predictionscoreminor_subclass_ to predictionscoreminorsubclass_
# Integrating dataset 2 with reference dataset
# Finding integration vectors
# Integrating data
# Computing nearest neighbors
# Running UMAP projection
# 18:19:01 Read 70494 rows
# 18:19:01 Processing block 1 of 1
# 18:19:01 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
# 18:19:02 Initializing by weighted average of neighbor coordinates using 1 thread
# 18:19:02 Commencing optimization for 67 epochs, with 2114812 positive edges
# Using method 'umap'
# 18:19:27 Finished

# save the result
saveRDS(query, paste0(work_dir, "/", subset, "/data/mathys2019_bakken2020training.rds"))
write.csv(query@meta.data, paste0(work_dir, "/", subset, "/data/seurat_metadata.csv"))

p1 <- DimPlot(object = reference,
              reduction = "umap",
              group.by = "cluster_label",
              label = FALSE,
              label.size = 3,
              repel = TRUE) + NoLegend()
p2 <- DimPlot(object = query,
              reduction = "ref.umap",
              group.by = "predicted.cluster",
              label = FALSE,
              label.size = 3,
              repel = TRUE) + NoLegend()
ggsave('mathys2019_bakken2020training_ref_and_projected.pdf',
       plot = p1 + p2,
       path = paste0(work_dir, "/", subset, "/plots/dimplots"),
       width = 20,
       height = 10)

p1 <- DimPlot(object = query,
              reduction = "umap",
              group.by = "seurat_clusters",
              label = TRUE,
              label.size = 3,
              repel = TRUE) + NoLegend()
p2 <- DimPlot(object = query,
              reduction = "umap",
              group.by = "predicted.minor_subclass",
              label = TRUE,
              label.size = 3,
              repel = TRUE)
ggsave('mathys2019_bakken2020training_predicted_vs_clusters.pdf',
       plot = p1 + p2,
       path = paste0(work_dir, "/", subset, "/plots/dimplots"),
       width = 20,
       height = 10)

p1 <- DimPlot(object = query,
              reduction = "umap",
              group.by = "predicted.major_subclass",
              cols = reference_major_subclass_colors,
              label = TRUE,
              label.size = 3,
              repel = TRUE)
p2 <- DimPlot(object = query,
              reduction = "umap",
              group.by = "broad.cell.type",
              cols = query_broad_cell_type_colors,
              label = TRUE,
              label.size = 3,
              repel = TRUE)
ggsave('mathys2019_bakken2020training_mathys_vs_azimuth.pdf',
       plot = p1 + p2,
       path = paste0(work_dir, "/", subset, "/plots/dimplots"),
       width = 20,
       height = 10)

# plot the scores
create_features_per_celltype <- function(seurat_object, celltype_column, features, plot_out_loc, reduction = 'umap', min.cutoff = NA, max.cutoff = NA){
  # do each cell type
  for(celltype in unique(seurat_object@meta.data[[celltype_column]])){
    # paste the output location together
    features_pasted <- paste(features, collapse = '_')
    full_file_loc <- paste(celltype, features_pasted, sep = '_')
    full_file_loc <- paste(full_file_loc, '.pdf', sep = '')
    full_file_loc <- gsub('/', '_', full_file_loc)
    # make the plot
    plot_to_save <- FeaturePlot(seurat_object[, !is.na(seurat_object@meta.data[[celltype_column]]) &
                                                seurat_object@meta.data[[celltype_column]] == celltype],
                                features = features,
                                reduction = reduction,
                                min.cutoff = min.cutoff,
                                max.cutoff = max.cutoff)
    # save the plot
    ggsave(full_file_loc,
           plot = plot_to_save,
           path = plot_out_loc,
           width = 10,
           height = 10)
  }
}
create_features_per_celltype(seurat_object = query,
                             celltype_column = "predicted.major_subclass",
                             features = c("predicted.major_subclass.score"),
                             plot_out_loc = paste0(work_dir, "/", subset, "/plots/prediction_scores/feature_plots"))
create_features_per_celltype(seurat_object = query,
                             celltype_column = "predicted.minor_subclass",
                             features = c("predicted.minor_subclass.score"),
                             plot_out_loc = paste0(work_dir, "/", subset, "/plots/prediction_scores/feature_plots"))

# as violins as well
create_violins <- function(seurat_object, split_column, plot_column, fill_colors = NULL){
  # create the table with data we need
  table_to_use <- seurat_object@meta.data[, c(split_column, plot_column)]
  violin_plot <- ggplot(data = NULL,
                        mapping=aes(x = table_to_use[[split_column]],
                                    y = table_to_use[[plot_column]],
                                    fill = table_to_use[[split_column]])
  ) +
    geom_violin() +
    geom_jitter(aes(colour = table_to_use[[split_column]]),
                width = 0.00005,
                alpha = 0.05) +
    ylab("prediction score") +
    xlab("cell type")
  if(!is.null(fill_colors)){
    paste0("bla")
    violin_plot <- violin_plot + scale_fill_manual(values = fill_colors) + scale_color_manual(values = fill_colors)
  }
  return(violin_plot)
}
scores_violins_major <- create_violins(seurat_object = query,
                                       split_column = 'predicted.major_subclass',
                                       plot_column = c('predicted.major_subclass.score'),
                                       fill_colors = reference_major_subclass_colors)
ggsave('mathys2019_bakken2020training_score_violin_major.pdf',
       plot = scores_violins_major,
       path = paste0(work_dir, "/", subset, "/plots/prediction_scores/violin_plots"),
       width = 20,
       height = 10)
scores_violins_minor <- create_violins(seurat_object = query,
                                       split_column = 'predicted.minor_subclass',
                                       plot_column = c('predicted.minor_subclass.score'))
ggsave('mathys2019_bakken2020training_score_violin_minor.pdf',
       plot = scores_violins_minor,
       path = paste0(work_dir, "/", subset, "/plots/prediction_scores/violin_plots"),
       width = 20,
       height = 10)

cell_types_to_barplot <- function(metadata_table, celltype_column, to_fraction = F, fill_colors = NULL, ylim = NULL, pointless = T, legendless = F){
  celltype_table <- data.frame(table(metadata_table[[celltype_column]]))
  celltype_table[['cell_type']] <- celltype_table[['Var1']]
  ylabel <- 'number of cells'
  xlabel <- 'cell type'
  if(to_fraction){
    celltype_table$Freq <- celltype_table$Freq / sum(celltype_table$Freq)
    ylabel <- 'fraction of cells'
  }
  p <- ggplot(data = celltype_table,
              aes(x = cell_type,
                  y = Freq,
                  fill = cell_type)) +
    geom_bar(stat = 'identity') +
    ylab(ylabel) +
    xlab(xlabel)
  if(!is.null(fill_colors)){
    p <- p + scale_fill_manual(values = fill_colors)
  }
  if(!is.null(ylim)){
    p <- p + ylim(ylim)
  }
  if(pointless){
    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks = element_blank())
  }
  if(legendless){
    p <- p + theme(legend.position = 'none')
  }
  return(p)
}
cell_types_to_barplot(metadata_table = query@meta.data,
                      celltype_column = 'predicted.major_subclass',
                      fill_colors = reference_major_subclass_colors)
ggsave(paste0(work_dir, "/", subset, '/plots/barplots/mathys2019_bakken2020training_predicted_cell_numbers.pdf'),
       width = 20,
       height = 10)
cell_types_to_barplot(metadata_table = query@meta.data,
                      celltype_column = 'predicted.major_subclass',
                      to_fraction = T,
                      fill_colors = reference_major_subclass_colors)
ggsave(paste0(work_dir, "/", subset, '/plots/barplots/mathys2019_bakken2020training_predicted_cell_fractions.pdf'),
       width = 20,
       height = 10)
cell_types_to_barplot(metadata_table = query@meta.data,
                      celltype_column = 'broad.cell.type',
                      fill_colors = query_broad_cell_type_colors)
ggsave(paste0(work_dir, "/", subset, '/plots/barplots/mathys2019_broad_cell_type_numbers.pdf'),
       width = 20,
       height = 10)
cell_types_to_barplot(metadata_table = query@meta.data,
                      celltype_column = 'broad.cell.type',
                      to_fraction = T,
                      fill_colors = query_broad_cell_type_colors)
ggsave(paste0(work_dir, "/", subset, '/plots/barplots/mathys2019_broad_cell_type_fractions.pdf'),
       width = 20,
       height = 10)

# get the overlap of clusters with a cell type
create_confusion_table <- function(assignment_table, truth_column, prediction_column){
  confusion_table <- NULL
  for(truth in unique(assignment_table[[truth_column]])){
    truth_rows <- assignment_table[assignment_table[[truth_column]] == truth, ]
    this_truth_number <- nrow(truth_rows)
    for(prediction in unique(assignment_table[[prediction_column]])){
      this_prediction_number <- nrow(truth_rows[truth_rows[[prediction_column]] == prediction, ])
      fraction <- NULL
      if (this_prediction_number > 0) {
        fraction <- this_prediction_number / this_truth_number
      } else{
        fraction <- 0
      }
      this_row <- data.frame(truth = c(truth),
                             prediction = c(prediction),
                             freq = c(fraction),
                             stringsAsFactors = F)
      if (is.null(confusion_table)) {
        confusion_table <- this_row
      } else {
        confusion_table <- rbind(confusion_table, this_row)
      }
    }
  }
  return(confusion_table)
}
clusters_celltype_overlap <- create_confusion_table(assignment_table = query@meta.data,
                                                    truth_column = 'broad.cell.type',
                                                    prediction_column = 'predicted.major_subclass')
clusters_celltype_overlap <- clusters_celltype_overlap[clusters_celltype_overlap[['freq']] != 0, ]
clusters_celltype_overlap <- clusters_celltype_overlap[order(clusters_celltype_overlap[['truth']], clusters_celltype_overlap[['freq']], decreasing = T), ]
colnames(clusters_celltype_overlap) <- c('broad.cell.type', 'predicted.major_subclass', 'fraction')
write.table(clusters_celltype_overlap,
            paste0(work_dir, "/", subset, "/tables/broad_cell_type_to_bakken_predic_overlap.tsv"),
            sep = '\t',
            row.names = F,
            col.names = T)

clusters_celltype_overlap <- create_confusion_table(assignment_table = query@meta.data,
                                                    truth_column = 'seurat_clusters',
                                                    prediction_column = 'predicted.major_subclass')
clusters_celltype_overlap <- clusters_celltype_overlap[clusters_celltype_overlap[['freq']] != 0, ]
clusters_celltype_overlap <- clusters_celltype_overlap[order(clusters_celltype_overlap[['truth']], clusters_celltype_overlap[['freq']], decreasing = T), ]
colnames(clusters_celltype_overlap) <- c('seurat_clusters', 'predicted.major_subclass', 'fraction')
write.table(clusters_celltype_overlap,
            paste0(work_dir, "/", subset, "/tables/clusters_to_bakken_predic_overlap.tsv"),
            sep = '\t',
            row.names = F,
            col.names = T)

# do a cell to cluster approach
add_imputed_meta_data <- function(seurat_object, column_to_transform, column_to_reference, column_to_create){
  seurat_object@meta.data[column_to_create] <- NA
  for(group in unique(seurat_object@meta.data[[column_to_transform]])){
    seurat_group <- seurat_object[, seurat_object@meta.data[[column_to_transform]] == group]
    best_group <- 'unknown'
    best_number <- 0
    for(reference in unique(seurat_group@meta.data[[column_to_reference]])){
      if(is.na(reference) == F){
        number_of_reference_in_group <- nrow(seurat_group@meta.data[seurat_group@meta.data[[column_to_reference]] == reference & is.na(seurat_group@meta.data[[column_to_reference]]) == F,])
        correctpercent <- number_of_reference_in_group/ncol(seurat_group)
        print(paste(group, "matches", reference, correctpercent, sep = " "))
        if(number_of_reference_in_group > best_number){
          best_number <- number_of_reference_in_group
          best_group <- reference
        }
      }
    }
    print(paste("setting ident:", best_group ,"for group", group, sep=" "))
    seurat_object@meta.data[seurat_object@meta.data[[column_to_transform]] == group,][column_to_create] <- best_group
    rm(seurat_group)
  }
  return(seurat_object)
}
query <- add_imputed_meta_data(seurat_object = query,
                               column_to_transform = 'seurat_clusters',
                               column_to_reference = 'predicted.major_subclass',
                               column_to_create = 'predicted.major_subclass.cttoclus')
clusters_cttoclus_overlap <- create_confusion_table(assignment_table = query@meta.data,
                                                    truth_column = 'predicted.major_subclass.cttoclus',
                                                    prediction_column = 'predicted.major_subclass')
clusters_cttoclus_overlap <- clusters_cttoclus_overlap[clusters_cttoclus_overlap[['freq']] != 0, ]
clusters_cttoclus_overlap <- clusters_cttoclus_overlap[order(clusters_cttoclus_overlap[['truth']], clusters_cttoclus_overlap[['freq']], decreasing = T), ]
write.table(clusters_cttoclus_overlap,
            paste0(work_dir, "/", subset, "/tables/clusters_to_bakken_predic_cttoclus_overlap.tsv"),
            sep = '\t',
            row.names = F,
            col.names = T)

query@meta.data[["broad.cell.type.nafilled"]] <- query@meta.data[["broad.cell.type"]]
query@meta.data[is.na(query@meta.data[["broad.cell.type.nafilled"]]), "broad.cell.type.nafilled"] <- "None"
query <- add_imputed_meta_data(seurat_object = query,
                               column_to_transform = 'broad.cell.type.nafilled',
                               column_to_reference = 'predicted.major_subclass',
                               column_to_create = 'predicted.major_subclass.cttoclus')
broad_cell_type_cttoclus_overlap <- create_confusion_table(assignment_table = query@meta.data,
                                                           truth_column = 'predicted.major_subclass.cttoclus',
                                                           prediction_column = 'predicted.major_subclass')
broad_cell_type_cttoclus_overlap <- broad_cell_type_cttoclus_overlap[broad_cell_type_cttoclus_overlap[['freq']] != 0, ]
broad_cell_type_cttoclus_overlap <- broad_cell_type_cttoclus_overlap[order(broad_cell_type_cttoclus_overlap[['truth']], broad_cell_type_cttoclus_overlap[['freq']], decreasing = T), ]
write.table(broad_cell_type_cttoclus_overlap,
            paste0(work_dir, "/", subset, "/tables/broad_cell_type_to_bakken_predic_cttoclus_overlap.tsv"),
            sep = '\t',
            row.names = F,
            col.names = T)