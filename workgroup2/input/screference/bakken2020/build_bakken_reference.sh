#!/bin/bash

mkdir data
mkdir logs
mkdir seurat_objects
mkdir reference

wget http://data.nemoarchive.org/publication_release/Lein_2020_M1_study_analysis/Transcriptomics/sncell/10X/human/processed/counts/counts/M1/Human_M1_10xV3_Matrix.RDS -P data
wget http://data.nemoarchive.org/publication_release/Lein_2020_M1_study_analysis/Transcriptomics/sncell/10X/human/processed/counts/counts/M1/Human_M1_10xV3_Metadata.feather -P data
echo "Human motorcortex data downloaded on: $(date)" > logs/download_data.log

Rscript scripts/integrate.R data/Human_M1_10xV3_Matrix.RDS data/Human_M1_10xV3_Metadata.feather seurat_objects/human_m1_integrated_full_ref.rds > logs/integrate.Rout 2>&1

for i in "all" 1000 2000 4000 8000 16000
do
  mkdir reference/downsample$i
	Rscript scripts/build_reference.R seurat_objects/human_m1_integrated_full_ref.rds $i reference/downsample$i > logs/build_reference_$i.Rout 2>&1
done