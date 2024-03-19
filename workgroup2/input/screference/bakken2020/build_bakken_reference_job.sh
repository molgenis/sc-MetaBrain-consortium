#!/bin/bash
#SBATCH --job-name=build_bakken_reference
#SBATCH --output=/groups/umcg-biogen/tmp01/input/processeddata/single-cell/screference/2023-12-19-Bakken2020/log/build_bakken_reference.out
#SBATCH --error=/groups/umcg-biogen/tmp01/input/processeddata/single-cell/screference/2023-12-19-Bakken2020/log/build_bakken_reference.out
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load R

Rscript /groups/umcg-biogen/tmp01/input/processeddata/single-cell/screference/2023-12-19-Bakken2020/integrate.R /groups/umcg-biogen/tmp01/input/rawdata/single-cell/Bakken2020/data/Human_M1_10xV3_Matrix.RDS /groups/umcg-biogen/tmp01/input/rawdata/single-cell/Bakken2020/data/Human_M1_10xV3_Metadata.feather /groups/umcg-biogen/tmp01/input/processeddata/single-cell/screference/2023-12-19-Bakken2020/seurat_object/human_m1_integrated_full_ref.rds

Rscript /groups/umcg-biogen/tmp01/input/processeddata/single-cell/screference/2023-12-19-Bakken2020/build_reference.R /groups/umcg-biogen/tmp01/input/processeddata/single-cell/screference/2023-12-19-Bakken2020/seurat_object/human_m1_integrated_full_ref.rds all /groups/umcg-biogen/tmp01/input/processeddata/single-cell/screference/2023-12-19-Bakken2020/reference/Bakken2020.rds
