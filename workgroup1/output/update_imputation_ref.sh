#!/bin/bash
for i in {1..22}
do
   echo "minimac4 --update-m3vcf ../imputation_m3vcf/chr$i.m3vcf.gz > chr$i.msav"
   singularity exec --bind /groups/umcg-biogen/tmp01/,/groups/umcg-biogen/tmp01/umcg-mvochteloo/simulated_home:/home/umcg-mvochteloo /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-09-06-scMetaBrain-WorkGroup1QC/20231006-0-WG1-pipeline-QC.sif minimac4 --update-m3vcf ../imputation_m3vcf/chr$i.m3vcf.gz > chr$i.msav
done
