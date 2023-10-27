#!/bin/bash

#wget https://www.dropbox.com/s/xso2vt3p9h2rh8m/1000G.tar.gz
#tar -xzf 1000G.tar.gz
#mv 1000G 1000G_b37

#wget https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz
#gunzip GRCh37_to_GRCh38.chain.gz
#rm GRCh37_to_GRCh38.chain.gz

#mkdir 1000G_b38
#cd 1000G_b38

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt

./1000G_b38_psam.py

singularity exec --bind /groups/umcg-biogen/tmp01/,/groups/umcg-biogen/tmp01/umcg-mvochteloo/simulated_home:/home/umcg-mvochteloo /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-09-06-scMetaBrain-WorkGroup1QC/20231006-0-WG1-pipeline-QC.sif awk 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$5}}' ../1000G_b37/all_phase3_filtered.pvar > all_phase3_filtered_b37.bed

singularity exec --bind /groups/umcg-biogen/tmp01/,/groups/umcg-biogen/tmp01/umcg-mvochteloo/simulated_home:/home/umcg-mvochteloo /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-09-06-scMetaBrain-WorkGroup1QC/20231006-0-WG1-pipeline-QC.sif CrossMap.py bed ../GRCh37_to_GRCh38.chain all_phase3_filtered_b37.bed all_phase3_filtered_b38.bed

singularity exec --bind /groups/umcg-biogen/tmp01/,/groups/umcg-biogen/tmp01/umcg-mvochteloo/simulated_home:/home/umcg-mvochteloo /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-09-06-scMetaBrain-WorkGroup1QC/20231006-0-WG1-pipeline-QC.sif awk 'BEGIN{{FS=OFS="\t"}}NR>1{{print $1":"$2":"$5"_"$6}}' all_phase3_filtered_b38.bed > 1000g_b37_variants.txt

rm pmerge_list.txt
for i in {1..22}
do
  # wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz

  singularity exec --bind /groups/umcg-biogen/tmp01/,/groups/umcg-biogen/tmp01/umcg-mvochteloo/simulated_home:/home/umcg-mvochteloo /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-09-06-scMetaBrain-WorkGroup1QC/20231006-0-WG1-pipeline-QC.sif plink2 \
      --threads 4 \
      --vcf /groups/umcg-biogen/tmp01/annotation/GenomeReference/1kg/1kg-20220422_3202_phased_SNV_INDEL_SV-b38/1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
      --psam 20130606_g1k_3202_samples_ped_population.psam \
      --keep 1000g_b37_samples.psam \
      --extract 1000g_b37_variants.txt \
      --new-id-max-allele-len 500 \
      --max-alleles 2 \
      --set-all-var-ids @:#:\$r_\$a \
      --make-pgen \
      --out 1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel

  echo 1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel >> pmerge_list.txt
done

# wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz

singularity exec --bind /groups/umcg-biogen/tmp01/,/groups/umcg-biogen/tmp01/umcg-mvochteloo/simulated_home:/home/umcg-mvochteloo /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-09-06-scMetaBrain-WorkGroup1QC/20231006-0-WG1-pipeline-QC.sif plink2 \
    --threads 4 \
    --vcf /groups/umcg-biogen/tmp01/annotation/GenomeReference/1kg/1kg-20220422_3202_phased_SNV_INDEL_SV-b38/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz \
    --psam 20130606_g1k_3202_samples_ped_population.psam \
    --split-par b38 \
    --keep 1000g_b37_samples.psam \
    --extract 1000g_b37_variants.txt \
    --max-alleles 2 \
    --new-id-max-allele-len 500 \
    --set-all-var-ids @:#:\$r_\$a \
    --make-pgen \
    --out 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2

echo 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2 >> pmerge_list.txt

singularity exec --bind /groups/umcg-biogen/tmp01/,/groups/umcg-biogen/tmp01/umcg-mvochteloo/simulated_home:/home/umcg-mvochteloo /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-09-06-scMetaBrain-WorkGroup1QC/20231006-0-WG1-pipeline-QC.sif plink2 \
    --threads 4 \
    --pmerge-list pmerge_list.txt pfile \
    --out 20220422_3202_phased_SNV_INDEL_SV_b38

singularity exec --bind /groups/umcg-biogen/tmp01/,/groups/umcg-biogen/tmp01/umcg-mvochteloo/simulated_home:/home/umcg-mvochteloo /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-09-06-scMetaBrain-WorkGroup1QC/20231006-0-WG1-pipeline-QC.sif cp 20220422_3202_phased_SNV_INDEL_SV_b38.pvar 20220422_3202_phased_SNV_INDEL_SV_b38_original.pvar

./1000G_b38_pvar.py
