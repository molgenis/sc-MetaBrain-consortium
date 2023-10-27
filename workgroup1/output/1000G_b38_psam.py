#!/usr/bin/env python3
import gzip

keep = []
with open("../1000G_b37/all_phase3_filtered.psam", "r") as f:
    for i, line in enumerate(f):
        if i == 0:
            continue
        keep.append(line.strip("\n").split("\t")[0])
f.close()

info_dict = {}
with open("20130606_g1k_3202_samples_ped_population.txt", "r") as f:
    for i, line in enumerate(f):
        if i == 0:
            continue
        family_id, sample_id, father_id, mother_id, sex, population, super_population = line.strip("\n").split(" ")
        info_dict[sample_id] = [family_id, sample_id, father_id, mother_id, sex, population, super_population]
f.close()

header = None
with gzip.open("/groups/umcg-biogen/tmp01/annotation/GenomeReference/1kg/1kg-20220422_3202_phased_SNV_INDEL_SV-b38/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz", "rt") as f:
    for line in f:
        if line.startswith("#CHROM"):
            header = line.strip("\n").split("\t")[9:]
            break
f.close()

with open("20130606_g1k_3202_samples_ped_population.psam", "w") as f:
    f.write("#FID\tIID\tPAT\tMAT\tSEX\tPopulation\tSuperPop\n")
    for sample_id in header:
        f.write("\t".join(info_dict[sample_id]) + "\n")
f.close()

with open("1000g_b37_samples.psam", "w") as f:
    f.write("#FID\tIID\tPAT\tMAT\tSEX\tPopulation\tSuperPop\n")
    for sample_id in header:
        if sample_id in keep:
            f.write("\t".join(info_dict[sample_id]) + "\n")
f.close()