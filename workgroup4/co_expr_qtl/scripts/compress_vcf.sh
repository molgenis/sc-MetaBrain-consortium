#!/bin/bash

ml HTSlib/1.19.1-GCCcore-11.3.0 

# Input parameters
input_vcf=$1
output_vcf=$2

# Compress the VCF file
zcat $input_vcf | bgzip -c > $output_vcf

# Index the compressed VCF file
tabix -p vcf $output_vcf