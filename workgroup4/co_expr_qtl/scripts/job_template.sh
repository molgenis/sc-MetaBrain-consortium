#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=3g
#SBATCH --cpus-per-task=2
#SBATCH -J JOBNAME
#SBATCH -o LOGPREFIX.log
#SBATCH -e LOGPREFIX.err

set -e
set -u

#ml Java/11-LTS
# ml Java/11.0.2
# ml Java/11.0.16
ml Java/11-LTS
# CHROM, BATCHFILE, OUTPREFIX
# EXP, GTE, GENOTYPE

threads=2
nice -n20 java -Xmx1500m \
        -Djava.util.concurrent.ForkJoinPool.common.parallelism=$threads \
        -Dmaximum.threads=$threads -Dthread.pool.size=$threads \
        -jar /groups/umcg-biogen/tmp02/tools/MbQTL-1.5.1-SNAPSHOT-jar-with-dependencies.jar \
        -m mbqtl \
        -a ANNOTATION \
        -e EXPRESSION \
        -g GTE \
        -v GENOTYPE \
        --chr CHROM \
        --perm 10 \
        --minobservations 15 \
        -gl BATCHFILE \
        -o OUTPREFIX \
        --expgroups GROUPS \
        --mingenotypecount 2 \
        --fisherzmeta \
        -sgl SNPGENELIMIT