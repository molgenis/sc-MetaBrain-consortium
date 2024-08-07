# mbQTL - snakemake

Snakemake pipeline build around mbQTL.

## Info

This implements preprocessing steps to create the correct input files as well as perform covariate and/or expression PCs correction. The mbQTL top effects files are automatically combined, multiple testing correction using [qvalue](https://github.com/StoreyLab/qvalue) is applied, and the number of eQTLs are counted.

## Installing

In order to run the pipeline you require snakemake. The pipeline was developed using snakemake version `5.26.1=0` as part of the [sc-eQTLgen WG1](https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/tree/scMetaBrain) conda environment ([snakemake.yaml](https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/blob/master/Demultiplexing/snakemake.yaml)). In order to install this you require [miniconda3](https://repo.anaconda.com/miniconda/) and subsequently run:
```
conda env create -f snakemake.yaml -n wg1_snakemake
``` 
Then, to activate the environment, run:

```
conda activate wg1_snakemake
``` 

## Arguments

See the [manual](https://github.com/molgenis/systemsgenetics/tree/master/mbQTL) of mbQTL for information on mbQTL arguments.

 * `bind_path`: directories that the singularity should have access to (comma seperated)
 * `singularity_image`: path to the singularity image containing the software. I used the [sc-MetaBrain WG3](https://github.com/sc-eQTLgen-consortium/WG3-pipeline-QTL/tree/scMetaBrain) singularity file from this created from this [Dockerfile](https://github.com/sc-eQTLgen-consortium/WG3-pipeline-QTL/blob/scMetaBrain/Dockerfile).
 * `repo_dir`: path to the base directory where the scripts are stored (i.e. parent of `scripts`)
 * `mbqtl_jar`: mbQTL jar [download](https://jenkins.harmjanwestra.nl/job/systemsgenetics_hjw/lastStableBuild/nl.systemsgenetics$MbQTL/)
 * `annotation`: `genes.gtf` of your alignment reference. If the file does not end with `.gtf` it assumes it is a mbQTL annotation file.
 * `gte`: genotype-expression-dataset mapping file. If the file endswith `.smf` it assumes it is a genotype-expression mapping file and will add a dataset column to it in order to be compatable with mbQTL.
 * `cov`: covariate file
 * `output_dir`: the output directory
 * `output_prefix`: the output filename prefix
 * `n_pcs`: how many PCs should be removed. If `cov` is also used these PCs are added to those covariates.
 * `n_genes`: how many genes should be tested per chunk. If empty, all genes are tested in 1 chunk.
 * `use_snpannotation`: whether or not the snpannotation option should be used. Default False.
 * `filter_vcf`: whether or not the input VCF should be filtered on variants / samples of interest before running QTL analysis. This can add substantial speed improvements. Only avaialble if `snplimit` or `snpgenelimit` is used. Default False.

### Additional to mbQTL manual:
 * In this implementation `genelimit`, `snplimit` and `snpgenelimit` cannot be combined.
 * The `--out` mbQTL argument is created by `output_dir` + `output_prefix`

## Usage  

#### Dry run:
This script show what rules will be executed.
```{R}
snakemake \
  --snakefile Snakefile \
  --configfile mbQTL.yaml \
  --dryrun \
  --cores 1 \
  --reason
```  
Be aware that `run_qtl` chunks are determined on the fly; dry-run will not show how many chunks will run. If you are using `n_genes` and observe 1 `run_qtl`; this is normal.

#### Run (local):
This script runs the rules in your current session.
```{R}
snakemake \
  --snakefile Snakefile \
  --configfile mbQTL.yaml \
  --cores 1
```

#### Run (SLURM):
This script runs the rules by submitting them to the SLURM queue.
```{R}
LOGDIR=TODO
mkdir -p $LOGDIR
mkdir -p $LOGDIR/log
mkdir -p $LOGDIR/slurm_log

nohup \
snakemake \
    --snakefile Snakefile \
    --configfile mbQTL.yaml \
    --rerun-incomplete \
    --jobs 10000 \
    --use-singularity \
    --restart-times 0 \
    --keep-going \
    --latency-wait 60 \
    --cluster \
       "sbatch \
       --job-name=snakemake_{rule}_{wildcards} \
       --nodes=1 \
       --cpus-per-task={threads} \
       --mem=\$(({resources.mem_per_thread_gb} * {threads}))G \
       --time={resources.time} \
       --output=$LOGDIR/slurm_log/{rule}_{wildcards}.out \
       --error=$LOGDIR/slurm_log/{rule}_{wildcards}.out \
       --export=NONE \
       --qos=regular \
       --parsable" \
    > $LOGDIR/log/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &

echo "Check status of command with:" ps -p $! -u
```  

#### Unlock:
This script unlocks the working directory if for some reason the manager process got killed.
```{R}
snakemake \
  --snakefile Snakefile \
  --configfile mbQTL.yaml \
  --unlock
```  

## Important

Note that expression PCs calculation is only performed over the samples that overlap between the expression identifiers in `gte` and the columns in `exp`.

## Author  

Martijn Vochteloo (m.vochteloo@umcg.nl) *(1)*

1. Department of Genetics, University Medical Center Groningen, University of Groningen, Hanzeplein 1, Groningen, The Netherlands

## License  

This project is licensed under the BSD 3-Clause "New" or "Revised" License - see the [LICENSE](LICENSE.txt) file for details
