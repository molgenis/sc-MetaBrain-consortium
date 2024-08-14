# mbQTL - snakemake

Snakemake pipeline build around [mbQTL](https://github.com/molgenis/systemsgenetics/tree/master/mbQTL).

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

See the README of [mbQTL](https://github.com/molgenis/systemsgenetics/tree/master/mbQTL) for information on the arguments.

**Pipeline specific inputs**:
 * `bind_path`: directories that the singularity should have access to (comma seperated).
 * `singularity_image`: path to the singularity image containing the software. I used the [sc-MetaBrain WG3](https://github.com/sc-eQTLgen-consortium/WG3-pipeline-QTL/tree/scMetaBrain) singularity file from this created from this [Dockerfile](https://github.com/sc-eQTLgen-consortium/WG3-pipeline-QTL/blob/scMetaBrain/Dockerfile).
 * `repo_dir`: path to the base directory where the scripts are stored (i.e. parent of `scripts`).
 * `mbqtl_jar`: mbQTL jar [download](https://jenkins.harmjanwestra.nl/job/systemsgenetics_hjw/lastBuild/nl.systemsgenetics$MbQTL/).
 * `output_dir`: the output directory.
 * `output_prefix`: the output filename prefix.
 * `cov`: covariate file. This is a tab-separated file the value for each covariate per sample. The first column contains the covariate names. The rest of the columns are the sample names. Non-numerical covariates are automatically one-hot encoded where the most abundant category is excluded. Expression PCs can be automatically added as covariates by using `n_pcs`.

**Pipeline specific settings**:
 * `include_modes`: which modes to run (options: `all`, `default`, `cov`, `covXPcs`, `XPcs`). For more info, see modes.
 * `n_pcs`: how many PCs should be removed from the expression matrix (e.g. `[0, 5, 10]`). If `cov` is also used these PCs are added to those covariates.
 * `n_genes`: how many genes should be tested per chunk. If empty, all genes are tested in 1 chunk.
 * `use_snpannotation`: whether or not the `snpannotation` option should be used. Default `False`. Automatically set to `False` if `snplimit` or `snpgenelimit` is used since it is faster without `snpannotation` then.
 * `filter_vcf`: whether or not the input VCF should be filtered on variants / samples of interest before running QTL analysis. This adds some pre-processing time but can add substantial speed improvements in the QTL analysis. Only available if `snplimit` or `snpgenelimit` is used. Note that it also filters on the samples in the `gte` file. Default `False`.

**mbQTL standard inputs**:
 * `annotation`: `genes.gtf` of your alignment reference. If the file does not end with `.gtf` it assumes it is a mbQTL annotation file.
 * `exp`: expression file. This is a tab-separated file containing the expression of each feature per sample. The first column contains the features names. The rest of the columns are the sample names.
 * `gte`: genotype-expression-dataset mapping file. If the file endswith `.smf` it assumes it is a genotype-expression mapping file and a dataset column will be added to it in order to be compatable with mbQTL.
 * `vcf`: bgzipped genotype VCF file. If no tabix index (`.tbi`) is present it will be generated by the pipeline.

**Additional to mbQTL manual**:
 * In this implementation `genelimit`, `snplimit` and `snpgenelimit` cannot be combined.
 * The `--out` mbQTL argument is created by `output_dir` + `output_prefix`
 * The `--outputall` mbQTL argument is automatically set to `True` if `snpgenelimit` is used.
 * The `--perm` mbQTL argument is automatically set to `0` if `snpgenelimit` is used.

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
Be aware that `run_qtl` chunks are determined on the fly when `n_genes` is used. As a result, the shown snakemake rules and job counts might not properly the rules that will be executed in reality.

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
The log output will be written to files in `log/` with `log/nohup_*` being the general pipeline log output.
The `slurm_log/` contains log output of SLURM executing your rule as a jobfile (e.g. runtime, memory usage etc.).

#### Unlock:
This script unlocks the working directory if for some reason the manager process got killed.
```{R}
snakemake \
  --snakefile Snakefile \
  --configfile mbQTL.yaml \
  --unlock
```  

#### Modes:
There are multiple 'modes' that can be used:
 * `default`: eQTL mapping on `exp` input without correction for `cov` or `Pcs`
 * `cov`: eQTL mapping on `cov` corrected `exp` input
 * `XPcs`: eQTL mapping on `XPcs` corrected `exp` input
 * `covXPcs`: eQTL mapping on `cov` and `XPcs` corrected `exp` input

The mode is determined automatically based on your input. For example: if you leave `cov` and `n_pcs` empty the mode will be `default`. If you supply a `cov` the mode will be `cov`. etc.

You can choose to run different combinations of modes at the same time but using `include_modes`. For example, if you supply a `cov` but use `include_modes: all` there will also be results generated where `cov` are not corrected for (i.e. mode `default`).

## Important

Note that expression PCs calculation is only performed over the samples that overlap between the expression identifiers in `gte` and the columns in `exp`.

## Output

Each eQTL run is outputted in a seperate folder: e.g. no covariate or PCs (`default`), cov (`cov`), 5 Pcs (`5Pcs`), or cov + 5 Pcs (`cov5Pcs`) all get their own folder in `output` containing the default mbQTL output files. In addition, the following extra files are created:
 * a `*.Pcs.png` and `*.Scree.png` figure containing visualisations of the expression matrix PCA. If covariates are corrected a plot after correction is created as well.
 * a `-TopEffectsWithqval.txt` file with 2 columns added:
   * `PvalueNominalThreshold`: nominal p-value thresholds based on the permutation beta distribution (`BetaDistAlpha` and `BetaDistBeta`).
   * `qval`: based on the nominal p-values (`MetaP`) if no permutation are run (`perm: 0`) or the permutation p-values (`BetaAdjustedMetaP`) if permutations are run (`perm: >0`).
 * a `-results.txt` file with the number of effects with nominal p-value (`MetaP`), permuted p-value (`BetaAdjustedMetaP`) and qvalue (`qval`) below significance threshold (`<0.05`) per QTL run (e.g. 0Pcs removed, 5Pcs removed, etc.).

## Author  

Martijn Vochteloo (m.vochteloo@umcg.nl) *(1)*

1. Department of Genetics, University Medical Center Groningen, University of Groningen, Hanzeplein 1, Groningen, The Netherlands

## License  

This project is licensed under the BSD 3-Clause "New" or "Revised" License - see the [LICENSE](LICENSE.txt) file for details