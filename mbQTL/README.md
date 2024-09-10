# mbQTL - snakemake

Snakemake pipeline build around [mbQTL](https://github.com/molgenis/systemsgenetics/tree/master/mbQTL).

## Info

The snakemake implements the following additions to the standard mbQTL software:
 * pre-process input files:
   * create gene annotation file from a `gtf` file
   * filter the input VCF on variants of interest to reduce disk IO; great speed increase if you use `--snplimit` or `--snpgenelimit`
 * calculate and visualise expression PCs (possibly per dataset in parallel)
 * correct expression input for covariates and / or N expression PCs (possibly per dataset in parallel)
 * QTL analysis in arbitrary number of batches on the HPC in parallel (including merging of output files)
 * perform multiple testing correction over top effects (Bonferroni, Benjamini-Hochberg, [qvalue](https://github.com/StoreyLab/qvalue))
 * summarise number of eQTLs detected (possibly per N expression PCs removed)

Furthermore, quality of life settings and setting warnings are added to make the use of mbQTL even easier.

## Installing

In order to run the pipeline you require snakemake. The pipeline was developed using snakemake version `5.26.1=0` in a conda [environment](snakemake.yaml)). In order to install this you require [miniconda3](https://repo.anaconda.com/miniconda/) and subsequently run:
```console
conda env create -f snakemake.yaml -n mbqtl
``` 
Then, to activate the environment, run:

```console
conda activate mbqtl
``` 

## Arguments

See the README of [mbQTL](https://github.com/molgenis/systemsgenetics/tree/master/mbQTL) for information on the arguments.

**Pipeline specific inputs**:
 * `bind_path`: directories that the singularity should have access to (comma seperated).
 * `singularity_image`: path to the singularity image containing the required software. Can be created from this [Dockerfile](Dockerfile).
 * `repo_dir`: path to the base directory where the scripts are stored (i.e. parent of `scripts`).
 * `output_dir`: the output directory.
 * `output_prefix`: the output filename prefix.
 * `cov`: covariate file. Tab-separated file containing the covariate values of each feature per sample. The first column contains the covariates. The rest of the columns are the sample names. Non-numerical covariates are automatically one-hot encoded where the most abundant category is excluded. If a `gte` file is given the dataset column will be used os covariate. Expression PCs can be automatically added as covariates by using `n_pcs`.

**Pipeline specific settings**:
 * `preflight_checks`: perform pre-flight checks such as checking the samples overlap between the input files. No other results will be generated. Default `False`.
 * `plot_pca`: whether or not to PCA visualise the expression matrix. Default `False`. 
 * `map_qtls`: whether or not to eQTLs should be mapped, e.g. if you wish to inspect the PCA plots first. Default `True`.
 * `include_modes`: which modes to run (options: `all`, `default`, `cov`, `covXPcs`, `XPcs`). For more info, see modes. Default `null`.
 * `force_mega`: force the covariate correction and / or eQTL mapping to be done over all samples at once (options: `all`, `cov`, `qtl`, `none`). Default: `null`.
 * `n_pcs`: how many PCs should be removed from the expression matrix (e.g. `[0, 5, 10]`). If `cov` is also used these PCs are added to those covariates. Default `null`.
 * `filter_vcf`: whether or not the input VCF should be filtered on variants / samples of interest before running QTL analysis. This adds some pre-processing time but can add substantial speed improvements in the QTL analysis. Only available if `snplimit` or `snpgenelimit` is used. Note that it also filters on the samples in the `gte` file. Default `False`.
 * `n_genes`: how many genes should be tested per chunk. If empty, all genes are tested in 1 chunk. Default `100`.
 * `use_snpannotation`: whether or not the `snpannotation` option should be used. Automatically set to `False` if `snplimit` or `snpgenelimit` is used since it is faster without `snpannotation` then. Default `False`.
 * `alpha`: QTL significance threshold. Default `0.05`.
 * `feature_name`: the info column that in your gtf that describes the names of the genes in your expression matrix. Default `gene_name`.
 * `autosomes_only`: whether or not to only include autosomal chromosomes. Default `True`.
 * `eval_n_pcs`: how many PCs to evaluate for outlier detection; same number of Pcs is plotted. Default `3`. 
 * `sample_outlier_zscore`: max absolute zscore calculated over `eval_n_pcs` number of PCs, higher than this number is defined as outlier. Note that outlier samples are not automatically excluded from the analysis. Default `3`.
 * `java_memory_buffer`: memory buffer in Gb to request in addition to what is set for `-Xmx` and `-Xms` to prevent out of memory isues in Java. Default `1`.
 * `force`: prevent snakemake from updating input settings that are unlogical. Use with caution. Default `False`.
 * `debug`: set logger to level DEBUG printing additional information. Default `False`.
 * `mbqtl_jar`: use this mbQTL jar instead of the one in the singularity image. Default: `null`.

**mbQTL standard inputs**:
 * `annotation`: `genes.gtf` of your alignment reference. If the file does not end with `.gtf` it assumes it is a mbQTL annotation file: tab seperated with gene ID, chromosome, gene start, gene end, and strand.
 * `exp`: expression file. This is a tab-separated file containing the expression of each feature per sample. The first column contains the features names. The rest of the columns are the sample names.
 * `gte`: genotype to expression coupling file. This is a tab-separated file with three columns; genotype ID, expression ID, dataset (name). There is no header. If the file endswith `.smf` it assumes it is a genotype-expression mapping file and a dataset column will be added to it in order to be compatable with mbQTL.
 * `vcf`: bgzipped genotype VCF file. If no tabix index (`.tbi`) is present it will be generated by the pipeline.

**Additional to mbQTL manual**:
 * In this implementation `snpgenelimit` has priority over `genelimit`, `snplimit`. As a result, if `snpgenelimit` is given and `filter_vcf` is True, any effects that are in the `snplimit` file but not in `snpgenelimit` file will not be tested. Similarly, if `ngenes`>1 the batches will be created over the genes in `snpgenelimit` and extra genes in `genelimit` will not be tested.
 * The `--out` mbQTL argument is created by `output_dir` + `output_prefix`
 * The `--outputall` mbQTL argument is automatically set to `True` if `snpgenelimit` is used.
 * The `--perm` mbQTL argument is automatically set to `0` if `snpgenelimit` is used.

The following arguments of mbQTL are not (yet) implemented: `--expgroups`, `--nriters`, `--sortbyz`, and `--testnonparseablechr`.

## Usage  

Before running the pipeline it is advised to set `preflight_checks` to True and run a `Run (local)` to check if all input and settings are valid and that the sample overlap is as expected. Be sure to check the top of the output as import warnings and info are printed.

#### Visualise pipeline:
This script that will generate a `dag.svg` file that shows the rules that will be executed and in what order.
```console
snakemake \
  --snakefile Snakefile \
  --configfile mbQTL.yaml \
  --dag | dot -Tsvg > dag.svg
```  

#### Dry run:
This script show what rules will be executed.
```console
snakemake \
  --snakefile Snakefile \
  --configfile mbQTL.yaml \
  --dryrun \
  --cores 1 \
  --reason
```  
Be aware that some variables are determined at runtime and as a result, the shown snakemake rules and job counts might not properly the rules that will be executed in reality.


#### Run (local):
This script runs the rules in your current session.
```console
snakemake \
  --snakefile Snakefile \
  --configfile mbQTL.yaml \
  --cores 1
```

#### Run (SLURM):
This script runs the rules by submitting them to the SLURM queue.
```console
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

#### Report:
This script creates a HTML report of the figures.
```console
snakemake \
  --snakefile Snakefile \
  --configfile mbQTL.yaml \
  --report mbQTL.html
```

#### Unlock:
This script unlocks the working directory if for some reason the manager process got killed.
```console
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

Please keep in mind that:
 * Expression PCs calculation is per dataset over the samples that overlap between the expression identifiers in `gte` and the columns in `exp`. 
 * Expression PCs are calculated without first removing covariates in `cov`.
 * Snakemake checks if the output files exist but not with what settings they were generated. Therefore, if you wish to rerun the eQTL mapping with different settings while keeping the old files it is recommended to change the complete `output_dir` rather than just changing the`output_prefix` since wrongly pre-processed input files might be used otherwise.

## Output

Each eQTL run is outputted in a seperate folder: e.g. no covariate or PCs (`default`), cov (`cov`), 5 Pcs (`5Pcs`), or cov + 5 Pcs (`cov5Pcs`) all get their own folder in `output` containing the default mbQTL output files. In addition, the following extra files are created:
 * if `plot_pca` is True, a `*.Pcs.png` and `*.Scree.png` figure containing visualisations of the expression matrix PCA. If covariates are corrected a plot after correction is created as well.
 * a `-TopEffectsWithqval.txt` file with a couple of columns added:
   * `BonfAdjustedMetaP`: nominal p-value (`MetaP`) multiplied by the number of tests performed for that gene (`NrTestedSNPs`).
   * `BonfBHAdjustedMetaP`: boneferroni corrected p-value (`BonfAdjustedMetaP`) converted to Benjamini-Hochberg FDR values.
   * `PvalueNominalThreshold`: nominal p-value thresholds based on the permutation beta distribution (`BetaDistAlpha` and `BetaDistBeta`, only if `perm: >0`).
   * `bh_fdr`: Benjamini-Hochberg FDR values calculated over the p-values (using `BetaAdjustedMetaP` if `perm: >0` else `MetaP`).
   * `qval`: [qvalue](https://github.com/StoreyLab/qvalue) q-values calculated over the p-values (using `BetaAdjustedMetaP` if `perm: >0` else `MetaP`).
 * a `-results.txt` file with the number of effects, tests, and significance values (e.g. `MetaP   BetaAdjustedMetaP`, `BonfAdjustedMetaP`, `BonfBHAdjustedMetaP`, `bh_fdr`, and `qval`) below the significance threshold (`<0.05`) per QTL run (e.g. 0Pcs removed, 5Pcs removed, etc.). If a meta-analysis was performed the per dataset, the number of nominal significant effects (based on `DatasetZScores`) are also counted.

## Author  

Martijn Vochteloo (m.vochteloo@umcg.nl) *(1)*

1. Department of Genetics, University Medical Center Groningen, University of Groningen, Hanzeplein 1, Groningen, The Netherlands

## License  

This project is licensed under the BSD 3-Clause "New" or "Revised" License - see the [LICENSE](LICENSE.txt) file for details
