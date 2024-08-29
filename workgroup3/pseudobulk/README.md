# pseudobulk

Snakemake pipeline to create pseudobulk from CelRanger output.

## Info

## Installing

In order to run the pipeline you require snakemake. The pipeline was developed using snakemake version `5.26.1=0` as part of the [sc-eQTLgen WG1](https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/tree/scMetaBrain) conda environment ([snakemake.yaml](https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/blob/master/Demultiplexing/snakemake.yaml)). In order to install this you require [miniconda3](https://repo.anaconda.com/miniconda/) and subsequently run:
```console
conda env create -f snakemake.yaml -n wg1_snakemake
``` 
Then, to activate the environment, run:

```console
conda activate wg1_snakemake
``` 

## Arguments

**Pipeline specific inputs**:
 * `bind_path`: directories that the singularity should have access to (comma seperated).
 * `singularity_image`: path to the singularity image containing the software. Can be created from this [Dockerfile](Dockerfile).
 * `repo_dir`: path to the base directory where the scripts are stored (i.e. parent of `scripts`).
 * `poolsheet`:
 * `psam`:
 * `cell_annotation`:
 * `droplet_type_annotation`:
 * `cell_type_annotation`:
 * `rb_genes`:
 * `mt_genes`:
 * `cell_type_pairing`:
 * `output_dir`: the output directory.
 * `palette`:

**Pipeline specific settings**:
 * `ancestry`:
 * `cell_level`:
 * `cell_type`:
 * `ncount_rna`:
 * `nfeature_rna`:
 * `percent_rb`:
 * `percent_mt`:
 * `malat1`:
 * `min_cells`:
 * `min_ind_expr`:
 * `min_cpm`:
 * `feature_name`:
 * `sample_aggregate`:
 * `aggregate_fun`:
 * `debug`: set logger to level DEBUG printing additional information. Default `False`.

## Usage  

Before running the pipeline it is adviced to perform a `dryrun` to check if all input and settings are valid. Be sure to check the top of the output as import warnings and info are printed.

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

#### Unlock:
This script unlocks the working directory if for some reason the manager process got killed.
```console
snakemake \
  --snakefile Snakefile \
  --configfile mbQTL.yaml \
  --unlock
```  

## Output

## Author  

Martijn Vochteloo (m.vochteloo@umcg.nl) *(1)*

1. Department of Genetics, University Medical Center Groningen, University of Groningen, Hanzeplein 1, Groningen, The Netherlands

## License  

This project is licensed under the BSD 3-Clause "New" or "Revised" License - see the [LICENSE](LICENSE.txt) file for details
