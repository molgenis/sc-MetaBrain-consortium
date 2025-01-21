# Replication

This code offers functionality to perform replication analyses between sets of expression quantitative trait loci (eQTL) summary statistics.

## How it works

The main code ([replication.py](scripts/replication.py)) offers several implementations of common eQTL summary statistics. The following are implemented:

 * [LIMIX](https://github.com/sc-eQTLgen-consortium/limix_qtl): eQTL mapping software
 * [mbQTL](https://github.com/molgenis/systemsgenetics/tree/master/mbQTL): eQTL mapping software
 * [eQTLMappingPipeline](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline): eQTL mapping software
 * [DeconQTL](https://github.com/molgenis/systemsgenetics/tree/master/Decon2/Decon-eQTL): interaction eQTL mapping software
 * [PICALO](https://github.com/molgenis/PICALO): interaction eQTL mapping software

But also certain published summary statistics are supported:
 * eQTLgenPhase2: Unpublished
 * Bryois: [Bryois et al. Nature Neuroscience 2022](https://doi.org/10.1038/s41593-022-01128-z)
 * Bryois_REDUCED: Bryois et al. bioRxiv 2021, personal use
 * Fujita: [Fujita et al. Nature Genetics 2024](https://doi.org/10.1038/s41588-024-01685-y)

Each type of summary statistics has a class that loads the data in one of three ways: 1) all effects, 2) top variant per gene, and 3) specific gene - variant combinations.

The steps of the program are as follows:
1. Load in the top effects from the discovery summary statistics
2. Look for these specific gene - variant combinations in the replication summary statistics
3. Overlap and harmonise the data
   * Replication significance is calculated using the **--fdr_calc_method** method only over the discovery significant effects (< **--alpha**).
4. Calculate concordance
   * AC: allelic concordance indicating the proportion of effects with the shared direction of effect
   * coef: Pearson correlation coefficient between eQTL effects
   * pi1: proportion of true positive within the replication data
   * Rb: correlation of effect sizes correcting for standard error of the effect sizes
5Visualise the replication

### Implementing custom summary statistics (no coding required)

New summary statistics can easily be added by following these steps:
1. Copy the example custom class JSON [custom.json](custom.json).
2. Fill in the top fields. Note that these are optional and can also be set using **--[discovery/replication]_all_filename** and **--[discovery/replication]_top_filename**.
   * **class_name**: the name of your custom class.
   * **n**: default sample size equal for all effects. Note: leave empty and use `columns - N` if your files have a sample size per effect.
   * **all_effects_filename**: the default filename for the file containing all gene - variant combinations.
   * **top_effects_filename**: same as **all_effects_filename** but then the file containing the top variant per gene.
4. Fill in the columns info denoting which column in your **[all/top]_effects_filename** represents which type of information. Some examples on how to define this are:
 * Option 0: data does not exist `[(null, null, null)]` or delete the line
 * Option 1: data is a full singular column `"A"` or `[("A", null, null)]`
 * Option 2: data is part of a singular column `[("A", "(a-zA-Z]+)_", null)]`, make sure that group 1 of the regex captures the info of interest
 * Option 3: data is two or more full columns `[("A", null, null), ("B", null, null)]`
 * Option 4: data is two or more full columns with a symbol inbetween `[("A", null, ":"), ("B", null, null)]`
 * Option 5: data is a combination of option 2 and option 3 `[("A", "(a-zA-Z]+)_", null), ("B", null, null)]`
 * Option 6: data needs info from other file `[("A", {"A": "a"}, null)]`, prepare info as a translate dictionary
4. Run [replication.py](scripts/replication.py) with your custom class using **--[discovery/replication]_class_settings**


Note that if you do not have a `top_effects_path` file you can leave this empty and the program will automatically select the top effects from the `all_effects_path` file. Similarly, if you do not have a `all_effects_path` you can leave this empty and the `top_effects_path` will be used instead; note that you might get very low overlap between two datasets because of this. Furthermore, you can use `<.+>` fields to use dynamic input for the **--[discovery/replication]_[all/top]_filename**:
* the field `<CT>` will be replaced with **--[discovery/replication]_cell_type**
* The field `<CHR>` will be replaced with chromosome 1 to 22, X, Y, and MT
* The field `<BATCH>` will be replaced with 0 to 1000
* Other `<.+>` fields can be used to define wildcards, the program will use `glob` to find a singular file that matches. If nothing matches, the wildcard is kept.

### Implementing custom summary statistics (coding required)

In some cases you might need to overwrite default code that loads or processes your input files. In order to do so you need to follow these steps:

New summary statistics can easily be added by following these steps:
1. Update the `METHODS` global variable to include your new method
2. Update the `get_method_class` function to return your new class if selected, e.g. `elif method == "METHOD": return METHOD`
3. Implement your new class:

```  
class METHOD(Dataset):
    def __init__(self, *args, **kwargs):
        super(METHOD, self).__init__(*args, **kwargs)
        self.class_name = "METHOD"
        
         # Update filenames. Only applies if --[all/top]_filename is not set.
        self.update_all_filename(filename="all_results.txt")
        self.update_top_filename(filename="top_results.txt")

        # Set file paths.
        self.set_all_effects_path()
        self.set_top_effects_path()

        # Columns that are in the original file.
        self.columns.update({
            "gene_hgnc": [(None, None, None)], # HGNC symbol of the gene
            "gene_ensembl": [(None, None, None)], # ENSEMBL id of the gene (no version)
            "SNP_rsid": [(None, None, None)], # RS id of the variant
            "SNP_chr:pos": [(None, None, None)], # chr:pos of the variant (no chr prefix)
            "alleles": [(None, None, None)], # alleles of the variant in format alleleA/alleleB
            "EA": [(None, None, None)], # the allele indicating the direction of the beta
            "OA": [(None, None, None)], # the other alelle
            "beta": [(None, None, None)], # the (interaction) eQTL beta
            "beta_se": [(None, None, None)], # the (interaction) eQTL standard error
            "n_tests": [(None, None, None)], # the number of tests (e.g. for top effects the variants had the lowest p-value out of N total variants)
            "nominal_pvalue": [(None, None, None)], # the nominal (interaction) eQTL p-value
            "permuted_pvalue": [(None, None, None)], # the gene level permutation based (interaction) eQTL p-value
            "bonferroni_pvalue": [(None, None, None)], # the bonferroni corrected p-value correcting the nominal p-value for the number of variants considered
            "zscore": [(None, None, None)], # the z-score of the nominal p-value
            "FDR": [(None, None, None)], # the global FDR of the eQTL effect
            "N": [(None, None, None)], # the number of samples
            "AF": [(None, None, None)], # the allele frequency
            "MAF": [(None, None, None)] # the minor allele frequency
        })
```

4. Overwrite functions inherited from the default `Dataset` class or add new function as necessary.


The `self.columns` variable defines how the standard columns are named in your specific summary statistic dataset. Try to fill in as many of them as are available to you. The format should always be `[(column name, regex / dict, suffix)]`, some examples:
 * Option 0: data does not exist `[(None, None, None)]`
 * Option 1: data is a full singular column `[("A", None, None)]`
 * Option 2: data is part of a singular column `[("A", "(a-zA-Z]+)_", None)]`, make sure that group 1 of the regex captures the info of interest
 * Option 3: data is two or more full columns `[("A", None, None), ("B", None, None)]`
 * Option 4: data is two or more full columns with a symbol inbetween `[("A", null, ":"), ("B", null, null)]`
 * Option 5: data is a combination of option 2 and option 3 `[("A", "(a-zA-Z]+)_", None), ("B", None, None)]`
 * Option 6: data needs info from other file `[("A", {"A": "a"}, None)]`, prepare info as a translate dictionary

### Required and autofilled columns

Required columns are:
 * `gene_hgnc` or `gene_ensembl` depending on **--gene** setting.
 * `SNP_rsid` or `SNP_chr:pos` depending on **--snp** setting.
 * `beta` or `zscore` depending on **--effect** setting.
 * `EA`
 * discovery: `nominal_pvalue`, `permuted_pvalue` or `bonferroni_pvalue` depending on **--pvalue** setting.
 * replication: `nominal_pvalue`

Certain columns will be automatically decuced from other columns if they are unavailable:
 * `OA` can be deduced if `EA` and `alleles` are available
 * `bonferroni_pvalue` can be deduced if `nominal_pvalue` and `n_tests` are available
 * `zscore` can be deduced if `beta` and `nominal_pvalue` are available
 * `nominal_pvalue` can be deduced if `zscore` is available
 * `MAF` can be deduced if `AF` is available

Furthermore, if `N` is unavailable in the summary statistics file you can define this using `self.n = N`.

Moreover, some columns can be approximated if allowed (**--allow_infer**, default `False`):
 * `zscore` can be approximated if `beta`, `MAF`, and `N` are available
 * `beta_se` can be approximated if `zscore`, `MAF`, and `N` are available

## Prerequisites  

This program is written in Python v3.7. The program requires the following packages to be installed:  

 * numpy (v1.21.6)
 * pandas (v1.3.5) 
 * matplotlib (v3.5.3)
 * h5py (v3.8.0)
 * natsort (v8.4.0)
 * scipy (v1.7.3)
 * statsmodels (v0.13.5)
 * seaborn (v0.12.2)
 * adjustText (v1.1.1)
 * rpy2 (v3.5.16)

Furthermore, it uses two R v4.1 scripts [Rb.R](scripts/Rb.R) and [qvalue_truncp.R](scripts/qvalue_truncp.R) that require the following packages to be installed:
 * BiocManager (v3.14)
 * ggplot2
 * reshape2
 * qvalue

See 'Installing' on how to install these packages. Note that the performance / reliability of this program on other versions is not guaranteed.

Note that a [Dockerfile](Dockerfile) is available will all required software.

## Usage  

```  
./replication.py -h
```  

 * **-h**, **--help**: show this help message and exit
 * **-v**, **--version**: show program's version number and exit

### Discovery arguments:
 * **--discovery_method**: Method of the discovery summary statistics.
 * **--discovery_path**: Basedir of the discovery summary statistics.
 * **--discovery_all_filename**: Update the default all effects filename of the replication class. Default: class default.
 * **--discovery_top_filename**: Update the default top effects filename of the replication class. Default: class default.
 * **--discovery_name**: Name of the discovery summary statistics.
 * **--discovery_cell_type**: Cell type of the discovery summary statistics. Default: None.
 * **--discovery_class_settings**: JSON file containing custom discovery class settings.
 * **--discovery_rm_dupl**: How to deal with duplicates in the discovery summary statistics. Options: none) throw error and exit, all) remove all duplicates, allbutfirst) removed all except the first effect. Default: 'none'.


### Replication arguments:
 * **--replication_[method/path/all_filename/top_filename/name/cell_type/class_settings]**: same as discovery but then for the replication summary statistics
 * **--replication_rm_dupl**: How to deal with duplicates in the replication summary statistics. Options: **none**) throw error and exit, **all**) remove all duplicates, **mismatches**) removed duplicates for which the effect allele does not match. Default: 'none'.

### Data extraction arguments:
 * **--gene**: Which gene format to select on. Options: hgnc, ensembl. Default: 'ensembl'.
 * **--snp**: Which variant format to select on. Options: chr:pos, rsid. Default: 'chr:pos'.
 * **--pvalue**: Which pvalue to use. Options: permuted, bonferroni, nominal. Default: 'permuted'.
 * **--effect**: What to consider as the effect column. Options: beta, zscore. Default: 'zscore'.
 * **--allow_infer**: Allow for inferring summary stats information. Note that inferred info are approximations and not exact. Default: False.

### Overlap arguments:
 * **--alpha**: The significance threshold to use. Default: 0.05.
 * **--fdr_calc_method**: The multiple testing correction method to use. Default: 'qvalues'.

### Visualise arguments:
 * **--log_modulus**: Transfer the effect column into log space while maintaining effect direction. Default: False.
 * **--cell_type_names**: A JSON file containing method specific to normalised cell type names. Default: None.
 * **--palette**: A color palette file. Default: None.
 * **--extensions**: The figure file extension. Default: 'png'.

### General arguments:
 * **--outdir**: The output directory. Default: current work directory.
 * **--force**: Whether to ignore previously loaded summary statistics. Default: False.
 * **--save**: Whether to store loaded summary statistics. Default: False.
 * **--verbose**: Print additional log messages. Default: False.
 * **--qvalue_truncp**: The path to the qvalues script. Default: 'qvalue_truncp.R'.
 * **--rb**: The path to the Rb script. Default: 'Rb.R'.

## Output

if **--save** is selected the individual and merged summary statistics will be saved in `replication_data`. A replication figure will be created in `replication_plot` containing three panels: left) all overlapping effects, middle) discovery significant effects, and right) discovery and replication significant effects.

## Snakemake automatisation usage

If you wish to create multiple replication plots this can be done using the supplied snakemake pipeline.

### Installing

In order to run the pipeline you require snakemake. The pipeline was developed using snakemake version `5.26.1=0` as part of the [sc-eQTLgen WG1](https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/tree/scMetaBrain) conda environment ([snakemake.yaml](https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/blob/master/Demultiplexing/snakemake.yaml)). In order to install this you require [miniconda3](https://repo.anaconda.com/miniconda/) and subsequently run:
```console
conda env create -f snakemake.yaml -n wg1_snakemake
``` 
Then, to activate the environment, run:

```console
conda activate wg1_snakemake
``` 

### Arguments

The different input can be defined using `disc_inputs` in a dictionary format where the tag will be replaced with a supplied values. For example:
```
disc_inputs: {
   "<TAG>": ["foo", "bar"]
}
```  
will result in a run where `<TAG>` in the `discovery_path` will be replaced with `foo` in the first run but with `bar` in the second run. Multiple tags can be supplied. Furthermore, all values for `disc_cell_types` will also be seperate runs where `<CT>` will be filled in similarly. The same applied for the `replication` input and settings.

### Usage

Before running the pipeline it is advised to perform a `dryrun` to check if all input and settings are valid. Be sure to check the top of the output as import warnings and info are printed.

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

### Output

Each replication run is outputted in a seperate folder. In addition, the following extra files are created:
 * a `_ReplicationStats.txt.gz` containing the replication stats from over all runs.
 * a `Repl_ReplicationStats.{extension}` file with a heatmap visualisation of the replication statistics for all comparisons.

## Author  

Martijn Vochteloo (m.vochteloo@umcg.nl) *(1)*

1. Department of Genetics, University Medical Center Groningen, University of Groningen, Hanzeplein 1, Groningen, The Netherlands

## License  

This project is licensed under the BSD 3-Clause "New" or "Revised" License - see the [LICENSE](LICENSE.txt) file for details
