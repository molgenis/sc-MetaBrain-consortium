# Replication

This code offers functionality to perform replication analyses between sets of expression quantitative trait loci (eQTL) summary statistics.

## How it works

The main code ([replication.py](replication.py)) offers several implementations of common eQTL summary statistics. The following are implemented:

 * [LIMIX](https://github.com/sc-eQTLgen-consortium/limix_qtl): eQTL mapping software
 * LIMIX_REDUCED: adapted from LIMIX but only expecting the [minimal_postporcessing](https://github.com/sc-eQTLgen-consortium/limix_qtl/blob/master/Limix_QTL/post_processing/minimal_postprocess.py) output
 * [mbQTL](https://github.com/molgenis/systemsgenetics/tree/master/mbQTL): eQTL mapping software
 * mbQTL_MetaBrain: adapted from mbQTL but expecting different formatting of gene and SNPs, personal use
 * [eQTLMappingPipeline](https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline): eQTL mapping software
 * [DeconQTL](https://github.com/molgenis/systemsgenetics/tree/master/Decon2/Decon-eQTL): interaction eQTL mapping software
 * [PICALO](https://github.com/molgenis/PICALO): interaction eQTL mapping software

But also certain published summary statistics are supported:
 * eQTLgenPhase2: Unpublished
 * Bryois: [Bryois  et al. Nature Neuroscience 2022](https://doi.org/10.1038/s41593-022-01128-z)
 * Bryois_REDUCED: adapted from Bryois, personal use
 * Fujita: [Fujita et al. Nature Genetics 2024](https://doi.org/10.1038/s41588-024-01685-y)

Each type of summary statistics has a class that loads the data in one of three ways: 1) all effects, 2) top variant per gene, and 3) specific gene - variant combinations. All these functions can be inherited from the default `Dataset` class and do not need to be re-implemented in most cases.

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
4. Visualise the replication

### Implementing custom summary statistics

New summary statistics can easily be added by following these steps:
1. Update the `METHODS` global variable to include your new method
2. Update the `get_method_class` function to return your new class if selected, e.g. `elif method == "METHOD": return METHOD`
3. Implement your new class:

```  
class METHOD(Dataset):
    def __init__(self, *args, **kwargs):
        super(METHOD, self).__init__(*args, **kwargs)
        self.class_name = "METHOD"

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

        # File paths.
        self.set_all_effects_path(effects_path=os.path.join(self.path, "celltype-eqtl-sumstats." + self.cell_type + ".tsv"))
        self.set_top_effects_path(effects_path=os.path.join(self.path, "celltype-eqtl-sumstats." + self.cell_type + ".tsv"))
```

The `self.columns` variable defines how the standard columns are named in your specific summary statistic dataset. Try to fill in as many of them as are available to you using the following formatting:
 * Option 0: data does not exist `[(None, None, sep)]`
 * Option 1: data is a full singular column `[("A", None, sep)]`
 * Option 2: data is part of a singular column `[("A", "(a-zA-Z]+)_", sep)]`, make sure that group 1 of the regex captures the info of interest
 * Option 3: data is two or more full columns `[("A", None, sep), ("B", None, sep)]`
 * Option 4: data is a combination of option 2 and option 3 `[("A", "(a-zA-Z]+)_", sep), ("B", None, sep)]`
 * Option 5: data needs info from other file `[("A", {"A": "a"}, sep)]`, prepare info as a translate dictionary

Certain columns will be automatically decuced from other columns if they are unavailable:
 * `OA` can be deduced if `EA` and `alleles` are available
 * `bonferroni_pvalue` can be deduced if `nominal_pvalue` and `n_tests` are available
 * `zscore` can be deduced if `beta` and `nominal_pvalue` are available
 * `MAF` can be deduced if `AF` is available

Furthermore, if `N` is unavailable in the summary statistics file you can define this using `self.n = N`.

Moreover, some columns can be approximated if allowed (**--allow_infer**, default `False`):
 * `zscore` can be approximated if `beta`, `MAF`, and `N` are available
 * `beta_se` can be approximated if `zscore`, `MAF`, and `N` are available

Next, you need to define the path where all gene - variant combinations are stored in `all_effects_path` and where the top variant per gene are stored in `top_effects_path`. Note that if you do not have a `top_effects_path` file you can leave this empty and the program will automatically select the top effects from the `all_effects_path` file. Similarly, if you do not have a `all_effects_path` you can leave this empty and the `top_effects_path` will be used instead; note that you might get very low overlap between two datasets because of this. Finally, if you are using single-cell eQTL summary stats you can use the `self.cell_type` argument to differentiate between these.

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

Furthermore, it uses two R v4.1 scripts [Rb.R](Rb.R) and [qvalue_truncp.R](qvalue_truncp.R) that require the following packages to be installed:
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
 * **--discovery_name**: Name of the discovery summary statistics.
 * **--discovery_cell_type**: Cell type of the discovery summary statistics. Default: None.

### Replication arguments:
 * **--replication_[method/path/name/cell_type]**: same as discovery but then for the replication summary statistics

### Data extraction arguments:
 * **--gene**: Which gene format to select on. Options: hgnc, ensembl. Default: ensembl.
 * **--snp**: Which variant format to select on. Options: chr:pos, rsid. Default: chr:pos.
 * **--pvalue**: Which pvalue to use. Options: permuted, bonferroni, nominal. Default: permuted.
 * **--effect**: What to consider as the effect column. Options: beta, zscore. Default: zscore.
 * **--allow_infer**: Allow for inferring summary stats information. Note that inferred info are approximations and not exact. Default: False.

### Overlap arguments:
 * **--rm_dupl**: How to deal with duplicates in the replication summary statistics. Options: none) throw error and exit, all) remove all duplicates, mismatches) removed duplicates for which the effect allele does not match. Default: none.
 * **--alpha**: The significance threshold to use. Default: 0.05.
 * **--fdr_calc_method**: The multiple testing correction method to use. Default: 'qvalues'.

### Visualise arguments:
 * **--log_modulus**: Transfer the effect column into log space while maintaining effect direction. Default: False.
 * **--palette**: A color palette file. Default: None.
 * **--extensions**: The figure file extension. Default: 'png'.

### General arguments:
 * **--outdir**: The output directory. Default: current work directory.
 * **--force**: Whether to ignore previously loaded summary statistics. Default: False.
 * **--save**: Whether to store loaded summary statistics. Default: False.

 * **--qvalue_truncp**: The path to the qvalues script. Default: qvalue_truncp.R.
 * **--rb**: The path to the Rb script. Default: Rb.R.

## Output

if **--save** is selected the individual and merged summary statistics will be saved in `replication_data`. A replication figure will be created in `replication_plot` containing three panels: left) all overlapping effects, middle) discovery significant effects, and right) discovery and replication significant effects.

## Author  

Martijn Vochteloo (m.vochteloo@umcg.nl) *(1)*

1. Department of Genetics, University Medical Center Groningen, University of Groningen, Hanzeplein 1, Groningen, The Netherlands

## License  

This project is licensed under the BSD 3-Clause "New" or "Revised" License - see the [LICENSE](LICENSE.txt) file for details
