#!/usr/bin/env python


rule export_vcf_samples:
    input:
        wg1_metadata = FILES["WG1-metadata"],
        vcf = lambda wildcards: VCF_INPUT_FILES[wildcards.md5_hash]
    output:
        samples = config["outputs"]["output_dir"] + "FilterVCF/{md5_hash}_samples.txt.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_input"]["export_vcf_samples_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_input"]["export_vcf_samples_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["combine_input"]["export_vcf_samples_time"]]
    threads: config["combine_input"]["export_vcf_samples_threads"]
    params:
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"],
        script = config["inputs"]["repo_dir"] + "scripts/export_vcf_samples.py",
    log: config["outputs"]["output_dir"] + "log/export_vcf_samples.{md5_hash}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --wg1_metadata {input.wg1_metadata} \
            --vcf {input.vcf} \
            --outfile {output.samples}
        """


rule filter_vcf:
    input:
        vcf = lambda wildcards: VCF_INPUT_FILES[wildcards.md5_hash],
        variants = lambda wildcards: config["inputs"]["variants_path"] if (config["inputs"]["variants_path"] is not None) and (wildcards.md5_hash in [md5_hash(vcf) for vcf in FILES["WG1-ImputedVCF"]]) else [],
        samples = config["outputs"]["output_dir"] + "FilterVCF/{md5_hash}_samples.txt.gz" if config["settings"]["filter_samples"] else ""
    output:
        vcf_gz = config["outputs"]["output_dir"] + "FilterVCF/{md5_hash}.vcf.gz",
        vcf_bgz = config["outputs"]["output_dir"] + "FilterVCF/{md5_hash}.vcf.bgz",
        index = config["outputs"]["output_dir"] + "FilterVCF/{md5_hash}.vcf.bgz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_input"]["filter_vcf_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_input"]["filter_vcf_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["combine_input"]["filter_vcf_time"]]
    threads: config["combine_input"]["filter_vcf_threads"]
    params:
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"],
        script = config["inputs"]["repo_dir"] + "scripts/filter_vcf.py",
        variants = lambda wildcards: "--variants " + config["inputs"]["variants_path"] if (config["inputs"]["variants_path"] is not None) and (wildcards.md5_hash in [md5_hash(vcf) for vcf in FILES["WG1-ImputedVCF"]]) else "",
        samples = "--samples " + config["outputs"]["output_dir"] + "FilterVCF/{md5_hash}_samples.txt.gz" if config["settings"]["filter_samples"] else "",
    log: config["outputs"]["output_dir"] + "log/filter_vcf.{md5_hash}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --vcf {input.vcf} \
            {params.variants} \
            {params.samples} \
            --outfile {output.vcf_gz}
        singularity exec --bind {params.bind} {params.sif} gunzip -c {output.vcf_gz} | \
            singularity exec --bind {params.bind} {params.sif} bgzip > {output.vcf_bgz}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf_bgz}
        """


def get_combine_vcf_input(wildcards):
    vcf_files = []
    if config["inputs"]["variants_path"] is None and not config["settings"]["filter_samples"]:
        if wildcards.vcf == "imputed_hg38":
            vcf_files = FILES["WG1-ImputedVCF"]
        elif wildcards.vcf == "het_filtered":
            vcf_files = FILES["WG1-HetVCF"]
    else:
        if wildcards.vcf == "imputed_hg38":
            vcf_files = [config["outputs"]["output_dir"] + "FilterVCF/{md5_hash}.vcf.bgz".format(md5_hash=md5_hash(vcf)) for vcf in FILES["WG1-ImputedVCF"]]
        elif wildcards.vcf == "het_filtered":
            vcf_files = [config["outputs"]["output_dir"] + "FilterVCF/{md5_hash}.vcf.bgz".format(md5_hash=md5_hash(vcf)) for vcf in FILES["WG1-HetVCF"]]

    return vcf_files


rule combine_vcf:
    input:
        vcf = get_combine_vcf_input
    output:
        vcf = config["outputs"]["output_dir"] + "CombineVCF/{vcf}.vcf.gz",
        index1 = config["outputs"]["output_dir"] + "CombineVCF/{vcf}.vcf.gz.csi",
        index2 = config["outputs"]["output_dir"] + "CombineVCF/{vcf}.vcf.gz.tbi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_input"]["combine_vcf_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_input"]["combine_vcf_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["combine_input"]["combine_vcf_time"]]
    threads: config["combine_input"]["combine_vcf_threads"]
    params:
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"]
    log: config["outputs"]["output_dir"] + "log/combine_vcf.{vcf}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools merge \
            --force-samples \
            --threads {threads} \
            {input.vcf} \
            --write-index \
            -Oz \
            --output {output.vcf}
        singularity exec --bind {params.bind} {params.sif} tabix -p vcf {output.vcf}
        """


def get_combine_files_input(wildcards):
    if wildcards.outfile == "update_sex.psam":
        return FILES["WG1-PSAM"]
    elif wildcards.outfile == "Final_Assignments_demultiplexing_doublets.tsv.gz":
        return FILES["WG1-metadata"]
    elif wildcards.outfile == "azimuth_all.metadata.tsv.gz":
        return FILES["WG2-metadata"]
    elif wildcards.outfile == "cell_annotation.tsv.gz":
        return FILES["CellAnnotation"]
    elif wildcards.outfile == "wg0_file_directories.tsv":
        return FILES["PoolSheet"]
    else:
        return []


rule combine_files:
    input:
        metadata = get_combine_files_input
    output:
        metadata = config["outputs"]["output_dir"] + "CombineFiles/{outfile}"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_input"]["combine_files_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_input"]["combine_files_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["combine_input"]["combine_files_time"]]
    threads: config["combine_input"]["combine_files_threads"]
    params:
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"],
        script = config["inputs"]["repo_dir"] + "scripts/combine_files.py",
    log: config["outputs"]["output_dir"] + "log/combine_files.{outfile}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --files {input.metadata} \
            --outfile {output.metadata}
        """