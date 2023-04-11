#!/usr/bin/env python3

"""
File:         initialize_wg3_files.py
Created:      2023/04/11
Last Changed:
Author:       M.Vochteloo

Copyright (C) 2022 M.Vochteloo
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""

# Standard imports.
from __future__ import print_function
from datetime import datetime
import argparse
import os
import re

# Third party imports.

# Local application imports.

# Metadata
__program__ = "Initialize Workgroup3 Files"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@rug.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)

"""
Syntax: 

### Mathys2019 ###
./initialize_wg3_files.py \
    --work_dir /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-04-11-WorkGroup3eQTLAndDEA \
    --dataset_outdir Mathys2019 \
    --image_folder /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-04-11-WorkGroup3eQTLAndDEA/WG3-pipeline-QTL/ \
    --wg1_folder /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2023-02-02-AMP_AD/Step1-Imputation/ \
    --wg2_folder /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-02-02-WorkGroup2CellType/2023-04-05-Mathys2019/all/ \
    --unimputed_folder /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2022-10-07-WorkGroup1QC/2023-02-02-AMP_AD/Step1-Imputation/harmonize_hg38/EUR_harmonised_hg38 \
    --ancestry_split_vcf vcf_merged_by_ancestries/ \
    --wg1_psam update_sex_ancestry/update_sex.psam \
    --gtf_annotation_file /groups/umcg-biogen/tmp01/input/processeddata/single-cell/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
    --cell_annotation /groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-02-02-WorkGroup2CellType/2023-04-05-Mathys2019/all/data/seurat_metadata.csv \
    --wg1_singlets_assigned ../Step2-DemultiplexingAndDoubletRemoval/QC_figures/seurat_object_all_pools_singlet_barcodes_final_assignments.rds
    

"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.step = getattr(arguments, 'step')
        self.work_dir = getattr(arguments, 'work_dir')
        dataset_outdir = getattr(arguments, 'dataset_outdir')
        self.image_folder = getattr(arguments, 'image_folder')
        self.top_dir = ",".join(getattr(arguments, 'top_dir'))

        # Pre-process the dataset output directory.
        date_str = datetime.now().strftime("%Y-%m-%d")
        if not re.match("^[0-9]{4}-(0[1-9]|1[0-2])-(0[1-9]|[1-2][0-9]|3[0-1])", dataset_outdir):
            dataset_outdir = "{}-{}".format(date_str, dataset_outdir)
        self.dataset_outdir = dataset_outdir

        # Step 1 arguments.
        self.data_prep_config = getattr(arguments, 'data_prep_config')
        self.wg1_folder = getattr(arguments, 'wg1_folder')
        self.wg2_folder = getattr(arguments, 'wg2_folder')
        self.unimputed_folder = getattr(arguments, 'unimputed_folder')
        self.ancestry_split_vcf = getattr(arguments, 'ancestry_split_vcf')
        self.wg1_psam = getattr(arguments, 'wg1_psam')
        self.gtf_annotation_file = getattr(arguments, 'gtf_annotation_file')
        self.cell_annotation = getattr(arguments, 'cell_annotation')
        self.wg1_singlets_assigned = getattr(arguments, 'wg1_singlets_assigned')
        self.genotype_dir = getattr(arguments, 'genotype_dir')
        self.wg1_wg2_qc_taged = getattr(arguments, 'wg1_wg2_qc_taged')
        self.number_of_genes_in_chunk = getattr(arguments, 'number_of_genes_in_chunk')

        if self.step is None or self.step == 1:
            for label, value in [("--data_prep_config", self.data_prep_config),
                                 ("--wg1_folder", self.wg1_folder),
                                 ("--wg2_folder", self.wg2_folder),
                                 ("--unimputed_folder", self.unimputed_folder),
                                 ("--ancestry_split_vcf", self.ancestry_split_vcf),
                                 ("--wg1_psam", self.wg1_psam),
                                 ("--gtf_annotation_file", self.gtf_annotation_file),
                                 ("--cell_annotation", self.cell_annotation),
                                 ("--wg1_singlets_assigned", self.wg1_singlets_assigned),
                                 ("--wg1_wg2_qc_taged", self.wg1_wg2_qc_taged),
                                 ("--number_of_genes_in_chunk", self.number_of_genes_in_chunk)]:
                if value is None:
                    print("Argument {} is required when --step equals {}.".format(label, self.step))
                    exit()

        # Step 2 arguments.
        self.qtl_config = getattr(arguments, 'qtl_config')

        if self.step is None or self.step == 2:
            for label, value in [("--qtl_config", self.qtl_config)]:
                if value is None:
                    print("Argument {} is required when --step equals {}.".format(label, self.step))
                    exit()

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # general arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit")
        parser.add_argument("--step",
                            type=int,
                            choices=[1, 2],
                            default=None,
                            help="Which step to prepare. Default: None.")
        parser.add_argument("--work_dir",
                            type=str,
                            required=True,
                            default=None,
                            help="The working directory.")
        parser.add_argument("--dataset_outdir",
                            type=str,
                            required=True,
                            default=None,
                            help="The name of the output directory where you "
                                 "would like all outputs/results saved.")
        parser.add_argument("--image_folder",
                            type=str,
                            required=True,
                            default=None,
                            help="The path to the folder within which the "
                                 "simgularity image & 'wp3.simg' "
                                 "'limixDec22.simg' is located")
        parser.add_argument("--top_dir",
                            nargs="*",
                            type=str,
                            default=["/groups/umcg-biogen/tmp01/",
                                     "/groups/umcg-biogen/tmp01/umcg-mvochteloo/simulated_home:/home/umcg-mvochteloo"],
                            help="List of paths to bind to Singularity. "
                                 "Default: ['/groups/umcg-biogen/tmp01/']")

        # arguments for step 1.
        parser.add_argument("--data_prep_config",
                            type=str,
                            default="wp3_input_create.yaml",
                            help="")
        parser.add_argument("--wg1_folder",
                            type=str,
                            default=None,
                            help="The path to the outputs from WG1")
        parser.add_argument("--wg2_folder",
                            type=str,
                            default=None,
                            help="The path to the outputs from WG2")
        parser.add_argument("--unimputed_folder",
                            type=str,
                            default=None,
                            help="The path, and base filename to the "
                                 "unimputed genotype file. These are created "
                                 "in WG1 and are in plink2 formats.")
        parser.add_argument("--ancestry_split_vcf",
                            type=str,
                            default=None,
                            help="Relative path, relative to WG1_folder, "
                                 "where the ancestry split imputation "
                                 "results are stored. ")
        parser.add_argument("--wg1_psam",
                            type=str,
                            default=None,
                            help="Relative path, relative to WG1 folder, "
                                 "psam file with donor information including "
                                 "both the user provided information and "
                                 "updates during the WG1 pipeline. "
                                 "(https://wg1-pipeline-qc.readthedocs.io/en/latest/Imputation/Imputation_Required_Input.html#plink2-reference-snp-genotype-pfiles)")
        parser.add_argument("--gtf_annotation_file",
                            type=str,
                            default=None,
                            help="GTF file path used when quantifying "
                                 "gene-expression. ")
        parser.add_argument("--cell_annotation",
                            type=str,
                            default=None,
                            help="Relative path, relative to WG1 folder, "
                                 "psam file with donor information including "
                                 "both the user provided information and "
                                 "updates during the WG1 pipeline. "
                                 "(https://wg1-pipeline-qc.readthedocs.io/en/latest/Imputation/Imputation_Required_Input.html#plink2-reference-snp-genotype-pfiles)")
        parser.add_argument("--wg1_singlets_assigned",
                            type=str,
                            default="demultiplexing/output/QC_figures/seurat_object_all_pools_singlet_barcodes_final_assignments.rds",
                            help="")
        parser.add_argument("--genotype_dir",
                            type=str,
                            default="genotype_input/",
                            help="")
        parser.add_argument("--wg1_wg2_qc_taged",
                            type=str,
                            default="WG1_WG2_summary/qc_tag.rds",
                            help="")
        parser.add_argument("--number_of_genes_in_chunk",
                            type=int,
                            default=500,
                            help="")

        # argument for step 2.
        parser.add_argument("--qtl_config",
                            type=str,
                            default="Qtl_wp3.yaml",
                            help="")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        dataset_work_dir = os.path.join(self.work_dir, self.dataset_outdir)
        if not os.path.exists(dataset_work_dir):
            os.makedirs(dataset_work_dir)

        if self.step is None or self.step == 1:
            print("Generating step 1 files")

            output_dir = os.path.join(dataset_work_dir, "Step1-DataPreparation")
            log_dir = os.path.join(output_dir, "log")
            for outdir in [output_dir, log_dir]:
                if not os.path.exists(outdir):
                    os.makedirs(outdir)

            # Fix some dumb shit
            if not self.wg1_folder.endswith("/"):
                self.wg1_folder += "/"
            if not self.wg2_folder.endswith("/"):
                self.wg2_folder += "/"
            if not dataset_work_dir.endswith("/"):
                dataset_work_dir += "/"
            if self.ancestry_split_vcf.startswith("/"):
                self.ancestry_split_vcf = self.ancestry_split_vcf[1:]
            if not self.ancestry_split_vcf.endswith("/"):
                self.ancestry_split_vcf += "/"
            if self.wg1_singlets_assigned.startswith("/"):
                self.wg1_singlets_assigned = self.wg1_singlets_assigned[1:]
            if self.genotype_dir.startswith("/"):
                self.genotype_dir = self.genotype_dir[1:]
            if not self.genotype_dir.endswith("/"):
                self.genotype_dir += "/"
            if self.wg1_wg2_qc_taged.startswith("/"):
                self.wg1_wg2_qc_taged = self.wg1_wg2_qc_taged[1:]

            self.generate_step1_files(
                wg3_folder=dataset_work_dir,
                output_dir=output_dir,
                log_dir=log_dir,
                snakefile=os.path.join(self.image_folder, "Create_wp3_inputs.smk"),
                configtemplate=os.path.join(self.image_folder, self.data_prep_config),
                configfile=os.path.join(output_dir, self.data_prep_config)
            )

            print("")

        if self.step is None or self.step == 2:
            print("Generating step 2 files")

            output_dir = os.path.join(dataset_work_dir, "Step2-QTLMapping")
            log_dir = os.path.join(output_dir, "log")
            for outdir in [output_dir, log_dir]:
                if not os.path.exists(outdir):
                    os.makedirs(outdir)

            self.generate_step2_files(
                wg3_folder=dataset_work_dir,
                output_dir=output_dir,
                log_dir=log_dir,
                snakefile=os.path.join(self.image_folder, "Qtl_Snakefile"),
                configtemplate=os.path.join(self.image_folder, self.qtl_config),
                configfile=os.path.join(output_dir, self.qtl_config)
            )

            print("")

    def generate_step1_files(self, wg3_folder, output_dir, log_dir, snakefile,
                             configtemplate, configfile):
        """
        """
        config_arguments = (
            ("image_folder", self.image_folder),
            ("top_dir", self.top_dir),
            ("WG1_folder", self.wg1_folder),
            ("WG2_folder", self.wg2_folder),
            ("WG3_folder", wg3_folder),
            ("unimputed_folder", self.unimputed_folder),
            ("ancestry_split_vcf", self.ancestry_split_vcf),
            ("wg1_psam", self.wg1_psam),
            ("gtf_annotation_file", self.gtf_annotation_file),
            ("cellAnnotation", self.cell_annotation),
            ("wg1SingletsAssigned", self.wg1_singlets_assigned),
            ("genotype_dir", self.genotype_dir),
            ("wg1_wg2_qc_taged", self.wg1_wg2_qc_taged),
            ("number_of_genes_in_chunk", self.number_of_genes_in_chunk)
        )
        self.write_configfile(
            template=configtemplate,
            arguments=config_arguments,
            outpath=configfile
        )

        self.write_dry_run_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir
        )

        self.write_unlock_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir
        )

        self.write_build_dag_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir,
            outfile="dag1"
        )

    def generate_step2_files(self, wg3_folder, output_dir, log_dir, snakefile,
                             configtemplate, configfile):
        """
        """
        config_arguments = (
            ("image_folder", self.image_folder),
            ("top_dir", self.top_dir),
            ("WG3_folder", wg3_folder),
            ("out_folder", os.path.join(wg3_folder, "output"))
        )
        self.write_configfile(
            template=configtemplate,
            arguments=config_arguments,
            outpath=configfile
        )

        self.write_dry_run_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir
        )

        self.write_unlock_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir
        )

        self.write_build_dag_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir,
            outfile="dag1"
        )

    def write_configfile(self, template, arguments, outpath):
        yaml_lines = []
        for line in open(template, 'r'):
            for label, argument in arguments:
                if line.startswith(label):
                    line = "{}: {}\n".format(label, argument)
                line = line.replace("\n", "")
            yaml_lines.append(line)

        self.write_lines_to_file(
            lines=yaml_lines,
            path=outpath
        )

    def write_build_dag_script(self, snakefile, configfile, output_dir,
                               outfile="dag"):
        lines = [
            '#!/bin/bash',
            '',
            'snakemake \\',
            '    --snakefile {} \\'.format(snakefile),
            '    --configfile {} \\'.format(configfile),
            '    --dag | \\',
            '    dot -Tsvg \\',
            '        > {}.svg'.format(outfile)
        ]
        self.write_lines_to_file(
            lines=lines,
            path=os.path.join(output_dir, "build_{}.sh".format(outfile))
        )

    def write_dry_run_script(self, snakefile, configfile, output_dir, cores=1):
        lines = [
            '#!/bin/bash',
            '',
            'snakemake \\',
            '    --snakefile {} \\'.format(snakefile),
            '    --configfile {} \\'.format(configfile),
            '    --dryrun \\',
            '    --cores {} \\'.format(cores),
            '    --reason'
        ]
        self.write_lines_to_file(
            lines=lines,
            path=os.path.join(output_dir, "dry_run.sh")
        )

    def write_unlock_script(self, snakefile, configfile, output_dir):
        lines = [
            '#!/bin/bash',
            '',
            'snakemake \\',
            '    --snakefile {} \\'.format(snakefile),
            '    --configfile {} \\'.format(configfile),
            '    --unlock'
        ]
        self.write_lines_to_file(
            lines=lines,
            path=os.path.join(output_dir, "unlock.sh")
        )

    def write_run_script(self, snakefile, configfile, output_dir, log_dir,
                         until=None, jobs=1, restart_times=2, latency_wait=30,
                         nodes=1, time="short", outfile="run"):
        time_dict = {
            "short": "05:59:00",
            "medium": "23:59:00",
            "long": "6-23:59:00"
        }

        if time not in time_dict.keys():
            print("Error, time not recognized.")
            exit()

        lines = [
            '#!/bin/bash',
            '',
            'nohup \\',
            '  snakemake \\',
            '    --snakefile {} \\'.format(snakefile),
            '    --configfile {} \\'.format(configfile),
            '    --rerun-incomplete \\',
            '    --jobs {} \\'.format(jobs),
            '    --use-singularity \\',
            '    --restart-times {} \\'.format(restart_times),
            '    --keep-going \\',
            '    --latency-wait {} \\'.format(latency_wait),
            '    --cluster \\',
            '       "sbatch \\',
            '       --nodes={} \\'.format(nodes),
            '       --cpus-per-task={threads} \\',
            '       --mem=\$(({resources.mem_per_thread_gb} * {threads}))G \\',
            '       --tmp=\$(({resources.disk_per_thread_gb} * {threads}))G \\',
            '       --time={} \\'.format(time_dict[time]),
            '       --output={} \\'.format(os.path.join(log_dir, '{rule}.out')),
            '       --error={} \\'.format(os.path.join(log_dir, '{rule}.out')),
            '       --export=ALL \\',
            '       --qos=regular" \\',
            '    > {}/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &'.format(log_dir),
            '',
            'echo "Check status of command with:" ps -p $! -u'
        ]

        if until is not None:
            lines.insert(6, '    --until {} \\'.format(until))

        self.write_lines_to_file(
            lines=lines,
            path=os.path.join(output_dir, "{}.sh".format(outfile))
        )

    @staticmethod
    def write_lines_to_file(lines, path):
        with open(path, "w") as f:
            for line in lines:
                f.write(line + "\n")
        f.close()
        print("\tSaved file: {}".format(os.path.basename(path)))

    def print_arguments(self):
        print("General arguments")
        print("  > Step:                              {}".format("1 & 2" if self.step is None else str(self.step)))
        print("  > Working directory:                 {}".format(self.work_dir))
        print("  > Dataset output directory:          {}".format(self.dataset_outdir))
        print("  > Image folder:                      {}".format(self.image_folder))
        print("  > Top dir:                           {}".format(self.top_dir))
        print("")

        if self.step is None or self.step == 1:
            print("[Step 1] Data preparation arguments:")
            print("  > Configuration file:                {}".format(self.data_prep_config))
            print("  > Workgroup 1 folder:                {}".format(self.wg1_folder))
            print("  > Workgroup 2 folder:                {}".format(self.wg2_folder))
            print("  > Unimputed genotype folder:         {}".format(self.unimputed_folder))
            print("  > Ancestry splitted VCF:             {}".format(self.ancestry_split_vcf))
            print("  > Workgroup 1 PSAM:                  {}".format(self.wg1_psam))
            print("  > GTF annotation file:               {}".format(self.gtf_annotation_file))
            print("  > Cell annotation:                   {}".format(self.cell_annotation))
            print("  > Workgroup 1 signlets assigned:     {}".format(self.wg1_singlets_assigned))
            print("  > Genotype folder:                   {}".format(self.genotype_dir))
            print("  > Workgroup1 Workgroup2 QC taged:    {}".format(self.wg1_wg2_qc_taged))
            print("  > Number of genes in chunk:          {}".format(self.number_of_genes_in_chunk))
            print("")

        if self.step is None or self.step == 2:
            print("[Step2] QTL arguments:")
            print("  > Configuration file:                {}".format(self.qtl_config))
            print("")


if __name__ == '__main__':
    m = main()
    m.start()