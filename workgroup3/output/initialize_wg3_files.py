#!/usr/bin/env python3

"""
File:         initialize_wg3_files.py
Created:      2023/04/11
Last Changed: 2023/12/06
Author:       M.Vochteloo

Copyright (C) 2022 University Medical Center Groningen.
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
./initialize_wg3_files.py -h
"""

class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.work_dir = getattr(arguments, 'work_dir')
        dataset_outdir = getattr(arguments, 'dataset_outdir')
        self.down_analyses_dir = getattr(arguments, 'down_analyses_dir')
        self.down_analyses_config = getattr(arguments, 'down_analyses_config')
        self.bind_path = getattr(arguments, 'bind_path')
        self.sif_path = getattr(arguments, 'sif_path')
        self.limix_sif_path = getattr(arguments, 'limix_sif_path')
        self.repo_dir = getattr(arguments, 'repo_dir')
        self.poolsheet_path = getattr(arguments, 'poolsheet_path')
        self.wg1_genotype_folder = getattr(arguments, 'wg1_genotype_folder')
        self.wg1_demultiplex_folder = getattr(arguments, 'wg1_demultiplex_folder')
        self.wg2_folder = getattr(arguments, 'wg2_folder')
        self.cell_annotation_file = getattr(arguments, 'cell_annotation_file')
        self.ref_dir = getattr(arguments, 'ref_dir')
        self.gtf_annotation_file = getattr(arguments, 'gtf_annotation_file')
        self.ancestry = getattr(arguments, 'ancestry')
        self.cell_level = getattr(arguments, 'cell_level')
        self.cell_types = getattr(arguments, 'cell_types')
        self.genome_build = getattr(arguments, 'genome_build')
        self.filter_samples = getattr(arguments, 'filter_samples')
        self.save_all_samples = getattr(arguments, 'no_save_all_samples')
        self.calculate_qtl = getattr(arguments, 'no_qtl')
        self.output_flat_qtl_results = getattr(arguments, 'output_flat_qtl_results')
        self.compress_qtl = getattr(arguments, 'no_qtl_compression')
        self.calculate_ld = getattr(arguments, 'no_ld')
        self.compress_ld = getattr(arguments, 'no_ld_compression')
        self.relative_wg1_imputed_genotype_vcf = getattr(arguments, 'relative_wg1_imputed_genotype_vcf')
        self.relative_wg1_psam = getattr(arguments, 'relative_wg1_psam')
        self.wg2_pairing = getattr(arguments, 'wg2_pairing')
        self.eqtl_chunks_n_genes = getattr(arguments, 'eqtl_chunks_n_genes')
        self.n_expression_pcs = "[" + ".join(getattr(arguments, 'n_expression_pcs')" + "]"
        self.exclude_temp_in_sbatch = getattr(arguments, 'exclude_temp_in_sbatch')

        # Pre-process the dataset output directory.
        date_str = datetime.now().strftime("%Y-%m-%d")
        if not re.match("^[0-9]{4}-(0[1-9]|1[0-2])-(0[1-9]|[1-2][0-9]|3[0-1])", dataset_outdir):
            dataset_outdir = "{}-{}".format(date_str, dataset_outdir)
        self.dataset_outdir = dataset_outdir

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
        parser.add_argument("--down_analyses_dir",
                            type=str,
                            required=True,
                            default=None,
                            help=".")
        parser.add_argument("--down_analyses_config",
                            type=str,
                            default="sceQTL-Gen_DA.yaml",
                            help="")
        parser.add_argument("--bind_path",
                            type=str,
                            required=True,
                            default=None,
                            help="List of paths to bind to Singularity.")
        parser.add_argument("--sif_path",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--limix_sif_path",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--repo_dir",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--poolsheet_path",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--wg1_genotype_folder",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--wg1_demultiplex_folder",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--wg2_folder",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--cell_annotation_file",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--ref_dir",
                            type=str,
                            required=True,
                            default=None,
                            help="This is the path to the directory containing "
                                 "the references files provided for the "
                                 "downstream analyses")
        parser.add_argument("--gtf_annotation_file",
                            type=str,
                            required=True,
                            default=None,
                            help="")
        parser.add_argument("--ancestry",
                            type=str,
                            required=False,
                            default="EUR",
                            help="")
        parser.add_argument("--cell_level",
                            type=str,
                            required=False,
                            default="L1",
                            help="")
        parser.add_argument("--cell_types",
                            nargs="*",
                            type=str,
                            required=False,
                            default=["AST", "END", "EX", "IN", "MIC", "OLI", "OPC", "PER"],
                            help="")
        parser.add_argument("--genome_build",
                            type=str,
                            required=False,
                            default="hg19",
                            help="")
        parser.add_argument("--filter_samples",
                            action='store_false',
                            help="")
        parser.add_argument("--no_save_all_samples",
                            action='store_false',
                            help="")
        parser.add_argument("--no_qtl",
                            action='store_false',
                            help="")
        parser.add_argument("--output_flat_qtl_results",
                            action='store_true',
                            help="")
        parser.add_argument("--no_qtl_compression",
                            action='store_false',
                            help="")
        parser.add_argument("--no_ld",
                            action='store_false',
                            help="")
        parser.add_argument("--no_ld_compression",
                            action='store_false',
                            help="")
        parser.add_argument("--relative_wg1_imputed_genotype_vcf",
                            type=str,
                            required=False,
                            default="vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz",
                            help="")
        parser.add_argument("--relative_wg1_psam",
                            type=str,
                            required=False,
                            default="update_sex_ancestry/update_sex.psam",
                            help="")
        parser.add_argument("--wg2_pairing",
                            type=str,
                            required=False,
                            default="azimuth_l1_brain.csv",
                            help="")
        parser.add_argument("--eqtl_chunks_n_genes",
                            type=int,
                            required=False,
                            default=150,
                            help="")
        parser.add_argument("--n_expression_pcs",
                            nargs="*",
                            type=int,
                            required=False,
                            default=[0, 2, 4, 6, 8, 10],
                            help="")
        parser.add_argument("--exclude_temp_in_sbatch",
                            action='store_true',
                            help="")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Generating files")

        output_dir = os.path.join(self.work_dir, self.dataset_outdir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        log_dir = os.path.join(output_dir, "log")
        slurm_log_dir = os.path.join(output_dir, "slurm_log")
        for outdir in [output_dir, log_dir, slurm_log_dir]:
            if not os.path.exists(outdir):
                os.makedirs(outdir)
        snakefile = os.path.join(self.down_analyses_dir, "Snakefile")
        configfile = os.path.join(output_dir, self.down_analyses_config)

        config_arguments = (
            ("bind_path", self.bind_path),
            ("singularity_image", self.sif_path),
            ("repo_dir", os.path.join(self.down_analyses_dir)),
            ("singularity_image", self.sif_path),
            ("limix_singularity_image", self.limix_sif_path),
            ("repo_dir", self.repo_dir),
            ("poolsheet_path", os.path.join(self.poolsheet_path)),
            ("wg1_genotype_folder", self.wg1_genotype_folder),
            ("wg1_demultiplex_folder", self.wg1_demultiplex_folder),
            ("wg2_folder", self.wg2_folder),
            ("cell_annotation_file", self.cell_annotation_file),
            ("ref_dir", self.ref_dir),
            ("gtf_annotation_file", self.gtf_annotation_file),
            ("output_dir", output_dir),
            ("ancestry", self.ancestry),
            ("cell_level", self.cell_level),
            ("cell_types", self.cell_types),
            ("genome_build", self.genome_build),
            ("filter_samples", self.filter_samples),
            ("save_all_samples", self.save_all_samples),
            ("calculate_qtl", self.calculate_qtl),
            ("output_flat_qtl_results", self.output_flat_qtl_results),
            ("compress_qtl", self.compress_qtl),
            ("calculate_ld", self.calculate_ld),
            ("compress_ld", self.compress_ld),
            ("relative_wg1_imputed_genotype_vcf", self.relative_wg1_imputed_genotype_vcf),
            ("relative_wg1_psam", self.relative_wg1_psam),
            ("wg2_pairing", self.wg2_pairing),
            ("eqtl_chunks_n_genes", self.eqtl_chunks_n_genes),
            ("n_expression_pcs", self.n_expression_pcs)
        )
        self.write_configfile(
            template=os.path.join(self.down_analyses_dir, self.down_analyses_config),
            arguments=config_arguments,
            outpath=configfile
        )

        cluster_status_script = self.write_cluster_status_script(
            output_dir=output_dir
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

        for i in range(1, 3):
            self.write_build_dag_script(
                snakefile=snakefile,
                configfile=configfile,
                output_dir=output_dir,
                outfile="dag{}".format(i)
            )

        for jobs in [1, 2, 22, 100, 1000, 10000]:
            self.write_run_script(
                snakefile=snakefile,
                configfile=configfile,
                output_dir=output_dir,
                log_dir=log_dir,
                slurm_log_dir=slurm_log_dir,
                cluster_status_script=cluster_status_script,
                jobs=jobs,
                outfile="run_{}jobs".format(jobs)
            )

        self.write_run_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir,
            log_dir=log_dir,
            slurm_log_dir=slurm_log_dir,
            cluster_status_script=cluster_status_script,
            jobs=3,
            restart_times=0,
            qos="priority",
            outfile="run_dev"
        )

        self.write_run_local_script(
            snakefile=snakefile,
            configfile=configfile,
            output_dir=output_dir
        )

        self.write_report_script(
            snakefile=snakefile,
            configfile=configfile,
            report=os.path.join(output_dir, "downstream_analyses_report.html"),
            output_dir=output_dir
        )

    def write_configfile(self, template, arguments, outpath):
        yaml_lines = []
        for line in open(template, 'r'):
            line = line.replace("\n", "")
            for label, argument in arguments:
                if line.startswith("  {}:".format(label)):
                    line = "  {}: {}".format(label, argument)
            yaml_lines.append(line)

        self.write_lines_to_file(
            lines=yaml_lines,
            path=outpath
        )

    def write_cluster_status_script(self, output_dir):
        """
        https://www.embl.org/groups/bioinformatics-rome/blog/2022/05/snakemake-profile-5-handling-memory-and-timeout-errors/
        """
        outfile = os.path.join(output_dir, "status-sacct.sh")

        lines = [
            '#!/usr/bin/env bash',
            '# Check status of Slurm job',
            'jobid="$1"',
            'if [[ "$jobid" == Submitted ]]',
            'then',
            '  echo smk-simple-slurm: Invalid job ID: "$jobid" >&2',
            '  echo smk-simple-slurm: Did you remember to add the flag --parsable to your sbatch call? >&2',
            '  exit 1',
            'fi',
            'output=`sacct -j "$jobid" --format State --noheader | head -n 1 | awk \'{print $1}\'`',
            'if [[ $output =~ ^(COMPLETED).* ]]',
            'then',
            '  echo success',
            'elif [[ $output =~ ^(RUNNING|PENDING|COMPLETING|CONFIGURING|SUSPENDED).* ]]',
            'then',
            '  echo running',
            'else',
            '  echo failed',
            'fi']

        self.write_lines_to_file(
            lines=lines,
            path=outfile
        )

        return outfile

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
            '        > {}.svg'.format(os.path.join(output_dir, outfile))
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
                         slurm_log_dir, cluster_status_script, nodes=1,
                         jobs=1, restart_times=0, latency_wait=60,
                         qos="regular", until=None, outfile="run"):
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
            '       --job-name=snakemake_{rule}_{wildcards} \\',
            '       --nodes={} \\'.format(nodes),
            '       --cpus-per-task={threads} \\',
            '       --mem=\$(({resources.mem_per_thread_gb} * {threads}))G \\',
            '       --tmp=\$(({resources.disk_per_thread_gb} * {threads}))G \\',
            '       --time={resources.time} \\',
            '       --output={}'.format(slurm_log_dir) + '/{rule}_{wildcards}.out \\',
            '       --error={}'.format(slurm_log_dir) + '/{rule}_{wildcards}.out \\',
            '       --export=NONE \\',
            '       --qos={} \\'.format(qos),
            '       --parsable" \\',
            '    --cluster-status {} \\'.format(cluster_status_script),
            '    > {}/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &'.format(log_dir),
            '',
            'echo "Check status of command with:" ps -p $! -u'
        ]

        if self.exclude_temp_in_sbatch:
            lines.remove('       --tmp=\$(({resources.disk_per_thread_gb} * {threads}))G \\')

        if until is not None:
            lines.insert(6, '    --until {} \\'.format(until))

        self.write_lines_to_file(
            lines=lines,
            path=os.path.join(output_dir, "{}.sh".format(outfile))
        )

    def write_run_local_script(self, snakefile, configfile, output_dir):
        lines = [
            '#!/bin/bash',
            '',
            'snakemake \\',
            '  --snakefile {} \\'.format(snakefile),
            '  --configfile {} \\'.format(configfile),
            '  --cores 1'
        ]
        self.write_lines_to_file(
            lines=lines,
            path=os.path.join(output_dir, "run_local.sh")
        )

    def write_report_script(self, snakefile, configfile, report, output_dir):
        lines = [
            '#!/bin/bash',
            '',
            'snakemake \\',
            '  --snakefile {} \\'.format(snakefile),
            '  --configfile {} \\'.format(configfile),
            '  --report {}'.format(report)
        ]
        self.write_lines_to_file(
            lines=lines,
            path=os.path.join(output_dir, "report.sh")
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
        print("  > Working directory:              {}".format(self.work_dir))
        print("  > Dataset output directory:       {}".format(self.dataset_outdir))
        print("  > Down. Analyses directory:       {}".format(self.down_analyses_dir))
        print("  > Down. Analyses configuration:   {}".format(self.down_analyses_config))
        print("  > Bind path:                      {}".format(self.bind_path))
        print("  > Singularity path:               {}".format(self.sif_path))
        print("  > Limix singularity path:         {}".format(self.limix_sif_path))
        print("  > Repository folder:              {}".format(self.repo_dir))
        print("  > Poolsheet path:                 {}".format(self.poolsheet_path))
        print("  > Workgroup 1 genotype folder:    {}".format(self.wg1_genotype_folder))
        print("  > Workgroup 1 demultiplex folder: {}".format(self.wg1_demultiplex_folder))
        print("  > Workgroup 2 folder:             {}".format(self.wg2_folder))
        print("  > Cell annot. file:               {}".format(self.cell_annotation_file))
        print("  > Reference directory:            {}".format(self.ref_dir))
        print("  > GTF annot. file:                {}".format(self.gtf_annotation_file))
        print("  > Ancestry:                       {}".format(self.ancestry))
        print("  > Cell level:                     {}".format(self.cell_level))
        print("  > Cell types:                     {}".format(", ".join(self.cell_types)))
        print("  > Genome build:                   {}".format(self.genome_build))
        print("  > Filter samples:                 {}".format(self.filter_samples))
        print("  > Save all samples:               {}".format(self.save_all_samples))
        print("  > Calculate QTL:                  {}".format(self.calculate_qtl))
        print("  > Output flat QTL results:        {}".format(self.output_flat_qtl_results))
        print("  > Compress QTL:                   {}".format(self.compress_qtl))
        print("  > Calculate LD:                   {}".format(self.calculate_ld))
        print("  > Compress LD:                    {}".format(self.compress_ld))
        print("  > Relative WG1 imp. geno. VCF:    {}".format(self.relative_wg1_imputed_genotype_vcf))
        print("  > Relative WG1 PSAM:              {}".format(self.relative_wg1_psam))
        print("  > WG2 pairing:                    {}".format(self.wg2_pairing))
        print("  > eQTL chunks N genes:            {}".format(self.eqtl_chunks_n_genes))
        print("  > N expression PCs:               {}".format(self.n_expression_pcs))
        print("  > Exclude TEMP in SBATCH:         {}".format(self.exclude_temp_in_sbatch))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()