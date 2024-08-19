#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse
import logging
import subprocess
import os
import glob
from datetime import datetime


"""
Syntax: 
./bam_to_cram.py -h
"""

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--bind", type=str, required=True, help="")
parser.add_argument("--singularity", type=str, required=True, help="")
parser.add_argument("--work_dir", type=str, required=True, help="")
parser.add_argument("--pattern", type=str, required=False, default="CellRanger/*/outs/possorted_genome_bam.bam", help="")
parser.add_argument("--reference", type=str, required=True, help="")
parser.add_argument("--dry_run", action='store_true', help="Print all info.")
args = parser.parse_args()

if args.pattern.startswith("/"):
    args.pattern = args.pattern[1:]

def get_logger():
    handlers = []

    # Construct stream handler.
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(logging.Formatter('%(asctime)s  [%(levelname)-4.4s]  %(message)s', "%H:%M:%S"))
    handlers.append(stream_handler)

    if not args.dry_run:
        # Construct file handler.
        file_handler = logging.FileHandler(os.path.join(args.work_dir, "{}_bam_to_cram.log".format(datetime.now().strftime("%Y-%m-%d"))))
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter('%(asctime)s  %(module)-16s %(levelname)-8s %(message)s', "%Y-%m-%d %H:%M:%S"))
        handlers.append(file_handler)

    # Construct logger object.
    logging.basicConfig(
        level=logging.DEBUG,
        handlers=handlers)

    return logging.getLogger(__name__)

LOGGER = get_logger()
LOGGER.info("bam_to_cram.py")
LOGGER.info("Starting program at: {}".format(datetime.now().strftime("%d-%m-%Y, %H:%M:%S")))

LOGGER.info("Options in effect:")
for arg in vars(args):
    LOGGER.info("  --{} {}".format(arg, getattr(args, arg)))
LOGGER.info("")

for subdir in args.bind.split(","):
    if not os.path.exists(subdir):
        LOGGER.error("Error, --bind '{}' does not exist.".format(subdir))
        exit()

if not os.path.exists(args.singularity):
    LOGGER.error("Error, --singularity does not exist.")
    exit()

if not os.path.exists(args.work_dir):
    LOGGER.error("Error, --work_dir does not exist.")
    exit()

if not os.path.exists(args.reference):
    LOGGER.error("Error, --reference does not exist.")
    exit()

def run_command(command):
    LOGGER.info("\tRunning command: {}".format(" ".join(command)))
    if args.dry_run:
        return True

    try:
        output = subprocess.check_output(command, stderr=subprocess.STDOUT).decode()
        LOGGER.info("\toutput: {}".format(output.rstrip("\n")))
    except subprocess.CalledProcessError as e:
        LOGGER.warning("  Subprocess cannot read some directories.")
        LOGGER.warning("\toutput: {}".format(e.output.decode().rstrip("\n")))
    except Exception as e:
        LOGGER.error("  Subprocess failed.")
        LOGGER.error("\toutput: {}".format(str(e).rstrip("\n")))

LOGGER.info("Saving md5sum of reference")
run_command(["md5sum", args.reference])
LOGGER.info("")

LOGGER.info("Searching bam files")
bam_fpaths = glob.glob(os.path.join(args.work_dir, args.pattern))

for i, bam_fpath in enumerate(bam_fpaths):
    LOGGER.info("{:,} / {:,} Processing: {}".format(i, len(bam_fpaths), bam_fpath))
    if not bam_fpath.endswith(".bam"):
        LOGGER.error("Expected file to end with '.bam'. Exciting.")
        exit()

    cram_fpath = bam_fpath[:-4] + ".cram"
    if os.path.exists(cram_fpath):
        LOGGER.warning("\tCram output file already exists. Skipping file.")
        LOGGER.info("")
        continue

    run_command(["md5sum", bam_fpath])
    run_command(["singularity", "exec", "--bind", args.bind, args.singularity, "samtools", "view", "-T", args.reference, "-C", "-o", cram_fpath, bam_fpath])
    run_command(["md5sum", cram_fpath])
    LOGGER.info("")

LOGGER.info("Finished program at: {}".format(datetime.now().strftime("%d-%m-%Y, %H:%M:%S")))