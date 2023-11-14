#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--indir", type=str, required=True, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import glob
import os
import numpy as np
import pandas as pd

def parse_value(line, pos, fun=None):
    str = line[pos[0]:pos[1]].replace(" ", "")
    if str == "":
        return np.nan
    if fun is not None:
        return fun(str)
    return str

def time_to_min(value):
    hour, min, sec = value.split(":")
    return (int(hour) * 60) + int(min) + (int(sec) / 60)

def size_to_gb(value):
    if "K" in value:
        return float(value.replace("K", "")) / 1e6
    elif "M" in value:
        return float(value.replace("M", "")) / 1e3
    else:
        return float(value)

data = []
for fpath in glob.glob(os.path.join(args.indir, "Step1-Imputation", "slurm_log", "*")):
    filename = os.path.basename(fpath).split(".")[0]
    chr = None
    if "chr" in filename:
        filename, chr = filename.split("_chr=")

    flag = False
    indices = {}
    with open(fpath, 'r') as f:
        for line in f:
            if line.startswith("JobID"):
                indices = {}

                start_of_word = 0
                prev_char = None
                for pos, char in enumerate(line):
                    if prev_char == " " and char != " ":
                        # end of word.
                        word = line[start_of_word:pos]
                        indices[word.replace(" ", "")] = [start_of_word, pos]
                        start_of_word = pos

                    prev_char = char

                flag = True
                continue

            if flag:
                job_id = parse_value(line, indices["JobID"])
                elapsed = parse_value(line, indices["Elapsed"], time_to_min)
                allo_cpus = parse_value(line, indices["AllocCPUS"], int)
                ave_cpu = parse_value(line, indices["AveCPU"], time_to_min)
                reg_mem = parse_value(line, indices["ReqMem"], float)
                max_vm_size = parse_value(line, indices["MaxVMSize"], size_to_gb)
                max_rss = parse_value(line, indices["MaxRSS"], size_to_gb)
                max_disk_read = parse_value(line, indices["MaxDiskRead"], size_to_gb)
                max_disk_write = parse_value(line, indices["MaxDiskWrite"], size_to_gb)

                data.append([filename, chr, job_id, elapsed, allo_cpus, ave_cpu, reg_mem, max_vm_size, max_rss, max_disk_read, max_disk_write])
                flag = False
    f.close()

df = pd.DataFrame(data, columns=["Filename", "CHR", "JobID", "Elapsed", "AllocCPUS", "AveCPU", "ReqMem", "MaxVMSize", "MaxRSS", "MaxDiskRead", "MaxDiskWrite"])
print(df)
print("")

df.drop(["CHR", "JobID"], axis=1, inplace=True)

print("Max per rule:")
print(df.groupby("Filename").max())
print("")

# print("Min per rule:")
# print(df.groupby("Filename").min())
# print("")
#
# print("Mean per rule:")
# print(df.groupby("Filename").mean())
# print("")