#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--indir", type=str, required=True, help="")
parser.add_argument("--poolsheet", type=str, required=False, default=None, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import glob
import os
import numpy as np
import pandas as pd

ncells = {}
if args.poolsheet is not None:
    print("Counting cells")
    poolsheet = pd.read_csv(args.poolsheet, sep="\t")
    for index, row in poolsheet.iterrows():
        metrics_summary = row["Bam"].replace("/possorted_genome_bam.bam", "/metrics_summary.csv")
        metrics_df = pd.read_csv(metrics_summary, sep=",", header=0, index_col=None)
        ncells[row["Pool"]] = int(str(metrics_df["Estimated Number of Cells"][0]).replace(",", ""))


def parse_value(line, pos, fun=None):
    str = line[pos[0]:pos[1]].replace(" ", "")
    if str == "":
        return np.nan
    if fun is not None:
        return fun(str)
    return str


def time_to_min(value):
    days = 0
    time = value
    if "-" in value:
        days, time = value.split("-")
    hour, min, sec = time.split(":")
    return (int(days) * 24 * 60) + (int(hour) * 60) + int(min) + (int(sec) / 60)


def size_to_gb(value):
    if "K" in value:
        return float(value.replace("K", "")) / 1e6
    elif "M" in value:
        return float(value.replace("M", "")) / 1e3
    elif "G" in value:
        return float(value.replace("G", ""))
    else:
        return float(value)

print("Parsing slurm log files")
data = []
for fpath in glob.glob(os.path.join(args.indir, "slurm_log", "*")):
    flag = False
    indices = {}
    with open(fpath, 'r') as f:
        for line in f:
            if line.startswith("rule "):
                rule = line.strip("\n").replace("rule ", "").rstrip(":")
            if line.startswith("    wildcards: "):
                wildcards = line.strip("\n").replace("    wildcards: ", "")

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

                data.append([rule, wildcards, job_id, elapsed, allo_cpus, ave_cpu, reg_mem, max_vm_size, max_rss, max_disk_read, max_disk_write])
                flag = False
    f.close()

print("Memory usage: ")
df = pd.DataFrame(data, columns=["Rule", "Wildcards", "JobID", "Elapsed", "AllocCPUS", "AveCPU", "ReqMem", "MaxVMSize", "MaxRSS", "MaxDiskRead", "MaxDiskWrite"])
df.sort_values(by="Elapsed", ascending=False, inplace=True)
print(df.iloc[0, 0], df.iloc[0, 1], df.iloc[0, 3])
if args.poolsheet is not None:
    ncells_values = []
    for _, row in df.iterrows():
        for key, value in ncells.items():
            if key in row["Wildcards"]:
                ncells_values.append(value)
    df["NCells"] = ncells_values
print(df)
# df.to_excel("memory_usage.xlsx")
# print("")

print(df["Rule"].value_counts())

df.drop(["Wildcards", "JobID"], axis=1, inplace=True)

print("Max per rule:")
print(df.groupby("Rule").max())
print("")

print("Min per rule:")
print(df.groupby("Rule").min())
print("")

print("Mean per rule:")
print(df.groupby("Rule").mean())
print("")

print("Sum per rule:")
print(df.groupby("Rule").sum())
print("")
