#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--logfile", type=str, required=True, help="")
parser.add_argument("--slurm_log", type=str, required=True, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import numpy as np
import glob
import os

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

slurm_info = {}
for filename in glob.glob(os.path.join(args.slurm_log, "*")):
    oom = False
    oot = False
    incomplete_files = False
    missing_output = False

    flag = False
    indices = {}
    with open(filename, 'r') as f:
        for line in f:
            if "oom-kill" in line:
                oom = True
            elif "DUE TO TIME LIMIT" in line:
                oot = True
            elif "IncompleteFilesException" in line:
                incomplete_files = True
            elif "MissingOutputException" in line:
                missing_output = True
            elif line.startswith("JobID"):
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
                job_id = int(parse_value(line, indices["JobID"]).replace(".batch", ""))
                elapsed = parse_value(line, indices["Elapsed"], time_to_min)
                allo_cpus = parse_value(line, indices["AllocCPUS"], int)
                ave_cpu = parse_value(line, indices["AveCPU"], time_to_min)
                reg_mem = parse_value(line, indices["ReqMem"], float)
                max_vm_size = parse_value(line, indices["MaxVMSize"], size_to_gb)
                max_rss = parse_value(line, indices["MaxRSS"], size_to_gb)
                max_disk_read = parse_value(line, indices["MaxDiskRead"], size_to_gb)
                max_disk_write = parse_value(line, indices["MaxDiskWrite"], size_to_gb)

                slurm_info[job_id] = {"filename": filename, "job_id": job_id, "elapsed": elapsed, "allo_cpus": allo_cpus, "ave_cpu": ave_cpu, "reg_mem": reg_mem, "max_vm_size": max_vm_size, "max_rss": max_rss, "max_disk_read": max_disk_read, "max_disk_write": max_disk_write, "OOM": oom, "OOT": oot, "IncompleteFiles": incomplete_files, "MissingOutput": missing_output}
                flag = False

failed_jobs = {}
with open(args.logfile, 'r') as f:
    error_message = False
    is_shell = False
    save = False

    for line in f:
        line = line.strip("\n")
        if line.replace(" ", "") == "":
            continue

        if line.startswith("Error in rule"):
            rule = line.replace("Error in rule ", "").replace(":", "")
            jobid = np.nan
            output = ""
            log = ""
            shell = []
            cluster_jobid = np.nan
            error_message = True

        if error_message:
            if is_shell:
                if line.startswith("        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)"):
                    is_shell = False
                    continue
                shell.append(line)

            if line.startswith("    jobid: "):
                jobid = int(line.replace("    jobid: ", ""))
            elif line.startswith("    output: "):
                output = line.replace("    output: ", "")
            elif line.startswith("    log: "):
                log = line.replace("    log: ", "").replace(" (check log file(s) for error message)", "")
            elif line.startswith("    shell:"):
                is_shell = True
            elif line.startswith("    cluster_jobid: "):
                cluster_jobid = int(line.replace("    cluster_jobid: ", ""))
                error_message = False
                save = True

        if save:
            slurm_job_info = slurm_info[cluster_jobid]
            slurm_job_info.update({"rule": rule, "jobid": jobid, "output": output, "log": log, "shell": shell, "cluster_jobid":  cluster_jobid})
            failed_jobs[jobid] = slurm_job_info
            save = False
            # exit()

rule_stats = {}
for job_id, failed_job in failed_jobs.items():
    if failed_job["rule"] not in rule_stats:
        rule_stats[failed_job["rule"]] = {"OOT": 0, "OOM": 0, "IncompleteFiles": 0, "MissingOutput": 0, "Other": 0}
    if failed_job["OOT"]:
        rule_stats[failed_job["rule"]]["OOT"] += 1
    if failed_job["OOM"]:
        rule_stats[failed_job["rule"]]["OOM"] += 1
    if failed_job["IncompleteFiles"]:
        rule_stats[failed_job["rule"]]["IncompleteFiles"] += 1
    if failed_job["MissingOutput"]:
        rule_stats[failed_job["rule"]]["MissingOutput"] += 1
    if not failed_job["OOT"] and not failed_job["OOM"] and not failed_job["IncompleteFiles"] and not failed_job["MissingOutput"]:
        rule_stats[failed_job["rule"]]["Other"] += 1

    print(job_id)
    for key, value in failed_job.items():
        print("\t{} = {}".format(key, value))
    print("")

print(rule_stats)