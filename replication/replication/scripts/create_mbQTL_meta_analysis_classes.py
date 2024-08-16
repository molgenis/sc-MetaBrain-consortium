#!/usr/bin/env python
# Author: M. Vochteloo

import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--top_effects", required=True, type=str, help="")
parser.add_argument("--template", required=True, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

import gzip
import json
import re
import copy

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)

def get_datasets_from_column(column):
    return re.search(r'\((.*?)\)', column).group(1).split(";")

print("Loading mbQTL TopEffects file header...")
dataset_columns = {"beta": "DatasetCorrelationCoefficients", "zscore": "DatasetZScores", "N": "DatasetSampleSizes"}
replication_columns = {}
overlapping_datasets = None
fhin = gzopen(args.top_effects, mode='r')
for i, line in enumerate(fhin):
    if i > 0:
        break
    values = line.rstrip("\n").split("\t")
    for value in values:
        for column, colname in dataset_columns.items():
            if value.startswith(colname):
                datasets = set(get_datasets_from_column(value))
                if overlapping_datasets is None:
                    overlapping_datasets = datasets
                overlapping_datasets = overlapping_datasets.intersection(datasets)

                replication_columns[column] = value
fhin.close()

print("\tFound {} datasets: {}".format(len(overlapping_datasets), ", ".join(overlapping_datasets)))
if len(overlapping_datasets) == 0:
    exit()

print("\nLoading class template...")
default_class_settings = json.load(gzopen(args.template))

print("\nCreating dataset class templates...")
numeric_regex = "[0-9.-]+"
for dataset in overlapping_datasets:
    dataset_class_settings = copy.deepcopy(default_class_settings)
    dataset_class_settings["class_name"] = dataset_class_settings["class_name"] + "_" + dataset
    for column, colname in replication_columns.items():
        dataset_order = {dataset: i for i, dataset in enumerate(get_datasets_from_column(colname))}
        regex_list = [numeric_regex] * len(dataset_order)
        regex_list[dataset_order[dataset]] = "(" + numeric_regex + ")"
        regex = ";".join(regex_list)
        dataset_class_settings["columns"][column] = [[colname, regex, None]]

    print("\tSaving dataset {}".format(dataset))
    with gzopen(args.out + '_' + dataset + '.json', 'w') as f:
        json.dump(dataset_class_settings, f, indent=4)
    f.close()

print("\nDone")