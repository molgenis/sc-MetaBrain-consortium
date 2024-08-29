#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import pandas as pd
import glob
import re
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--workdir", required=True, type=str, help="")
args = parser.parse_args()

os.makedirs(args.workdir, exist_ok=True)

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

def get_key_value(setting):
    value = re.match("([0-9]+)", setting).group(1)
    key = setting.lstrip(value)
    return key, value


print("\nLoading ncells stats ...")
data = {}
order = ["Cell type"]
for i, fpath in enumerate(glob.glob(os.path.join(args.workdir, "*/*/*.pseudobulk.stats.tsv"))):
    fparts = fpath.split(os.sep)

    match = re.match("([A-Za-z]+).pseudobulk.stats.tsv", os.path.basename(fpath))
    cell_type = match.group(1)

    settings = {"Cell type": cell_type}
    for fpart in fparts[-3:-1]:
        for setting in fpart.split("_"):
            key, value = get_key_value(setting)
            settings[key] = value
            if key not in order:
                order.append(key)

    stats = pd.read_csv(fpath, sep="\t", header=0, index_col=0)
    if "Dataset" in stats.columns:
        settings["ncells"] = stats.loc[cell_type, "ncells"]
        if "ncells" not in order:
            order.append("ncells")
    else:
        for dataset in stats.columns:
            settings[dataset] = stats.loc[cell_type, dataset]
            if dataset not in order:
                order.append(dataset)
        settings["ncells"] = stats.loc[cell_type, :].sum()
        if "ncells" not in order:
            order.append("ncells")
    data[i] = settings
if len(data) == 0:
    exit()
stats = pd.DataFrame(data).T.loc[:, order]
print("\tLoaded cell stats with shape: {}".format(stats.shape))

print("\nLoaded stats:")
stats.sort_values(by=["Cell type", "ncells"], ascending=[True, False], inplace=True)
print(stats)

print("\nSaving files")
stats.to_csv(os.path.join(args.workdir, "ncells.stats.tsv"), sep="\t", header=True, index=True)

print("Done")