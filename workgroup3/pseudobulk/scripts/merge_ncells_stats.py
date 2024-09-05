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
    if setting[0].isdigit():
        value = re.match("([0-9]+)", setting).group(1)
    elif setting.startswith("None"):
        value = "None"
    elif setting.startswith("True"):
        value = "True"
    elif setting.startswith("False"):
        value = "False"
    elif setting.startswith("Bryois"):
        value = "Bryois"
    elif setting.startswith("Fujita"):
        value = "Fujita"
    else:
        print("Error in get_key_value for setting: {}".format(setting))
        exit()
    key = setting.lstrip(value)
    return key, value


print("\nLoading ncells stats ...")
data = {}
setting_keys = []
datasets = set()
for i, fpath in enumerate(glob.glob(os.path.join(args.workdir, "*/*/*.pseudobulk.stats.tsv"))):
    fparts = fpath.split(os.sep)

    match = re.match("([A-Za-z]+).pseudobulk.stats.tsv", os.path.basename(fpath))
    cell_type = match.group(1)

    settings = {"Cell type": cell_type}
    for fpart in fparts[-3:-1]:
        for setting in fpart.split("_"):
            key, value = get_key_value(setting)
            settings[key] = value
            if key not in setting_keys:
                setting_keys.append(key)

    stats = pd.read_csv(fpath, sep="\t", header=0, index_col=0)
    settings.update(stats.groupby("Dataset").sum()["ncells"].to_dict())
    settings["ncells"] = stats["ncells"].sum()
    datasets.update(stats["Dataset"].unique())
    data[i] = settings
if len(data) == 0:
    exit()

stats = pd.DataFrame(data).T
print("\tLoaded cell stats with shape: {}".format(stats.shape))

# Post-processing.
datasets = list(datasets)
datasets.sort()
order = ["Cell type"] + datasets + ["ncells"]
stats = stats.loc[:, order]
stats.sort_values(by=["Cell type", "ncells"], ascending=[True, False], inplace=True)

print("\nLoaded stats:")
print(stats)

print("\nSaving files")
stats.to_csv(os.path.join(args.workdir, "ncells.stats.tsv"), sep="\t", header=True, index=True)

print("Done")