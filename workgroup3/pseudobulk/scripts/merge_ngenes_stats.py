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


print("\nLoading ngenes stats ...")
data = {}
order = ["Cell type"]
for i, fpath in enumerate(glob.glob(os.path.join(args.workdir, "*/*/*/*.pseudobulk.*.stats.tsv"))):
    fparts = fpath.split(os.sep)

    match = re.match("([A-Za-z]+).pseudobulk.([0-9A-Za-z_]+).stats.tsv", os.path.basename(fpath))
    cell_type = match.group(1)
    norm = match.group(2)

    settings = {"Cell type": cell_type}
    for fpart in fparts[-4:-1]:
        for setting in fpart.split("_"):
            key, value = get_key_value(setting)
            settings[key] = value
            if key not in order:
                order.append(key)

    stats = pd.read_csv(fpath, sep="\t", header=0, index_col=None)
    for index, row in stats.iterrows():
        if index == stats.index[0]:
            key = row["filter"]
            settings[key] = row["ngenes"]
            if key not in order:
                order.append(key)
        else:
            key = str(index) + "_" + row["filter"]
            settings[key] = row["ngenes"]
            if key not in order:
                order.append(key)

        if index == stats.index[-1]:
            settings["ngenes"] = row["ngenes"]
            if "ngenes" not in order:
                order.append("ngenes")
    data[i] = settings
if len(data) == 0:
    exit()
order.remove("ngenes")
order.append("ngenes")
stats = pd.DataFrame(data).T.loc[:, order]
print("\tLoaded cell stats with shape: {}".format(stats.shape))

print("\nLoaded stats:")
stats.sort_values(by=["Cell type", "Norm", "ngenes"], ascending=[True, True, False], inplace=True)
print(stats)

print("\nSaving files")
stats.to_csv(os.path.join(args.workdir, "ngenes.stats.tsv"), sep="\t", header=True, index=True)

print("Done")