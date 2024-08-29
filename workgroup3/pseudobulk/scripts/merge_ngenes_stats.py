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
    key = setting.replace(value, "")
    return key, value


print("\nLoading ngenes stats ...")
data = {}
for i, fpath in enumerate(glob.glob(os.path.join(args.workdir, "*/*/*.pseudobulk.*.*.stats.tsv"))):
    print(fpath)
    fparts = fpath.split(os.sep)

    match = re.match("([A-Za-z]+).pseudobulk.([0-9A-Za-z_]+).([0-9A-Za-z_]+).stats.tsv", os.path.basename(fpath))
    cell_type = match.group(1)
    filter_fpart = match.group(2)
    norm = match.group(3)

    if norm == "TMM":
        method = "Bryois"
    elif norm == "log2CPM_QN":
        method = "Fujita"
    else:
        print("Error, unexpected norm value '{]'".format(norm))
        exit()

    settings = {"Cell type": cell_type, "method": method}
    for fpart in fparts[-3:-1] + [filter_fpart]:
        for setting in fpart.split("_"):
            key, value = get_key_value(setting)
            settings[key] = value

    stats = pd.read_csv(fpath, sep="\t", header=0, index_col=None)
    for index, row in stats.iterrows():
        if index == stats.index[0]:
            settings[row["filter"]] = row["ngenes"]
        else:
            settings[str(index) + "_" + row["filter"]] = row["ngenes"]

        if index == stats.index[-1]:
            settings["ngenes"] = row["ngenes"]
    data[i] = settings
if len(data) == 0:
    exit()
stats = pd.DataFrame(data).T
print("\tLoaded cell stats with shape: {}".format(stats.shape))

print("\nLoaded stats:")
order = [column for column in stats.columns if column != "ngenes"] + ["ngenes"]
stats = stats.loc[:, order]
stats.sort_values(by=["Cell type", "method", "ngenes"], ascending=[True, True, False], inplace=True)

print("\nSaving files")
stats.to_csv(os.path.join(args.workdir, "ngenes.stats.tsv"), sep="\t", header=True, index=True)

print("Done")