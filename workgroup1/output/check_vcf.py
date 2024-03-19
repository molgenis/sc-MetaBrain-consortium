#!/usr/bin/env python3

import gzip
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import os

ANCESTRY = "EUR"
CHR = "X"
nrows = 1000


# old_indir = ""
# new_indir = ""
# info = []

def get_columns(inpath):
    header = None
    skiprows = 0

    if inpath.endswith(".gz"):
        f = gzip.open(inpath, 'rt')
    else:
        f = open(inpath, 'r')

    for i, line in enumerate(f):
        if line.startswith("##"):
            skiprows += 1
            continue

        header = line.strip("\n").split("\t")
        break

    f.close()

    return skiprows, header

# def parse_file(inpath, nrows):
#     header = None
#     data = []
#     with gzip.open(inpath, 'rt') as f:
#         for i, line in enumerate(f):
#             if len(data) == nrows:
#                 break
#             if line.startswith("##"):
#                 continue
#
#             values = line.strip("\n").split("\t")
#
#             if header is None:
#                 header = values
#             else:
#                 if values[1] != CHR:
#                     continue
#                 data.append(values)
#     f.close()
#
#     return pd.DataFrame(data, colums=header)

def split_info(s):
    data = []
    for _, value in s.iteritems():
        fields = value.split(";")

        row_data = {}
        for field in fields:
            splitted_fields = field.split("=")
            if len(splitted_fields) == 2:
                try:
                    row_data[splitted_fields[0]] = float(splitted_fields[1])
                except ValueError:
                    row_data[splitted_fields[0]] = splitted_fields[1]
            else:
                row_data[splitted_fields[0]] = ""

        data.append(pd.Series(row_data))

    return pd.DataFrame(data)

def split_format(df):
    data = []
    sample_index = df.columns.get_loc("FORMAT") + 1
    for _, row in df.iterrows():
        format = row["FORMAT"].split(":")

        row_data = {}
        for sample, value in zip(df.columns[sample_index:], row[sample_index:]):
            values = value.split(":")
            for key, value in zip(format, values):
                if key == "GP":
                    for i, sub_value in enumerate(value.split(",")):
                        row_data[sample + "-" + key + "[" + str(i) + "]"] = float(sub_value)
                elif key == "DS":
                    row_data[sample + "-" + key] = float(value)
                else:
                    row_data[sample + "-" + key] = value

        data.append(pd.Series(row_data))

    return pd.DataFrame(data)


def visualise_data(plot_data, title="", filename="plot"):
    keys = list(plot_data.keys())
    keys.sort()

    nplots = len(keys)
    ncols = math.ceil(np.sqrt(nplots))
    nrows = math.ceil(nplots / ncols)

    sns.set_style("ticks")
    fig, axes = plt.subplots(nrows=nrows,
                             ncols=ncols,
                             sharex='none',
                             sharey='none',
                             figsize=(6 * ncols, 6 * nrows))
    sns.set(color_codes=True)

    row_index = 0
    col_index = 0
    for i in range(ncols * nrows):
        if nrows == 1 and ncols == 1:
            ax = axes
        elif nrows == 1 and ncols > 1:
            ax = axes[col_index]
        elif nrows > 1 and ncols == 1:
            ax = axes[row_index]
        else:
            ax = axes[row_index, col_index]

        if i < nplots:
            data = plot_data[keys[i]]
            if data.shape[0] == 0:
                continue

            sns.despine(fig=fig, ax=ax)

            sns.regplot(x="x", y="y", data=data, ci=95,
                        scatter_kws={'facecolors': "#808080",
                                     'edgecolors': "#808080"},
                        line_kws={"color": "#b22222"},
                        ax=ax
                        )

            ax.axline((0, 0), slope=1, ls='--', color="#000000", alpha=0.15, zorder=-1)

            # Set annotation.
            pearson_coef, _ = stats.pearsonr(data["y"], data["x"])
            ax.annotate(
                'total N = {:,}'.format(data.shape[0]),
                xy=(0.03, 0.94),
                xycoords=ax.transAxes,
                color="#000000",
                fontsize=14,
                fontweight='bold')
            ax.annotate(
                'total r = {:.4f}'.format(pearson_coef),
                xy=(0.03, 0.88),
                xycoords=ax.transAxes,
                color="#000000",
                fontsize=14,
                fontweight='bold')

            ax.set_xlabel("Drew",
                          fontsize=10,
                          fontweight='bold')
            ax.set_ylabel("Martijn",
                          fontsize=10,
                          fontweight='bold')
            ax.set_title(keys[i],
                         fontsize=14,
                         fontweight='bold')
        else:
            ax.set_axis_off()

        col_index += 1
        if col_index > (ncols - 1):
            col_index = 0
            row_index += 1

    fig.suptitle(title,
                 fontsize=16,
                 color="#000000",
                 weight='bold')

    if not os.path.exists("plot"):
        os.makedirs("plot")

    plt.tight_layout()
    outpath = os.path.join("plot/{}.png".format(filename))
    fig.savefig(outpath)
    print("\tSaved figure: {}".format(os.path.basename(outpath)))
    plt.close()

for title, subdir1, subdir2 in info:
    if subdir2 is None:
        subdir2 = subdir1
    print(title)

    inpath1 = os.path.join(old_indir, subdir1)
    print(inpath1)
    skiprows1, usecols1 = get_columns(inpath1)
    df1 = pd.read_csv(inpath1, sep="\t", skiprows=skiprows1, nrows=nrows, usecols=usecols1[:11])
    info_df1 = split_info(df1["INFO"])
    sample_df1 = split_format(df1)
    df1 = pd.concat([df1.iloc[:, :7], info_df1, sample_df1], axis=1)
    df1.index = df1["#CHROM"].astype(str) + ":" + df1["POS"].astype(str) + ":" + df1["REF"] + "/" + df1["ALT"]
    print(df1)
    del info_df1, sample_df1

    inpath2 = os.path.join(new_indir, subdir2)
    print(inpath2)
    skiprows2, usecols2 = get_columns(inpath2)
    df2 = pd.read_csv(inpath2, sep="\t", skiprows=skiprows2, nrows=nrows, usecols=usecols2[:11])
    info_df2 = split_info(df2["INFO"])
    sample_df2 = split_format(df2)
    df2.drop(["INFO"], axis=1, inplace=True)
    df2 = pd.concat([df2.iloc[:, :7], info_df2, sample_df2], axis=1)
    df2.index = df2["#CHROM"].astype(str) + ":" + df2["POS"].astype(str) + ":" + df2["REF"] + "/" + df2["ALT"]
    print(df2)
    del info_df2, sample_df2

    # r2_df = df1[["R2"]].merge(df2[["R2"]], left_index=True, right_index=True)
    # r2_df["max"] = r2_df.max(axis=1)
    # mask = (r2_df["max"] >= 0.4).to_numpy()
    # df1 = df1.loc[mask, :]
    # df2 = df2.loc[mask, :]
    #
    # print(df1)
    # print(df2)

    print("\nComparing:")
    print("  Fully identical:", df1.equals(df2))
    print("")

    if df1.equals(df2):
        continue

    if df1.index.tolist() != df2.index.tolist():
        print("  Indices differ, filtering on overlapping")
        overlap = [index for index in df1.index if index in df2.index]
        df1 = df1.loc[overlap, :]
        df2 = df2.loc[overlap, :]
        print("  Data frames identical after row filter:", df1.equals(df2))
    else:
        print("  Indices are identical")
    print("")

    if df1.columns.tolist() != df2.columns.tolist():
        print("  Columns differ, filtering on overlapping")
        overlap = [column for column in df1.columns if column in df2.columns]
        df1 = df1.loc[:, overlap]
        df2 = df2.loc[:, overlap]
        print("  Data frames identical after col filter:", df1.equals(df2))
    else:
        print("  Columns are identical")
    print("")

    print(df1)
    print(df2)


    if not df1.equals(df2):
        plot_data = {}
        for col in df1.columns:
            if col in df2.columns:
                if df1[col].equals(df2[col]):
                    print("  ", col, "identical")
                elif df1[col].values.tolist() == df2[col].values.tolist():
                    print("  ",col, "identical values, different indices")
                elif pd.api.types.is_numeric_dtype(df1[col]) and pd.api.types.is_numeric_dtype(df1[col]):
                    corr_df = df1[[col]].merge(df2[[col]], left_index=True, right_index=True)
                    corr_df.columns = ["x", "y"]
                    corr_df.dropna(inplace=True)
                    pearson_coef, pearson_p = stats.pearsonr(corr_df["x"], corr_df["y"])
                    print("  ", col, "different", pearson_coef, pearson_p)

                    plot_data[col] = corr_df
                else:
                    print("  ",col, "different values, different indices")

        visualise_data(plot_data=plot_data, title=title, filename=title.lower())

    print("\n----------------------------------------\n")

    del df1, df2
