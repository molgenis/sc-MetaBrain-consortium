#!/usr/bin/env python
# Author: M. Vochteloo

import warnings
warnings.simplefilter("ignore", UserWarning)

import numpy as np
import pandas as pd
import scanpy
import h5py
import argparse
import gzip
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--poolsheet", required=True, type=str, help="")
parser.add_argument("--pool", required=True, type=str, help="")
parser.add_argument("--psam", required=True, type=str, help="")
parser.add_argument("--cell_annotation", required=True, type=str, help="")
parser.add_argument("--droplet_type_annotation", required=True, type=str, help="")
parser.add_argument("--cell_type_annotation", required=True, type=str, help="")
parser.add_argument("--ct_pairing", required=False, type=str, default=None, help="")
parser.add_argument("--rb_genes", required=True, type=str, help="")
parser.add_argument("--mt_genes", required=True, type=str, help="")
parser.add_argument("--feature_name", required=False, type=str, default="HGNC", choices=["HGNC", "ENSG", "HGNC_ENSG"], help="")
parser.add_argument("--sample_aggregate", required=False, type=str, default="Assignment", help="")
parser.add_argument("--ancestry", required=False, type=str, default="EUR", help="")
parser.add_argument("--cell_level", required=False, type=str, default="L1", help="")
parser.add_argument("--ncount_rna", required=False, type=int, default=500, help="")
parser.add_argument("--nfeature_rna", required=False, type=int, default=0, help="")
parser.add_argument("--complexity", required=False, type=int, default=100, help="")
parser.add_argument("--percent_rb", required=False, type=int, default=100, help="")
parser.add_argument("--percent_mt", required=False, type=int, default=5, help="")
parser.add_argument("--malat1", required=False, type=int, default=0, help="")
parser.add_argument("--cap_barcodes", required=False, type=int, default=None, help="")
parser.add_argument("--cr_barcodes", action="store_true", default=False, help="")
parser.add_argument("--min_cells", required=False, type=int, default=5, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

os.makedirs(args.out, exist_ok=True)
data_out = os.path.join(args.out, "data")
os.makedirs(data_out, exist_ok=True)

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)

def load_file(fpath, sep="\t", header=0, index_col=None, must_contain=None):
    data = []
    columns = None
    index = 0
    with gzopen(fpath, mode="r") as fh:
        for index, line in enumerate(fh):
            if index == 0 or (index + 1) % 1e5 == 0:
                print("\tParsed {:,} lines".format(index + 1), end='\r')

            if index != header and must_contain is not None and must_contain not in line:
                continue

            values = line.rstrip("\n").split(sep)
            if index == header:
                columns = values
                continue
            data.append(values)
    print("\tParsed {:,} lines".format(index + 1), end='\n')

    if len(data) == 0:
        print("\tError, failed to load {}: no data".format(os.path.basename(fpath)))
        exit()

    df = pd.DataFrame(data, columns=columns)
    if index_col is not None:
        if len(set(df.columns)) != df.shape[1]:
            print("Error, not all columns are unique.")
            exit()

        column = df.columns[index_col]
        df.index = df[column]
        df.drop(column, axis=1, inplace=True)

    print("\tLoaded dataframe {} with shape: {}".format(os.path.basename(fpath), df.shape))
    return df

def load_metadata(barcodes_fpath):
    print("  Loading Barcodes")
    barcodes = load_file(barcodes_fpath, header=None).rename(columns={0: "Barcode"})
    barcodes.reset_index(drop=False, inplace=True)
    barcodes["Barcode"] = barcodes["Barcode"].str.split("-", n=1, expand=True)[0] + "_" + str(args.pool)
    barcodes["Pool"] = args.pool

    print("  Loading cell annotations")
    cell_annot = load_file(args.cell_annotation, must_contain=args.pool)
    cr_barcodes = None
    if "CellBender" in cell_annot:
        cr_barcodes = [barcode.split("_")[0] + "-1" for barcode in cell_annot.loc[cell_annot["CellBender"] == "False", "Barcode"].values]
    elif "CellRanger" in cell_annot:
        cr_barcodes = [barcode.split("_")[0] + "-1" for barcode in cell_annot.loc[cell_annot["CellRanger"] == "True", "Barcode"].values]
    else:
        print("  Warning, could not add CellRanger column.")

    print("  Loading droplet type annotations")
    droplet_annot = load_file(args.droplet_type_annotation, must_contain=args.pool)

    print("  Loading cell type annotations")
    type_annot = load_file(args.cell_type_annotation, must_contain=args.pool)

    print("  Loading PSAM")
    psam = load_file(args.psam).rename(columns={"#FID": "FID"})

    ct_pairs = None
    if args.ct_pairing is not None:
        print("  Loading cell type pairing")
        ct_pairs = load_file(args.ct_pairing, sep=";")

    merge_incl = []
    if "Dataset" in cell_annot.columns and "Dataset" in droplet_annot.columns and "Dataset" in type_annot.columns and "Dataset" in psam.columns:
        print("  Warning, found 'Dataset' column. Will include this as key in merge.")
        merge_incl = ["Dataset"]

    print("\nMerging metadata ...")
    metadata = barcodes.merge(cell_annot, on="Barcode", how="inner").merge(droplet_annot, on=["Pool", "Barcode"] + merge_incl, how="inner").merge(type_annot, on=["Pool", "Barcode"] + merge_incl, how="left").merge(psam, left_on=["Assignment"] + merge_incl, right_on=["IID"] + merge_incl, how="left")
    if ct_pairs is not None:
        metadata = metadata.merge(ct_pairs, on=ct_pairs.columns[1], how="left")

    # Adding extra columns.
    metadata["Assignment_Run_Lane"] = metadata["Assignment"].astype(str) + ";;" + metadata["sequencing_run"].astype(str) + "_" + metadata["sequencing_lane"].astype(str)

    # Reorder back into original order.
    metadata = metadata.sort_values(by="index", ascending=True).drop(["index"], axis=1)

    # Check if all barcodes are unique.
    if len(metadata["Barcode"].unique()) != metadata.shape[0]:
        print("Error, not all barcodes are unique in the metadata.")
        exit()

    print("\tLoaded metadata with shape: {}".format(metadata.shape))
    return metadata, cr_barcodes

def load_counts(counts_fpath):
    adata = scanpy.read_10x_h5(counts_fpath)
    count_matrix = adata.X
    n_barcodes, n_genes = count_matrix.shape
    genes = adata.var_names.to_numpy()
    barcodes = adata.obs_names.to_numpy()
    if np.size(genes) != n_genes:
        print("Error, matrix size does not match gene annotation.")
        exit()
    if np.size(barcodes) != n_barcodes:
        print("Error, matrix size does not match barcodes annotation.")
        exit()
    print("\tLoaded matrix with {:,} barcodes and {:,} genes.".format(n_barcodes, n_genes))

    # Check for unique barcodes.
    if len(np.unique(barcodes)) != np.size(barcodes):
        print("Error, not all barcodes are unique.")
        exit()

    # Parse the gene info
    gene_ids = adata.var['gene_ids']
    gene_ids = gene_ids.reset_index(drop=False)
    gene_ids = gene_ids.to_numpy()
    if not (genes == gene_ids[:, 0]).all():
        print("Error, 'count_data.var_names' are expected to have the same order as 'count_data.var['gene_ids'].")
        exit()

    # Set gene names.
    features = get_features(m=gene_ids, indices={"HGNC": 0, "ENSG": 1})

    # Remove duplicates.
    u, c = np.unique(features, return_counts=True)
    dup = u[c > 1]
    features_mask = np.in1d(features, dup)
    if np.sum(features_mask) > 0:
        print("\tRemoving {:,} duplicate genes.".format(np.sum(features_mask)))
        adata = adata[:, ~features_mask]
        adata.var_names = features[~features_mask]

    # Parse the barcode info.
    # barcode_ids = count_data.var['feature_types'].to_numpy()

    print("\tLoaded expression with shape: {}".format(count_matrix.shape))
    return adata

def get_features(m, indices):
    if args.feature_name == "HGNC":
        return m[:, indices["HGNC"]]
    elif args.feature_name == "ENSG":
        return m[:, indices["ENSG"]]
    elif args.feature_name == "HGNC_ENSG":
        return np.char.add(np.char.add(m[:, indices["HGNC"]], "_"), m[:, indices["ENSG"]])
    else:
        print("Unexpected feature name '{}'".format(args.feature_name))
        exit()

def load_features(fpath):
    gene_ids = load_file(fpath=fpath).to_numpy()
    features = get_features(m=gene_ids, indices={"HGNC": 1, "ENSG": 0})
    return features

def get_malat1_feature():
    features = get_features(m=np.array([["ENSG00000251562", "MALAT1"]]), indices={"HGNC": 1, "ENSG": 0})
    return features[0]

def calc_barcode_qc(adata, cr_barcodes):
    print("  Adding UMI and feature counts")
    info = pd.DataFrame(np.hstack((np.sum(adata.X, axis=1), np.sum(adata.X != 0, axis=1))),
                        columns=["nCount_RNA", "nFeature_RNA"], index=adata.obs_names)

    print("  Adding complexity")
    info["complexity"] = np.nan
    mask = info["nCount_RNA"] > 0
    info.loc[mask, "complexity"] = (info.loc[mask, "nFeature_RNA"] / info.loc[mask, "nCount_RNA"]) * 100
    del mask

    print("  Adding ribosomal %")
    rb_features = load_features(args.rb_genes)
    rb_mask = np.in1d(adata.var_names, rb_features)
    info["rb"] = np.sum(adata.X[:, rb_mask], axis=1)

    info["percent.rb"] = np.nan
    mask = info["nCount_RNA"] > 0
    info.loc[mask, "percent.rb"] = (info.loc[mask, "rb"] / info.loc[mask, "nCount_RNA"]) * 100
    del rb_features, rb_mask, mask

    print("  Adding mitochondiral %")
    mt_features = load_features(args.mt_genes)
    mt_mask = np.in1d(adata.var_names, mt_features)
    info["mt"] = np.sum(adata.X[:, mt_mask], axis=1)

    info["percent.mt"] = np.nan
    mask = info["nCount_RNA"] > 0
    info.loc[mask, "percent.mt"] = (info.loc[mask, "mt"] / info.loc[mask, "nCount_RNA"]) * 100
    del mt_features, mt_mask, mask

    print("  Adding MALAT1")
    info["MALAT1"] = 0
    malat1_feature = get_malat1_feature()
    malat1_mask = adata.var_names == malat1_feature
    info["MALAT1"] = adata.X[:, malat1_mask].toarray()
    del malat1_feature, malat1_mask

    print("  Adding nCount_RNA index")
    info["nCount_RNAIndex"] = info[["nCount_RNA"]].rank(method="first", ascending=False)

    print("  Adding CellRanger flag")
    if cr_barcodes is None:
        info["CellRanger"] = True
    else:
        info["CellRanger"] = False
        info.loc[info.index.isin(cr_barcodes), "CellRanger"] = True

    return info

def apply_barcode_qc(barcode_qc):
    print(barcode_qc)

    n_barcodes = barcode_qc.shape[0]
    print("  Input: {:,} barcodes".format(n_barcodes))

    ncount_rna_mask = barcode_qc["nCount_RNA"] >= args.ncount_rna
    print("\tnCount_RNA >= {} yields {:,} [{:.2f}%] barcodes".format(args.ncount_rna, np.sum(ncount_rna_mask), (100 / n_barcodes) * np.sum(ncount_rna_mask)))

    nfeature_rna_mask = barcode_qc["nFeature_RNA"] >= args.nfeature_rna
    print("\tnFeature_RNA >= {} yields {:,} [{:.2f}%] barcodes".format(args.nfeature_rna, np.sum(nfeature_rna_mask), (100 / n_barcodes) * np.sum(nfeature_rna_mask)))

    complexity_mask = barcode_qc["complexity"] <= args.complexity
    print("\tcomplexity <= {} yields {:,} [{:.2f}%] barcodes".format(args.complexity, np.sum(complexity_mask), (100 / n_barcodes) * np.sum(complexity_mask)))

    percent_rb_mask = barcode_qc["percent.rb"] <= args.percent_rb
    print("\tpercent_rb <= {} yields {:,} [{:.2f}%] barcodes".format(args.percent_rb, np.sum(percent_rb_mask), (100 / n_barcodes) * np.sum(percent_rb_mask)))

    percent_mt_mask = barcode_qc["percent.mt"] <= args.percent_mt
    print("\tpercent_mt <= {} yields {:,} [{:.2f}%] barcodes".format(args.percent_mt, np.sum(percent_mt_mask), (100 / n_barcodes) * np.sum(percent_mt_mask)))

    malat1_mask = barcode_qc["MALAT1"] >= args.malat1
    print("\tMALAT1 >= {} yields {:,} [{:.2f}%] barcodes".format(args.malat1, np.sum(malat1_mask), (100 / n_barcodes) * np.sum(malat1_mask)))

    cap_barcodes = np.inf
    if args.cap_barcodes is not None:
        cap_barcodes = args.cap_barcodes
    cap_barcodes_mask = barcode_qc["nCount_RNAIndex"] <= cap_barcodes
    print("\tmax cells <= {} yields {:,} [{:.2f}%] barcodes".format(cap_barcodes, np.sum(cap_barcodes_mask), (100 / n_barcodes) * np.sum(cap_barcodes_mask)))

    cr_barcodes_mask = pd.Series(True, index=barcode_qc.index, name="index")
    if args.cr_barcodes:
        cr_barcodes_mask = barcode_qc["CellRanger"]
    print("\tCellRanger barcodes yields {:,} [{:.2f}%] barcodes".format(np.sum(cr_barcodes_mask), (100 / n_barcodes) * np.sum(cr_barcodes_mask)))

    mask = ncount_rna_mask & nfeature_rna_mask & complexity_mask & percent_rb_mask & percent_mt_mask & malat1_mask & cap_barcodes_mask & cr_barcodes_mask

    print("  Pass-QC: {:,} [{:.2f}%] barcodes".format(np.sum(mask), (100 / n_barcodes) * np.sum(mask)))

    barcode_qc["tag"] = "Outlier"
    barcode_qc.loc[mask, "tag"] = "NotOutlier"

    return barcode_qc

def filter_split(adata, metadata):
    ncells = adata.X.shape[0]

    stats_data = []
    samples = metadata[args.sample_aggregate].unique()
    n_samples = 0
    for sample in samples:
        if sample == "doublet" or sample == "unassigned" or sample == "":
            continue

        # Filter on the cells from the current samples.
        sample_mask = metadata[args.sample_aggregate] == sample
        n_sample_cells = np.sum(sample_mask)
        print("  {} has {:,} input cells [{:.2f}%]".format(sample, n_sample_cells, (100 / ncells) * n_sample_cells))

        if n_sample_cells == 0:
            stats_data.append([sample, np.nan, ncells, n_sample_cells, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
            continue

        # Loop over the different cell types that this sample has.
        cell_types = metadata.loc[sample_mask, args.cell_level].unique()
        for cell_type in cell_types:
            if pd.isnull(cell_type):
                continue

            # Filter on both sample and cel type.
            sample_ct_mask = (metadata[args.sample_aggregate] == sample) & (metadata[args.cell_level] == cell_type)
            n_sample_ct_cells = np.sum(sample_ct_mask)
            print("\t{} has {:,} input cells".format(cell_type, n_sample_ct_cells))

            if n_sample_ct_cells == 0:
                stats_data.append([sample, cell_type, ncells, n_sample_cells, n_sample_ct_cells, np.nan, np.nan, np.nan, np.nan, np.nan])
                continue

            # Filter on the other characteristics.
            mask1 = sample_ct_mask & (metadata["DropletType"] == "singlet")
            mask2 = sample_ct_mask & (metadata["tag"] == "NotOutlier")
            mask3 = sample_ct_mask & (metadata["Provided_Ancestry"] == args.ancestry)
            mask4 = sample_ct_mask & (metadata["cell_treatment"] == "UT")
            mask = sample_ct_mask & mask1 & mask2 & mask3 & mask4

            print("\t  DropletType - Singlet: N = {:,} [{:.2f}%]".format(np.sum(mask1), (100 / n_sample_ct_cells) * np.sum(mask1)))
            print("\t  tag - NotOutlier: N = {:,} [{:.2f}%]".format(np.sum(mask2), (100 / n_sample_ct_cells) * np.sum(mask2)))
            print("\t  Provided_Ancestry - {}: N = {:,} [{:.2f}%]".format(args.ancestry, np.sum(mask3), (100 / n_sample_ct_cells) * np.sum(mask3)))
            print("\t  cell_treatment - UT: N = {:,} [{:.2f}%]".format(np.sum(mask4), (100 / n_sample_ct_cells) * np.sum(mask4)))
            print("\t{} has {:,} cells pass QC [{:.2f}%]\n".format(cell_type, np.sum(mask), (100 / n_sample_ct_cells) * np.sum(mask)))

            n_sample_ct_pass_cells = np.sum(mask)
            stats_data.append([sample, cell_type, ncells, n_sample_cells, n_sample_ct_cells, np.sum(mask1), np.sum(mask2), np.sum(mask3), np.sum(mask4), n_sample_ct_pass_cells])
            
            # Stop if there are not enough cells.
            if n_sample_ct_pass_cells < args.min_cells:
                continue

            # Aggregate the cells.
            save_file(df=pd.DataFrame(np.sum(adata.X[mask, :], axis=1), index=adata.obs_names[mask], columns=[sample]).T, outpath=os.path.join(data_out, f"{args.pool}.{sample}.{cell_type}.raw.weights.txt.gz"))

            # Save as h5.
            save_filtered_counts_h5(fpath=os.path.join(data_out, f"{args.pool}.{sample}.{cell_type}.raw.counts.h5"), adata=adata[mask, :])
            n_samples += 1

    stats_df = None
    if len(stats_data) > 0:
        stats_df = pd.DataFrame(stats_data, columns=["sample", "cell type", "ncells_pool", "ncells_sample", "ncells_sample_cell_type", "ncells_singlet", "ncells_pass_qc", "ncells_ancestry", "ncells_treatment", "ncells"])
        print("\tCombined {:,} cells over {:,} samples.".format(stats_df["ncells"].sum(), n_samples))

    return stats_df


def save_file(df, outpath, header=True, index=False, sep="\t"):
    compression = 'infer'
    if outpath.endswith('.gz'):
        compression = 'gzip'

    df.to_csv(outpath, sep=sep, index=index, header=header,
              compression=compression)


def save_filtered_counts_h5(fpath, adata):
    """
    Write hdf5 file from Cell Ranger v4 or later versions.
    Source: https://github.com/scverse/anndata/issues/595
    """
    if os.path.exists(fpath):
        raise FileExistsError(f"There already is a file `{fpath}`.")

    int_max = lambda x: int(max(np.floor(len(str(int(max(x)))) / 4), 1) * 4)
    str_max = lambda x: max([len(i) for i in x])

    # Save.
    w = h5py.File(fpath, 'w')
    grp = w.create_group("matrix")
    grp.create_dataset("barcodes", data=np.array(adata.obs_names, dtype=f'|S{str_max(adata.obs_names)}'))
    grp.create_dataset("data", data=np.array(adata.X.data, dtype=f'<i{int_max(adata.X.data)}'))
    ftrs = grp.create_group("features")
    # this group will lack the following keys:
    # '_all_tag_keys', 'feature_type', 'genome', 'id', 'name', 'pattern', 'read', 'sequence'
    ftrs.create_dataset("feature_type", data=np.array(adata.var.feature_types, dtype=f'|S{str_max(adata.var.feature_types)}'))
    ftrs.create_dataset("genome", data=np.array(adata.var.genome, dtype=f'|S{str_max(adata.var.genome)}'))
    ftrs.create_dataset("id", data=np.array(adata.var.gene_ids, dtype=f'|S{str_max(adata.var.gene_ids)}'))
    ftrs.create_dataset("name", data=np.array(adata.var.index, dtype=f'|S{str_max(adata.var.index)}'))
    grp.create_dataset("indices", data=np.array(adata.X.indices, dtype=f'<i{int_max(adata.X.indices)}'))
    grp.create_dataset("indptr", data=np.array(adata.X.indptr, dtype=f'<i{int_max(adata.X.indptr)}'))
    grp.create_dataset("shape", data=np.array(list(adata.X.shape)[::-1], dtype=f'<i{int_max(adata.X.shape)}'))


#############################################

print("\nLoading poolsheet ...")
fpaths_df = load_file(args.poolsheet, must_contain=args.pool)
if fpaths_df.shape[0] > 1:
    fpaths_df = fpaths_df.loc[fpaths_df["Pool"] == args.pool, :].reset_index()
fpaths = fpaths_df.to_dict("index")
del fpaths_df

if len(fpaths) != 1:
    print("Error, pool is not unique.")
    exit()
if fpaths[0]["Pool"] != args.pool:
    print("Error, unexpected Pool value.")
    exit()

counts_fpath = fpaths[0]["Counts"]
barcodes_fpath = fpaths[0]["Barcodes"]
print("  Counts fpath: " + counts_fpath)
print("  Barcodes fpath: " + barcodes_fpath)

print("\nCreating metadata ...")
metadata, cr_barcodes = load_metadata(barcodes_fpath=barcodes_fpath)

print("\nLoading counts matrix ...")
adata = load_counts(counts_fpath=counts_fpath)

print("\nCalculating barcode QC stats ...")
barcode_qc = calc_barcode_qc(adata=adata, cr_barcodes=cr_barcodes)

print("\nApplying barcode QC ...")
barcode_qc = apply_barcode_qc(barcode_qc=barcode_qc)
if barcode_qc.shape[0] == 0:
    print("Error, 0 barcodes passed QC.")
    exit()

print("\nMerging metadata and barcode QC ...")
# Add the Barcode_Pool column.
barcode_qc.reset_index(drop=False, inplace=True)
barcode_qc = barcode_qc.rename(columns={"index": "barcode"})
barcode_qc["Barcode"] = barcode_qc["barcode"].str.split("-", n=1, expand=True)[0] + "_" + str(args.pool)

# Merge on Barcode and reorder by the count matrix (barcode column).
metadata = metadata.merge(barcode_qc, on="Barcode", how="inner").set_index("barcode").loc[adata.obs_names, :]

print("\tSaving file")
metadata.to_csv(os.path.join(args.out, str(args.pool) + ".full.metadata.tsv.gz"), sep="\t", header=True, index=False, compression="gzip")

print("\nPseudobulking per cell type ...")
stats = filter_split(adata=adata, metadata=metadata)

print("\nBarcode selection stats output:")
print(stats)

if stats is not None:
    print("\nBarcode selection summary stats:")
    sumstats = stats.loc[:, [col for col in stats.columns if col not in ["sample", "ncells_pool", "ncells_sample"]]].groupby("cell type").sum()
    sumstats.sort_values(by="ncells", ascending=False, inplace=True)
    sumstats = sumstats.T
    sumstats["Total"] = sumstats.sum(axis=1)
    print(sumstats)
    del sumstats

print("Done")