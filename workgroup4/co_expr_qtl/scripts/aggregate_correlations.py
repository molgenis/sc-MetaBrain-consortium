import argparse
import gzip
import os
import glob
import time
import struct
import sys
import math

import cProfile
import pstats
profiler = cProfile.Profile()
profiler.enable()

parser = argparse.ArgumentParser(description="")
parser.add_argument("--indir", required=True, type=str,  help="Input directory")
parser.add_argument("--geneannotation", required=False, type=str, help="Gene chromosome annotation")
parser.add_argument("--feature_name", required=False, type=str, default="HGNC", choices=["HGNC", "ENSG", "HGNC_ENSG"], help="")
parser.add_argument("--chr", required=False, type=str, help="Limit correlation calculation to gene pairs where the first gene is on specified chromosome")
parser.add_argument("--egenelist", required=False, type=str, help="List of all valid eQTL genes")
parser.add_argument("--coegenelist", required=False, type=str, help="List of all valid co-eQTL genes")
parser.add_argument("--binary_in", action="store_true", default=False, help="Input in binary format")
parser.add_argument("--binary_out", action="store_true", default=False, help="Output in binary format")
parser.add_argument("--out", required=True, type=str, help="Output file name")
args = parser.parse_args()

if os.path.dirname(args.out) != "":
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

if args.geneannotation is not None and args.chr is None:
    print("Warning, --geneannotation is ignored if --chr is not given.")
if args.chr is not None and args.geneannotation is None:
    print("Error, --geneannotation must be given in if --chr is given.")
    exit()
    
N_BYTES_DOUBLE = 8


def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        if mode == "w":
            return gzip.open(file, mode + 't', 4)
        else:
            return gzip.open(file, mode + 't')
    else:
        return open(file, mode)

def load_annotation(annotation, chr=None):
    print("Loading gene annotation...")
    annot_genes = set()
    pos = {}
    with gzopen(annotation, mode='r') as f:
        for i, line in enumerate(f):
            values = line.rstrip().split("\t")
            if i == 0:
                pos = {value: index for index, value in enumerate(values)}
                continue
            if chr is not None and values[pos["Chr"]] != chr:
                continue

            if args.feature_name == "HGNC":
                annot_genes.add(values[pos["GeneSymbol"]])
            elif args.feature_name == "ENSG":
                annot_genes.add(values[pos["Gene"]])
            elif args.feature_name == "HGNC_ENSG":
                annot_genes.add(values[pos["Gene"] + "_" + values[pos["GeneSymbol"]]])
            else:
                print("Unexpected feature name '{}'".format(args.feature_name))
                exit()
    f.close()
    print("  Loaded {:,} genes{}".format(len(annot_genes), " for chr '{}'".format(chr) if chr is not None else ""))

    return annot_genes

def read_as_set(fpath):
    fh = gzopen(fpath,'r')
    data = set()
    for line in fh:
        data.add(line.strip())
    return data

def parse_binary_header(fh):
    n_egenes, n_coegenes, n_correlations, n_feature_chars = struct.unpack(">4i", fh.read(16))

    feature_str = struct.unpack(f'>{n_feature_chars}s', fh.read(n_feature_chars))[0].decode('ascii')
    features = feature_str.split(",")
    del feature_str

    n_features = len(features)
    is_egene = struct.unpack(f'>{n_features}?', fh.read(n_features))

    return n_egenes, n_coegenes, n_correlations, n_feature_chars, features, is_egene


def parse_text_header(fh):
    n_egenes = int(fh.readline().rstrip("\n"))
    n_coegenes = int(fh.readline().rstrip("\n"))
    n_correlations = int(fh.readline().rstrip("\n"))
    n_feature_chars = int(fh.readline().rstrip("\n"))

    feature_str = fh.readline().rstrip("\n")
    features = feature_str.split(",")
    del feature_str

    is_egene = [bool(val) for val in fh.readline().rstrip("\n").split(",")]

    return n_egenes, n_coegenes, n_correlations, n_feature_chars, features, is_egene

def clear(list, binary):
    if binary:
        nanfloatbytes = struct.pack('>d',float('nan'))
        for i in range(0, len(list), N_BYTES_DOUBLE):
            end = i + N_BYTES_DOUBLE
            list[i:end] = nanfloatbytes[0:8]
    else:
        nan = "NaN"
        for i in range(len(list)):
            list[i] = nan

def to_time_str(nanoseconds):
    seconds = nanoseconds / 1e9
    minutes = (seconds - (seconds % 60)) / 60
    hours = (minutes - (minutes % 60 )) / 60
    remainingmins = minutes % 60
    remainingseconds = seconds % 60
    return f"{hours:.0f}h, {remainingmins:.0f}m, {remainingseconds:.0f}s"

####################################################################################

sys.stdout.flush()

# find matching files
filename_suffix = ".corr." + ("dat" if args.binary_in else "txt.gz")

if args.chr is not None:
    filename_suffix = ".chr" + args.chr + ".corr." + ("dat" if args.binary_in else "txt.gz")

print(f"Looking for correlation files ending with *{filename_suffix} in {args.indir}")
files = glob.glob(args.indir + f"/*{filename_suffix}*")
n_files = len(files)
if n_files == 0:
    print("No matching files found.")
    exit()
files.sort() # make sure sample order is equal between chromosomes
print(f"\t{n_files:,} files found with such pattern")

# extract the genes.
print(f"Extracting eGenes and co-eGenes from input file")
featureset = set()
egeneset = set()
for file in files:
    if file.endswith(".dat"):
        fh = gzopen(file, mode='rb')
        _, _, _, _, file_features, file_is_egene = parse_binary_header(fh=fh)
        fh.close()
    else:
        fh = gzopen(file, mode='r')
        _, _, _, _, file_features, file_is_egene = parse_text_header(fh=fh)
        fh.close()

    file_egenes = [feature for feature, is_egene in zip(file_features, file_is_egene) if is_egene]
    featureset.update(file_features)
    egeneset.update(file_egenes)
print(f"\t{len(featureset):,} unique features found")
print(f"\t{len(egeneset):,} unique eGenes found")

print(f"Filtering eGenes and co-eGenes")
# read a set of genes that we want to limit to (if any)
if args.egenelist is not None:
    egeneset = egeneset.intersection(read_as_set(fpath=args.egenelist))
    print(f"\t{len(egeneset):,} eGenes after filtering on --egenelist input")

coegeneset = set(featureset)
if args.coegenelist is not None:
    coegeneset = coegeneset.intersection(read_as_set(fpath=args.coegenelist))
    print(f"\t{len(coegeneset):,} co-eGenes after filtering on --coegenelist input")

if args.geneannotation is not None and args.chr is not None:
    # limit the eGene set to only genes on the chromosome of interest.
    # Load the gene to chromosome annotation
    chr_genes = load_annotation(annotation=args.geneannotation, chr=args.chr)
    egeneset = egeneset.intersection(chr_genes)
    print(f"\t{len(egeneset):,} eGenes after filtering on chr '{args.chr}'")
    if len(egeneset) == 0:
        print("Specified to run on chromosome {}, but no genes in the data match to this string.".format(args.chr))
        exit()
    del chr_genes

# Sort the egenes and co-egenes.
egenes = sorted(list(egeneset))
egenes_indices = {egene: index for index, egene in enumerate(egenes)}
n_egenes = len(egenes)
coegenes = sorted(list(coegeneset))
coegenes_indices = {coegene: index for index, coegene in enumerate(coegenes)}
n_coegenes = len(coegenes)
del featureset, egeneset, coegeneset
print(f"\toutput will have {n_egenes:,} eGenes and {n_coegenes:,} co-eGenes.")

# Write the header to the output file.
fho = None
fho_cols = None
if args.binary_out:
    fho = open(args.out + ".dat", 'wb')
    fho.write(struct.pack('>i', n_files)) # first 'header number' in DoubleMatrixDataset format

    fho_cols = open(args.out + ".cols.txt", 'wt')
else:
    fho = gzopen(args.out,'w')
    fho.write("Sample") # first part of the header

sys.stdout.flush()

# # Write the gene pairs to the output file.
# print("Making gene pair index")
# i = 0
# gene_pairs_map = {}
# gene_pairs_index = 0
# for i, egene in enumerate(egenes):
#     for j, coegene in enumerate(coegenes):
#         if i == j:
#             continue
#         gene_pair = f"{egene}_{coegene}"
#         gene_pairs_map[gene_pair] = gene_pairs_index
#         if args.binary_out:
#             fho_cols.write(gene_pair + "\n")
#         else:
#             fho.write("\t" + gene_pair)  # add more to the header
#         gene_pairs_index += 1
#     if (i + 1) % 100 == 0:
#         print(f"\tMapped {i + 1:,}/{n_egenes:,} - {gene_pairs_index:,} pairs ...",end='\r')
# print(f"\tMapped {i + 1:,}/{n_egenes:,} - {gene_pairs_index:,} pairs found", end='\n')
# # this is a bit sketchy since this index can overflow the max index of an array/list...
# nrcols = gene_pairs_index
# del gene_pairs_index

# Write the gene pairs to the output file. We do not have to store a
# position dict here since the eGene and co-eGene lists are sorted
# and therefore we can calculate the position of the output value
# on the fly later.
print("Making gene pair index")
for egene in egenes:
    for coegene in coegenes:
        gene_pair = f"{egene}_{coegene}"
        if args.binary_out:
            fho_cols.write(gene_pair + "\n")
        else:
            fho.write("\t" + gene_pair)  # add more to the header
nrcols = n_egenes * n_coegenes
print(f"{nrcols} gene pairs")

if args.binary_out:
    fho_cols.close()
    fho.write(struct.pack('>i', nrcols)) # second 'header number' in DoubleMatrixDataset format
else:
    fho.write("\n") # write end of the header

sys.stdout.flush()

# iterate files
if args.binary_out:
    byteout = bytearray(N_BYTES_DOUBLE * nrcols) # initialize byte buffer once
    clear(byteout, binary=True) # replace 0 with nan
else:
    outln = ["NaN"] * nrcols # initialize string buffer once

filectr = 0
startTime = time.time_ns()

if args.binary_out:
    fho_rows = open(args.out + ".rows.txt", 'wt')

for file in files:
    # Get the sample ID.
    filename_prefix = os.path.basename(file).replace(filename_suffix, "")
    pool,sample = filename_prefix.split(".")[0:2]
    samplename = pool + "." + sample

    # Parse the file completely.
    fh = None
    header_info = None
    if file.endswith(".dat"):
        # Open the file and read the header.
        fh = gzopen(file, mode='rb')
        header_info = parse_binary_header(fh=fh)
    else:
        fh = gzopen(file, mode='r')
        header_info = parse_text_header(fh=fh)

    file_n_egenes, file_n_coegenes, file_n_correlations, _, file_features, _ = header_info
    print(f"parsing sample {samplename} with {file_n_egenes:,} eGenes, {file_n_coegenes:,} co-eGenes, and {len(file_features):,} features.")

    # Loop over the file.
    corrctr = 0
    writtenctr = 0
    for corrctr in range(file_n_correlations):
        if file.endswith(".dat"):
            value_index, beta = struct.unpack(">id", fh.read(12))
        else:
            values = fh.readline().rstrip("\n").split("\t")
            value_index = int(values[0])
            beta = float(values[1])
            del values
            
        # Translate the value index to i, j positions and find the corresponding gene names.
        i = math.floor(value_index / file_n_coegenes)
        j = value_index if i == 0 else value_index % file_n_coegenes
        feature_i = file_features[i]
        feature_j = file_features[j]
            
        # The correlation files only have gene_pairs and not their complements. We can check if we need to
        # save both by looking them both op on the gene_pairs_map. If this is the case, we can
        # store the same value twice.
        # out_indices = [gene_pairs_map.get(f"{feature_i}_{feature_j}"), gene_pairs_map.get(f"{feature_j}_{feature_i}")]

        # We look op the position of the genes in the output matrix. Correlation are only stored
        # one way (geneA-geneB = geneB-geneA). We therefore need to look up both combinations
        # of eGene and co-eGene.
        out_indices = []
        egene_i = egenes_indices.get(feature_i)
        coegene_j = coegenes_indices.get(feature_j)
        if egene_i is not None and coegene_j is not None:
            out_indices.append((egene_i * n_coegenes) + coegene_j)

        egene_j = egenes_indices.get(feature_j)
        coegene_i = coegenes_indices.get(feature_i)
        if egene_j is not None and coegene_i is not None:
            out_indices.append((egene_j * n_coegenes) + coegene_i)

        for out_index in out_indices:
            if out_index is None:
                continue
            print(f"{feature_i} [{i}]\t{feature_j} [{j}]\tvalue_index={value_index}\tbeta={beta}\tout={out_index}")

            # Store the correlation value.
            if args.binary_out:
                byteidx = out_index * N_BYTES_DOUBLE
                byteidxend = byteidx + N_BYTES_DOUBLE
                v = struct.pack('>d', beta) # TODO: I am unpacking and packing the same number here
                byteout[byteidx:byteidxend] = v[0:N_BYTES_DOUBLE]
            else:
                outln[out_index] = str(beta)

            writtenctr += 1
        if corrctr % 100000 == 0:
            print(f"\t{corrctr + 1:,} correlations parsed {writtenctr:,} written, file {filectr + 1:,}/{n_files:,}", end='\r')
    print(f"\t{corrctr + 1:,} correlations parsed {writtenctr:,} written, file {filectr + 1:,}/{n_files:,}", end='\n')

    fh.close()

    print(f"\tWriting sample")
    if args.binary_out:
        fho.write(byteout)
        fho_rows.write(samplename + "\n")
    else:
        fho.write(samplename+"\t")
        fho.write("\t".join(outln))
        fho.write("\n")
        fho.flush()

    filectr += 1

    current_time = time.time_ns()
    delta_time = current_time - startTime
    perc_done = (filectr / n_files) * 100

    time_per_sample = delta_time / filectr
    time_left = time_per_sample * (n_files - filectr)

    print(f"\n{perc_done:.2f} % done. Time spent: {to_time_str(delta_time)}. Time left: {to_time_str(time_left)}.")
    sys.stdout.flush()

if args.binary_out:
    fho_rows.close()
    
fho.close()
print("Done.")

profiler.disable()
stats = pstats.Stats(profiler).sort_stats('cumtime')
stats.print_stats(40)