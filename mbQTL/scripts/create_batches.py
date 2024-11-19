#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import os
import gzip

parser = argparse.ArgumentParser(description="")
parser.add_argument("--genes", required=True, type=str, help="")
parser.add_argument("--header_index", required=False, type=int, default=None, help="")
parser.add_argument("--gene_index", required=False, type=int, default=0, help="")
parser.add_argument("--annotation", required=False, type=str, default=None, help="")
parser.add_argument("--expgroups", required=False, type=str, help="")
parser.add_argument("--n_genes", required=False, type=int, default=200, help="")
parser.add_argument("--chunk_overflow", required=False, type=float, default=0.1, help="% of genes above --n_genes that is allowed")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

os.makedirs(os.path.dirname(args.out), exist_ok=True)

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


print("Loading genes...")
genes = set()
n_duplicated = 0
with gzopen(args.genes, mode='r') as fhin:
    for i, line in enumerate(fhin):
        values = line.rstrip("\n").split("\t")
        if args.header_index is not None and i == args.header_index:
            continue

        gene = values[args.gene_index]
        if gene in genes:
            print("\tWarning, gene '{}' is duplicated. Skipping gene.".format(gene))
            n_duplicated += 1
            continue

        genes.add(gene)
fhin.close()
print("  loaded {:,} genes, {:,} genes duplicates found".format(len(genes), n_duplicated))
del n_duplicated

groups = None
if args.expgroups is not None:
    print("Loading groups...")
    n_skipped = 0
    group_genes = set()
    groups = {}
    with gzopen(args.expgroups, mode='r') as f:
        for i, line in enumerate(f):
            gene, group = line.rstrip("\n").split("\t")
            if gene in group_genes:
                print("\tError, gene '{}' is duplicated.".format(gene))
                exit()

            if gene not in genes:
                n_skipped += 1
                continue

            if group not in groups:
                groups[group] = set()
            groups[group].add(gene)

            group_genes.add(gene)
    f.close()
    print("  Loaded {:,} groups with {:,} genes, {:,} genes were skipped".format(len(groups), len(group_genes), n_skipped))
    if group_genes != genes:
        print("  Warning, could not find group for {:,} / {:,} genes".format(len(genes) - len(group_genes), len(genes)))
        genes = group_genes
    del n_skipped, group_genes
else:
    groups = {gene: {gene} for gene in genes}

ordered_groups = None
group_order = None
if args.annotation is None:
    ordered_groups = {group_id: sorted(list(group_genes)) for group_id, group_genes in groups.items()}
    group_order = sorted(list(ordered_groups.keys()))
else:
    print("Loading gene annotation...")
    n_skipped = 0
    gene_to_pos = {}
    numeric_chrs = set()
    str_chrs = set()
    with gzopen(args.annotation, mode='r') as f:
        for i, line in enumerate(f):
            values = line.rstrip("\n").split("\t")
            if i == 0:
                pos = {value: index for index, value in enumerate(values)}
                continue

            gene = values[pos["Gene"]]
            if gene in gene_to_pos:
                print("\tError, gene '{}' is duplicated.".format(gene))
                exit()

            if gene not in genes:
                n_skipped += 1
                continue

            chr = values[pos["Chr"]]
            if chr.isnumeric():
                numeric_chrs.add(int(chr))
            else:
                str_chrs.add(chr)

            chrstart = values[pos["ChrStart"]]
            gene_to_pos[gene] = (chr, chrstart)
    f.close()
    print("  Loaded {:,} genes, skipped {:,} genes".format(len(gene_to_pos), n_skipped))
    if gene_to_pos.keys() != genes:
        print("  Warning, could not find annotation for {:,} / {:,} genes".format(len(genes) - len(gene_to_pos), len(genes)))

    chrs = [str(x) for x in sorted(list(numeric_chrs))] + sorted(list(str_chrs))
    n_chrs = len(chrs)
    chr_order = {chr: index for index, chr in enumerate(chrs)}
    print("  Loaded {:,} chromosomes: {}".format(len(chrs), ", ".join(chrs)))
    del n_skipped, chrs

    print("Order groups on chromosome - position")
    # Here we order the genes within the group as well as the groups. We are
    # not assuming a group is just on one chromosome though so we sort on the
    # first gene per group.
    ordered_groups = {}
    group_order = []
    n_genes = 0
    for group_id, group_genes in groups.items():
        ordered_group_genes = []
        for gene in group_genes:
            gene_pos = gene_to_pos.get(gene)
            if gene_pos is None:
                continue
            ordered_group_genes.append((gene, chr_order[gene_pos[0]], gene_pos[1]))
            n_genes += 1
        if len(ordered_group_genes) == 0:
            continue
        ordered_group_genes.sort(key=lambda x: (x[1], x[2], x[0]))
        ordered_groups[group_id] = [group_gene[0] for group_gene in ordered_group_genes]
        group_order.append((group_id, ordered_group_genes[0][1], ordered_group_genes[0][2]))
    group_order.sort(key=lambda x: (x[1], x[2], x[0]))
    group_order = [group[0] for group in group_order]
    print("  Sorted {:,} groups with {:,} genes in total".format(len(group_order), n_genes))
    del chr_order, n_genes

del groups, genes

print("Creating chunks ...")
max_chunk_size = args.n_genes + (0 if args.expgroups is None else (args.n_genes * args.chunk_overflow))
print("  maximum chunk size: {:,}".format(max_chunk_size))

chunk_id = 0
total_n_genes = 0
chunk = []
for group_id in group_order:
    group_genes = ordered_groups[group_id]
    n_genes = len(group_genes)
    if n_genes > args.n_genes:
        print("  Warning, group {} has {:,} genes with exceeds the intended chunk size of {:,}. This group will be saved as a single chunk.".format(group_id, n_genes, args.n_genes))

    # Check if adding this group would exceed the chunk size. It might have been
    # that even without this group we exceeded the chunk size but we want to keep
    # groups within one chunk so whatever.
    if len(chunk) + n_genes > max_chunk_size:
        print("  Writing chunk {} with {:,} genes".format(chunk_id, len(chunk)))
        fhout = gzopen(os.path.join(args.out + "-chunk" + str(chunk_id) + "-genes.txt"), mode='w')
        for gene in chunk:
            fhout.write(gene + "\n")
        fhout.close()
        chunk = []
        chunk_id += 1

    # Add the current group to the chunk.
    chunk.extend(group_genes)

    # Save the total number of genes.
    total_n_genes += n_genes

if len(chunk) > 0:
    print("  Writing chunk {} with {:,} genes".format(chunk_id, len(chunk)))
    fhout = gzopen(os.path.join(args.out + "-chunk" + str(chunk_id) + "-genes.txt"), mode='w')
    for gene in chunk:
        fhout.write(gene + "\n")
    fhout.close()
    chunk = []
    chunk_id += 1

print("  Saved {:,} genes into {:,} chunks".format(total_n_genes, chunk_id))

if chunk_id == 0:
    print("Error, no chunks created.")
    exit(1)

print("Done")