#!/usr/bin/env Rscript
# Author: M. Vochteloo

# options parser
.libPaths("/usr/local/lib/R/site-library")
shhh <- suppressPackageStartupMessages
shhh(library(argparse))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--exp", required=TRUE, type="character", help="")
parser$add_argument("--min_cpm", required=FALSE, type="numeric", default=2., help="")
parser$add_argument("--out", required=TRUE, type="character", help="")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

print("Options in effect:")
for (name in names(args)) {
	print(paste0("  --", name, " ", args[[name]]))
}
print("")

shhh(library(readr))
shhh(library(edgeR))
shhh(library(limma))

print("Loading expression data")
expression <- read.table(file = args$exp, row.names = 1, check.names = FALSE)
print(paste0("  Loaded expression with shape (", nrow(expression), ", ", ncol(expression), ")"))
n_genes_input = dim(expression)[1]

#################################
############# cpm.R #############
#################################

# filter out low-expression genes
y <- DGEList(counts = expression)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=F]
n_filterbyexpr_removed <- dim(y)[1]
print(paste0("  Removed ", n_genes_input - n_filterbyexpr_removed, " genes by filterByExpr()"))

# TMM normalization
y <- calcNormFactors(y, method = "TMM")

# voom transform
v <- voom(y, plot=F)
logcpm <- v$E

# remove genes if mean log2CPM < min_cpm
mean_logcpm <- apply(logcpm, 1, mean)
logcpm <- logcpm[mean_logcpm > args$min_cpm,]
n_cpm_removed <- dim(logcpm)[1]
print(paste0("  Removed ", n_filterbyexpr_removed - n_cpm_removed, " genes due to >", args$min_cpm, " average CPM filter"))

##########################################
############# adjust-batch.R #############
##########################################

# # read batch info
# b <- read.table(batch_file, sep="\t", header=T)
# batches <- b$batch
# names(batches) <- b$donor
#
# # reorder batches to have the same order with log2cpm
# batches <- batches[colnames(logcpm)]
# stopifnot(names(batches) == colnames(logcpm))
#
# # adjust for batch effect
# if(all(table(batches) == 1)){
#   res <- logcpm  # if there is no meaningful batch
# } else {
#   res <- ComBat(dat = logcpm, batch = batches)
# }

########################################
############# quant-norm.R #############
########################################

# quantile normalizarion
logcpm <- t(apply(logcpm, 1, rank, ties.method = "average"))
logcpm <- qnorm(logcpm / (ncol(logcpm) + 1))

########################################

print("Writing output")
write.table(logcpm, gzfile(paste0(args$out, "log2CPM_QN.tsv.gz")), quote = F, sep = "\t", row.names=TRUE, col.names=NA)
print(paste0("  Written expression with shape (", nrow(expression), ", ", ncol(expression), ")"))

filter.stats <- data.frame(c("Input", "FilterByExp", "CPMFilter"), c(n_genes_input, n_filterbyexpr_removed, n_cpm_removed))
colnames(filter.stats) <- c("filter", "ngenes")
write.table(filter.stats, paste0(args$out, "log2CPM_QN.stats.tsv"), quote = F, sep = "\t", row.names=TRUE, col.names=NA)

print("Done")