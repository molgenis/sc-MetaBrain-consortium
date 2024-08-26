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

#################################
############# cpm.R #############
#################################

# filter out low-expression genes
y <- DGEList(counts = expression)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=F]

# TMM normalization
y <- calcNormFactors(y, method = "TMM")

# voom transform
v <- voom(y, plot=F)
logcpm <- v$E

# remove genes if mean log2CPM < min_cpm
mean_logcpm <- apply(logcpm, 1, mean)
logcpm <- logcpm[mean_logcpm > args$min_cpm,]

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
write.table(logcpm, gzfile(args$out), quote = F, sep = "\t", row.names=TRUE, col.names=NA)
print(paste0("  Written expression with shape (", nrow(expression), ", ", ncol(expression), ")"))

print("Done")