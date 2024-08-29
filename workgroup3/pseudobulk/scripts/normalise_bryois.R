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
parser$add_argument("--min_obs", required=FALSE, type="numeric", default=10, help="")
parser$add_argument("--min_cpm", required=FALSE, type="numeric", default=1, help="")
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
shhh(library(matrixStats))
shhh(library(dplyr))
shhh(library(edgeR))

print("Loading expression data")
expression <- read.table(file = args$exp, row.names = 1, check.names = FALSE)
print(paste0("  Loaded expression with shape (", nrow(expression), ", ", ncol(expression), ")"))
n_genes_input = dim(expression)[1]

##Filter based on gene variance
print("Filtering genes")
rowVarInfo <- rowVars(as.matrix(expression))
expression <- expression[which(rowVarInfo != 0),]
n_var_removed <- dim(expression)[1]
print(paste0("  Removed ", n_genes_input - n_var_removed, " genes due to no variance"))

##Filter on expressed in 10 individuals
expression <- expression[which(rowSums(expression > 0) > args$min_obs), ]
n_obs_removed <- dim(expression)[1]
print(paste0("  Removed ", n_var_removed - n_obs_removed, " genes due to >0 counts in >", args$min_obs, " individuals filter"))

##Filter on CPM
expression <- expression[which(rowMeans(expression * 10^6 / colSums(expression)) > args$min_cpm), ]
n_cpm_removed <- dim(expression)[1]
print(paste0("  Removed ", n_obs_removed - n_cpm_removed, " genes due to >", args$min_cpm, " average CPM filter"))

print("Normalise expression")
expression <- DGEList(counts = expression) %>% #TMM normalize the data using edgeR
    edgeR::calcNormFactors(.) %>%
    edgeR::cpm(.) %>%
    as.data.frame()

print("Writing output")
write.table(expression, gzfile(paste0(args$out, "TMM.tsv.gz")), quote = F, sep = "\t", row.names=TRUE, col.names=NA)
print(paste0("  Written expression with shape (", nrow(expression), ", ", ncol(expression), ")"))

filter.stats <- data.frame(c("Input", "Varfilter", "MinObsFilter", "CPMFilter"), c(n_genes_input, n_var_removed, n_obs_removed, n_cpm_removed))
colnames(filter.stats) <- c("filter", "ngenes")
write.table(filter.stats, paste0(args$out, "TMM.stats.tsv"), quote = F, sep = "\t", row.names=TRUE, col.names=NA)

print("Done")