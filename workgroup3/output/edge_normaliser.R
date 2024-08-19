#!/usr/bin/env Rscript

# options parser
.libPaths("/usr/local/lib/R/site-library")
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list <- list(
  make_option(c("--countsFile"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--filterFile"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--out"), action="store", default=NA, type='character',
              help="")
)
opt <- parse_args(OptionParser(option_list=option_list))

print("Options in effect:")
for (name in names(opt)) {
	print(paste0("  --", name, " ", opt[[name]]))
}
print("")

shhh(library(data.table))
shhh(library(dplyr))
shhh(library(edgeR))

print("Loading data")
counts <- read.table(opt$countsFile, check.names=FALSE)

if (!is.na(opt$filterFile)) {
    print(dim(counts))
    filter <- read.table(opt$filterFile, sep="\t")
    genes <- data.frame(do.call('rbind', strsplit(as.character(filter$V1),'_',fixed=TRUE)))
    counts <- counts[rownames(counts) %in% genes$X1, ]
    print(dim(counts))
}

print(paste0('TMM/PM normalising expression matrix'))
normalised_counts <- counts %>%
    edgeR::DGEList(counts = .) %>% #TMM normalize the data using edgeR
    edgeR::calcNormFactors(.) %>%
    edgeR::cpm(.) %>%
    as.data.frame()

print("Writing normalised pseudo bulk data")
write.table(normalised_counts, file = paste0(opt$out, "TMM.txt"), sep = "\t", row.names=TRUE, quote=FALSE, col.names = NA)