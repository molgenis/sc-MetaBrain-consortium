#!/usr/bin/env Rscript
# Author: M. Vochteloo

# options parser
.libPaths("/usr/local/lib/R/site-library")
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list <- list(
  make_option(c("--data1"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--name1"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--data2"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--name2"), action="store", default=NA, type='character',
              help=""))
opt <- parse_args(OptionParser(option_list=option_list))

print("Options in effect:")
for (name in names(opt)) {
	print(paste0("  --", name, " ", opt[[name]]))
}
print("")

shhh(library(Seurat))

seurat1 <- readRDS(opt$data1)
print(str(seurat1))
counts1 <- GetAssayData(object = seurat1, assay = "RNA", slot = "counts") # 32547 21283
data1 <- GetAssayData(object = seurat1, assay = "data", slot = "data") # 32547 21283

head(counts1)

seurat2 <- readRDS(opt$data2)
print(str(seurat2))
counts2 <- GetAssayData(object = seurat2, assay = "RNA", slot = "counts") # 32531 18553
data2 <- GetAssayData(object = seurat2, assay = "data", slot = "data") # 32531 18553

head(counts2)

col.overlap <- intersect(colnames(counts1), colnames(counts2)) #17664
row.overlap <- intersect(rownames(counts1), rownames(counts2)) #32470

overlap1 <- counts1[row.overlap, col.overlap]
overlap2 <- counts2[row.overlap, col.overlap]
identical <- overlap1 == overlap2
pcnt.identical <- sum(identical@x) / length(identical@x)
pcnt.identical

different1 <- overlap1[!identical@x]
different1 <- overlap1[!identical@x]