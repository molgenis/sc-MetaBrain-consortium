#!/usr/bin/env Rscript
# Author: M. Vochteloo

# options parser
.libPaths("/usr/local/lib/R/site-library")
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list <- list(
  make_option(c("--input"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--beta_dist_a"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--beta_dist_b"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--nom_threshold"), action="store", default="signif_threshold", type='character',
              help=""),
  make_option(c("--pvalue"), action="store", default="pvalue", type='character',
              help=""),
  make_option(c("--qvalue"), action="store", default="qvalue", type='character',
              help=""),
  make_option(c("--alpha"), action="store", default=0.05, type='numeric',
              help=""),
  make_option(c("--data_out"), action="store", default=NA, type='character',
              help="Output main directory"),
  make_option(c("--plot_out"), action="store", default=NA, type='character',
              help="Output plot directory"),
  make_option(c("--suffix"), action="store", default="_qvalues", type='character',
              help="Output suffix")
)
opt <- parse_args(OptionParser(option_list=option_list))

print("Options in effect:")
for (name in names(opt)) {
	print(paste0("  --", name, " ", opt[[name]]))
}
print("")

# check date is provided
if (any(unlist(lapply(c(opt$input, opt$data_out, opt$plot_out), is.na)))) {
  stop("required parameters (input, data_out & plot_out) must be provided.")
}

shhh(library(readr))
shhh(library(ggplot2))
shhh(library(qvalue))

print("Loading data")
data <- as.data.frame(read_delim(opt$input, delim = "\t"))
data <- data[order(data[[opt$pvalue]]), ]

cols <- colnames(data)
data <- data[!is.na(data[[opt$pvalue]]), ]
print(paste0("  ", nrow(data), " genes loaded"))
print(paste0("  significant at pvalue<", opt$alpha, ": ", sum(data[[opt$pvalue]] < opt$alpha), " genes"))

print("Calculating qvalues")
qobj <- qvalue(p = data[[opt$pvalue]])
data[[opt$qvalue]] <- qobj$qvalues
print(paste0("  significant at qvalue<", opt$alpha, ": ", sum(data[[opt$qvalue]] < opt$alpha), " genes"))

# got this from: https://github.com/francois-a/fastqtl/blob/master/R/calculateSignificanceFastQTL.R
# determine global min(p) significance threshold and calculate nominal p-value threshold for each gene
ub <- sort(data[data[[opt$qvalue]] > opt$alpha, opt$pvalue])[1]  # smallest p-value above FDR
lb <- -sort(-data[data[[opt$qvalue]] <= opt$alpha, opt$pvalue])[1]  # largest p-value below FDR
pthreshold <- (lb + ub) / 2
print(paste0("  * min p-value threshold @ FDR ", opt$alpha, ": ", round(pthreshold, 2)))

if (!is.na(opt$beta_dist_a) & !is.na(opt$beta_dist_b)) {
    print("Calculating nominal threshold")
    data[[opt$nom_threshold]] <- signif(qbeta(pthreshold, data[[opt$beta_dist_a]], data[[opt$beta_dist_b]], ncp=0, lower.tail=TRUE, log.p=FALSE), 6)
    cols <- c(cols, opt$nom_threshold, opt$qvalue)
} else {
    cols <- c(cols, opt$qvalue)
}
data <- data[order(data[[opt$qvalue]]), cols]

print("Saving data")
write.table(data, paste0(opt$data_out, opt$suffix, ".txt"), quote = F, sep = "\t", row.names = FALSE)

dir.create(dirname(opt$plot_out), recursive = TRUE, showWarnings = FALSE)
setwd(dirname(opt$plot_out))

print("Creating figures")
plot1 <- plot(qobj)
plot2 <- hist(qobj)
ggsave(plot2, filename = paste0(opt$plot_out, "_hist.png"), width = 29.7, height = 21 ,units = c("cm"))

print("Done")