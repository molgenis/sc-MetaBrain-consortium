#!/usr/bin/env Rscript
# Author: M. Vochteloo

# options parser
.libPaths("/usr/local/lib/R/site-library")
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list <- list(
  make_option(c("--input"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--nom_pvalue"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--n_tests"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--beta_adj_pvalue"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--beta_dist_a"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--beta_dist_b"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--bonf_pvalue_col"), action="store", default="bonf_pvalue", type='character',
              help=""),
  make_option(c("--bonf_bh_fdr_col"), action="store", default="bonf_bh_fdr", type='character',
              help=""),
  make_option(c("--bh_fdr_col"), action="store", default="bh_fdr", type='character',
              help=""),
  make_option(c("--qvalue_col"), action="store", default="qvalue", type='character',
              help=""),
  make_option(c("--nom_thresh_col"), action="store", default="signif_threshold", type='character',
              help=""),
  make_option(c("--alpha"), action="store", default=0.05, type='numeric',
              help=""),
  make_option(c("--data_out"), action="store", default=NA, type='character',
              help="Output main directory"),
  make_option(c("--plot_out"), action="store", default=NA, type='character',
              help="Output plot directory"),
  make_option(c("--suffix"), action="store", default="WithMulTest", type='character',
              help="Output suffix")
)
opt <- parse_args(OptionParser(option_list=option_list))

print("Options in effect:")
for (name in names(opt)) {
	print(paste0("  --", name, " ", opt[[name]]))
}
print("")

# check date is provided
if (any(unlist(lapply(c(opt$input, opt$nom_pvalue, opt$data_out, opt$plot_out), is.na)))) {
  stop("required parameters (input, nom_pvalue, data_out & plot_out) must be provided.")
}

shhh(library(readr))
shhh(library(ggplot2))
shhh(library(qvalue))

print("Loading data")
data <- as.data.frame(read_delim(opt$input, delim = "\t"))

cols <- colnames(data)
data <- data[!is.na(data[[opt$nom_pvalue]]), ]
print(paste0("  ", nrow(data), " genes loaded"))
print(paste0("  significant at nominal pvalue<", opt$alpha, ": ", sum(data[[opt$nom_pvalue]] < opt$alpha), " genes"))

data <- data[order(data[[opt$nom_pvalue]]), ]

pvalue <- NULL
if (!is.na(opt$beta_adj_pvalue) & (opt$beta_adj_pvalue %in% colnames(data))) {
    print(paste0("  significant at beta adujusted pvalue<", opt$alpha, ": ", sum(data[[opt$beta_adj_pvalue]] < opt$alpha), " genes"))
    pvalue <- opt$beta_adj_pvalue
} else if (!is.na(opt$nom_pvalue) & (opt$nom_pvalue %in% colnames(data))) {
    pvalue <- opt$nom_pvalue
} else {
    print("Error, both nom_pvalue and beta_adj_pvalue are unavailable.")
    quit()
}
print(paste0("Using ", pvalue, " as p-value columns for multiple testing correction."))


if (!is.na(opt$n_tests) & (opt$n_tests %in% colnames(data))) {
    print("Calculating bonferroni pvalues")
    data[[opt$bonf_pvalue]] <- data[[opt$nom_pvalue]] * data[[opt$n_tests]]
    data[data$bonf_pvalue > 1] <- 1
    print(paste0("  significant at bonferroni p-value<", opt$alpha, ": ", sum(data[[opt$bonf_pvalue]] < opt$alpha), " genes"))
    cols <- c(cols, opt$bonf_pvalue)

    print("Calculating two-step FDR")
    data[[opt$bonf_bh_fdr_col]] <- p.adjust(data[[opt$bonf_pvalue]], method = 'hochberg', n = length(data[[opt$bonf_pvalue]]))
    print(paste0("  significant at two step FDR<", opt$alpha, ": ", sum(data[[opt$bonf_bh_fdr_col]] < opt$alpha), " genes"))
    cols <- c(cols, opt$bonf_bh_fdr_col)
}

print("Calculating Benjamini-Hochberg FDR")
data[[opt$bh_fdr_col]] <- p.adjust(data[[pvalue]], method = 'hochberg', n = length(data[[pvalue]]))
print(paste0("  significant at BH-FDR<", opt$alpha, ": ", sum(data[[opt$bh_fdr_col]] < opt$alpha), " genes"))
cols <- c(cols, opt$bh_fdr_col)

print("Calculating qvalues")
qobj <- qvalue(p = data[[pvalue]])
data[[opt$qvalue_col]] <- qobj$qvalues
print(paste0("  significant at qvalue<", opt$alpha, ": ", sum(data[[opt$qvalue_col]] < opt$alpha), " genes"))
cols <- c(cols, opt$qvalue_col)

# got this from: https://github.com/francois-a/fastqtl/blob/master/R/calculateSignificanceFastQTL.R
# determine global min(p) significance threshold and calculate nominal p-value threshold for each gene
ub <- sort(data[data[[opt$qvalue_col]] > opt$alpha, pvalue])[1]  # smallest p-value above FDR
lb <- -sort(-data[data[[opt$qvalue_col]] <= opt$alpha, pvalue])[1]  # largest p-value below FDR
pthreshold <- (lb + ub) / 2
print(paste0("  * min p-value threshold @ qvalue<", opt$alpha, ": ", round(pthreshold, 2)))

if (!is.na(opt$beta_dist_a) & !is.na(opt$beta_dist_b) & (opt$beta_dist_a %in% colnames(data)) & (opt$beta_dist_b %in% colnames(data))) {
    print("Calculating nominal threshold")
    data[[opt$nom_thresh_col]] <- signif(qbeta(pthreshold, data[[opt$beta_dist_a]], data[[opt$beta_dist_b]], ncp=0, lower.tail=TRUE, log.p=FALSE), 6)
    cols <- c(cols, opt$nom_thresh_col)
}
data <- data[order(data[[opt$qvalue]]), cols]

print("Saving data")
write.table(data, paste0(opt$data_out, opt$suffix, ".txt"), quote = F, sep = "\t", row.names = FALSE)

dir.create(dirname(opt$plot_out), recursive = TRUE, showWarnings = FALSE)
setwd(dirname(opt$plot_out))

print("Creating figures")
png(paste0(opt$plot_out, "_overview.png"))
plot1 <- plot(qobj)
dev.off()

plot2 <- hist(qobj)
ggsave(plot2, filename = paste0(opt$plot_out, "_hist.png"), width = 29.7, height = 21 ,units = c("cm"))

print("Done")