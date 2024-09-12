qvalue: Q-value estimation for false discovery rate control.

First examine that the distribution of the p-values and confirmed they are well-behaved. An important assumption behind the estimation performed in this package is that null p-values follow a Uniform(0,1) distribution. Therefore, the p-values should be relatively flat at the right tail of the histogram, an “U-shaped” p-value histogram is a red flag. U-shaped p-value histograms can indicate that a one-sided test was performed on data where there is signal in both directions, or it can indicate that there is dependence among the variables in the data.

source: Bioconductor’s qvalue package Version 2.36.0. John D. Storey and Andrew J.. Bass Princeton University.