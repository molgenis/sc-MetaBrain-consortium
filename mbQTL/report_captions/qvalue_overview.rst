qvalue: Q-value estimation for false discovery rate control.

The main purpose of the upper-left plot is to gauge the reliability of the π0 estimate, where the estimated π0 is plotted versus the tuning parameter λ. The variable λ is called lambda in the package; it can be fixed or automatically handled. As λ gets larger, the bias of the estimate decreases, yet the variance increases.Comparing your final estimate of π0 to this plot gives a good sense as to its quality. The remaining plots show how many tests are significant, as well as how many false positives to expect for each q-value cut-off.

source: Bioconductor’s qvalue package Version 2.36.0. John D. Storey and Andrew J.. Bass Princeton University.