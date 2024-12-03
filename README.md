# SIGNAL (SIngle-cell Group techNical vAriations Learning)

`SIGNAL` is a R package for data integration by learning group technical variations.
<img src="https://github.com/yzhou1999/SIGNAL/blob/main/docs/logo.jpg" width="120" style="display: inline;">

# Installation
To run `SIGNAL`, install from GitHub through ``devtools`` directly:
```R
install.packages('devtools')
library(devtools)
devtools::install_github("yzhou1999/SIGNAL")
```

# Usage

For usage examples, please see the `vignettes` directory of the repository.

* [Supervised integration using SIGNAL](https://yzhou1999.github.io/SIGNAL/articles/Supervised_integration.html)
* [SIGNAL integration for unlabeled datasets](https://yzhou1999.github.io/SIGNAL/articles/Unsupervised_integration.html)
* [Multi-scale analysis using SIGNAL](https://yzhou1999.github.io/SIGNAL/articles/Multiscale_analysis.html)
* [Knowledge transfer from reference to query using SIGNAL](https://yzhou1999.github.io/SIGNAL/articles/Knowledge_transfer.html)


# Dependencies
SIGNAL has been successfully installed and used on Windows, Linux and Mac OS (R version >= 4.0.2). The dependencies including: base, stats, Matrix, matrixStats, sparseMatrixStats, mclust, Rcpp, RcppEigen, irlba, BiocNeighbors, bigstatsr, dplyr, RSpectra.