<!-- badges: start -->
[![R-CMD-check](https://github.com/mengxu98/monocle/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mengxu98/monocle/actions/workflows/R-CMD-check.yaml)
[![version](https://img.shields.io/github/r-package/v/mengxu98/monocle?label=version&color=blue)](https://github.com/mengxu98/monocle/blob/master/DESCRIPTION)
<!-- badges: end -->

# monocle

Monocle is an analysis toolkit for single-cell RNA-Seq experiments. This repository is the [mengxu98/monocle](https://github.com/mengxu98/monocle) fork of [cole-trapnell-lab/monocle-release](https://github.com/cole-trapnell-lab/monocle-release), with performance-oriented maintenance (including Rcpp/OpenMP fast paths for negative-binomial testing and geometry).

To use this package you will need the R statistical computing environment and several Bioconductor/CRAN packages.

## Installation

Install the development version from GitHub with [pak](https://github.com/r-lib/pak):

```r
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}
pak::pak("mengxu98/monocle")
```

The Bioconductor release of Monocle remains available at:

https://bioconductor.org/packages/release/bioc/html/monocle.html

## News

Upstream project news historically appeared at:

http://cole-trapnell-lab.github.io/monocle-release/
