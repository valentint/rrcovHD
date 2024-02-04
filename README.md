
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `rrcovHD`: Robust Multivariate Methods for High Dimensional Data

<!-- badges: start -->

[![CRAN
version](https://www.r-pkg.org/badges/version/rrcovHD)](https://cran.r-project.org/package=rrcovHD)
[![R-CMD-check](https://github.com/valentint/rrcovHD/workflows/R-CMD-check/badge.svg)](https://github.com/valentint/rrcovHD/actions)
[![downloads](https://cranlogs.r-pkg.org/badges/rrcovHD)](https://cran.r-project.org/package=rrcovHD)
[![downloads](https://cranlogs.r-pkg.org/badges/grand-total/rrcovHD)](https://cran.r-project.org/package=rrcovHD)
[![license](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![Codecov test
coverage](https://codecov.io/gh/valentint/rrcovHD/branch/master/graph/badge.svg)](https://app.codecov.io/gh/valentint/rrcovHD?branch=master)
<!-- badges: end -->

The package `rrcovHD` provides robust multivariate methods for high
dimensional data including outlier detection (Filzmoser and Todorov
(2013) <doi:10.1016/j.ins.2012.10.017>), robust sparse PCA (Croux et
al. (2013) <doi:10.1080/00401706.2012.727746>, Todorov and Filzmoser
(2013) <doi:10.1007/978-3-642-33042-1_31>), robust PLS (Todorov and
Filzmoser (2014) <doi:10.17713/ajs.v43i4.44>), and robust sparse
classification (Ortner et al. (2020) <doi:10.1007/s10618-019-00666-8>).

## Installation

The `rrcovHD` package is on CRAN (The Comprehensive R Archive Network)
and the latest release can be easily installed using the command

    install.packages("rrcovHD")
    library(rrcovNA)

## Building from source

To install the latest stable development version from GitHub, you can
pull this repository and install it using

    ## install.packages("remotes")
    remotes::install_github("valentint/rrcovHD")

Of course, if you have already installed `remotes`, you can skip the
first line (I have commented it out).

## Example

This is a basic example which shows you if the package is properly
installed:

``` r

library(rrcovHD)
#> Loading required package: rrcov
#> Loading required package: robustbase
#> Scalable Robust Estimators with High Breakdown Point (version 1.7-5)
#> Robust Multivariate Methods for High Dimensional Data (version 0.2-7)

data(pottery)
dim(pottery)        # 27 observations in 2 classes, 6 variables
#> [1] 27  7
head(pottery)
#>     SI   AL   FE  MG   CA   TI origin
#> 1 55.8 14.0 10.2 4.9  5.0 0.88  Attic
#> 2 51.2 12.5 10.1 4.4  4.8 0.86  Attic
#> 3 57.1 14.0  8.3 6.4 11.2 0.75  Attic
#> 4 53.8 13.1  9.3 4.9  6.6 0.81  Attic
#> 5 59.4 14.8  9.8 5.5  5.4 0.89  Attic
#> 6 56.2 14.0  9.9 4.9  5.4 0.89  Attic

## Build the SIMCA model. Use RSimca for a robust version
rs <- RSimca(origin~., data=pottery)
rs
#> Call:
#> RSimca(origin ~ ., data = pottery)
#> 
#> Prior Probabilities of Groups:
#>     Attic  Eritrean 
#> 0.4814815 0.5185185 
#> 
#> Pca objects for Groups:
#> 
#> Call:
#> PcaHubert(x = class, k = k[i], kmax = kmax[i], trace = trace)
#> Importance of components:
#>                           PC1    PC2
#> Standard deviation     5.2672 0.8564
#> Proportion of Variance 0.7186 0.1804
#> Cumulative Proportion  0.7186 0.8990
#> 
#> Call:
#> PcaHubert(x = class, k = k[i], kmax = kmax[i], trace = trace)
#> Importance of components:
#>                           PC1
#> Standard deviation     3.2934
#> Proportion of Variance 0.8102
#> Cumulative Proportion  0.8102
summary(rs)
#> 
#> Call:
#> RSimca(formula = origin ~ ., data = pottery)
#> 
#> Prior Probabilities of Groups:
#>     Attic  Eritrean 
#> 0.4814815 0.5185185 
#> 
#> Pca objects for Groups:
#> 
#> Call:
#> PcaHubert(x = class, k = k[i], kmax = kmax[i], trace = trace)
#> Importance of components:
#>                           PC1    PC2
#> Standard deviation     5.2672 0.8564
#> Proportion of Variance 0.7186 0.1804
#> Cumulative Proportion  0.7186 0.8990
#> 
#> Call:
#> PcaHubert(x = class, k = k[i], kmax = kmax[i], trace = trace)
#> Importance of components:
#>                           PC1
#> Standard deviation     3.2934
#> Proportion of Variance 0.8102
#> Cumulative Proportion  0.8102
```

## Community guidelines

### Report issues and request features

If you experience any bugs or issues or if you have any suggestions for
additional features, please submit an issue via the
[*Issues*](https://github.com/valentint/rrcovHD/issues) tab of this
repository. Please have a look at existing issues first to see if your
problem or feature request has already been discussed.

### Contribute to the package

If you want to contribute to the package, you can fork this repository
and create a pull request after implementing the desired functionality.

### Ask for help

If you need help using the package, or if you are interested in
collaborations related to this project, please get in touch with the
package maintainer.
