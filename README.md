
<!-- README.md is generated from README.Rmd. Please edit that file -->

# simstandard <img src="man/figures/logo.png" align="right" height="140/"/>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/simstandard)](https://cran.r-project.org/package=simstandard)
[![Rdoc](https://www.rdocumentation.org/badges/version/simstandard)](https://www.rdocumentation.org/packages/simstandard)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![Travis build
status](https://travis-ci.org/wjschne/simstandard.svg?branch=master)](https://travis-ci.org/wjschne/simstandard)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/wjschne/simstandard?branch=master&svg=true)](https://ci.appveyor.com/project/wjschne/simstandard)
[![Coverage
status](https://codecov.io/gh/wjschne/simstandard/branch/master/graph/badge.svg)](https://codecov.io/github/wjschne/simstandard?branch=master)
[![R-CMD-check](https://github.com/wjschne/simstandard/workflows/R-CMD-check/badge.svg)](https://github.com/wjschne/simstandard/actions)
<!-- badges: end -->

Sometimes you have a structural model with standardized path
coefficients, structural coefficients, and correlations, but you do not
know the error and disturbance variances. The purpose of `simstandard`
is to calculate these variances and then simulate multivariate normal
data based on your model.

## Installation

You can either install simstandard from CRAN or install the development
version of simstandard from github.

### Option 1: Install the most recent stable release from CRAN

You can install simstandard from CRAN by running this code:

``` r
install.packages("simstandard")
```

### Option 2: Install the development version from GitHub

To install the development version of simstandard, you need to check if
devtools is installed. If not, run this:

``` r
install.packages("devtools")
```

Once you are sure you have devtools installed, you can install the
development version of simstandard from GitHub by running this code:

``` r
devtools::install_github("wjschne/simstandard")
```

## Example

The `simstandard` package uses [lavaan
syntax](https://lavaan.ugent.be/tutorial/syntax1.html) to specify
models.

``` r
library(simstandard)
model <- "
A =~ 0.5 * A1 + 0.8 * A2
B =~ 0.6 * B1 + 0.7 * B2
B ~ 0.8 * A
C ~~ 0.5 * A
"
data <- sim_standardized(m = model, n = 500)

knitr::kable(head(data), digits = 2)
```

|     A1 |     A2 |     B1 |     B2 |      C |      A |      B |  e\_A1 |  e\_A2 |  e\_B1 |  e\_B2 |   d\_B |
| -----: | -----: | -----: | -----: | -----: | -----: | -----: | -----: | -----: | -----: | -----: | -----: |
|   1.18 | \-0.56 | \-0.65 | \-0.22 |   0.35 | \-0.22 | \-0.08 |   1.29 | \-0.38 | \-0.60 | \-0.16 |   0.09 |
|   0.36 | \-0.42 | \-1.73 | \-0.38 | \-0.43 |   0.51 |   0.00 |   0.10 | \-0.83 | \-1.74 | \-0.38 | \-0.41 |
| \-0.19 |   0.24 | \-0.59 | \-0.59 | \-0.12 |   0.51 | \-0.76 | \-0.44 | \-0.17 | \-0.13 | \-0.06 | \-1.17 |
|   1.93 |   0.09 |   0.73 |   1.51 | \-0.33 |   1.42 |   1.98 |   1.22 | \-1.04 | \-0.45 |   0.13 |   0.84 |
|   0.76 | \-0.24 |   0.05 |   0.22 | \-1.04 |   0.44 |   0.48 |   0.54 | \-0.59 | \-0.24 | \-0.12 |   0.13 |
|   0.54 | \-0.69 | \-0.19 | \-1.04 | \-0.90 | \-0.60 | \-1.32 |   0.84 | \-0.21 |   0.60 | \-0.12 | \-0.84 |

See more in the [tutorial for this
package](https://wjschne.github.io/simstandard/articles/simstandard_tutorial.html).
