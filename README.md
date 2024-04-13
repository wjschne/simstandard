
<!-- README.md is generated from README.Rmd. Please edit that file -->

# simstandard <a href="https://wjschne.github.io/simstandard/"><img src="man/figures/logo.png" align="right" height="120" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/wjschne/simstandard/workflows/R-CMD-check/badge.svg)](https://github.com/wjschne/simstandard/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/simstandard)](https://CRAN.R-project.org/package=simstandard)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Codecov test
coverage](https://codecov.io/gh/wjschne/simstandard/branch/main/graph/badge.svg)](https://app.codecov.io/gh/wjschne/simstandard?branch=main)
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
the remotes packages is installed. If not, run this:

``` r
install.packages("remotes")
```

Once you are sure you have the remotes package installed, you can
install the development version of simstandard from GitHub by running
this code:

``` r
remotes::install_github("wjschne/simstandard")
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

|    A1 |   A2 |    B1 |    B2 |     C |     A |     B |  e_A1 |  e_A2 |  e_B1 |  e_B2 |   d_B |
|------:|-----:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|
| -0.72 | 0.70 |  0.53 | -0.81 | -1.05 |  0.65 | -0.22 | -1.04 |  0.18 |  0.66 | -0.66 | -0.74 |
| -0.13 | 0.26 |  0.48 |  0.66 |  1.51 |  1.51 |  1.89 | -0.89 | -0.95 | -0.65 | -0.66 |  0.68 |
|  0.67 | 0.39 | -0.53 |  0.28 |  0.56 |  0.21 | -0.37 |  0.56 |  0.22 | -0.31 |  0.54 | -0.53 |
| -0.55 | 0.15 |  0.82 | -1.60 |  0.71 | -0.27 | -0.46 | -0.42 |  0.36 |  1.10 | -1.28 | -0.25 |
|  0.62 | 0.10 | -1.36 | -0.78 | -0.87 | -0.46 |  0.11 |  0.86 |  0.47 | -1.43 | -0.86 |  0.48 |
|  1.56 | 1.85 |  1.01 | -0.03 |  0.43 |  1.35 |  1.29 |  0.89 |  0.77 |  0.24 | -0.93 |  0.21 |

See more in the [tutorial for this
package](https://wjschne.github.io/simstandard//articles/simstandard_tutorial.html).
