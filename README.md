
<!-- README.md is generated from README.Rmd. Please edit that file -->

# simstandard <a href="https://wjschne.github.io/simstandard/"><img src="man/figures/logo.png" align="right" height="120" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/wjschne/simstandard/workflows/R-CMD-check/badge.svg)](https://github.com/wjschne/simstandard/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/simstandard)](https://CRAN.R-project.org/package=simstandard)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
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

|    A1 |    A2 |    B1 |    B2 |     C |     A |     B |  e_A1 |  e_A2 |  e_B1 |  e_B2 |   d_B |
|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|
| -0.04 |  1.33 | -0.71 |  1.41 |  1.92 |  1.02 |  1.52 | -0.55 |  0.51 | -1.62 |  0.35 |  0.70 |
|  0.84 | -0.06 | -0.70 | -0.79 | -1.33 | -0.79 | -0.83 |  1.23 |  0.57 | -0.20 | -0.21 | -0.19 |
|  1.84 |  0.45 | -0.14 |  0.19 |  0.59 |  1.38 |  0.93 |  1.16 | -0.65 | -0.70 | -0.46 | -0.17 |
| -1.88 | -0.73 | -1.73 | -0.54 | -1.00 | -1.15 | -0.62 | -1.31 |  0.19 | -1.36 | -0.11 |  0.30 |
|  0.24 |  0.01 |  1.10 | -0.01 |  1.07 |  0.51 | -0.03 | -0.01 | -0.40 |  1.11 |  0.01 | -0.44 |
|  0.51 |  0.75 |  1.65 |  1.48 |  2.45 |  1.62 |  1.50 | -0.29 | -0.54 |  0.75 |  0.44 |  0.20 |

See more in the [tutorial for this
package](https://wjschne.github.io/simstandard//articles/simstandard_tutorial.html).
