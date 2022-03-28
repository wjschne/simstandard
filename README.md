
<!-- README.md is generated from README.Rmd. Please edit that file -->

# simstandard [<img src="man/figures/logo.png" style="float: right;" width="140" alt="logo" />](https://wjschne.github.io/simstandard/)

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
| -0.18 | -0.36 | -1.73 |  0.06 |  0.91 | -0.54 | -1.24 |  0.10 |  0.08 | -0.99 |  0.93 | -0.80 |
| -1.94 | -1.27 | -0.11 | -1.07 | -1.68 | -2.67 | -1.98 | -0.60 |  0.87 |  1.08 |  0.31 |  0.16 |
| -2.67 | -1.10 | -0.87 | -0.76 | -1.22 | -1.76 | -1.34 | -1.79 |  0.30 | -0.07 |  0.17 |  0.07 |
|  0.60 |  1.15 |  1.55 |  1.30 | -0.44 |  1.45 |  1.48 | -0.12 |  0.00 |  0.66 |  0.26 |  0.32 |
| -0.58 | -0.18 | -0.37 | -1.33 |  0.01 | -1.22 | -0.93 |  0.03 |  0.79 |  0.19 | -0.68 |  0.05 |
|  1.35 |  1.09 |  1.07 |  0.64 |  1.59 |  1.93 |  1.30 |  0.39 | -0.45 |  0.30 | -0.27 | -0.25 |

See more in the [tutorial for this
package](https://wjschne.github.io/simstandard//articles/simstandard_tutorial.html).
