
<!-- README.md is generated from README.Rmd. Please edit that file -->

# simstandard <img src="man/figures/logo.png" align="right" height="140/"/>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/simstandard)](https://cran.r-project.org/package=simstandard)
[![Rdoc](https://www.rdocumentation.org/badges/version/simstandard)](https://www.rdocumentation.org/packages/simstandard)  
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
| \-1.56 | \-1.44 | \-0.63 |   0.00 | \-0.78 | \-1.45 | \-0.32 | \-0.84 | \-0.28 | \-0.44 |   0.22 |   0.84 |
|   1.56 | \-1.03 |   0.67 | \-0.74 |   0.17 | \-1.25 | \-0.39 |   2.19 | \-0.02 |   0.91 | \-0.46 |   0.61 |
| \-0.22 | \-1.48 | \-1.19 | \-0.98 | \-1.50 | \-0.83 | \-1.62 |   0.19 | \-0.81 | \-0.22 |   0.15 | \-0.95 |
|   2.41 |   2.41 |   0.58 |   0.05 |   1.92 |   2.63 |   1.53 |   1.10 |   0.31 | \-0.34 | \-1.03 | \-0.57 |
| \-0.35 |   0.22 |   0.60 |   0.59 | \-0.10 |   0.38 | \-0.70 | \-0.54 | \-0.08 |   1.02 |   1.08 | \-1.01 |
|   1.48 |   1.32 |   1.03 |   0.55 |   0.97 |   1.18 |   0.72 |   0.89 |   0.38 |   0.60 |   0.05 | \-0.22 |

See more in the [tutorial for this
package](https://wjschne.github.io/simstandard/articles/simstandard_tutorial.html).
