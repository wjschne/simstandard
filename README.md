
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

|    A1 |    A2 |    B1 |    B2 |     C |     A |     B | e\_A1 | e\_A2 | e\_B1 | e\_B2 |  d\_B |
|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|
| -1.31 | -1.03 | -0.80 | -0.97 |  1.08 | -0.11 |  0.55 | -1.26 | -0.94 | -1.13 | -1.35 |  0.63 |
|  0.38 |  1.24 |  0.64 |  0.75 | -0.73 |  0.00 |  1.20 |  0.38 |  1.24 | -0.08 | -0.08 |  1.20 |
| -0.83 | -0.55 |  2.12 |  1.21 | -0.25 |  0.12 |  1.02 | -0.88 | -0.65 |  1.51 |  0.50 |  0.93 |
|  1.13 | -0.23 |  1.51 |  1.40 |  0.92 |  0.40 |  1.31 |  0.93 | -0.55 |  0.73 |  0.48 |  0.99 |
|  0.90 | -0.62 | -0.85 | -0.28 | -0.30 | -0.19 | -0.43 |  0.99 | -0.47 | -0.60 |  0.02 | -0.28 |
|  0.05 |  0.63 |  0.90 |  1.14 |  1.03 |  0.46 |  0.76 | -0.18 |  0.26 |  0.45 |  0.61 |  0.39 |

See more in the [tutorial for this
package](https://wjschne.github.io/simstandard/articles/simstandard_tutorial.html).
