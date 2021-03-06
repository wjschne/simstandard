
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
| -1.55 | -0.64 | -0.44 | -1.74 |  1.26 | -1.78 | -1.02 | -0.66 |  0.78 |  0.17 | -1.03 |  0.41 |
| -1.83 | -1.05 | -0.19 |  0.54 |  0.11 | -0.79 | -0.62 | -1.43 | -0.42 |  0.18 |  0.97 |  0.01 |
|  0.77 | -0.75 | -0.92 | -0.19 | -0.99 | -0.70 | -0.43 |  1.12 | -0.19 | -0.66 |  0.11 |  0.12 |
|  1.06 |  1.43 | -0.10 |  0.61 |  2.30 |  1.76 |  1.13 |  0.18 |  0.03 | -0.77 | -0.17 | -0.28 |
| -0.68 |  0.37 |  0.61 | -0.16 |  0.32 | -0.51 | -0.26 | -0.42 |  0.77 |  0.77 |  0.02 |  0.14 |
|  0.65 | -0.31 | -0.45 |  0.99 | -1.48 | -0.68 | -0.10 |  0.99 |  0.24 | -0.39 |  1.06 |  0.45 |

See more in the [tutorial for this
package](https://wjschne.github.io/simstandard//articles/simstandard_tutorial.html).
