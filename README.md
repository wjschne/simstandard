
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
|  1.70 | -0.25 |  1.11 |  1.30 | -0.72 |  1.28 |  1.93 |  1.06 | -1.27 | -0.04 | -0.05 |  0.91 |
| -0.18 | -0.92 |  0.63 | -0.67 | -0.56 | -0.72 | -0.17 |  0.18 | -0.35 |  0.74 | -0.55 |  0.40 |
|  0.80 |  0.04 | -0.52 |  0.49 | -0.24 |  0.38 | -0.57 |  0.61 | -0.26 | -0.18 |  0.88 | -0.87 |
| -0.66 | -0.56 | -0.99 | -0.08 | -0.19 | -0.73 | -0.62 | -0.30 |  0.02 | -0.62 |  0.35 | -0.03 |
|  1.14 |  0.50 |  0.03 | -0.39 |  1.43 |  1.04 |  0.56 |  0.62 | -0.33 | -0.31 | -0.79 | -0.27 |
| -0.30 | -0.27 |  1.12 |  0.97 |  0.34 | -0.08 |  1.03 | -0.27 | -0.21 |  0.50 |  0.25 |  1.09 |

See more in the [tutorial for this
package](https://wjschne.github.io/simstandard/articles/simstandard_tutorial.html).
