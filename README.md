
<!-- README.md is generated from README.Rmd. Please edit that file -->

# simstandard <img src="man/figures/logo.png" align="right" height=140/>

[![CRAN
status](https://www.r-pkg.org/badges/version/simstandard)](https://cran.r-project.org/package=simstandard)
[![Rdoc](http://www.rdocumentation.org/badges/version/simstandard)](http://www.rdocumentation.org/packages/simstandard)
[![develVersion](https://img.shields.io/badge/devel%20version-0.2.0-blue.svg?style=flat)](https://github.com/wjschne/simstandard)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Travis build
status](https://travis-ci.org/wjschne/simstandard.svg?branch=master)](https://travis-ci.org/wjschne/simstandard)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/wjschne/simstandard?branch=master&svg=true)](https://ci.appveyor.com/project/wjschne/simstandard)
[![Coverage
status](https://codecov.io/gh/wjschne/simstandard/branch/master/graph/badge.svg)](https://codecov.io/github/wjschne/simstandard?branch=master)

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
syntax](http://lavaan.ugent.be/tutorial/syntax1.html) to specify models.

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
| \-1.22 |   0.58 |   1.08 |   1.43 |   0.48 | \-0.49 |   0.12 | \-0.97 |   0.97 |   1.00 |   1.35 |   0.51 |
| \-0.53 |   0.52 |   0.66 |   0.40 | \-0.76 |   1.77 |   0.41 | \-1.41 | \-0.90 |   0.41 |   0.11 | \-1.00 |
| \-0.08 |   0.35 |   1.21 |   2.25 |   0.27 |   0.44 |   0.65 | \-0.30 |   0.00 |   0.82 |   1.80 |   0.30 |
|   1.23 |   0.06 | \-2.34 | \-1.45 |   0.16 |   0.39 | \-0.63 |   1.04 | \-0.26 | \-1.96 | \-1.00 | \-0.95 |
| \-0.38 | \-0.86 | \-0.66 | \-2.23 | \-0.96 | \-1.69 | \-1.48 |   0.47 |   0.50 |   0.23 | \-1.19 | \-0.12 |
|   0.63 |   0.55 | \-0.65 | \-0.56 | \-0.02 |   0.22 | \-0.12 |   0.51 |   0.38 | \-0.58 | \-0.47 | \-0.30 |

See more in the [tutorial for this
package](https://wjschne.github.io/simstandard/articles/simstandard_tutorial.html).
