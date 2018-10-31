
<!-- README.md is generated from README.Rmd. Please edit that file -->

# simstandard

[![CRAN
status](https://www.r-pkg.org/badges/version/simstandard)](https://cran.r-project.org/package=simstandard)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Travis build
status](https://travis-ci.org/wjschne/simstandard.svg?branch=master)](https://travis-ci.org/wjschne/simstandard)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/wjschne/simstandard?branch=master&svg=true)](https://ci.appveyor.com/project/wjschne/simstandard)

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
| \-0.13 |   0.16 |   0.66 | \-0.59 |   1.52 | \-0.09 | \-0.07 | \-0.09 |   0.23 |   0.70 | \-0.54 |   0.00 |
| \-0.81 | \-0.47 | \-0.29 |   1.17 | \-0.79 | \-0.87 |   0.22 | \-0.37 |   0.23 | \-0.42 |   1.02 |   0.91 |
|   0.54 | \-1.87 | \-2.08 | \-1.37 |   0.82 | \-1.28 | \-1.20 |   1.18 | \-0.85 | \-1.36 | \-0.53 | \-0.18 |
|   0.50 | \-0.53 |   0.24 |   0.73 | \-0.07 |   0.27 |   0.48 |   0.36 | \-0.75 | \-0.05 |   0.39 |   0.26 |
| \-0.09 | \-0.47 |   0.79 |   0.13 |   0.65 |   0.42 | \-0.09 | \-0.29 | \-0.80 |   0.84 |   0.19 | \-0.42 |
| \-0.54 | \-1.06 | \-0.49 | \-0.20 | \-0.26 | \-1.78 | \-1.18 |   0.34 |   0.36 |   0.22 |   0.63 |   0.24 |

See more in the [tutorial for this
package](https://wjschne.github.io/simstandard/articles/simstandard_tutorial.html).
