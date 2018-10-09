
<!-- README.md is generated from README.Rmd. Please edit that file -->

# simstandard

Sometimes you have a structural model with standardized path
coefficients, structural coefficients, and correlations, but you do not
know the error and disturbance variances. The purpose of `simstandard`
is to calculate these variances and then simulate multivariate normal
data based on your model.

## Installation

You can install `simstandard` from github with:

``` r
# install.packages("devtools")
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
|   0.05 | \-1.18 | \-1.00 |   1.02 | \-0.12 |   0.02 | \-0.22 |   0.04 | \-1.19 | \-0.87 |   1.17 | \-0.23 |
| \-0.56 |   0.07 |   0.09 | \-0.11 |   0.38 | \-0.69 | \-0.68 | \-0.22 |   0.63 |   0.50 |   0.36 | \-0.13 |
| \-1.05 |   0.67 | \-1.13 |   0.17 | \-0.02 |   0.08 | \-0.15 | \-1.09 |   0.61 | \-1.04 |   0.27 | \-0.21 |
| \-1.40 | \-0.95 |   0.67 | \-0.32 | \-0.21 |   0.34 | \-0.32 | \-1.57 | \-1.23 |   0.86 | \-0.10 | \-0.60 |
| \-0.38 | \-1.95 | \-2.62 | \-0.67 | \-1.20 | \-1.93 | \-1.04 |   0.58 | \-0.41 | \-2.00 |   0.06 |   0.50 |
|   0.59 | \-0.04 | \-0.86 |   0.97 |   2.10 |   0.66 |   0.91 |   0.26 | \-0.57 | \-1.41 |   0.33 |   0.38 |

See more in the [tutorial for this
package](https://wjschne.github.io/simstandard/articles/simstandard_tutorial.html).
