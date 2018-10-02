
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
C ~ 0.5 * A1
"
data <- sim_standardized(m = model, n = 500)
```

See more in the [tutorial for this
package](https://wjschne.github.io/simstandard/articles/simstandard_tutorial.html).
