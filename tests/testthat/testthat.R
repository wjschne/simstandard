library(testthat)
library(simstandard)

test_that(
  "correlation model",
  code = {
    m <- "x ~~ 0.5 * y"
    expect_equal(sim_standardized_matrices(m)$Correlations$R[1, 2], 0.5)

    sim_standardized_matrices(m)
  }
)
