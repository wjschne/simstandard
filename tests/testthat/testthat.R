library(testthat)
library(simstandard)

test_that(
  "correlation model",
  code = {
    m <- "x ~~ 0.5 * y"
    expect_equal(sim_standardized_matrices(m)$Correlations$R[1, 2], 0.5)


  }
)

test_that(
  "composites",
  code = {
    m <- "A =~ 0.8 * A1 + 0.8 * A2 + 0.8 * A3"
    o <- sim_standardized_matrices(m)
    m_c <- o$Coefficients$composite_score
    m_c_total <- sum(sign(m_c))
    expect_equal(m_c_total, 3)


  }
)
