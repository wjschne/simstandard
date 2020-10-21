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
  "composite threshold",
  code = {
    m <- "A =~ 0.8 * A1 + 0.8 * A2 + 0.2 * A3"
    expect_silent(sim_standardized_matrices(m, composite_threshold = 0.3))

  }
)





test_that(
  "no formative variables", {
    m <- "A <~ 0.8 * A1 + 0.8 * A2 + 0.8 * A3"
    expect_error(sim_standardized_matrices(m),
                 "Formative variables \\(defined with <~\\) are not allowed for this function.")
    }
  )


test_that(
  "user-set variances",
  code = {
    m <- "A =~ 0.8 * A1 + 0.8 * A2 + 0.8 * A3
         A ~~ 2 * A"
    expect_error(sim_standardized_matrices(m), "All variances are set automatically to create standardized data. You may not set variances manually. Remove the following parameter:\nA ~~ 2 \\* A")
  }
)

test_that(
  "unset paths",
  code = {
    m <- "A =~ A1 + 0.8 * A2 + 0.8 * A3"
    expect_warning(sim_standardized_matrices(m), "Because the following relationship was not set, it is assumed to be 0:
A =~ A1")
  }
)


test_that(
  "v_error_y",
  code = {
    m <- "
    A =~ 0.5 * A1 + 0.8 * A2 + 0.8 * A3
    B ~ 0.5 * A
    "
    expect_silent(sim_standardized_matrices(m))
  }
)

test_that(
  "v_latent_endogenous",
  code = {
    m <- "
    A =~ 0.5 * A1 + 0.8 * A2 + 0.8 * A3
    B =~ 0.5 * B1 + 0.8 * B2 + 0.8 * B3
    B ~ 0.5 * A
    "
    expect_silent(sim_standardized_matrices(m))
  }
)


test_that(
  "simstandardized",
  code = {
    m <- "
    A ~ 0.5 * B
    "
    expect_silent(sim_standardized(m))
  }
)


test_that(
  "model_complete",
  code = {
    m <- "
    A =~ 0.5 * A1 + 0.8 * A2 + 0.8 * A3
    B =~ 0.5 * B1 + 0.8 * B2 + 0.8 * B3
    B ~ 0.5 * A
    "
    expect_silent(model_complete(m))
  }
)


test_that(
  "add_factor_scores and add_composite_scores",
  code = {
    m <- "
    A =~ 0.5 * A1 + 0.8 * A2 + 0.8 * A3
    "
    d <- data.frame(A1 = 1, A2 = 2, A3 = 0)
    expect_silent(add_factor_scores(d, m, CI = T))
    expect_silent(add_composite_scores(d, m))
  }
)


test_that(
  "matrix2lavaan",
  code = {
    m <- data.frame(X = c(0.7,0.8,0,0),
                    Y = c(0,0,0.8,0.9))
    rownames(m) <- c("A", "B", "C", "D")
    s <- matrix(0.5, nrow = 1, ncol = 1, dimnames = list("Y", "X"))
    Sigma <- matrix(c(1, 0.3,
                      0.3, 1),
                    nrow = 2,
                    ncol = 2,
                    dimnames = list(c("B","C"),
                                    c("B","C")))
    expect_silent(matrix2lavaan(measurement_model = m,
                                structural_model = s,
                                covariances = Sigma))
  }
)


test_that(
  "lav2ram",
  {
    m <- "A =~ 0.7 * A1 + 0.8 * A2 + 0.9 * A3 + 0.3 * B1
    B =~ 0.7 * B1 + 0.8 * B2 + 0.9 * B3
    B ~ 0.6 * A"
    # Make model m free
    m_free <- fixed2free(m)
    # Generate data based on model m
    d <- sim_standardized(
      m,
      n = 100000,
      latent = FALSE,
      errors = FALSE)

    # Evaluate the fit of model m_free on data d
    library(lavaan)
    lav_results <- sem(
      model = m_free,
      data = d)
    expect_silent(lav2ram(lav_results))
  })


test_that(
  "add_factor_scores and add_composite_scores",
  code = {
    m <- "
    A =~ 0.5 * A1 + 0.8 * A2 + 0.8 * A3
    "
    d <- data.frame(A1 = 1, A2 = 2, A3 = 0)

    expect_silent(add_composite_scores(d, m))
  }
)


test_that(
  "check_matrix2lavaan",
  code = {
    m <- data.frame(X = c(0.7,0.8,0,0),
                    Y = c(0,0,0.8,0.9))
    # rownames(m) <- c("A", "B", "C", "D")
    s <- matrix(0.5, nrow = 1, ncol = 1, dimnames = list("Y", "X"))
    Sigma <- matrix(c(1, 0.3,
                      0.3, 1),
                    nrow = 2,
                    ncol = 2,
                    dimnames = list(c("B","C"),
                                    c("B","C")))
    expect_error(matrix2lavaan(measurement_model = m,
                                structural_model = s,
                                covariances = Sigma))

    rownames(m) <- c("A", "B", "C", "D")
    m <- as.matrix(m)
    m_no_colnames <- m
    colnames(m_no_colnames) <- NULL
    m_no_rownames <- m
    rownames(m_no_rownames) <- NULL
    expect_error(matrix2lavaan(measurement_model = m_no_colnames,
                               structural_model = s,
                               covariances = Sigma))
    expect_error(matrix2lavaan(measurement_model = m_no_rownames,
                               structural_model = s,
                               covariances = Sigma))

    m_as_list <- as.list(m)
    expect_error(matrix2lavaan(measurement_model = m_as_list,
                               structural_model = s,
                               covariances = Sigma))


  }
)


test_that(
  "add_factor_scores and add_composite_scores",
  code = {
    m <- "
    A =~ 0.5 * A1 + 0.8 * A2 + 0.8 * A3
    "
    expect_silent(sim_standardized(m,
                                   factor_scores = TRUE,
                                   composites = TRUE,
                                   matrices = TRUE))
  }
)

test_that(
  "max_iterations",
  code = {
    m <- "
    A =~ 0.5 * A1 + 0.8 * A2 + 0.8 * A3
    "
    expect_error(sim_standardized_matrices(m, max_iterations = 0))
  }
)
