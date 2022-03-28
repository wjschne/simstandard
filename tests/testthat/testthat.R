library(testthat)
library(simstandard)


test_that("path greater than 1",
          code = {
            m <- "A ~~ -0.99 * B
            C ~ 1.01 * A + 0.1 * B"
            expect_warning(
              sim_standardized_matrices(m))

          })
test_that("correlation model",
          code = {
            m <- "x ~~ 0.5 * y"
            expect_equal(sim_standardized_matrices(m)$Correlations$R[1, 2], 0.5)


          })

test_that("composite threshold",
          code = {
            m <- "A =~ 0.8 * A1 + 0.8 * A2 + 0.2 * A3"
            expect_silent(
              sim_standardized_matrices(m, composite_threshold = 0.3))

          })

test_that("no formative variables", {
  m <- "A <~ 0.8 * A1 + 0.8 * A2 + 0.8 * A3"
  expect_error(
    sim_standardized_matrices(m),
    paste0("Formative variables \\(defined with <~\\) ",
           "are not allowed for this function.")
  )
})

test_that("user-set variances",
          code = {
            m <- "A =~ 0.8 * A1 + 0.8 * A2 + 0.8 * A3
         A ~~ 2 * A"
            expect_error(
              sim_standardized_matrices(m),
              paste0("All variances are set automatically to create ",
                     "standardized data. You may not set variances manually. ",
                     "Remove the following parameter:\nA ~~ 2 \\* A"
            ))
          })

test_that(
  "unset paths",
  code = {
    m <- "A =~ A1 + 0.8 * A2 + 0.8 * A3"
    expect_warning(
      sim_standardized_matrices(m),
      "Because the following relationship was not set, it is assumed to be 0:
A =~ A1")
          })


test_that("v_error_y",
          code = {
            m <- "
    A =~ 0.5 * A1 + 0.8 * A2 + 0.8 * A3
    B ~ 0.5 * A
    "
            expect_silent(sim_standardized_matrices(m))
          })

test_that("v_latent_endogenous",
          code = {
            m <- "
    A =~ 0.5 * A1 + 0.8 * A2 + 0.8 * A3
    B =~ 0.5 * B1 + 0.8 * B2 + 0.8 * B3
    B ~ 0.5 * A
    "
            expect_silent(sim_standardized_matrices(m))
          })


test_that("simstandardized",
          code = {
            m <- "
    A ~ 0.5 * B
    "
            expect_silent(sim_standardized(m))
          })


test_that("model_complete",
          code = {
            m <- "
    A =~ 0.5 * A1 + 0.8 * A2 + 0.8 * A3
    B =~ 0.5 * B1 + 0.8 * B2 + 0.8 * B3
    B ~ 0.5 * A
    "
            expect_silent(model_complete(m))
          })


test_that("add_factor_scores and add_composite_scores",
          code = {
            m <- "
    A =~ 0.5 * A1 + 0.8 * A2 + 0.8 * A3
    "
            d <- data.frame(A1 = 1, A2 = 2, A3 = 0)
            expect_silent(add_factor_scores(d, m, CI = T))
            expect_silent(add_composite_scores(d, m))
          })


test_that("matrix2lavaan",
          code = {
            m <- data.frame(X = c(0.7, 0.8, 0, 0),
                            Y = c(0, 0, 0.8, 0.9))
            rownames(m) <- c("A", "B", "C", "D")
            s <- matrix(0.5,
                        nrow = 1,
                        ncol = 1,
                        dimnames = list("Y", "X"))
            sigma <- matrix(
              c(1, 0.3,
                0.3, 1),
              nrow = 2,
              ncol = 2,
              dimnames = list(c("B", "C"),
                              c("B", "C"))
            )
            expect_silent(matrix2lavaan(
              measurement_model = m,
              structural_model = s,
              covariances = sigma
            ))
          })


test_that("lav2ram", code = {
  m <- "
          A =~ 0.7 * A1 + 0.8 * A2 + 0.9 * A3 + 0.3 * B1
          B =~ 0.7 * B1 + 0.8 * B2 + 0.9 * B3
          B ~ 0.6 * A"
  # Make model m free
  m_free <- fixed2free(m)
  # Generate data based on model m
  d <- sim_standardized(m,
                        n = 100000,
                        latent = FALSE,
                        errors = FALSE)

  # Evaluate the fit of model m_free on data d
  library(lavaan)
  lav_results <- sem(model = m_free,
                     data = d)
  expect_silent(lav2ram(lav_results))
})


test_that("add_factor_scores and add_composite_scores",
          code = {
            m <- "
    A =~ 0.5 * A1 + 0.8 * A2 + 0.8 * A3
    "
            d <- data.frame(A1 = 1, A2 = 2, A3 = 0)

            expect_silent(add_composite_scores(d, m))
          })


test_that("check_matrix2lavaan",
          code = {
            m <- data.frame(X = c(0.7, 0.8, 0, 0),
                            Y = c(0, 0, 0.8, 0.9))
            s <- matrix(0.5,
                        nrow = 1,
                        ncol = 1,
                        dimnames = list("Y", "X"))
            sigma <- matrix(
              c(1, 0.3,
                0.3, 1),
              nrow = 2,
              ncol = 2,
              dimnames = list(c("B", "C"),
                              c("B", "C"))
            )
            expect_error(matrix2lavaan(
              measurement_model = m,
              structural_model = s,
              covariances = sigma
            ))

            rownames(m) <- c("A", "B", "C", "D")
            m <- as.matrix(m)
            m_no_colnames <- m
            colnames(m_no_colnames) <- NULL
            m_no_rownames <- m
            rownames(m_no_rownames) <- NULL
            expect_error(
              matrix2lavaan(
                measurement_model = m_no_colnames,
                structural_model = s,
                covariances = sigma
              )
            )
            expect_error(
              matrix2lavaan(
                measurement_model = m_no_rownames,
                structural_model = s,
                covariances = sigma
              )
            )

            m_as_list <- as.list(m)
            expect_error(
              matrix2lavaan(
                measurement_model = m_as_list,
                structural_model = s,
                covariances = sigma
              )
            )


          })


test_that("add_factor_scores and add_composite_scores",
          code = {
            m <- "
    A =~ 0.5 * A1 + 0.8 * A2 + 0.8 * A3
    "
            expect_silent(sim_standardized(
              m,
              factor_scores = TRUE,
              composites = TRUE,
              matrices = TRUE
            ))
          })

test_that("max_iterations",
          code = {
            m <- "
    A =~ 0.5 * A1 + 0.8 * A2 + 0.8 * A3
    "
            expect_error(sim_standardized_matrices(m, max_iterations = 0))
          })


test_that(
  "missing composite data",
  code = {
    m <- "A =~ 0.5 * A1 + 0.8 * A2 + 0.8 * A3"
    d <- data.frame(A1 = 1, A2 = 2)
    expect_error(
      add_composite_scores(d, m),
      paste0("Some observed variables specified in the model",
             " are missing from the data\\.")
      )
    })

test_that("composite with mu and sigma specified",
          code = {
            m <- "
    A =~ 0.8 * A1 + 0.8 * A2 + 0.8 * A3
    "
            d_ss <- data.frame(A1 = 115, A2 = 115, A3 = 115) %>%
              add_composite_scores(m, mu = 100, sigma = 15) %>%
              add_factor_scores(m, mu = 100, sigma = 15)
            d_z <- data.frame(A1 = 1, A2 = 1, A3 = 1) %>%
              add_composite_scores(m) %>%
              add_factor_scores(m)
            expect_equal(d_ss$A_Composite, d_z$A_Composite * 15 + 100)
            expect_equal(d_ss$A_FS, d_z$A_FS * 15 + 100)
          })


test_that("composite and factor scores with names suffixes",
          code = {
            m <- "
    A =~ 0.8 * A1 + 0.8 * A2 + 0.8 * A3
    "
            d <- data.frame(A1 = 115, A2 = 115, A3 = 115)
            v_composite <-
              colnames(add_composite_scores(d,
                                            m,
                                            names_suffix = "_cs"))[4]
            v_fs <-
              colnames(add_factor_scores(d,
                                         m,
                                         names_suffix = "_factor_score"))[4]
            expect_equal(v_composite,  "A_cs")
            expect_equal(v_fs,  "A_factor_score")
          })


test_that("get model-implied matrix from sim_standardized_matrices output.",
          code = {
            m <- "
    A =~ 0.8 * A1 + 0.8 * A2 + 0.8 * A3
    "
            fit <- sim_standardized_matrices(m)
            expect_equal(get_model_implied_correlations(m),
                         get_model_implied_correlations(fit))
          })


test_that("get model-implied matrix from sim_standardized_matrices output.",
          code = {
            m <- "
    A =~ 0.5 * A1 + 0.8 * A2 + 0.8 * A3
    B =~ 0.5 * B1 + 0.8 * B2 + 0.8 * B3
    B ~ 0.5 * A
    "
            fit <- sim_standardized_matrices(m)
            get_factor_score_coefficients(m, latent = F, errors = T)
            expect_equal(get_model_implied_correlations(m),
                         get_model_implied_correlations(fit))
          })


# Model of Reading
m_reading <- "
Ga =~ 0.83 * Ga1 + 0.92 * Ga2 + 0.95 * Ga3
Gc =~ 0.88 * Gc1 + 0.71 * Gc2 + 0.85 * Gc3
RD =~ 0.93 * RD1 + 0.87 * RD2 + 0.85 * RD3
RC =~ 0.91 * RC1 + 0.86 * RC2 + 0.90 * RC3
Ga ~~ 0.68 * Gc
RD ~  0.57 * Ga + 0.33 * Gc
RC ~  0.05 * Ga + 0.40 * Gc  + 0.43 * RD
"
test_that("Number of factor scores",
          code = {
            m_names <- get_model_names(m_reading)
            expect_equal(m_names$v_factor_score,
                         paste0(c("Ga", "Gc", "RD", "RC"), "_FS")
                         )
          })

