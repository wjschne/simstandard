#' Return model characteristics
#'
#' Function that takes a lavaan model with standardized parameters and returns a list with model characteristics
#'
#' This function supports the `~` operator for regressions, the `~~` for covariances (but not variances), and the `=~` latent variable loadings. It does not support intercepts (e.g,. `y ~ 1`), thresholds, scaling factors, formative factors, or equality constraints.
#' @export
#' @param m Structural model represented by lavaan syntax
#' @param max_iterations Maximum number of iterations before the algorithm fails
#' @param composite_threshold Loadings with absolute values less than this threshold will not be counted as composite indicators
#' @return list of path and covariance coefficients
#' @examples
#' library(simstandard)
#' # lavaan model
#' m = "Latent_1 =~ 0.8 * Ob_1 + 0.7 * Ob_2 + 0.4 * Ob_3"
#'
#' sim_standardized_matrices(m)
sim_standardized_matrices <- function(m,
                                      max_iterations = 100,
                                      composite_threshold = NULL) {

  # Parameter Table
  pt <- lavaan::lavParTable(m, fixed.x = F)

  # Checks----

  # Check for formative variables
  if (any(pt$op == "<~")) stop(
    "Formative variables (defined with <~) are not allowed for this function.")

  # Check for user-set variances
  if (any((pt$user != 0) & (pt$lhs == pt$rhs) & (pt$op == "~~"))) {
    pt_manual_var <- pt[(pt$user != 0) & (pt$lhs == pt$rhs) & (pt$op == "~~"), ]
    pt_manual_rows <- paste0(pt_manual_var$lhs,
      " ", pt_manual_var$op,
      " ", ifelse(is.na(pt_manual_var$ustart),
        "",
        paste0(pt_manual_var$ustart, " * ")
      ),
      pt_manual_var$rhs,
      collapse = "\n"
    )
    stop(paste0(
      "All variances are set automatically to create standardized data.",
      " You may not set variances manually. ",
      ifelse(nrow(pt_manual_var) > 1,
        "Remove the following parameters:\n",
        "Remove the following parameter:\n"
      ),
      pt_manual_rows
    ))
  }


  # Check for unset paths and covariances
  # if (any(((pt$op == "~") | (pt$op == "=~")) & (pt$free == 1),na.rm = T)) {
  if (any(pt$free == 1, na.rm = T)) {
    pt_unset <- pt[pt$free == 1, ]
    pt_unset_rows <- paste0(pt_unset$lhs,
      " ", pt_unset$op,
      " ", pt_unset$rhs,
      collapse = "\n"
    )
    warning(paste0(
      ifelse(nrow(pt_unset) > 1,
        "Because the following relations were not set, they are assumed to be 0:\n",
        "Because the following relationship was not set, it is assumed to be 0:\n"
      ),
      pt_unset_rows
    ))
  }

  # Check for paths greater than 1
  if (any(abs(pt$ustart) > 1, na.rm = T)) {
    pt_greater <- pt[abs(pt$ustart) > 1, ]
    pt_greater_rows <- paste0(pt_greater$lhs,
      " ", pt_greater$op,
      " ", pt_greater$ustart,
      " * ", pt_greater$rhs,
      collapse = "\n"
    )
    warningmessage <- paste0(
      "Although it is sometimes possible to set standardized parameters greater than 1 or less than -1, it is rare to do so. More often than not, it causes model convergence problems. Check to make sure you set such a value on purpose. ",
      ifelse(nrow(pt_greater) > 1,
        "The following paths were set to values outside the range of -1 to 1:\n",
        "The following path was set to a value outside the range of -1 to 1:\n"
      ),
      pt_greater_rows
    )

    warning(warningmessage)
  }


  # Variable Names----
  v_all <- unique(c(pt$lhs, pt$rhs))
  v_latent <- unique(pt$lhs[pt$op == "=~"])
  v_observed <- v_all[!(v_all %in% v_latent)]
  v_indicator <- unique(pt$rhs[pt$op == "=~"])
  v_y <- unique(pt$lhs[pt$op == "~"])
  v_latent_endogenous <- v_latent[
    (v_latent %in% v_y) | (v_latent %in% v_indicator)]
  v_latent_exogenous <- v_latent[!(v_latent %in% v_latent_endogenous)]

  v_observed_endogenous <- v_observed[
    v_observed %in% v_y | v_observed %in% v_indicator]
  v_observed_exogenous <- v_observed[!v_observed %in% v_observed_endogenous]
  v_observed_indicator <- v_observed[v_observed %in% v_indicator]
  v_latent_indicator <- v_latent[v_latent %in% v_indicator]
  v_observed_y <- v_observed_endogenous[!(v_observed_endogenous %in% v_observed_indicator)]
  v_order <- c(v_observed, v_latent)

  if (length(v_observed_y) == 0) {
    v_error_y <- character(0)
  } else {
    v_error_y <- paste0("e_", v_observed_y)
  }

  if (length(v_latent_endogenous) > 0) {
    v_disturbance <- paste0("d_", v_latent_endogenous)
  } else {
    v_disturbance <- character(0)
  }

  if (length(v_observed_endogenous) > 0) {
    v_error <- paste0("e_", v_observed_endogenous)
  } else {
    v_error <- character(0)
  }

  v_exogenous <- c(v_latent_exogenous, v_observed_exogenous)
  v_endogenous <- c(v_latent_endogenous, v_observed_endogenous)
  v_residual <- c(v_disturbance, v_error)
  v_ellipse <- c(v_latent, v_residual)
  v_source <- c(v_exogenous, v_residual)
  v_factor_score <- v_ellipse[!(v_ellipse %in% v_error_y)]

  # Set unspecified parameters to 0
  pt[is.na(pt[, "ustart"]), "ustart"] <- 0

  # Make RAM matrices----

  # Names for A, S and new S matrices
  vS <- vA <- c(v_endogenous, v_exogenous)

  # Number of Variables
  k <- length(vA)

  # Initialize A matrix and exogenous correlation matrix
  exo_cor <- A <- matrix(0, k, k, dimnames = list(vA, vA))



  # Assign loadings to A
  for (i in pt[pt[, "op"] == "=~", "id"]) {
    A[pt$rhs[i], pt$lhs[i]] <- pt$ustart[i]
  }

  # Assign regressions to A
  for (i in pt[pt[, "op"] == "~", "id"]) {
    A[pt$lhs[i], pt$rhs[i]] <- pt$ustart[i]
  }

  # Assign correlations to exo_cor
  diag(exo_cor) <- 1
  for (i in pt[(pt[, "op"] == "~~") & pt$lhs != pt$rhs, "id"]) {
    exo_cor[pt$lhs[i], pt$rhs[i]] <- pt$ustart[i]
    exo_cor[pt$rhs[i], pt$lhs[i]] <- pt$ustart[i]
  }

  # Solving for error variances and correlation matrix----

  # Column of k ones
  v1 <- matrix(1, k)

  # Initial estimate of error variances
  varS <- as.vector(v1 - (A * A) %*% v1)
  S <- diag(varS) %*% exo_cor %*% diag(varS)

  # Initial estimate of the correlation matrix
  iA <- solve(diag(k) - A)
  R <- iA %*% S %*% t(iA)

  # Set interaction count at 0
  iterations <- 0

  # Find values for S matrix
  while ((round(sum(diag(R)), 10) != k) * (iterations < max_iterations)) {
    iA <- solve(diag(k) - A)
    R <- iA %*% S %*% t(iA)
    sdS <- diag(diag(S) ^ 0.5)
    S <- diag(diag(diag(k) - R)) + (sdS %*% exo_cor %*% sdS)
    diag(S)[diag(S) < 0] <- 0.00000001
    iterations <- iterations + 1
  }
  if (iterations == max_iterations) {
    stop(paste0(
      "Model did not converge after ",
      max_iterations,
      " iterations because at least one variable had a negative variance: ",
      paste0(colnames(A)[(diag(S) == 0.00000001)], collapse = ", ")
    ))
  }

  dimnames(S) <- dimnames(A)

  # Filter Matrix
  filter_matrix <- diag((vA %in% v_observed) * 1)
  dimnames(filter_matrix) <- dimnames(A)

  # Big Matrices----

  A_residual_diag <- sqrt(diag(S[v_endogenous, v_endogenous, drop = F]))
  if (length(A_residual_diag) > 1) {
    A_residual <- diag(A_residual_diag)
  } else {
    if (length(A_residual_diag) == 1) {
      A_residual <- matrix(A_residual_diag, nrow = 1, ncol = 1)
    } else {
      A_residual <- matrix(nrow = 0, ncol = 0)
    }
  }

  A_big <- rbind(
    cbind(A[v_endogenous, , drop = F], A_residual),
    matrix(0, nrow = nrow(A), ncol = nrow(A) + length(v_residual))
  )
  v_big <- c(vA, v_residual)
  dimnames(A_big) <- list(v_big, v_big)

  # Initialize S_big
  S_big <- matrix(0,
    nrow = length(v_big),
    ncol = length(v_big),
    dimnames = dimnames(A_big)
  )

  # Insert off-diagonal values of S into S_big
  S_big[
    v_source,
    v_source
  ] <- exo_cor[c(v_exogenous, v_endogenous),
    c(v_exogenous, v_endogenous),
    drop = F
  ]
  # Insert diagonal values of S_big
  diag(S_big) <- c(
    rep(0, length(v_endogenous)),
    rep(1, length(v_exogenous) + length(v_residual))
  )


  # Compute big correlation matrix of all variables
  iA_big <- solve(diag(nrow(A_big)) - A_big)
  R_big <- iA_big %*% S_big %*% t(iA_big)

  R_xx <- R_big[v_observed_indicator, v_observed_indicator, drop = F]
  R_xy <- R_big[v_observed_indicator, v_factor_score, drop = F]


  # Factor and composite scores ----

  if (length(v_observed_indicator) > 0) {
    i_Rxx <- solve(R_xx)

    A_factor_score <- i_Rxx %*% R_xy

    if (length(v_factor_score) > 0) {
      v_FS <- paste0(v_factor_score, "_FS")
    } else {
      v_FS <- character(0)
    }

    colnames(A_factor_score) <- v_factor_score

    if (length(v_latent) > 0) {
      v_composite_score <- paste0(v_latent, "_Composite")
    } else {
      v_composite_score <- character(0)
    }


    A_composite_direct <- sign(A_big[v_observed_indicator, v_latent, drop = F])

    if (!is.null(composite_threshold)) {
      A_composite_direct <- A_composite_direct *
        (abs(A_big[v_observed_indicator, v_latent, drop = F]) >
           composite_threshold)
    }

    # Has Direct Indicators
    Has_direct <- (colSums(abs(A_composite_direct)) > 0) * 1

    # Second-order factors
    A_composite_second_order <- sign(
      A_composite_direct %*%
        A[v_latent, v_latent, drop = F] %*%
        diag(1 - Has_direct,
             nrow = length(Has_direct)))

    # Third-order factors
    Has_direct_second <- (colSums(abs(A_composite_second_order)) > 0) * 1
    A_composite_third_order <- sign(
      A_composite_second_order %*%
        A[v_latent, v_latent, drop = F] %*%
        diag(1 - Has_direct_second,
             nrow = length(Has_direct_second)))

    # Fourth-order factors
    Has_direct_third <- (colSums(abs(A_composite_third_order)) > 0) * 1
    A_composite_fourth_order <- sign(
      A_composite_third_order %*%
        A[v_latent, v_latent, drop = F] %*%
        diag(1 - Has_direct_third,
             nrow = length(Has_direct_third)))




    A_composite <- A_composite_direct + A_composite_second_order + A_composite_third_order + A_composite_fourth_order



    CM_composite <- t(A_composite) %*%
      R[v_observed_indicator,
        v_observed_indicator,
        drop = FALSE] %*%
      A_composite

    A_composite_w <- A_composite %*%
      diag(diag(CM_composite) ^ -0.5, nrow = nrow(CM_composite))

    colnames(A_composite_w) <- v_composite_score

    colnames(A_factor_score) <- v_FS

    # Factor weights
    fw <- cbind(A_factor_score, A_composite_w)
    W <- matrix(0,
      nrow = nrow(A_big),
      ncol = nrow(A_big) + ncol(fw),
      dimnames = list(
        rownames(A_big),
        c(rownames(A_big), colnames(fw))
      )
    )

    diag(W) <- 1

    W[rownames(fw), colnames(fw)] <- fw

    # Grand correlation matrix of all variables
    R_all <- stats::cov2cor(t(W) %*% R_big %*% W)


    v_order_all <- c(v_order, v_residual, v_FS, v_composite_score)

    R_all <- R_all[v_order_all, v_order_all]
  } else {
    R_all <- R_big
    A_factor_score <- NULL
    A_composite_w <- NULL
    v_FS <- character(0)
    v_composite_score <- character(0)
  }

  # Make complete lavaan model syntax
  lavaan_variances <- paste0(
    v_order,
    " ~~ ",
    diag(S[v_order, v_order]),
    " * ",
    v_order,
    collapse = "\n")

  # Factor Score Validity
  factor_score_validity <- diag(R_all[v_FS,
                                      v_factor_score,
                                      drop = FALSE])
  names(factor_score_validity) <- v_FS
  factor_score_se <- sqrt(1 - factor_score_validity ^ 2)

  # Composite Score Validity
  composite_score_validity <- diag(R_all[v_latent,
                                         v_composite_score,
                                         drop = FALSE])
  names(composite_score_validity) <- v_composite_score


  # Return list ----

  l_names <- list(
    v_observed = v_observed,
    v_latent = v_latent,
    v_latent_exogenous = v_latent_exogenous,
    v_latent_endogenous = v_latent_endogenous,
    v_observed_exogenous = v_observed_exogenous,
    v_observed_endogenous = v_observed_endogenous,
    v_observed_indicator = v_observed_indicator,
    v_disturbance = v_disturbance,
    v_error = v_error,
    v_residual = v_residual,
    v_factor_score = v_FS,
    v_composite_score = v_composite_score
  )

  list(
    RAM_matrices = list(
      A = A[v_order, v_order],
      S = S[v_order, v_order],
      filter_matrix = filter_matrix[v_order, v_order],
      iA = iA[v_order, v_order]
    ),
    Correlations = list(
      R = R[v_order, v_order],
      R_all = R_all
    ),
    Coefficients = list(
      factor_score = A_factor_score,
      factor_score_validity = factor_score_validity,
      factor_score_se = factor_score_se,
      composite_score = A_composite_w,
      composite_score_validity = composite_score_validity
    ),
    lavaan_models = list(
      model_without_variances = m,
      model_with_variances = paste0(m, "\n# Variances\n", lavaan_variances),
      model_free = fixed2free(m)
    ),
    v_names = l_names,
    iterations = iterations
  )
}


#' Generates simulated data with standardized parameters.
#'
#' This function takes a lavaan model with standardized parameters and simulates latent scores, errors, disturbances, and observed scores.
#'
#' This function supports the `~` operator for regressions, the `~~` for covariances (but not variances), and the `=~` latent variable loadings. It does not support intercepts (e.g,. `y ~ 1`), thresholds, scaling factors, formative factors, or equality constraints.
#'
#' @export
#' @param m Structural model represented by lavaan syntax
#' @param n Number of simulated cases
#' @param observed Include observed variables
#' @param latent Include latent variables
#' @param errors Include observed error and latent disturbances variables
#' @param factor_scores Include factor score variables
#' @param composites Include composite variables
#' @param matrices Include matrices as attribute of tibble
#' @param ... Arguments passed to `simstandardized_matrices`
#' @return tibble with standardized data
#' @examples
#' library(simstandard)
#' # Lavaan model
#' m = "Latent_1 =~ 0.8 * Ob_1 + 0.7 * Ob_2 + 0.4 * Ob_3"
#'
#' # simulate 10 cases
#' sim_standardized(m, n = 10)
sim_standardized <- function(
  m,
  n = 1000,
  observed = TRUE,
  latent = TRUE,
  errors = TRUE,
  factor_scores = FALSE,
  composites = FALSE,
  matrices = FALSE,
  ...) {

  # Get main object from sim_standardized_matrices
  o <- sim_standardized_matrices(m,...)

  # Names of variables in S (Symmetric) matrix

  S_names <- c(
    o$v_names$v_observed_exogenous,
    o$v_names$v_observed_endogenous,
    o$v_names$v_latent_exogenous,
    o$v_names$v_latent_endogenous
  )

  # Simulate exogenous variables in S matricx
  u <- mvtnorm::rmvnorm(n = n, sigma = o$RAM_matrices$S[S_names, S_names, drop = F])
  colnames(u) <- c(
    o$v_names$v_observed_exogenous,
    o$v_names$v_error,
    o$v_names$v_latent_exogenous,
    o$v_names$v_disturbance
  )

  # Create all variables from exogenous variables
  v <- u %*% t(o$RAM_matrices$iA[S_names, S_names, drop = F])

  # Make blank matrix with n rows
  d_blank <- matrix(nrow = n, ncol = 0)

  # Extract observed indicators of latent varibles
  d_observed_indicators <- v[, o$v_names$v_observed_indicator, drop = F]

  # Calculate estimated factor scores
  if (length(o$v_names$v_observed_indicator) > 0) {
    d_factor_scores <- d_observed_indicators %*% o$Coefficients$factor_score
  } else {
    d_factor_scores <- d_blank
  }

  # Calculate composite scores
  if (length(o$v_names$v_observed_indicator) > 0) {
    d_composite_scores <- d_observed_indicators %*%
      o$Coefficients$composite_score
  } else {
    d_composite_scores <- d_blank
  }

  # Make data to be returned
  d <- tibble::as_tibble(
    cbind(
      v[, c(o$v_names$v_observed, o$v_names$v_latent), drop = F],
      u[, c(o$v_names$v_disturbance, o$v_names$v_error), drop = F],
      d_factor_scores,
      d_composite_scores
    )
  )



  # Decide which variables to return
  v_include <- character(0)
  if (observed) v_include <- c(v_include, o$v_names$v_observed)
  if (latent) v_include <- c(v_include, o$v_names$v_latent)
  if (errors) v_include <- c(v_include, o$v_names$v_error)
  if (errors) v_include <- c(v_include, o$v_names$v_disturbance)
  if (factor_scores) v_include <- c(v_include, o$v_names$v_factor_score)
  if (composites) v_include <- c(v_include, o$v_names$v_composite_score)

  d <- d[, v_include]

  # Attach metadata as attribute
  if (matrices) attr(d, "matrices") <- o

  # Return tibble
  return(d)
}

#' Remove fixed parameters from a lavaan model
#'
#' @export
#' @param m Structural model represented by lavaan syntax
#' @return character string representing lavaan model
#' @importFrom rlang .data
#' @examples
#' library(simstandard)
#' # lavaan model with fixed parameters
#' m = "
#' Latent_1 =~ 0.9 * Ob_11 + 0.8 * Ob_12 + 0.7 * Ob_13
#' Latent_2 =~ 0.9 * Ob_21 + 0.6 * Ob_22 + 0.4 * Ob_23
#' "
#' # Same model, but with fixed parameters removed.
#' m_free <- fixed2free(m)
#' cat(m_free)
fixed2free <- function(m){
  m %>%
    lavaan::lavaanify(fixed.x = FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs) %>%
    dplyr::group_by(.data$lhs, .data$op) %>%
    dplyr::summarise(rhs = paste(.data$rhs, collapse = " + ")) %>%
    dplyr::arrange(dplyr::desc(.data$op)) %>%
    tidyr::unite("l", .data$lhs, .data$op, .data$rhs, sep = " ") %>%
    dplyr::pull(.data$l) %>%
    paste(collapse = "\n")

}


#' Function that takes a lavaan model with standardized paths and loadings and returns a complete lavaan model syntax with standardized variances
#'
#' @export
#' @param m Structural model represented by lavaan syntax
#' @return character string representing lavaan model
#' @examples
#' library(simstandard)
#' # lavaan model
#' m = "
#' Latent_1 =~ 0.9 * Ob_11 + 0.8 * Ob_12 + 0.7 * Ob_13
#' Latent_2 =~ 0.9 * Ob_21 + 0.6 * Ob_22 + 0.4 * Ob_23
#' Latent_2 ~ 0.6 * Latent_1
#' "
#' # Same lavaan syntax, but with standardized variances
#' m_complete <- model_complete(m)
#' cat(m_complete)
model_complete <- function(m){
  sim_standardized_matrices(m)$lavaan_models$model_with_variances
}

#' Add factor scores to observed data
#'
#' @export
#' @param d A data.frame with observed data in standardized form (i.e, z-scores)
#' @param m A character string with lavaan model
#' @param CI Add confidence intervals? Defaults to `FALSE`. If `TRUE`, For each factor score, a lower and upper bound of the confidence interval is created. For example, the lower bound of factor score `X` is `X_LB`, and the upper bound is `X_UB`.
#' @param p confidence interval proportion. Defaults to 0.95
#' @param ... parameters passed to simstandardized_matrices
#' @return data.frame with observed data and estimated factor scores
#' @examples
#' library(simstandard)
#' # lavaan model
#' m = "
#' X =~ 0.9 * X1 + 0.8 * X2 + 0.7 * X3
#' "
#'
#' # Make data.frame for two cases
#' d <- data.frame(
#'   X1 = c(1.2, -1.2),
#'   X2 = c(1.5, -1.8),
#'   X3 = c(1.8, -1.1))
#'
#' # Compute factor scores for two cases
#' add_factor_scores(d, m)
add_factor_scores <- function(d, m, CI = FALSE, p = 0.95, ...) {
  sm <- sim_standardized_matrices(m, ...)

  # Coefficients for estimated factor scores
  v_FS <- paste0(sm$v_names$v_latent, "_FS")
  latent_factor_score <- sm$Coefficients$factor_score[, v_FS, drop = FALSE]

  # Remove _FS from factor score names
  colnames(latent_factor_score) <- stringr::str_remove_all(
    colnames(latent_factor_score), "_FS")

  # Get observed score names
  v_observed <- rownames(sm$Coefficients$factor_score)

  # Get observed data
  d_observed <- as.matrix(d[, v_observed, drop = FALSE])

  # Make factor scores
  d_factor_score <- d_observed %*% latent_factor_score

  # Bind factor scores to observed data
  d_all <- cbind(as.data.frame(d), as.data.frame(d_factor_score))

  if (CI) {
    # Make CI
    z <- -1 * stats::qnorm((1 - p) / 2)
    FS_se <- sm$Coefficients$factor_score_se[v_FS]
    d_lower_bound <- d_factor_score - z * FS_se
    colnames(d_lower_bound) <- paste0(colnames(d_factor_score), "_LB")
    d_upper_bound <- d_factor_score + z * FS_se
    colnames(d_upper_bound) <- paste0(colnames(d_factor_score), "_UB")
    d_all <- cbind(d_all, d_lower_bound, d_upper_bound)
  }
  d_all
}

#' Add composite scores to observed data
#'
#' @export
#' @param d A data.frame with observed data in standardized form (i.e, z-scores)
#' @param m A character string with lavaan model
#' @param ... parameters passed to simstandardized_matrices
#' @return data.frame with observed data and estimated factor scores
#' @examples
#' library(simstandard)
#' # lavaan model
#' m = "
#' X =~ 0.9 * X1 + 0.8 * X2 + 0.7 * X3
#' "
#'
#' # Make data.frame for two cases
#' d <- data.frame(
#'   X1 = c(1.2, -1.2),
#'   X2 = c(1.5, -1.8),
#'   X3 = c(1.8, -1.1))
#'
#' # Compute factor scores for two cases
#' add_composite_scores(d, m)
add_composite_scores <- function(d, m, ...) {
  sm <- sim_standardized_matrices(m, ...)


  # Get composite score names
  v_composite <- sm$v_names$v_composite_score

  # Coefficients for composite scores
  l_composite_score <- sm$Coefficients$composite_score

  # Get observed score names
  v_observed <- rownames(sm$Coefficients$composite_score)

  # Get observed data
  d_observed <- as.matrix(d[, v_observed, drop = FALSE])

  # Make factor scores
  d_composite_score <- d_observed %*% l_composite_score

  # Bind factor scores to observed data
  d_all <- cbind(as.data.frame(d), as.data.frame(d_composite_score))
  d_all
}

#' Create lavaan model syntax from matrix coefficients
#'
#' @export
#' @param measurement_model A matrix or data.frame with measurement model loadings. Column names are latent variables. Row names or the first column of a data.frame are indicator variables.
#' @param structural_model A matrix or data.frame with structural model coefficients (i.e., regressions). Column names are "causal" variables. Row names or the first column of a data.frame are "effect" variables.
#' @param covariances A matrix or data.frame with model covariances. Column names must match the row names. If a data.frame, row variable names can be specified in the first column.
#' @return a character string with lavaan syntax
#' @examples
#' library(simstandard)
#' # Specifying the measurement model:
#' # For a data.frame, the column names are latent variables,
#' # and the indicators can be specified as rownames.
#' m <- data.frame(X = c(0.7,0.8,0,0),
#'                 Y = c(0,0,0.8,0.9))
#' rownames(m) <- c("A", "B", "C", "D")
#' # Indicator variables can also be specified
#' # as the first column variable
#' # with subsequent column names as latent variables
#' m <- data.frame(Indicators = c("A", "B", "C", "D"),
#'                 X = c(0.7,0.8,0,0),
#'                 Y = c(0,0,0.8,0.9))
#' # Alternately, a matrix can be used:
#' m <- matrix(c(0.7,0.8,0,0,
#'               0,0,0.8,0.9),
#'               ncol = 2,
#'               dimnames = list(c("A", "B", "C", "D"),
#'                               c("X", "Y")))
#' # Specifying the structural coefficients:
#' # The regression coefficients of the structural model can be
#' # specified as either a data.frame or a matrix. Column names
#' # are the predictors and row names are the criterion variables.
#' # With a data.frame, criterion variables can alternataly be
#' # specified with as the first column.
#' s <- matrix(0.5, nrow = 1, ncol = 1, dimnames = list("Y", "X"))
#' # The covariance matrix must be symmetric. Can also be specified
#' # as a data. frame.
#' Sigma <- matrix(c(1, 0.3,
#'                   0.3, 1),
#'                 nrow = 2,
#'                 ncol = 2,
#'                 dimnames = list(c("B","C"),
#'                                 c("B","C")) )
#' model <- matrix2lavaan(measurement_model = m,
#'                        structural_model = s,
#'                        covariances = Sigma)
#' cat(model)
matrix2lavaan <- function(
  measurement_model = NULL,
  structural_model = NULL,
  covariances = NULL) {
  lav_m <- character(0)
  lav_s <- character(0)
  lav_c <- character(0)

  # Measurement model ----
  if (!is.null(measurement_model)) {


    measurement_model <- check_matrix2lavaan(
      m = measurement_model,
      mname = "Measurement model")

    testcol_m <- colnames(measurement_model)[1]
    lav_m <- measurement_model %>%
      dplyr::rename(Test = !!testcol_m) %>%
      tidyr::gather(key = "Construct",
                    value = "Loading",
                    -1,
                    factor_key = T) %>%
      dplyr::filter(.data$Loading != 0) %>%
      dplyr::mutate(
        model = paste0(.data$Loading, " * ", .data$Test)) %>%
      dplyr::group_by(.data$Construct) %>%
      dplyr::summarise(
        model = paste0(
          .data$model,
          collapse = " + ")) %>%
      dplyr::summarise(
        model = paste0(
          .data$Construct,
          " =~ " ,
          .data$model,
          collapse = "\n")) %>%
      dplyr::pull(.data$model)

  }

  # Structural model ----
  if (!is.null(structural_model)) {

    structural_model <- check_matrix2lavaan(
      m = structural_model,
      mname = "Structural model")

    testcol_s <- colnames(structural_model)[1]
    lav_s <- structural_model %>%
      dplyr::rename(Criterion = !!testcol_s) %>%
      tidyr::gather(key = "Predictor",
                    value = "Coefficient",
                    -1,
                    factor_key = T) %>%
      dplyr::filter(.data$Coefficient != 0) %>%
      dplyr::mutate(
        model = paste0(.data$Coefficient, " * ", .data$Predictor)) %>%
      dplyr::group_by(.data$Criterion) %>%
      dplyr::summarise(
        model = paste0(
          .data$model,
          collapse = " + ")) %>%
      dplyr::summarise(
        model = paste0(
          .data$Criterion,
          " ~ " ,
          .data$model,
          collapse = "\n")) %>%
      dplyr::pull(.data$model)

  }

  # Covariances ----
  if (!is.null(covariances)) {


    covariances <- check_matrix2lavaan(
      m = covariances,
      mname = "Covariances")

    mcovariances <- as.matrix(covariances[,-1, drop = FALSE])
    rownames(mcovariances) <- colnames(mcovariances)
    if (!isSymmetric(mcovariances)) stop("covariances must be symmetric.")

    for (j in 2:ncol(covariances)) {
      for (i in 1:nrow(covariances)) {
        if (i + 1 < j) covariances[i,j] <- NA
        if (i + 1 == j & covariances[i,j] == 1) covariances[i,j] <- NA
      }
    }

    testcol_covariances <- colnames(covariances)[1]
    lav_c <- covariances %>%
      dplyr::rename(Test = !!testcol_covariances) %>%
      tidyr::gather(key = "Construct",
                    value = "Coefficient",
                    -1,
                    factor_key = T) %>%
      dplyr::filter(!is.na(.data$Coefficient)) %>%
      dplyr::filter(.data$Coefficient != 0) %>%
      dplyr::mutate(
        model = paste0(.data$Coefficient, " * ", .data$Test)) %>%
      dplyr::group_by(.data$Construct) %>%
      dplyr::summarise(
        model = paste0(
          .data$model,
          collapse = " + ")) %>%
      dplyr::summarise(
        model = paste0(
          .data$Construct,
          " ~~ " ,
          .data$model,
          collapse = "\n")) %>%
      dplyr::pull(.data$model)



  }

  paste(lav_m, lav_s, lav_c, sep = "\n")

}

#' Extract standardized RAM matrices from a lavaan object
#'
#' @export
#' @param fit An object of class lavaan
#' @return list of RAM matrices A (asymmetric paths), S (symmetric paths), and F (filter matrix)
lav2ram <- function(fit) {
  pt <- lavaan::standardizedSolution(fit)
  pt$id <- 1:nrow(pt)
  v_all <- unique(c(pt$lhs, pt$rhs))
  v_latent <- unique(pt$lhs[pt$op == "=~"])
  v_observed <- v_all[!(v_all %in% v_latent)]

  k <- length(v_all)
  A <- matrix(0, nrow = k, ncol = k, dimnames = list(v_all, v_all))
  S <- matrix(0, nrow = k, ncol = k, dimnames = list(v_all, v_all))
  F <- A
  diag(F) <- 1
  F <- F[v_observed,]

  # Assign loadings to A
  for (i in pt[pt[, "op"] == "=~", "id"]) {
    A[pt$rhs[i], pt$lhs[i]] <- pt$est.std[i]
  }

  # Assign regressions to A
  for (i in pt[pt[, "op"] == "~", "id"]) {
    A[pt$lhs[i], pt$rhs[i]] <- pt$est.std[i]
  }

  # Assign correlations to exo_cor

  for (i in pt[(pt[, "op"] == "~~"), "id"]) {
    S[pt$lhs[i], pt$rhs[i]] <- pt$est.std[i]
    S[pt$rhs[i], pt$lhs[i]] <- pt$est.std[i]
  }

  list(A = A, S = S, F = F)

}

#' Checks matrices for matrix2lavaan function
#'
#' @param m matrix, data.frame or tibble
#' @param mname Name of m
#' @keywords internal
#' @usage NULL
check_matrix2lavaan <- function(m, mname) {

  if ("matrix" %in% class(m)) {
    if (is.null(rownames(m))) stop(paste(mname, "must have row names."))
    if (is.null(colnames(m))) stop(paste(mname, "must have column names."))
    m <- as.data.frame(m) %>%
      tibble::rownames_to_column(var = "Test")
  }

  if (!("data.frame" %in% class(m))) {
    stop(
      paste(mname, "must be a data.frame, tibble, or matrix."))
  }

  if (!(purrr::map_chr(m[, 1, drop = FALSE],class) %in% c("character", "factor"))) {
    if (any(rownames(m) == as.character(seq(1, length(
      rownames(m)
    ))))) {
      stop(
        paste(mname, "must either have indicator variable names in the first column or as rownames of the data.frame.")
      )
    } else
      m <-
        tibble::rownames_to_column(m, var = "Test")

  }

  # allnumeric_s <- purrr::map_lgl(m[, -1, drop = FALSE],is.numeric) %>%
  #   all
  #
  # if (!allnumeric_s) stop(paste("All columns of", tolower(mname), "must be numeric except for the first column."))

  m

}


#' Return model-implied correlation matrix
#'
#' Function that takes a lavaan model with standardized parameters and returns a model-implied correlation matrix
#' @export
#' @param m Structural model represented by lavaan syntax
#' @param observed Include observed variables
#' @param latent Include latent variables
#' @param errors Include observed error and latent disturbances variables
#' @param factor_scores Include factor score variables
#' @param composites Include composite variables
#' @param ... parameters passed to the `sim_standardized_matrices` function
#' @return correlation matrix
#' @examples
#' library(simstandard)
#' # lavaan model
#' m = "Latent_1 =~ 0.8 * Ob_1 + 0.7 * Ob_2 + 0.4 * Ob_3"
#'
#' get_model_implied_correlations(m)
get_model_implied_correlations <- function(m,
                                           observed = TRUE,
                                           latent = FALSE,
                                           errors = FALSE,
                                           factor_scores = FALSE,
                                           composites = FALSE,
                                           ...) {
  fit <- sim_standardized_matrices(m, ...)

  # Variable names
  v_names <- character(0)

  # Observed Variable Names
  if (observed) v_names <- c(v_names, fit$v_names$v_observed)

  # Latent Variable Names
  if (latent) v_names <- c(v_names, fit$v_names$v_latent)

  # Error Variable Names
  if (errors) v_names <- c(v_names, fit$v_names$v_residual)

  # Factor-Score Variable Names
  if (factor_scores) v_names <- c(v_names, paste0(fit$v_names$v_latent,"_FS"))

  # Composite Variable Names
  if (composites) v_names <- c(v_names, fit$v_names$v_composite_score)

  # Return correlation matrix
  fit$Correlations$R_all[v_names, v_names]
}
