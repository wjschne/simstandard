#' Function that takes a lavaan model with standardized parameters returns a list with model caracteristics
#'
#' @export
#' @param m Structural model represented by lavaan Syntax
#' @param max_iterations Maximum number of iterations before the algorithm fails
#' @importFrom lavaan lavParTable
#' @return list of path and covariance coefficients
#' @examples
#' # lavaan model
#' m = "Latent_1 =~ 0.8 * Ob_1 + 0.8 * Ob_2"
#'
#' sim_standardized_matrices(m)
sim_standardized_matrices <- function(m, max_iterations = 100) {

  # Parameter Table
  pt <- lavParTable(m, fixed.x = F)

  # Checks----

  # Check for formative variables
  if (any(pt$op == "<~")) stop("Formative variables (defined with <~) are not allowed for this function.")

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
  v_latent_endogenous <- v_latent[(v_latent %in% v_y) | (v_latent %in% v_indicator)]
  v_latent_exogenous <- v_latent[!(v_latent %in% v_latent_endogenous)]

  v_observed_endogenous <- v_observed[v_observed %in% v_y | v_observed %in% v_indicator]
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
    sdS <- diag(diag(S)^0.5)
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

  # Initialise S_big
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



    A_composite <- sign(A_big[v_observed_indicator, v_latent, drop = F]) + sign(A[v_observed_indicator, v_latent, drop = F] %*% A[v_latent, v_latent, drop = F])



    CM_composite <- t(A_composite) %*% R[v_observed_indicator,
      v_observed_indicator,
      drop = F
    ] %*% A_composite

    A_composite_w <- A_composite %*% diag(diag(CM_composite)^-0.5,
      nrow = nrow(CM_composite)
    )
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
      composite_score = A_composite_w
    ),
    v_names = l_names,
    iterations = iterations
  )
}


#' Function that takes a lavaan model with standardized parameters and simulates latent scores, errors, disturbances, and observed scores
#'
#' @export
#' @param m Structural model represented by lavaan Syntax
#' @param n Number of simulated cases
#' @param observed Include observed variables
#' @param latent Include latent variables
#' @param errors Include observed error variables
#' @param disturbances Include latent disturbances variables
#' @param factor_scores Include factor score variables
#' @param composites Include composite variables
#' @param matrices Include matrices as attribute
#' @return tibble with standardized data
#' @importFrom mvtnorm rmvnorm
#' @examples
#' # Lavaan model
#' m = "Latent_1 =~ 0.8 * Ob_1 + 0.8 * Ob_2"
#'
#' # simulate 10 cases
#' sim_standardized(m, n = 10)
sim_standardized <- function(
  m,
  n = 1000,
  observed = TRUE,
  latent = TRUE,
  errors = TRUE,
  disturbances = TRUE,
  factor_scores = FALSE,
  composites = FALSE,
  matrices = FALSE) {

  # Get main object
  o <- sim_standardized_matrices(m)
  # Simulate exogenous variables
  S_names <- c(
    o$v_names$v_observed_exogenous,
    o$v_names$v_observed_endogenous,
    o$v_names$v_latent_exogenous,
    o$v_names$v_latent_endogenous
  )

  u <- rmvnorm(n = n, sigma = o$RAM_matrices$S[S_names, S_names, drop = F])
  colnames(u) <- c(
    o$v_names$v_observed_exogenous,
    o$v_names$v_error,
    o$v_names$v_latent_exogenous,
    o$v_names$v_disturbance
  )
  v <- u %*% t(o$RAM_matrices$iA[S_names, S_names, drop = F])
  d_blank <- matrix(nrow = n, ncol = 0)
  # d_observed <- v[ , o$v_names$v_observed, drop = F]
  # d_latent <- v[ , o$v_names$v_latent, drop = F]
  # d_disturbance <- u[ , o$v_names$v_disturbance, drop = F]
  # d_errors <- u[ , o$v_names$v_error, drop = F]
  d_observed_indicators <- v[, o$v_names$v_observed_indicator, drop = F]

  if (length(o$v_names$v_observed_indicator) > 0) {
    d_factor_scores <- d_observed_indicators %*% o$Coefficients$factor_score
  } else {
    d_factor_scores <- d_blank
  }

  if (length(o$v_names$v_observed_indicator) > 0) {
    d_composite_scores <- d_observed_indicators %*% o$Coefficients$composite_score
  } else {
    d_composite_scores <- d_blank
  }
  d <- tibble::as_tibble(
    cbind(
      v[, c(o$v_names$v_observed, o$v_names$v_latent), drop = F],
      u[, c(o$v_names$v_disturbance, o$v_names$v_error), drop = F],
      d_factor_scores,
      d_composite_scores
    )
  )

  if (matrices) attr(d, "matrices") <- o

  v_include <- character(0)
  if (observed) v_include <- c(v_include, o$v_names$v_observed)
  if (latent) v_include <- c(v_include, o$v_names$v_latent)
  if (errors) v_include <- c(v_include, o$v_names$v_error)
  if (disturbances) v_include <- c(v_include, o$v_names$v_disturbance)
  if (factor_scores) v_include <- c(v_include, o$v_names$v_factor_score)
  if (composites) v_include <- c(v_include, o$v_names$v_composite_score)


  return(d[,v_include])
}