% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{get_model_implied_correlations}
\alias{get_model_implied_correlations}
\title{Return model-implied correlation matrix}
\usage{
get_model_implied_correlations(
  m,
  observed = TRUE,
  latent = FALSE,
  errors = FALSE,
  factor_scores = FALSE,
  composites = FALSE,
  ...
)
}
\arguments{
\item{m}{Structural model represented by lavaan syntax or output
of sim_standardized_matrices function.}

\item{observed}{Include observed variables}

\item{latent}{Include latent variables}

\item{errors}{Include observed error and latent disturbances variables}

\item{factor_scores}{Include factor score variables}

\item{composites}{Include composite variables}

\item{...}{parameters passed to the `sim_standardized_matrices` function}
}
\value{
A correlation matrix
}
\description{
Function that takes a lavaan model with standardized parameters and
returns a model-implied correlation matrix
}
\examples{
library(simstandard)
# lavaan model
m = "Latent_1 =~ 0.8 * Ob_1 + 0.7 * Ob_2 + 0.4 * Ob_3"

get_model_implied_correlations(m)
}
