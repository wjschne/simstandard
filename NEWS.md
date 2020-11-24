# simstandard 0.6.0
-   Added `get_model_implied_correlations` function, which returns the model-implied correlation matrix of observed variables, latent variables, error terms, factor scores, and composite variables.

# simstandard 0.5.0

-   Added `composite_score_validity` to list returned by `sim_standardized_matrices`
-   Added `add_composite_scores` function to add composite scores to new data.

# simstandard 0.4.0

-   Fixed bug that prevents computation of composite scores of third-order and fourth-order latent variables.

# simstandard 0.3.0

-   Added the `matrix2lavaan` function to provide a convenient method of creating lavaan syntax from matrices.
-   Added the `lav2ram` function to extract standardized RAM matrices from a lavaan object.
-   The `sim_standardized_matrices` function has a new argument, `composite_threshold`. If this argument is specified, variables with loadings below the threshold are not used as indicators of the composite scores.
-   Removed the semPlot package from suggests list

# simstandard 0.2.1

-   Fixed the method of finding indicators for composite variables. A composite now is only created from direct indicators unless the latent variable is a higher-order factor with no direct indicators.

# simstandard 0.2.0

-   Added the `fixed2free` function, which takes a `lavaan` syntax model with fixed parameters and returns a `lavaan` syntax model in which all parameters are free.
-   Added the `model_complete` function, which takes a `lavaan` syntax model with standardized loadings, structure coefficients, and covariances, and returns a `lavaan` syntax model with all standardized coefficients, including standardized variances.
-   Added the `add_factor_scores` function, which adds predicted factor scores to a data.frame.

# simstandard 0.1.0

-   Initial release
