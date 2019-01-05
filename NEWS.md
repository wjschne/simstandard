# simstandard 0.3.0

* Added the `matrix2lavaan` function to provide a convenient method of creating lavaan syntax from matrices.
* Added the `lav2ram` function to extract standardized RAM matrices from a lavaan object.
* The `sim_standardized_matrices` function has a new argument, `composite_threshold`. If this argument is specified, variables with loadings below the threshold are not used as indicators of the composite scores.
* Removed the semPlot package from suggests list

# simstandard 0.2.1

* Fixed the method of finding indicators for composite variables. A composite now is only created from direct indicators unless the latent variable is a higher-order factor with no direct indicators.

# simstandard 0.2.0

* Added the `fixed2free` function, which takes a `lavaan` syntax model with fixed parameters and returns a `lavaan` syntax model in which all parameters are free.
* Added the `model_complete` function, which takes a `lavaan` syntax model with standardized loadings, structure coefficients, and covariances, and returns a `lavaan` syntax model with all standardized coefficients, including standardized variances.
* Added the `add_factor_scores` function, which adds predicted factor scores to a data.frame.

# simstandard 0.1.0

* Initial release
