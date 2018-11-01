# simstandard 0.2.0

* Added the `fixed2free` function, which takes a `lavaan` syntax model with fixed parameters and returns a `lavaan` syntax model in which all parameters are free.
* Added the `model_complete` function, which takes a `lavaan` syntax model with standardized loadings, structure coefficients, and covariances, and returns a `lavaan` syntax model with all standardized coefficients, including standardized variances.
* Added the `add_factor_scores` function, which adds predicted factor scores to a data.frame.

# simstandard 0.1.0

* Initial release
