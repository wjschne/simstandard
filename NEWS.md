# simstandard 0.7.0 *development version*
- Changed Depends field from R (>= 3.4.0) to R (>= 3.5.0) because of dependency on mvtnorm package.
- Fixed bug in how composite score indicators are assigned for higher-order factors. Previous code did not distinguish between indicator paths and regression paths.

# simstandard 0.6.3 *2021-05-07*
- Added `get_factor_score_validity_se` function to return factor score standard errors.
- Added `get_model_names` function to return a list of model variable names.


# simstandard 0.6.2 *2021-01-21*
- Added `get_factor_score_coefficients` function to return factor score coefficients
- Added `get_factor_score_validity` function to return factor score validity coefficients
- Added `v_factor_score_disturbance` and `v_factor_score_residual` to `v_names` list returned by `sim_standardized_matrices`. 
- The `v_factor_score` list now only returns factor score names associated with the latent variables.

# simstandard 0.6.1 *2020-12-22*
-   Can specify a mean and standard deviation in the `add_composite_scores` and `add_factor_scores` functions.
-   The `add_factor_scores` function now appends `_FS` to the factor score names.


# simstandard 0.6.0 *2020-11-25*
-   Added `get_model_implied_correlations` function, which returns the model-implied correlation matrix of observed variables, latent variables, error terms, factor scores, and composite variables.

# simstandard 0.5.0 *2020-10-22*

-   Added `composite_score_validity` to list returned by `sim_standardized_matrices`
-   Added `add_composite_scores` function to add composite scores to new data.

# simstandard 0.4.0 *2019-01-08*

-   Fixed bug that prevents computation of composite scores of third-order and fourth-order latent variables.

# simstandard 0.3.0 *2019-01-07*

-   Added the `matrix2lavaan` function to provide a convenient method of creating lavaan syntax from matrices.
-   Added the `lav2ram` function to extract standardized RAM matrices from a lavaan object.
-   The `sim_standardized_matrices` function has a new argument, `composite_threshold`. If this argument is specified, variables with loadings below the threshold are not used as indicators of the composite scores.
-   Removed the semPlot package from suggests list

# simstandard 0.2.1 _2018-11-09_

-   Fixed the method of finding indicators for composite variables. A composite now is only created from direct indicators unless the latent variable is a higher-order factor with no direct indicators.

# simstandard 0.2.0 _2018-11-08_

-   Added the `fixed2free` function, which takes a `lavaan` syntax model with fixed parameters and returns a `lavaan` syntax model in which all parameters are free.
-   Added the `model_complete` function, which takes a `lavaan` syntax model with standardized loadings, structure coefficients, and covariances, and returns a `lavaan` syntax model with all standardized coefficients, including standardized variances.
-   Added the `add_factor_scores` function, which adds predicted factor scores to a data.frame.

# simstandard 0.1.0 _2018-10-06_

-   Initial release
