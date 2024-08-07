# rpmodel 1.2.3

* Corrected ORCID in DESCRIPTION

# rpmodel 1.2.2

* Changed naming of auxiliary functions (e.g., `calc_patm()` instead of `patm()`)
* Minor commenting changes

# rpmodel 1.2.1

* Some minor bugfixes handling unintended usage.
* chi for C4 photosynthesis

# rpmodel 1.2

* Implemented tests and continuous integration. Refactoring.

# rpmodel 1.1.3

* Corrected vectorization. This solves issue [#6](https://github.com/geco-bern/rpmodel/issues/6)

# rpmodel 1.1.2

* Set negative values of the temperature-dependent `kphio` to zero.

# rpmodel 1.1.1

* Corrected the argument description for `kphio` in the manual of `rpmodel()`.

# rpmodel 1.1

* Corrected vectorization. This fixes wrong treatment of vector arguments and solves [issue #3](https://github.com/geco-bern/rpmodel/issues/3).

# rpmodel 1.0.6

* Set negative VPD values to zero. Resolves issue [issue #2](https://github.com/geco-bern/rpmodel/issues/2).

# rpmodel 1.0.5

* Added `jmax` (maximum rate of RuBP regeneration) to the list of returned variables by the `rpmodel()` function.
* Made `calc_viscosity_h2o()` and `calc_density_h2o()` to public functions, exported as part of the rpmodel package.

# rpmodel 1.0.4

* Bugfix: vectorized all outputs from `rpmodel()`.

# rpmodel 1.0.3

* Version at final publication of Stocker et al. (2019) GMD.
    - Updated parameters `kphio`, `soilm_par_a`, and `soilm_par_b` after updated calibration (slightly changed calibration data set after applying updated data filtering criteria, see Stocker et al. (2019) GMD).

# rpmodel 1.0.2

* Preparations for submission to CRAN
    - Vectorized all if-statements using ifelse().
    - Several small fixes in documentation to make it build on several platforms, while resolving devtools::release().

# rpmodel 1.0.1

* First version for submission to Geoscientific Model Development. Package available through github (github.com/rpmodel)

# rpmodel 1.0

* First release
