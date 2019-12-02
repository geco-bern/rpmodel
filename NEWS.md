## rpmodel 1.0.1

* First version for submission to Geoscientific Model Development. Package available through github (github.com/rpmodel)

## rpmodel 1.0.2

* Preparations for submission to CRAN
    - Vectorized all if-statements using ifelse().
    - Several small fixes in documentation to make it build on several platforms, while resolving devtools::release().
    
## rpmodel 1.0.3

* Version at final publication of Stocker et al. (2019) GMD.
    - Updated parameters `kphio`, `soilm_par_a`, and `soilm_par_b` after updated calibration (slightly changed calibration data set after applying updated data filtering criteria, see Stocker et al. (2019) GMD).

## rpmodel 1.0.4

* Bugfix: vectorized all outputs from `rpmodel()`.
