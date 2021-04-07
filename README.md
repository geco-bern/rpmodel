## Purpose

`rpmodel` provides an implementation of the P-model (Prentice et al., 2014; Wang et al., 2017; Stocker et al., 2020) for predicting acclimated photosynthetic parameters, assimilation, and dark respiration rates as a function of the environment. The main function is `rpmodel()` which returns a list of variables that are mutually consistent within the theory of the P-model (see [Usage](./articles/usage.html) ). Further functions used within `rpmodel()` are also provided through the package.

**Important note:**

The P-model predicts how photosynthesis acclimates to a changing environment, coordinating stomatal conductance, Vcmax and Jmax. This yields a model that has the form of a light use efficiency model, where gross primary production scales linearly with absorbed light, as described in [Stocker et al. 2020](https://doi.org/10.5194/gmd-13-1545-2020). It is important to note that this implies that the P-model is valid only for simulating responses to the environment that evolve over the time scale at which the photosynthetic machinery (e.g., Rubisco) can be assumed to acclimate. Sensible choices are on the order of a couple of weeks to a month. In other words, the arguments (climatic forcing), provided to `rpmodel()` should represent typical daytime mean values, averaged across a couple of weeks. The output is then representative also for average values across the same time scale.


## Usage

This loads the `rpmodel` package and executes the `rpmodel()` function without $J_{\text{max}}$ limitation (argument `method_jmaxlim = "none"`), and with a temperature-independent quantum yield efficiency (argument `do_ftemp_kphio = FALSE`):
```r
library(rpmodel)
out_pmodel <- rpmodel( 
  tc             = 20,           # temperature, deg C
  vpd            = 1000,         # Pa,
  co2            = 400,          # ppm,
  fapar          = 1,            # fraction  ,
  ppfd           = 300,          # mol/m2/d,
  elv            = 0,            # m.a.s.l.,
  kphio          = 0.049977,     # quantum yield efficiency as calibrated for setup ORG by Stocker et al. 2020 GMD,
  beta           = 146,          # unit cost ratio a/b,
  c4             = FALSE,
  method_optci   = "prentice14",
  method_jmaxlim = "none",
  do_ftemp_kphio = FALSE,        # corresponding to setup ORG
  do_soilmstress = FALSE,        # corresponding to setup ORG
  verbose        = TRUE
  )
```

For more information and examples see [Usage](articles/usage.html).

## Installation

### Stable release
`rpmodel` is available on CRAN [here](https://CRAN.R-project.org/package=rpmodel). To install and load, run the following commands in your R terminal:
```r
install.packages("rpmodel")
library(rpmodel)
```

### Development release
To install and load the latest version of the rpmodel package (development release, not yet on CRAN) run the following command in your R terminal: 
```r
if(!require(devtools)){install.packages(devtools)}
devtools::install_github( "stineb/rpmodel", build_vignettes = TRUE )
library(rpmodel)
```

## Author and contact

Benjamin Stocker
benjamin.stocker@gmail.com

## References

Stocker, B. D., Wang, H., Smith, N. G., Harrison, S. P., Keenan, T. F., Sandoval, D., Davis, T., and Prentice, I. C.: P-model v1.0: an optimality-based light use efficiency model for simulating ecosystem gross primary production, Geosci. Model Dev., 13, 1545–1581, https://doi.org/10.5194/gmd-13-1545-2020, 2020.

Wang, H., Prentice, I. C., Keenan, T. F., Davis, T. W., Wright, I. J., Cornwell, W. K.,Evans, B. J., and Peng, C.:  Towards a universal model for carbon dioxide uptake by plants, Nat Plants, 3, 734–741, 2017.

Prentice,  I. C.,  Dong,  N.,  Gleason,  S. M.,  Maire,  V.,  and Wright,  I. J.:  Balancingthe costs of carbon gain and water transport:  testing a new theoretical frameworkfor  plant  functional  ecology, Ecology  Letters,  17,  82–91, 10.1111/ele.12211, 2014.

## Acknowledgement

This project was funded by Marie Sklodowska-Curie fellowship H2020-MSCA-IF-2015, project FIBER, grant number 701329.
