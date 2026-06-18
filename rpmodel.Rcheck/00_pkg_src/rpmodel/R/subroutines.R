#' Dampen inputs of rpmodel
#'
#' Applies an exponential dampening input time series with specified time scale.
#'
#' @param vec A numeric vector for the time series of a daily meteorological
#'  variable used as input for \code{rpmodel} (temperature, 
#'  vapour pressure deficit, CO2, or atmospheric pressure). The length of
#'  \code{x} must be at least 365, i.e., corresponding to one year.
#' @param tau The time scale of dampening (e-folding time scale of a 
#'  perturbation). Must be smaller or equal to 365 d.
#' @return A numeric vector of equal length as \code{x} with damped variation.
#'  The dampening is calculated as:
#'  \deqn{
#'   S(t+1) - S(t) = (X(t+1) - S(t)) / \tau
#'  }
#'  Where \eqn{X} is the daily varying time series given by argument \code{x},
#'  \eqn{S} is the dampened time returned by this function, and \eqn{\tau} is 
#'  the decay time scale of a perturbation, given by argument \code{tau}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' dampen_vec(
#'  vec = 20 * (sin(doy*pi/(365)))^2 + rnorm(365, mean = 0, sd = 5),
#'  tau = 40 
#'  )
#' }

dampen_vec <- function( vec, tau ){
  
  if (length(vec)<365){
    
    stop("dampen_vec(): Aborting.
         Length of argument 'vec' must be at least 365.")
    
  } else {
    
    ## Add one year's data to the head of 'vec' to avoid "boundary effects"
    vec <- c(vec[1:365], vec)
    
    ## apply dampening
    vec_damped <- rep(NA, length(vec))
    vec_damped[1] <- vec[1] # this is the unwanted "boundary effect"
    for (idx in 2:length(vec)){
      dvar <- (1.0/tau) * (vec[idx] - vec_damped[idx - 1])
      vec_damped[idx] <- vec_damped[idx - 1] + dvar
    }
    
    ## remove the first year again, that was added artifically before
    vec_damped <- vec_damped[366:length(vec_damped)]
    
  }
  
  return( vec_damped )
}

#' Calculates an empirical soil moisture stress factor
#'
#' Calculates an empirical soil moisture stress factor as a function of relative
#' soil moisture (fraction of field capacity).
#'
#' @param soilm Relative soil moisture as a fraction
#' of field capacity (unitless). Defaults to 1.0 (no soil moisture stress).
#' @param meanalpha Local annual mean ratio of actual over potential 
#'  evapotranspiration, measure for average aridity. Defaults to 1.0.
#' @param apar_soilm (Optional, used only if \code{do_soilmstress==TRUE})
#'  Parameter determining the sensitivity of the empirical soil moisture stress
#'  function. Defaults to 0.0, the empirically fitted value as presented in Stocker
#'  et al. (2019) Geosci. Model Dev. for model setup 'FULL' (corresponding to a setup
#'  with \code{method_jmaxlim="wang17", do_ftemp_kphio=TRUE, do_soilmstress=TRUE}).
#' @param bpar_soilm (Optional, used only if \code{do_soilmstress==TRUE}) Parameter
#'  determining the sensitivity of the empirical soil moisture stress function. Defaults
#'  to ~0.6, the empirically fitted value as presented in Stocker et al. (2019) Geosci. Model Dev.
#'  for model setup 'FULL' (corresponding to a setup with \code{method_jmaxlim="wang17",
#'  do_ftemp_kphio=TRUE, do_soilmstress=TRUE}).
#'
#' @details The soil moisture stress factor is calculated using a quadratic function that
#' is 1 above \code{soilm} = 0.6 and has a sensitivity, given by the y-axis cutoff,
#' (zero soil moisture), determined  by average aridity (argument \code{meanalpha}) as:
#' \deqn{
#'     \beta = q(\theta - \theta*)^2 + 1
#' }
#' for \eqn{\theta < \theta*} and \eqn{\beta = 1.0} otherwise. \eqn{\theta*} is fixed at 0.6.
#' \eqn{q} is the sensitivity parameter and is calculated as a linear function of average aridity,
#' quantified by the local annual mean ratio of actual over potential evapotranspiration, termed
#' \eqn{\alpha}:
#' \deqn{
#'      q=(\beta0-1)/(\theta*-\theta0)^2
#' }
#' \eqn{\theta0} is 0.0, and
#' \deqn{
#'      \beta0 = a + b \alpha
#' }
#' \eqn{a} is given by argument \code{apar}, \eqn{b} is given by argument \code{bpar}.
#'
#' @return A numeric value for \eqn{\beta}
#'
#' @examples
#' ## Relative reduction (%) in GPP due to soil moisture stress at
#' ## relative soil water content ('soilm') of 0.2:
#' print((calc_soilmstress(0.2)-1)*100 )
#'
#' @references  Stocker, B. et al. Geoscientific Model Development Discussions (in prep.)
#'
#' @export

calc_soilmstress <- function(
  soilm,
  meanalpha = 1.0,
  apar_soilm = 0.0,
  bpar_soilm = 0.685
){
  
  # Fixed parameters
  x0 <- 0.0
  x1 <- 0.6
  
  y0 <- (apar_soilm + bpar_soilm * meanalpha)
  
  beta <- (1.0 - y0) / (x0 - x1)^2
  outstress <- 1.0 - beta * ( soilm - x1 )^2
  
  ## bound between 0 and 1
  outstress <- pmin(pmax(outstress, 0), 1)
  # and set to 1.0 above soil moisture threshold x1.
  outstress[soilm > x1] <- 1.0

  return(outstress)
}


#' Calculates atmospheric pressure
#'
#' Calculates atmospheric pressure as a function of elevation, by default assuming 
#' standard atmosphere (101325 Pa at sea level)
#'
#' @param elv Elevation above sea-level (m.a.s.l.)
#' @param patm0 (Optional) Atmospheric pressure at sea level (Pa), defaults to 101325 Pa.
#'
#' @details The elevation-dependence of atmospheric pressure is computed by 
#' assuming a linear decrease in temperature with elevation and a mean 
#' adiabatic lapse rate (Berberan-Santos et al., 1997):
#' \deqn{
#'    p(z) = p0 ( 1 - Lz / TK0) ^ ( g M / (RL) )
#' }
#' where \eqn{z} is the elevation above mean sea level (m, argument \code{elv}), 
#' \eqn{g} is the gravity constant (9.80665 m s-2), \eqn{p0} is the atmospheric 
#' pressure at 0 m a.s.l. (argument \code{patm0}, defaults to 101325 Pa), 
#' \eqn{L} is the mean adiabatic lapse rate (0.0065 K m-2), 
#' \eqn{M} is the molecular weight for dry air (0.028963 kg mol-1), 
#' \eqn{R} is the universal gas constant (8.3145 J mol-1 K-1), and \eqn{TK0}
#' is the standard temperature (298.15 K, corresponds to 25 deg C).
#'
#' @return A numeric value for \eqn{p}
#'
#' @examples print("Standard atmospheric pressure, in Pa, corrected for 1000 m.a.s.l.:")
#' print(calc_patm(1000))
#' 
#' @references  Allen, R. G., Pereira, L. S., Raes, D., Smith, M.: 
#'              FAO Irrigation and Drainage Paper No. 56, Food and 
#'              Agriculture Organization of the United Nations, 1998
#'
#' @export
#' 
calc_patm <- function( elv, patm0 = 101325 ){
  
  # Define constants:
  kTo <- 298.15    # base temperature, K (Prentice, unpublished)
  kL  <- 0.0065    # adiabiatic temperature lapse rate, K/m (Allen, 1973)
  kG  <- 9.80665   # gravitational acceleration, m/s^2 (Allen, 1973)
  kR  <- 8.3145    # universal gas constant, J/mol/K (Allen, 1973)
  kMa <- 0.028963  # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
  
  # Convert elevation to pressure, Pa:
  out <- patm0*(1.0 - kL*elv/kTo)^(kG*kMa/(kR*kL))
  
  return(out)
}

#' Calculates the Michaelis Menten coefficient for Rubisco-limited photosynthesis
#'
#' Calculates the Michaelis Menten coefficient of Rubisco-limited assimilation
#' as a function of temperature and atmospheric pressure.
#'
#' @param tc Temperature, relevant for photosynthesis (deg C)
#' @param patm Atmospheric pressure (Pa)
#'
#' @details The Michaelis-Menten coefficient \eqn{K} of Rubisco-limited
#' photosynthesis is determined by the Michalis-Menten constants for
#' O2 and CO2 (Farquhar, 1980) according to:
#' \deqn{
#'   K = Kc ( 1 + pO2 / Ko)
#' }
#' where \eqn{Kc} is the Michaelis-Menten constant for CO2 (Pa), \eqn{Ko} is
#' the Michaelis-Menten constant for O2 (Pa), and \eqn{pO2} is the partial
#' pressure of oxygen (Pa), calculated as \eqn{0.209476 p}, where \eqn{p} is
#' given by argument \code{patm}.  \eqn{Kc} and \eqn{Ko} follow a temperature
#' dependence, given by the Arrhenius Equation \eqn{f} (implemented by
#' \link{ftemp_arrh}):
#' \deqn{
#'    Kc = Kc25 f(T, \Delta Hkc)
#' }
#' \deqn{
#'    Ko = Ko25 f(T, \Delta Hko)
#' }
#' Values \eqn{\Delta Hkc} (79430 J mol-1), \eqn{\Delta Hko} (36380 J mol-1),
#' \eqn{Kc25} (39.97 Pa), and \eqn{Ko25} (27480 Pa) are taken from Bernacchi
#' et al. (2001) and have been converted from values given therein to units of Pa
#' by multiplication with the standard atmosphere (101325 Pa). \eqn{T} is given
#' by the argument \code{tc}.
#'
#' @references Farquhar,  G.  D.,  von  Caemmerer,  S.,  and  Berry,  J.  A.:
#'             A  biochemical  model  of photosynthetic CO2 assimilation in leaves of
#'             C 3 species, Planta, 149, 78–90, 1980.
#'
#'             Bernacchi,  C.  J.,  Singsaas,  E.  L.,  Pimentel,  C.,  Portis,  A.
#'             R.  J.,  and  Long,  S.  P.:Improved temperature response functions
#'             for models of Rubisco-limited photosyn-thesis, Plant, Cell and
#'             Environment, 24, 253–259, 2001
#'
#' @return A numeric value for \eqn{K} (in Pa)
#'
#' @examples print("Michaelis-Menten coefficient at 20 degrees Celsius and standard atmosphere (in Pa):")
#' print(calc_kmm(20, 101325))
#'
#' @export
#'
calc_kmm <- function( tc, patm ) {
  
  dhac   <- 79430      # (J/mol) Activation energy, Bernacchi et al. (2001)
  dhao   <- 36380      # (J/mol) Activation energy, Bernacchi et al. (2001)
  kco    <- 2.09476e5  # (ppm) O2 partial pressure, Standard Atmosphere
  
  ## k25 parameters are not dependent on atmospheric pressure
  kc25 <- 39.97   # Pa, value based on Bernacchi et al. (2001), converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.
  ko25 <- 27480   # Pa, value based on Bernacchi et al. (2001), converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.
  
  ## conversion to Kelvin
  tk <- tc + 273.15
  
  kc <- kc25 * ftemp_arrh( tk, dha=dhac )
  ko <- ko25 * ftemp_arrh( tk, dha=dhao )
  
  po  <- kco * (1e-6) * patm         # O2 partial pressure
  kmm <- kc * (1.0 + po/ko)
  
  return(kmm)
}


#' Calculates the CO2 compensation point
#'
#' Calculates the photorespiratory CO2 compensation point in absence of dark
#' respiration, \eqn{\Gamma*} (Farquhar, 1980).
#'
#' @param tc Temperature, relevant for photosynthesis (degrees Celsius)
#' @param patm Atmospheric pressure (Pa)
#'
#' @details The temperature and pressure-dependent photorespiratory
#' compensation point in absence of dark respiration \eqn{\Gamma* (T,p)}
#' is calculated from its value at standard temperature (\eqn{T0 = 25 }deg C)
#' and atmospheric pressure (\eqn{p0 = 101325} Pa), referred to as \eqn{\Gamma*0},
#' quantified by Bernacchi et al. (2001) to 4.332 Pa (their value in molar
#' concentration units is multiplied here with 101325 Pa to yield 4.332 Pa).
#' \eqn{\Gamma*0} is modified by temperature following an Arrhenius-type temperature
#' response function \eqn{f(T, \Delta Ha)} (implemented by \link{ftemp_arrh})
#' with activation energy \eqn{\Delta Ha = 37830} J mol-1  and is corrected for
#' atmospheric pressure \eqn{p(z)} (see \link{calc_patm}) at elevation \eqn{z}.
#' \deqn{
#'       \Gamma* = \Gamma*0 f(T, \Delta Ha) p(z) / p_0
#' }
#' \eqn{p(z)} is given by argument \code{patm}.
#'
#' @references Farquhar,  G.  D.,  von  Caemmerer,  S.,  and  Berry,  J.  A.:
#'             A  biochemical  model  of photosynthetic CO2 assimilation in leaves of
#'             C 3 species, Planta, 149, 78–90, 1980.
#'
#'             Bernacchi,  C.  J.,  Singsaas,  E.  L.,  Pimentel,  C.,  Portis,  A.
#'             R.  J.,  and  Long,  S.  P.:Improved temperature response functions
#'             for models of Rubisco-limited photosyn-thesis, Plant, Cell and
#'             Environment, 24, 253–259, 2001
#'
#' @return A numeric value for \eqn{\Gamma*} (in Pa)
#'
#' @examples print("CO2 compensation point at 20 degrees Celsius and standard atmosphere (in Pa):")
#' print(calc_gammastar(20, 101325))
#'
#' @export
#'
calc_gammastar <- function( tc, patm ) {
  
  # (J/mol) Activation energy, Bernacchi et al. (2001)
  dha    <- 37830
  
  # Pa, value based on Bernacchi et al. (2001), 
  # converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.
  gs25_0 <- 4.332
  
  gammastar <- gs25_0 * patm / calc_patm(0.0) * ftemp_arrh( (tc + 273.15), dha=dha )
  
  return( gammastar )
}


#' Calculates the temperature dependence of the quantum yield efficiency
#'
#' Calculates the temperature dependence of the quantum yield efficiency
#' following the temperature dependence of the maximum quantum yield of photosystem II
#' in light-adapted tobacco leaves, determined by Bernacchi et al. (2003)
#'
#' @param tc Temperature, relevant for photosynthesis (degrees Celsius)
#' @param c4 Boolean specifying whether fitted temperature response for C4 plants
#' is used. Defaults to \code{FALSE} (C3 photoynthesis temperature resposne following
#' Bernacchi et al., 2003 is used).
#'
#' @details The temperature factor for C3 photosynthesis (argument \code{c4 = FALSE}) is calculated
#' based on Bernacchi et al. (2003) as
#' 			\deqn{
#' 				\phi(T) = 0.352 + 0.022 T - 0.00034 T^2
#'       }
#'
#' The temperature factor for C4 (argument \code{c4 = TRUE}) photosynthesis is calculated based on
#' pers. comm. by David Orme, correcting values provided in Cai & Prentice (2020). Corrected 
#' parametrisation is:
#' 			\deqn{
#' 				\phi(T) = -0.064 + 0.03 T - 0.000464 T^2 
#'       }
#'
#' The factor \eqn{\phi(T)} is to be multiplied with leaf absorptance and the fraction
#' of absorbed light that reaches photosystem II. In the P-model these additional factors
#' are lumped into a single apparent quantum yield efficiency parameter (argument \code{kphio}
#' to function \link{rpmodel}).
#'
#' @return A numeric value for \eqn{\phi(T)}
#'
#' @examples
#' ## Relative change in the quantum yield efficiency
#' ## between 5 and 25 degrees celsius (percent change):
#' print(paste((ftemp_kphio(25.0)/ftemp_kphio(5.0)-1)*100 ))
#'
#' @references  
#' Bernacchi, C. J., Pimentel, C., and Long, S. P.:  In vivo temperature
#' 				response func-tions  of  parameters required  to  model  RuBP-limited
#' 				photosynthesis,  Plant  Cell Environ., 26, 1419–1430, 2003
#' Cai, W., and Prentice, I. C.: Recent trends in gross primary production 
#'        and their drivers: analysis and modelling at flux-site and global scales,
#'        Environ. Res. Lett. 15 124050 https://doi.org/10.1088/1748-9326/abc64e, 2020
#'
#' @export
#'
ftemp_kphio <- function( tc, c4 = FALSE ){
  
  if (c4){
    ftemp = -0.064 + 0.03 * tc - 0.000464 * tc^2     # correcting erroneous values provided in Cai & Prentice, 2020, according to D. Orme (issue #19) 
    # XXX THIS IS NOT CORRECT: ftemp = -0.008 + 0.00375 * tc - 0.58e-4 * tc^2   # Based on calibrated values by Shirley
  } else {
    ftemp <- 0.352 + 0.022 * tc - 3.4e-4 * tc^2
  }
  
  ## avoid negative values
  ftemp <- ifelse(ftemp < 0.0, 0.0, ftemp)
  
  return(ftemp)
}


#' Calculates the temperature response of dark respiration
#'
#' Given the dark respiration at the reference temperature 25 degress Celsius,
#' this function calculates its temperature-scaling factor following Heskel et al. 2016.
#'
#' @param tc Temperature (degrees Celsius)
#'
#' @details To correct for effects by temperature Heskel et al. 2016,
#' and given the reference temperature \eqn{T0 =} 25 deg C, this calculates the temperature
#' scaling factor to calculate dark respiration at temperature \eqn{T} (argument \code{tc}) as:
#' \deqn{
#'      fr = exp( 0.1012 (T0-T) - 0.0005(T0^2 - T^2) ) 
#' }
#' where \eqn{T} is given in degrees Celsius.
#' 
#' @return A numeric value for \eqn{fr}
#' 
#' @references  Heskel,  M.,  O’Sullivan,  O.,  Reich,  P.,  Tjoelker,  M.,  Weerasinghe,  L.,  Penillard,  A.,Egerton, J., 
#'              Creek, D., Bloomfield, K., Xiang, J., Sinca, F., Stangl, Z., Martinez-De La Torre, A., Griffin, K., 
#'              Huntingford, C., Hurry, V., Meir, P., Turnbull, M.,and Atkin, O.:  Convergence in the temperature response 
#'              of leaf respiration across biomes and plant functional types, Proceedings of the National Academy of Sciences,
#'              113,  3832–3837,  doi:10.1073/pnas.1520282113,2016.
#'
#' @examples 
#' ## Relative change in Rd going (instantaneously, i.e. not 
#' ## acclimatedly) from 10 to 25 degrees (percent change):
#' print( (ftemp_inst_rd(25)/ftemp_inst_rd(10)-1)*100 )
#'
#' @export
#'
ftemp_inst_rd <- function( tc ){
  
  # loal parameters
  apar <- 0.1012
  bpar <- 0.0005
  
  fr <- exp( apar * (tc - 25.0) - bpar * (tc^2 - 25.0^2) )
  
  return(fr)
}



#' Calculates the instantaneous temperature response of Vcmax
#'
#' Given Vcmax at a reference temperature (argument \code{tcref})
#' this function calculates its temperature-scaling factor following modified Arrhenius
#' kinetics based on Kattge & Knorr (2007). Calculates \eqn{f} for the conversion
#' \deqn{
#'    V = f Vref
#' }
#'
#' @param tcleaf Leaf temperature, or in general the temperature relevant for photosynthesis
#' (degrees Celsius)
#' @param tcgrowth (Optional) Growth temperature, in the P-model, taken to be equal to \code{tcleaf}
#' (in degrees Celsius). Defaults to \code{tcgrowth = tcleaf}.
#' @param tcref Reference temperature (in degrees Celsius)
#'
#' @details The function is given by Kattge & Knorr (2007) as
#' \deqn{
#' 		fv = f(T, \Delta Hv) A/B
#' }
#' where \eqn{f(T, \Delta Hv)} is a regular Arrhenius-type temperature response function (see
#' \link{ftemp_arrh}) with \eqn{Hv=71513} J mol-1,
#' \deqn{
#' 		A = 1 + exp( (T0 \Delta S - Hd) / (T0 R) )
#' }
#' and
#' \deqn{
#' 	    B = 1 + exp( (T \Delta S - Hd) / (TK R) )
#' }
#' Here, \eqn{T} is in Kelvin, \eqn{T0=293.15} K, \eqn{Hd = 200000} J mol-1 is the deactivation
#' energy and \eqn{R} is the universal gas constant and is 8.3145 J mol-1 K-1, and
#' \deqn{
#' 		\Delta S = aS - bS T
#' }
#' with \eqn{aS = 668.39} J mol-1 K-1, and \eqn{bS = 1.07} J mol-1 K-2, and \eqn{T} given in
#' degrees Celsius (!)
#'
#' @references Kattge, J. and Knorr, W.:  Temperature acclimation in a biochemical model of
#' photosynthesis: a reanalysis of data from 36 species, Plant, Cell and Environment, 30,1176–1190, 2007.
#'
#' @return A numeric value for \eqn{fv}
#'
#' @examples 
#' ## Relative change in Vcmax going (instantaneously, i.e. 
#' ## not acclimatedly) from 10 to 25 degrees (percent change):
#' print((ftemp_inst_vcmax(25)/ftemp_inst_vcmax(10)-1)*100 )
#'
#' @export

ftemp_inst_vcmax <- function(
  tcleaf,
  tcgrowth = tcleaf,
  tcref = 25.0
){
  
  # loal parameters
  Ha    <- 71513  # activation energy (J/mol)
  Hd    <- 200000 # deactivation energy (J/mol)
  Rgas  <- 8.3145 # universal gas constant (J/mol/K)
  a_ent <- 668.39 # offset of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K)
  b_ent <- 1.07   # slope of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K^2)
  
  tkref <- tcref + 273.15  # to Kelvin
  
  # conversion of temperature to Kelvin, tcleaf is the instantaneous leaf temperature in degrees C.
  tkleaf <- tcleaf + 273.15
  
  # calculate entropy following Kattge & Knorr (2007), negative slope and y-axis intersect is when expressed as a function of temperature in degrees Celsius, not Kelvin !!!
  dent <- a_ent - b_ent * tcgrowth   # 'tcgrowth' corresponds to 'tmean' in Nicks, 'tc25' is 'to' in Nick's
  fva <- ftemp_arrh( tkleaf, Ha, tkref = tkref )
  fvb <- (1 + exp( (tkref * dent - Hd)/(Rgas * tkref) ) ) / (1 + exp( (tkleaf * dent - Hd)/(Rgas * tkleaf) ) )
  fv  <- fva * fvb
  
  return( fv )
}


#' Calculates the instantaneous temperature response of Jmax
#'
#' Given Jmax at a reference temperature (argument \code{tcref})
#' this function calculates its temperature-scaling factor following
#' modified Arrhenius kinetics based on Kattge & Knorr (2007). 
#' Calculates \eqn{f} for the conversion
#' \deqn{
#'    V = f Vref
#' }
#'
#' @param tcleaf Leaf temperature, or in general the temperature relevant for 
#'  photosynthesis (degrees Celsius)
#' @param tcgrowth (Optional) Growth temperature, in the P-model, taken to be
#'  equal to \code{tcleaf} (in degrees Celsius). Defaults to 
#'  \code{tcgrowth = tcleaf}.
#' @param tcref Reference temperature (in degrees Celsius)
#'
#' @details The function is given by Kattge & Knorr (2007) as
#' \deqn{
#' 		fv = f(T, \Delta Hv) A/B
#' }
#' where \eqn{f(T, \Delta Hv)} is a regular Arrhenius-type temperature response
#'  function (see \link{ftemp_arrh}) with \eqn{Hv=49884} J mol-1,
#' \deqn{
#' 		A = 1 + exp( (T0 \Delta S - Hd) / (T0 R) )
#' }
#' and
#' \deqn{
#' 	    B = 1 + exp( (T \Delta S - Hd) / (TK R) )
#' }
#' Here, \eqn{T} is in Kelvin, \eqn{T0=293.15} K, \eqn{Hd = 200000} J mol-1 is
#' the deactivation energy and \eqn{R} is the universal gas constant and is 
#' 8.3145 J mol-1 K-1, and
#' \deqn{
#' 		\Delta S = aS - bS T
#' }
#' with \eqn{aS = 659.70} J mol-1 K-1, and \eqn{bS = 0.75} J mol-1 K-2, and
#'  \eqn{T} given in degrees Celsius (!)
#'
#' @references Kattge, J. and Knorr, W.:  Temperature acclimation in a 
#' biochemical model of photosynthesis: a reanalysis of data from 36 species,
#' Plant, Cell and Environment, 30,1176–1190, 2007.
#'
#' @return A numeric value for \eqn{fv}
#'
#' @examples
#' # Relative change in Jmax going (instantaneously, i.e.
#' # not acclimatedly) from 10 to 25 degrees (percent change):
#' print((ftemp_inst_jmax(25)/ftemp_inst_jmax(10)-1)*100 )
#'
#' @export

ftemp_inst_jmax <- function(
  tcleaf,
  tcgrowth = tcleaf,
  tcref = 25.0
){
  
  # local parameters
  Ha    <- 49884  # activation energy (J/mol)
  Hd    <- 200000 # deactivation energy (J/mol)
  Rgas  <- 8.3145 # universal gas constant (J/mol/K)
  
  # offset and slope (a / b) of entropy vs. temperature relationship 
  # from Kattge & Knorr (2007) (J/mol/K)
  a_ent <- 659.70 
  b_ent <- 0.75   
  
  # conversion to Kelvin
  tkref <- tcref + 273.15
  
  # conversion of temperature to Kelvin, 
  # tcleaf is the instantaneous leaf temperature in degrees C.
  tkleaf <- tcleaf + 273.15
  
  # calculate entropy following Kattge & Knorr (2007), negative slope and 
  # y-axis intersect is when expressed as a function of temperature in
  # degrees Celsius, not Kelvin !!!
  
  # 'tcgrowth' corresponds to 'tmean' in Nicks, 'tc25' is 'to' in Nick's
  dent <- a_ent - b_ent * tcgrowth   
  fva <- ftemp_arrh( tkleaf, Ha, tkref = tkref )
  fvb <- (1 + exp( (tkref * dent - Hd)/(Rgas * tkref) ) ) / 
    (1 + exp( (tkleaf * dent - Hd)/(Rgas * tkleaf) ) )
  fv  <- fva * fvb
  
  return( fv )
}


#' Viscosity of water
#'
#' Calculates the viscosity of water as a function of temperature and atmospheric
#' pressure.
#'
#' @param tc numeric, air temperature (tc), degrees C
#' @param p numeric, atmospheric pressure (p), Pa
#'
#' @return numeric, viscosity of water (mu), Pa s
#'
#' @examples print("Density of water at 20 degrees C and standard atmospheric pressure:")
#' print(density_h2o(20, 101325))
#'
#' @references  Huber, M. L., R. A. Perkins, A. Laesecke, D. G. Friend, J. V.
#' Sengers, M. J. Assael, ..., K. Miyagawa (2009) New
#' international formulation for the viscosity of H2O, J. Phys.
#' Chem. Ref. Data, Vol. 38(2), pp. 101-125.
#'
#' @export
#'
viscosity_h2o <- function(tc, p) {
  
  # Define reference temperature, density, and pressure values:
  tk_ast  <- 647.096    # Kelvin
  rho_ast <- 322.0      # kg/m^3
  mu_ast  <- 1e-6       # Pa s
  
  # Get the density of water, kg/m^3
  rho <- density_h2o(tc, p)
  
  # Calculate dimensionless parameters:
  tbar  <- (tc + 273.15)/tk_ast
  tbarx <- tbar^(0.5)
  tbar2 <- tbar^2
  tbar3 <- tbar^3
  rbar  <- rho/rho_ast
  
  # Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
  mu0 <- 1.67752 + 2.20462/tbar + 0.6366564/tbar2 - 0.241605/tbar3
  mu0 <- 1e2*tbarx/mu0
  
  # Create Table 3, Huber et al. (2009):
  h_array <- array(0.0, dim=c(7,6))
  h_array[1,] <- c(0.520094, 0.0850895, -1.08374, -0.289555, 0.0, 0.0)  # hj0
  h_array[2,] <- c(0.222531, 0.999115, 1.88797, 1.26613, 0.0, 0.120573) # hj1
  h_array[3,] <- c(-0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.0) # hj2
  h_array[4,] <- c(0.161913,  0.257399, 0.0, 0.0, 0.0, 0.0) # hj3
  h_array[5,] <- c(-0.0325372, 0.0, 0.0, 0.0698452, 0.0, 0.0) # hj4
  h_array[6,] <- c(0.0, 0.0, 0.0, 0.0, 0.00872102, 0.0) # hj5
  h_array[7,] <- c(0.0, 0.0, 0.0, -0.00435673, 0.0, -0.000593264) # hj6
  
  # Calculate mu1 (Eq. 12 & Table 3, Huber et al., 2009):
  mu1 <- 0.0
  ctbar <- (1.0/tbar) - 1.0
  # print(paste("ctbar",ctbar))
  # for i in xrange(6):
  for (i in 1:6){
    coef1 <- ctbar^(i-1)
    # print(paste("i, coef1", i, coef1))
    coef2 <- 0.0
    for (j in 1:7){
      coef2 <- coef2 + h_array[j,i] * (rbar - 1.0)^(j-1)
    }
    mu1 <- mu1 + coef1 * coef2
  }
  mu1 <- exp( rbar * mu1 )
  # print(paste("mu1",mu1))
  
  # Calculate mu_bar (Eq. 2, Huber et al., 2009)
  #   assumes mu2 = 1
  mu_bar <- mu0 * mu1
  
  # Calculate mu (Eq. 1, Huber et al., 2009)
  mu <- mu_bar * mu_ast    # Pa s
  
  return( mu )
}


#' Calculates the Arrhenius-type temperature response
#'
#' Given a kinetic rate at a reference temperature (argument \code{tkref})
#' this function calculates its temperature-scaling factor
#' following Arrhenius kinetics.
#'
#' @param tk Temperature (Kelvin)
#' @param dha Activation energy (J mol-1)
#' @param tkref Reference temperature (Kelvin)
#'
#' @details To correct for effects by temperature following Arrhenius kinetics,
#' and given a reference temperature \eqn{T_0}, \eqn{f} calculates the temperature
#' scaling. Arrhenius kinetics are described by an equation of form
#' \eqn{x(T)= exp(c - \Delta H_a / (T R))}. The temperature-correction
#' function \eqn{f(T, \Delta H_a)} is thus given by \eqn{f=x(T)/x(T_0)} which is:
#' \deqn{
#'      f = exp( \Delta H_a (T - T_0) / (T_0 R T_K) )
#' }
#' \eqn{\Delta H_a} is given by argument \code{dha}. \eqn{T} is given by argument
#' \code{tk} and has to be provided in Kelvin. \eqn{R} is the universal gas constant
#' and is 8.3145 J mol-1 K-1. Note that this is equivalent to 
#' \deqn{
#'      f = exp( (\Delta H_a/R) (1/T_0 - 1/T) )  
#' }
#'
#' @return A numeric value for \eqn{f}
#'
#' @examples 
#' # Relative rate change from 25 to 10 degrees Celsius (percent change)
#' print( (1.0-ftemp_arrh( 283.15, 100000, tkref = 298.15))*100 )
#'
#' @export

ftemp_arrh <- function(
  tk,
  dha,
  tkref = 298.15
){
  
  # Note that the following forms are equivalent:
  # ftemp = exp( dha * (tk - 298.15) / (298.15 * kR * tk) )
  # ftemp = exp( dha * (tc - 25.0)/(298.15 * kR * (tc + 273.15)) )
  # ftemp = exp( (dha/kR) * (1/298.15 - 1/tk) )
  
  kR   <- 8.3145     # Universal gas constant, J/mol/K
  ftemp <- exp( dha * (tk - tkref) / (tkref * kR * tk) )
  
  return(ftemp)
}


#' Density of water
#'
#' Calculates the density of water as a function of temperature and atmospheric
#' pressure, using the Tumlirz Equation.
#'
#' @param tc numeric, air temperature (tc), degrees C
#' @param p numeric, atmospheric pressure (p), Pa
#'
#' @return numeric, density of water, kg/m^3
#'
#' @examples 
#'  # Density of water at 20 degrees C and standard atmospheric pressure
#'  print(density_h2o(20, 101325))
#'
#' @references  F.H. Fisher and O.E Dial, Jr. (1975) Equation of state of
#' pure water and sea water, Tech. Rept., Marine Physical
#' Laboratory, San Diego, CA.
#'
#' @export
#'
density_h2o <- function(tc, p){
  
  # Calculate lambda, (bar cm^3)/g:
  my_lambda <- 1788.316 +
    21.55053*tc +
    -0.4695911*tc*tc +
    (3.096363e-3)*tc*tc*tc +
    -(7.341182e-6)*tc*tc*tc*tc
  
  # Calculate po, bar
  po <- 5918.499 +
    58.05267*tc +
    -1.1253317*tc*tc +
    (6.6123869e-3)*tc*tc*tc +
    -(1.4661625e-5)*tc*tc*tc*tc
  
  # Calculate vinf, cm^3/g
  vinf <- 0.6980547 +
    -(7.435626e-4)*tc +
    (3.704258e-5)*tc*tc +
    -(6.315724e-7)*tc*tc*tc +
    (9.829576e-9)*tc*tc*tc*tc +
    -(1.197269e-10)*tc*tc*tc*tc*tc +
    (1.005461e-12)*tc*tc*tc*tc*tc*tc +
    -(5.437898e-15)*tc*tc*tc*tc*tc*tc*tc +
    (1.69946e-17)*tc*tc*tc*tc*tc*tc*tc*tc +
    -(2.295063e-20)*tc*tc*tc*tc*tc*tc*tc*tc*tc
  
  # Convert pressure to bars (1 bar <- 100000 Pa)
  pbar <- (1e-5)*p
  
  # Calculate the specific volume (cm^3 g^-1):
  v <- vinf + my_lambda/(po + pbar)
  
  # Convert to density (g cm^-3) -> 1000 g/kg; 1000000 cm^3/m^3 -> kg/m^3:
  rho <- (1e3/v)
  
  return(rho)
}

#' CO2 partial pressure
#'
#' Calculates CO2 partial pressure from concentration in ppm.
#'
#' @param co2 Atmospheric CO2 concentration (ppm)
#' @param patm Atmospheric pressure (Pa).
#'
#'
#' @return CO2 partial pressure in Pa.
#'
#' @export
#'
co2_to_ca <- function(co2, patm){
  # Pa, atms. CO2
  ca   <- ( 1.0e-6 ) * co2 * patm
  return( ca )
}

optimal_chi <- function(kmm, gammastar, ns_star, ca, vpd, beta, c4){
  
  # Input:    - float, 'kmm' : Pa, Michaelis-Menten coeff.
  #           - float, 'ns_star'  : (unitless) viscosity correction factor for water
  #           - float, 'vpd' : Pa, vapor pressure deficit
  # Output:   float, ratio of ci/ca (chi)
  # Features: Returns an estimate of leaf internal to ambient CO2
  #           partial pressure following the "simple formulation".
  # Depends:  - kc
  #           - ns
  #           - vpd
  
  ## Avoid negative VPD (dew conditions), resolves issue #2 (https://github.com/stineb/rpmodel/issues/2)
  vpd <- ifelse(vpd < 0, 0, vpd)
  
  ## leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio
  xi  <- sqrt( (beta * ( kmm + gammastar ) ) / ( 1.6 * ns_star ) )
  chi <- gammastar / ca + ( 1.0 - gammastar / ca ) * xi / ( xi + sqrt(vpd) )
  
  if (c4){
    
    out <- list(
      xi = xi,
      chi = chi,
      mc = 1.0,
      mj = 1.0,
      mjoc = 1.0
    )    
    
  } else {

    ## alternative variables
    gamma <- gammastar / ca
    kappa <- kmm / ca
    
    ## use chi for calculating mj
    mj <- (chi - gamma) / (chi + 2.0 * gamma)
    
    ## mc
    mc <- (chi - gamma) / (chi + kappa)
    
    ## mj:mv
    mjoc <- (chi + kappa) / (chi + 2.0 * gamma)
    
    # format output list
    out <- list(
      xi = xi,
      chi = chi,
      mc = mc,
      mj = mj,
      mjoc = mjoc
    )
    
  }
  
  return(out)
}

# ## wrap if condition in a function to allow vectorization
# mj <- function(ns_star, vpd, vacg, vbkg, vdcg, gammastar){
#
#   vsr <- sqrt( 1.6 * ns_star * vpd / vbkg )
#
#   # Based on the mc' formulation (see Regressing_LUE.pdf)
#   mj <- vdcg / ( vacg + 3.0 * gammastar * vsr )
#
#   return(mj)
# }

lue_vcmax_wang17 <- function(out_optchi, kphio, c_molmass, soilmstress){
  
  ## Include effect of Jmax limitation
  len <- length(out_optchi[[1]])

  kc <- 0.41   # Jmax cost coefficient

  # ## Following eq. 17 in Stocker et al., 2020 GMD
  # tmp <- 1.0 - (kc / out_optchi$mj)^(2.0/3.0)
  # mprime <- ifelse(tmp > 0, out_optchi$mj * sqrt(tmp), NA)  # avoid square root of negative number

  ## original code - is equivalent to eq. 17-based formulation above
  tmp <- out_optchi$mj^2 - kc^(2.0/3.0) * (out_optchi$mj^(4.0/3.0))
  mprime <- ifelse(tmp > 0, sqrt(tmp), NA)  # avoid square root of negative number
  
  out <- list(
    
    mprime = mprime,
    
    ## Light use efficiency (gpp per unit absorbed light)
    lue = kphio * mprime * c_molmass * soilmstress,
    
    ## Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
    vcmax_unitiabs = soilmstress * kphio * mprime / out_optchi$mc,
    
    ## complement for non-smith19
    omega      = rep(NA, len),
    omega_star = rep(NA, len)
    
  )
  
  return(out)
}


lue_vcmax_smith19 <- function(out_optchi, kphio, c_molmass, soilmstress){
  
  len <- length(out_optchi[[1]])
  
  # Adopted from Nick Smith's code:
  # Calculate omega, see Smith et al., 2019 Ecology Letters
  calc_omega <- function( theta, c_cost, m ){
    
    cm <- 4 * c_cost / m                        # simplification term for omega calculation
    v  <- 1/(cm * (1 - theta * cm)) - 4 * theta # simplification term for omega calculation
    
    # account for non-linearities at low m values
    capP <- (((1/1.4) - 0.7)^2 / (1-theta)) + 3.4
    aquad <- -1
    bquad <- capP
    cquad <- -(capP * theta)
    m_star <- (4 * c_cost) / polyroot(c(aquad, bquad, cquad))
    
    out_omega <- ifelse(  m < Re(m_star[1]),
                      -( 1 - (2 * theta) ) - sqrt( (1 - theta) * v),
                      -( 1 - (2 * theta))  + sqrt( (1 - theta) * v)
    )
    return(out_omega)
  }
  
  ## constants
  theta <- 0.85    # should be calibratable?
  c_cost <- 0.05336251
  
  
  ## factors derived as in Smith et al., 2019
  omega <- calc_omega( theta = theta, c_cost = c_cost, m = out_optchi$mj )          # Eq. S4
  omega_star <- 1.0 + omega - sqrt( (1.0 + omega)^2 - (4.0 * theta * omega) )       # Eq. 18
  
  ## Effect of Jmax limitation
  mprime <- out_optchi$mj * omega_star / (8.0 * theta)
  
  ## Light use efficiency (gpp per unit absorbed light)
  lue <- kphio * mprime * c_molmass * soilmstress
  
  # calculate Vcmax per unit aborbed light
  vcmax_unitiabs  <- kphio * out_optchi$mjoc * omega_star / (8.0 * theta) * soilmstress   # Eq. 19
  
  out <- list(
    lue            = lue,
    vcmax_unitiabs = vcmax_unitiabs,
    omega          = omega,
    omega_star     = omega_star
  )
  
  return(out)
}


lue_vcmax_none <- function(out_optchi, kphio, c_molmass, soilmstress){
  ## Do not include effect of Jmax limitation
  len <- length(out_optchi[[1]])
  
  out <- list(
    
    ## Light use efficiency (gpp per unit absorbed light)
    lue = kphio * out_optchi$mj * c_molmass * soilmstress,
    
    ## Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
    vcmax_unitiabs = kphio * out_optchi$mjoc * soilmstress,
    
    ## complement for non-smith19
    omega               = rep(NA, len),
    omega_star          = rep(NA, len)
  )
  
  return(out)
}


lue_vcmax_c4 <- function( kphio, c_molmass, soilmstress ){
  
  len <- length(kphio)
  out <- list(
    ## Light use efficiency (gpp per unit absorbed light)
    lue = kphio * c_molmass * soilmstress,
    
    ## Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
    vcmax_unitiabs = kphio * soilmstress,
    
    ## complement for non-smith19
    omega               = rep(NA, len),
    omega_star          = rep(NA, len)
  )
  
  return(out)
}

## #' Larger quadratic root
## #' 
## #' Solves quadratic equation given by y = a*x^2 + bx + c
## #  Based on MAESTRA equivalent (B. Medlyn)
## #'
## #' @param A a parameter in the quadratic equation
## #' @param B b parameter in the quadratic equation
## #' @param C c parameter in the quadratic equation
## #'
## #' @return larger quardratic root
## #' @export

QUADP <- function(A,B,C){
  
  if (any(is.na(c(A,B,C)))){
    return(NA)
  } else {
    if((B^2 - 4*A*C) < 0){
      warning("IMAGINARY ROOTS IN QUADRATIC")
      return(0)
    }
    
    if(identical(A,0)){
      if(identical(B,0)){
        return(0)
      } else {
        return(-C/B)
      }
    } else {
      return((- B + sqrt(B^2 - 4*A*C)) / (2*A))
    }
  }
  
}

## #' Minor quadratic root
## #' 
## #' Solves quadratic equation given by y = a*x^2 + bx + c
## #' #  Based on MAESTRA equivalent (B. Medlyn)
## #' @param A a parameter in the quadratic equation
## #' @param B b parameter in the quadratic equation
## #' @param C c parameter in the quadratic equation
## #' 
## #' @return minor quardratic root
## #' @export

QUADM <- function(A,B,C){
  
  if (any(is.na(c(A,B,C)))){
    return(NA)
  } else {
    if((B^2 - 4*A*C) < 0){
      warning("IMAGINARY ROOTS IN QUADRATIC")
      return(0)
    }
    
    if(identical(A,0)){
      if(identical(B,0)){
        return(0)
      } else {
        return(-C/B)
      }
    } else {
      return((- B - sqrt(B^2 - 4*A*C)) / (2*A))
    }
  }
  
}

